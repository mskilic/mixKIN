from __future__ import annotations

from itertools import combinations
from typing import Dict, List, Tuple

import numpy as np
import pandas as pd


RELATION_THRESHOLDS = [
    (0.03125, "Unrelated"),
    (0.08838835, "Third"),
    (0.1767767, "Second"),
    (0.3535534, "First"),
]


def classify_relationship(theta: float) -> str:
    if np.isnan(theta):
        return "NA"
    if theta < RELATION_THRESHOLDS[0][0]:
        return "Unrelated"
    if theta < RELATION_THRESHOLDS[1][0]:
        return "Third"
    if theta < RELATION_THRESHOLDS[2][0]:
        return "Second"
    if theta < RELATION_THRESHOLDS[3][0]:
        return "First"
    return "MZT"


def _source_callrate(genotypes: np.ndarray, sample_indices: np.ndarray) -> np.ndarray:
    if sample_indices.size == 0:
        return np.zeros(genotypes.shape[0], dtype=float)
    subset = genotypes[:, sample_indices]
    return np.mean(~np.isnan(subset), axis=1)


def estimate_source_freqs(
    genotypes: np.ndarray,
    sample_ids: List[str],
    populations: List[str],
    source_pops: List[str],
) -> Tuple[pd.DataFrame, Dict[str, int], Dict[str, np.ndarray]]:
    pop_arr = np.array(populations, dtype=object)
    source_sizes: Dict[str, int] = {}
    callrates: Dict[str, np.ndarray] = {}
    freq_cols: Dict[str, np.ndarray] = {}

    for pop in source_pops:
        idx = np.where(pop_arr == pop)[0]
        if idx.size == 0:
            raise ValueError(f"No samples found for source population '{pop}'")
        source_sizes[pop] = int(idx.size)
        source_gt = genotypes[:, idx]
        freq_cols[pop] = np.nanmean(source_gt / 2.0, axis=1)
        callrates[pop] = _source_callrate(genotypes, idx)

    df = pd.DataFrame(freq_cols)
    return df, source_sizes, callrates


def _correction_factor(size: float, ploidy: str, avg_callrate: float) -> float:
    effective_n = max(size * avg_callrate, 1.01)
    ploidy_factor = 2.0 if str(ploidy).lower() == "diploid" else 1.0
    scale = ploidy_factor * effective_n
    if scale <= 1.0:
        return 1.0
    return scale / (scale - 1.0)


def run_kinship(
    target_ids: List[str],
    sample_ids: List[str],
    target_genotypes: np.ndarray,
    target_sample_indices: List[int],
    target_snp_ids: List[str],
    source_freqs: pd.DataFrame,
    source_pops: List[str],
    q_df: pd.DataFrame,
    source_meta: Dict[str, Dict[str, object]],
    source_sizes: Dict[str, int],
    source_callrates: Dict[str, np.ndarray],
    default_source_size: float,
    default_source_ploidy: str = "pseudohaploid",
    epmr_mode: str = "overlap",
) -> pd.DataFrame:
    source_mat = source_freqs[source_pops].to_numpy(dtype=float)
    valid_source = np.all(~np.isnan(source_mat), axis=1)

    gt = target_genotypes[valid_source, :]
    source_mat = source_mat[valid_source, :]

    if gt.shape[0] == 0:
        raise ValueError("No SNPs left after filtering for source-population frequency availability")

    q = q_df.set_index("ID").reindex(target_ids)
    if q[source_pops].isna().any().any():
        missing_ids = q[q[source_pops].isna().any(axis=1)].index.tolist()
        raise ValueError(f"Missing ancestry proportions for targets: {', '.join(missing_ids)}")

    q_mat = q[source_pops].to_numpy(dtype=float)

    avg_callrates = {}
    for p_idx, pop in enumerate(source_pops):
        if pop in source_callrates:
            avg_callrates[pop] = float(np.nanmean(source_callrates[pop][valid_source]))
        else:
            avg_callrates[pop] = float(np.mean(~np.isnan(source_mat[:, p_idx])))

    within = {}
    for p_idx, pop in enumerate(source_pops):
        p = source_mat[:, p_idx]
        size = float(source_meta.get(pop, {}).get("size", source_sizes.get(pop, default_source_size)))
        ploidy = str(source_meta.get(pop, {}).get("ploidy", default_source_ploidy))
        corr = _correction_factor(size=size, ploidy=ploidy, avg_callrate=avg_callrates[pop])
        within[pop] = corr * (2.0 * p * (1.0 - p))

    between = {}
    for a_idx, a in enumerate(source_pops):
        for b_idx, b in enumerate(source_pops):
            if a_idx == b_idx:
                continue
            pa = source_mat[:, a_idx]
            pb = source_mat[:, b_idx]
            between[(a, b)] = pa * (1.0 - pb) + (1.0 - pa) * pb

    rows = []
    sample_index = {sid: i for sid, i in zip(target_ids, target_sample_indices)}

    for id1, id2 in combinations(target_ids, 2):
        i1 = sample_index[id1]
        i2 = sample_index[id2]

        g1 = gt[:, i1]
        g2 = gt[:, i2]
        overlap = (~np.isnan(g1)) & (~np.isnan(g2))
        overlap_n = int(np.sum(overlap))
        if overlap_n == 0:
            rows.append(
                {
                    "ID1": id1,
                    "ID2": id2,
                    "Theta": np.nan,
                    "PMR": np.nan,
                    "EPMR": np.nan,
                    "OverlappingSNPs": 0,
                    "Relationship": "NA",
                }
            )
            continue

        g1_o = g1[overlap]
        g2_o = g2[overlap]

        mismatches = np.zeros_like(g1_o, dtype=float)
        mismatches[np.abs(g1_o - g2_o) == 2] = 1.0
        mismatches[(g1_o == 1) | (g2_o == 1)] = 0.5
        
        pmr = float(np.mean(mismatches))

        q1 = q.loc[id1, source_pops].to_numpy(dtype=float)
        q2 = q.loc[id2, source_pops].to_numpy(dtype=float)

        epmr_site = np.zeros(gt.shape[0], dtype=float)
        for a in source_pops:
            epmr_site += (q1[source_pops.index(a)] * q2[source_pops.index(a)]) * within[a]
        for a in source_pops:
            for b in source_pops:
                if a == b:
                    continue
                epmr_site += (q1[source_pops.index(a)] * q2[source_pops.index(b)]) * between[(a, b)]

        if epmr_mode == "all":
            epmr = float(np.mean(epmr_site))
        else:
            epmr = float(np.mean(epmr_site[overlap]))
            
        theta = 1.0 - (pmr / epmr) if epmr > 0 else np.nan

        rows.append(
            {
                "ID1": id1,
                "ID2": id2,
                "Theta": theta,
                "PMR": pmr,
                "EPMR": epmr,
                "OverlappingSNPs": overlap_n,
                "Relationship": classify_relationship(theta),
            }
        )

    return pd.DataFrame(rows)

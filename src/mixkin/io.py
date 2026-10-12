from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional

import numpy as np
import pandas as pd


@dataclass
class DatasetData:
    sample_ids: List[str]
    populations: List[str]
    snp_ids: List[str]
    genotypes: np.ndarray


def _read_eigenstrat(prefix: str) -> DatasetData:
    ind_path = Path(f"{prefix}.ind")
    snp_path = Path(f"{prefix}.snp")
    geno_path = Path(f"{prefix}.geno")

    ind_df = pd.read_csv(ind_path, sep=r"\s+", header=None, names=["id", "sex", "pop"])
    snp_df = pd.read_csv(
        snp_path,
        sep=r"\s+",
        header=None,
        names=["snp", "chr", "cM", "pos", "a1", "a2"],
    )

    sample_ids = ind_df["id"].astype(str).tolist()
    populations = ind_df["pop"].astype(str).tolist()
    snp_ids = snp_df["snp"].astype(str).tolist()

    n_ind = len(sample_ids)
    geno_rows: List[np.ndarray] = []
    with open(geno_path, "r", encoding="utf-8") as handle:
        for line in handle:
            row = line.strip()
            if not row:
                continue
            arr = np.frombuffer(row.encode("ascii"), dtype=np.uint8) - ord("0")
            if arr.size != n_ind:
                raise ValueError(f"EIGENSTRAT row length mismatch: expected {n_ind}, got {arr.size}")
            arr = arr.astype(float)
            arr[arr == 9] = np.nan
            arr = 2.0 - arr
            geno_rows.append(arr)

    genotypes = np.vstack(geno_rows) if geno_rows else np.empty((0, n_ind), dtype=float)
    return DatasetData(sample_ids=sample_ids, populations=populations, snp_ids=snp_ids, genotypes=genotypes)


def _read_plink(prefix: str) -> DatasetData:
    from pandas_plink import read_plink1_bin

    bed = Path(f"{prefix}.bed")
    bim = Path(f"{prefix}.bim")
    fam = Path(f"{prefix}.fam")

    xr = read_plink1_bin(str(bed), str(bim), str(fam), verbose=False)
    # Ensure variants are the first dimension
    if xr.dims[0] == "sample":
        xr = xr.transpose("variant", "sample")
    geno = xr.values.astype(float)

    sample_ids = [str(v) for v in xr.coords["iid"].values.tolist()] if "iid" in xr.coords else [str(v) for v in xr.coords["sample"].values.tolist()]
    snp_ids = [str(v) for v in xr.coords["snp"].values.tolist()] if "snp" in xr.coords else [str(v) for v in xr.coords["variant"].values.tolist()]

    fam_df = pd.read_csv(fam, sep=r"\s+", header=None)
    populations = fam_df.iloc[:, 0].astype(str).tolist()

    return DatasetData(sample_ids=sample_ids, populations=populations, snp_ids=snp_ids, genotypes=geno)


def read_dataset(prefix: str, dataset_type: str) -> DatasetData:
    if dataset_type == "eigenstrat":
        return _read_eigenstrat(prefix)
    if dataset_type == "plink":
        return _read_plink(prefix)
    raise ValueError(f"Unsupported dataset type: {dataset_type}")


def read_id_file(path: Optional[str]) -> Optional[List[str]]:
    if not path:
        return None
    values = pd.read_csv(path, header=None).iloc[:, 0].astype(str).tolist()
    return [v for v in values if v]


def read_q_file(path: str) -> tuple[pd.DataFrame, list[str]]:
    q = pd.read_csv(path, sep=None, engine="python")
    
    col_names = list(q.columns)
    id_col = col_names[0]
    rename = {id_col: "ID"}
    q = q.rename(columns=rename)
    
    source_pops = list(q.columns)[1:]

    q["ID"] = q["ID"].astype(str)
    raw_vals = q[source_pops].apply(pd.to_numeric, errors="coerce").fillna(0.0)

    out_of_bounds = (raw_vals < 0.0) | (raw_vals > 1.0)
    if out_of_bounds.any().any():
        num_oob = int(out_of_bounds.any(axis=1).sum())
        print(f"Warning: Found {num_oob} individual(s) with ancestry proportions outside [0, 1]. Values have been clipped to [0, 1].")

    clipped_vals = raw_vals.clip(0.0, 1.0)
    totals = clipped_vals.sum(axis=1)

    not_sum_one = ~np.isclose(totals, 1.0, atol=1e-3)
    if not_sum_one.any():
        num_not_one = int(not_sum_one.sum())
        print(f"Warning: Found {num_not_one} individual(s) whose ancestry proportions do not sum to 1. Proportions have been normalized to sum to 1.")

    valid = totals > 0
    clipped_vals.loc[valid] = clipped_vals.loc[valid].div(totals[valid], axis=0)
    q[source_pops] = clipped_vals
    return q, source_pops


def read_source_freq_file(path: str, source_pops: List[str]) -> pd.DataFrame:
    df = pd.read_csv(path, sep=None, engine="python")
    if "ID" not in df.columns:
        rename = {df.columns[0]: "ID"}
        df = df.rename(columns=rename)

    missing = [p for p in source_pops if p not in df.columns]
    if missing:
        raise ValueError(f"Source frequency file is missing columns: {', '.join(missing)}")

    out = df[["ID", *source_pops]].copy()
    out["ID"] = out["ID"].astype(str)
    for pop in source_pops:
        out[pop] = pd.to_numeric(out[pop], errors="coerce")
    return out



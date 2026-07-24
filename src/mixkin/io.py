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
    base = Path(prefix)
    ind_path = base.with_suffix(".ind")
    snp_path = base.with_suffix(".snp")
    geno_path = base.with_suffix(".geno")

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

    base = Path(prefix)
    bed = base.with_suffix(".bed")
    bim = base.with_suffix(".bim")
    fam = base.with_suffix(".fam")

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


def read_q_file(path: str, source_pops: List[str]) -> pd.DataFrame:
    q = pd.read_csv(path, sep=None, engine="python")
    if "ID" not in q.columns:
        rename = {q.columns[0]: "ID"}
        q = q.rename(columns=rename)

    missing = [p for p in source_pops if p not in q.columns]
    if missing:
        raise ValueError(f"Q file is missing ancestry columns for source populations: {', '.join(missing)}")

    q = q[["ID", *source_pops]].copy()
    q["ID"] = q["ID"].astype(str)
    for pop in source_pops:
        q[pop] = pd.to_numeric(q[pop], errors="coerce").fillna(0.0).clip(0.0, 1.0)

    total = q[source_pops].sum(axis=1)
    valid = total > 0
    q.loc[valid, source_pops] = q.loc[valid, source_pops].div(total[valid], axis=0)
    return q


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


def read_source_meta(path: Optional[str]) -> Dict[str, Dict[str, object]]:
    if not path:
        return {}
    df = pd.read_csv(path, sep=None, engine="python")
    required = {"population", "size"}
    if not required.issubset(df.columns):
        raise ValueError("source-meta must include columns: population,size and optional ploidy")

    meta: Dict[str, Dict[str, object]] = {}
    for _, row in df.iterrows():
        pop = str(row["population"])
        ploidy = str(row.get("ploidy", "pseudohaploid")).lower()
        size = float(row["size"])
        meta[pop] = {"size": size, "ploidy": ploidy}
    return meta

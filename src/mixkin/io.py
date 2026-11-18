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
    chromosomes: Optional[List[str]] = None
    positions: Optional[np.ndarray] = None


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
    chromosomes = snp_df["chr"].astype(str).tolist()
    positions = pd.to_numeric(snp_df["pos"], errors="coerce").fillna(0).to_numpy(dtype=np.int64)

    n_ind = len(sample_ids)
    n_snp = len(snp_ids)

    with open(geno_path, "rb") as handle:
        head_bytes = handle.read(100)

    head_str = head_bytes.split(b"\x00")[0].decode("ascii", errors="ignore").strip()
    head_parts = head_str.split()

    if head_parts and head_parts[0] == "TGENO":
        raise ValueError(
            "Transposed binary EIGENSTRAT (TGENO) format is not supported. "
            "Please convert to standard EIGENSTRAT, binary ANCESTRYMAP, or PLINK format."
        )

    if len(head_parts) >= 3 and head_parts[0] == "GENO":
        hdr_n_ind = int(head_parts[1])
        hdr_n_snp = int(head_parts[2])
        if hdr_n_ind != n_ind:
            raise ValueError(
                f"Binary ANCESTRYMAP header individual count ({hdr_n_ind}) does not match .ind file ({n_ind})"
            )
        if hdr_n_snp != n_snp:
            raise ValueError(
                f"Binary ANCESTRYMAP header SNP count ({hdr_n_snp}) does not match .snp file ({n_snp})"
            )

        row_len = max(48, (n_ind + 3) // 4)
        bytes_needed = (n_ind + 3) // 4
        with open(geno_path, "rb") as handle:
            handle.seek(row_len)
            raw_bytes = handle.read(n_snp * row_len)
            if len(raw_bytes) < n_snp * row_len:
                raise ValueError(
                    f"Binary ANCESTRYMAP file truncated: expected {n_snp * row_len} data bytes, got {len(raw_bytes)}"
                )
            raw = np.frombuffer(raw_bytes, dtype=np.uint8).reshape(n_snp, row_len)

        if n_snp > 0 and n_ind > 0:
            unpacked = np.unpackbits(raw[:, :bytes_needed], axis=1)[:, : (2 * n_ind)]
            geno_code = (2 * unpacked[:, ::2] + unpacked[:, 1::2]).astype(np.uint8)
            lut = np.array([2.0, 1.0, 0.0, np.nan], dtype=float)
            genotypes = lut[geno_code]
        else:
            genotypes = np.empty((n_snp, n_ind), dtype=float)
    else:
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
        if genotypes.shape[0] != n_snp:
            raise ValueError(
                f"EIGENSTRAT SNP count mismatch: .snp has {n_snp} SNPs, .geno has {genotypes.shape[0]} rows"
            )
    return DatasetData(
        sample_ids=sample_ids,
        populations=populations,
        snp_ids=snp_ids,
        genotypes=genotypes,
        chromosomes=chromosomes,
        positions=positions,
    )


def _read_plink(prefix: str) -> DatasetData:
    bed = Path(f"{prefix}.bed")
    bim = Path(f"{prefix}.bim")
    fam = Path(f"{prefix}.fam")

    fam_df = pd.read_csv(fam, sep=r"\s+", header=None)
    populations = fam_df.iloc[:, 0].astype(str).tolist()

    try:
        from pandas_plink import read_plink1_bin

        xr = read_plink1_bin(str(bed), str(bim), str(fam), verbose=False)
        # Ensure variants are the first dimension
        if xr.dims[0] == "sample":
            xr = xr.transpose("variant", "sample")
        geno = xr.values.astype(float)

        sample_ids = [str(v) for v in xr.coords["iid"].values.tolist()] if "iid" in xr.coords else [str(v) for v in xr.coords["sample"].values.tolist()]
        snp_ids = [str(v) for v in xr.coords["snp"].values.tolist()] if "snp" in xr.coords else [str(v) for v in xr.coords["variant"].values.tolist()]

        if "chrom" in xr.coords and "pos" in xr.coords:
            chromosomes = [str(v) for v in xr.coords["chrom"].values.tolist()]
            positions = pd.to_numeric(xr.coords["pos"].values, errors="coerce")
            positions = np.nan_to_num(positions, nan=0.0).astype(np.int64)
        else:
            bim_df = pd.read_csv(bim, sep=r"\s+", header=None, usecols=[0, 3], names=["chr", "pos"])
            chromosomes = bim_df["chr"].astype(str).tolist()
            positions = pd.to_numeric(bim_df["pos"], errors="coerce").fillna(0).to_numpy(dtype=np.int64)
    except ImportError:
        bim_df = pd.read_csv(
            bim,
            sep=r"\s+",
            header=None,
            usecols=[0, 1, 3],
            names=["chr", "snp", "pos"],
        )
        sample_ids = fam_df.iloc[:, 1].astype(str).tolist()
        snp_ids = bim_df["snp"].astype(str).tolist()
        chromosomes = bim_df["chr"].astype(str).tolist()
        positions = pd.to_numeric(bim_df["pos"], errors="coerce").fillna(0).to_numpy(dtype=np.int64)

        n_ind = len(sample_ids)
        n_snp = len(snp_ids)
        bytes_per_snp = (n_ind + 3) // 4
        with open(bed, "rb") as handle:
            magic = handle.read(3)
            if len(magic) < 3 or magic[0] != 0x6C or magic[1] != 0x1B or magic[2] != 0x01:
                raise ValueError(f"Invalid or unsupported PLINK .bed file: {bed}")
            raw = np.frombuffer(handle.read(), dtype=np.uint8)
        raw = raw.reshape(n_snp, bytes_per_snp)
        geno_code = np.empty((n_snp, bytes_per_snp * 4), dtype=np.uint8)
        for i in range(4):
            geno_code[:, i::4] = (raw >> (2 * i)) & 0b11
        geno_code = geno_code[:, :n_ind]
        lut = np.array([0.0, np.nan, 1.0, 2.0], dtype=float)
        geno = lut[geno_code]

    return DatasetData(
        sample_ids=sample_ids,
        populations=populations,
        snp_ids=snp_ids,
        genotypes=geno,
        chromosomes=chromosomes,
        positions=positions,
    )


def read_dataset(prefix: str, dataset_type: str) -> DatasetData:
    if dataset_type == "eigenstrat":
        return _read_eigenstrat(prefix)
    if dataset_type == "plink":
        return _read_plink(prefix)
    raise ValueError(f"Unsupported dataset type: {dataset_type}")


def _normalize_chrom(c: str) -> str:
    s = str(c).strip().lower()
    if s.startswith("chr"):
        s = s[3:]
    return s


def filter_dataset_chromosomes(
    ds: DatasetData,
    autosomes: int = 22,
    target_chr: Optional[str] = None,
) -> DatasetData:
    if ds.chromosomes is None:
        return ds

    norm_chroms = [_normalize_chrom(c) for c in ds.chromosomes]

    if target_chr:
        allowed = {_normalize_chrom(tok) for tok in target_chr.split(",") if tok.strip()}
        keep_mask = np.array([c in allowed for c in norm_chroms], dtype=bool)
        if not np.any(keep_mask) and autosomes == 22:
            fallback = set(allowed)
            if "x" in allowed or "23" in allowed:
                fallback.update({"x", "23"})
            if "y" in allowed or "24" in allowed:
                fallback.update({"y", "24"})
            keep_mask = np.array([c in fallback for c in norm_chroms], dtype=bool)
        if not np.any(keep_mask):
            raise ValueError(f"No SNPs found matching chromosome filter: '{target_chr}'")
    else:
        if autosomes < 1:
            raise ValueError("--autosomes must be a positive integer")
        allowed = {str(i) for i in range(1, autosomes + 1)}
        keep_mask = np.array([c in allowed for c in norm_chroms], dtype=bool)
        if not np.any(keep_mask):
            raise ValueError(f"No SNPs found on autosomes 1-{autosomes} (or chr1-chr{autosomes})")

    if np.all(keep_mask):
        return ds

    snp_ids = [s for s, k in zip(ds.snp_ids, keep_mask) if k]
    chromosomes = [c for c, k in zip(ds.chromosomes, keep_mask) if k]
    positions = ds.positions[keep_mask] if ds.positions is not None else None
    genotypes = ds.genotypes[keep_mask, :]

    return DatasetData(
        sample_ids=ds.sample_ids,
        populations=ds.populations,
        snp_ids=snp_ids,
        genotypes=genotypes,
        chromosomes=chromosomes,
        positions=positions,
    )


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



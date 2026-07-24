from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

from .io import (
    read_dataset,
    read_id_file,
    read_q_file,
    read_source_freq_file,
    read_source_meta,
)
from .model import estimate_source_freqs, run_kinship


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="MixKin: admixture-aware kinship estimation")

    parser.add_argument("--input-mode", choices=["whole", "target"], required=True)
    parser.add_argument("--dataset-type", choices=["plink", "eigenstrat"], required=True)
    parser.add_argument("--prefix", required=True, help="Dataset prefix")
    parser.add_argument("--source-pops", help="Comma-separated source population IDs")
    parser.add_argument("--q-file", required=True, help="Ancestry proportions file (ID + source population columns)")
    parser.add_argument("--target-ids", help="File with target IDs (one per line)")
    parser.add_argument("--target-pop", help="Population label in .ind/.fam to use as targets")

    parser.add_argument("--source-freq-file", help="Source frequencies table (required for --input-mode target)")

    parser.add_argument("--source-meta", help="CSV/TSV with columns: population,size,ploidy")
    parser.add_argument("--default-source-size", type=float, default=20.0)
    parser.add_argument("--default-source-ploidy", choices=["pseudohaploid", "diploid"], default="pseudohaploid")
    parser.add_argument("--epmr-mode", choices=["overlap", "all"], default="overlap", help="Calculate EPMR using only overlapping SNPs between targets or all available source SNPs")
    parser.add_argument("--out", default="mixkin_results.tsv")
    return parser.parse_args()


def _resolve_target_ids(all_ids: list[str], populations: list[str], target_ids_file: str | None, target_pop: str | None) -> list[str]:
    ids_from_file = read_id_file(target_ids_file)
    if ids_from_file:
        valid = [sid for sid in ids_from_file if sid in set(all_ids)]
        missing = sorted(set(ids_from_file) - set(valid))
        if missing:
            sys.exit(f"Error: Target IDs not found in dataset: {', '.join(missing)}")
        return valid

    if target_pop:
        if target_pop not in populations:
            sys.exit(f"Error: Target population '{target_pop}' not found in the dataset.")
        return [sid for sid, pop in zip(all_ids, populations) if pop == target_pop]

    return list(all_ids)


def check_file_exists(filepath: str | None, name: str):
    if filepath is not None:
        if not Path(filepath).is_file():
            sys.exit(f"Error: {name} file '{filepath}' does not exist.")

def check_prefix_exists(prefix: str, dataset_type: str):
    base = Path(prefix)
    exts = [".ind", ".snp", ".geno"] if dataset_type == "eigenstrat" else [".bed", ".bim", ".fam"]
    for ext in exts:
        if not base.with_suffix(ext).is_file():
            sys.exit(f"Error: Dataset file '{base.with_suffix(ext)}' does not exist.")

def print_banner():
    banner = """
======================================================================
  mixKIN version 1.0.0
  Developed by Muhammed Siddik Kilic
======================================================================"""
    print(banner, flush=True)

def main() -> int:
    print_banner()
    args = parse_args()
    
    check_prefix_exists(args.prefix, args.dataset_type)
    check_file_exists(args.q_file, "--q-file")
    check_file_exists(args.target_ids, "--target-ids")
    check_file_exists(args.source_freq_file, "--source-freq-file")
    check_file_exists(args.source_meta, "--source-meta")

    if args.input_mode == "target" and args.source_freq_file and not args.source_meta:
        sys.exit("Error: --source-meta parameter is required when providing an external allele frequency table (--source-freq-file).")

    if args.source_pops:
        source_pops = [v.strip() for v in args.source_pops.split(",") if v.strip()]
    elif args.input_mode == "target" and args.source_freq_file:
        df_cols = pd.read_csv(args.source_freq_file, sep=None, engine="python", nrows=0).columns.tolist()
        source_pops = [c for i, c in enumerate(df_cols) if i > 0 and c != "ID"]
    else:
        raise ValueError("--source-pops is required unless --input-mode target and --source-freq-file are provided.")

    source_meta = read_source_meta(args.source_meta)

    ds = read_dataset(args.prefix, args.dataset_type)
    target_ids = _resolve_target_ids(
        all_ids=ds.sample_ids,
        populations=ds.populations,
        target_ids_file=args.target_ids,
        target_pop=args.target_pop,
    )

    if not target_ids:
        raise ValueError("No target individuals found in dataset")
    if len(target_ids) < 2:
        raise ValueError("At least two target individuals are required")

    target_idx = [ds.sample_ids.index(x) for x in target_ids]

    q_df = read_q_file(args.q_file, source_pops)

    if args.input_mode == "whole":
        missing_pops = [p for p in source_pops if p not in ds.populations]
        if missing_pops:
            sys.exit(f"Error: Source populations not found in the dataset: {', '.join(missing_pops)}")

        if not args.target_pop and not args.target_ids:
            raise ValueError("--target-pop or --target-ids is required for --input-mode whole to identify target individuals")

        source_freqs, source_sizes, source_callrates = estimate_source_freqs(
            genotypes=ds.genotypes,
            sample_ids=ds.sample_ids,
            populations=ds.populations,
            source_pops=source_pops,
        )

        # In whole mode, SNPs naturally align. Ensure source_freqs aligns.
        src = pd.DataFrame({"ID": ds.snp_ids}).join(source_freqs)
        source_freqs = src[source_pops]

    else:
        if not args.source_freq_file:
            raise ValueError("--source-freq-file is required for --input-mode target")

        source_freq_df = read_source_freq_file(args.source_freq_file, source_pops)
        source_freqs = pd.DataFrame({"ID": ds.snp_ids}).merge(source_freq_df, on="ID", how="left")
        source_freqs = source_freqs[source_pops]
        source_sizes = {}
        source_callrates = {}

    results = run_kinship(
        target_ids=target_ids,
        sample_ids=ds.sample_ids,
        target_genotypes=ds.genotypes,
        target_sample_indices=target_idx,
        target_snp_ids=ds.snp_ids,
        source_freqs=source_freqs,
        source_pops=source_pops,
        q_df=q_df,
        source_meta=source_meta,
        source_sizes=source_sizes,
        source_callrates=source_callrates,
        default_source_size=args.default_source_size,
        default_source_ploidy=args.default_source_ploidy,
        epmr_mode=args.epmr_mode,
    )

    out_path = Path(args.out)
    results.to_csv(out_path, sep="\t", index=False)
    print(f"MixKin finished. Wrote {len(results)} pairs to {out_path}")
    return 0

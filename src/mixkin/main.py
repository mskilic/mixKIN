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
)
from .model import estimate_source_freqs, run_kinship


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="MixKin: admixture-aware kinship estimation")

    parser.add_argument("-i", "--input", required=True, help="Input file prefix")
    parser.add_argument("-q", "--ancestry-prop", required=False, help="Ancestry proportions file (Q-matrix)")
    parser.add_argument("-t", "--target-pop", help="Target population label (required if -q is omitted)")
    parser.add_argument("-s", "--source-pop", help="Source population label (defaults to target-pop if -q is omitted)")
    parser.add_argument("-e", "--epmr", help="EPMR calculation override ('median' or float). Only valid if -q is omitted or Q-file has one source population.")
    parser.add_argument("-p", "--source-freq", help="Source frequency file")
    parser.add_argument("-o", "--out", default="mixkin.out", help="Output file name")
    parser.add_argument("-y", "--ploidy", default="1", help="Ploidy values (1, 2, or comma-separated sequence)")
    parser.add_argument("-v", "--overlap", action="store_true", help="Calculate EPMR using only overlapping SNPs")

    parser.add_argument("-n", "--nind", help="Number of individuals for each source population (comma-separated)")
    parser.add_argument("-g", "--geno-rate", help="Average genotyping rate for each source population (comma-separated)")
    
    return parser.parse_args()


def _detect_dataset_type(prefix: str) -> str:
    base = Path(prefix)
    has_plink = base.with_suffix(".bed").is_file() and base.with_suffix(".bim").is_file() and base.with_suffix(".fam").is_file()
    has_eigenstrat = base.with_suffix(".geno").is_file() and base.with_suffix(".snp").is_file() and base.with_suffix(".ind").is_file()
    
    if has_plink and has_eigenstrat:
        print("Both PLINK and Eigenstrat formats detected. Proceeding with PLINK.")
        return "plink"
    elif has_plink:
        return "plink"
    elif has_eigenstrat:
        return "eigenstrat"
    else:
        sys.exit(f"Error: Could not find valid PLINK or Eigenstrat files for prefix '{prefix}'")

def check_file_exists(filepath: str | None, name: str):
    if filepath is not None:
        if not Path(filepath).is_file():
            sys.exit(f"Error: {name} file '{filepath}' does not exist.")

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
    
    dataset_type = _detect_dataset_type(args.input)
    ds = read_dataset(args.input, dataset_type)
    
    check_file_exists(args.source_freq, "--source-freq")
    external_freq = bool(args.source_freq)

    if args.ancestry_prop:
        check_file_exists(args.ancestry_prop, "--ancestry-prop")
        q_df, source_pops = read_q_file(args.ancestry_prop)
        query_ids_from_q = q_df["ID"].tolist()
    else:
        if not args.target_pop:
            sys.exit("Error: -t/--target-pop is required when -q/--ancestry-prop is omitted.")
        query_ids_from_q = [sid for sid, pop in zip(ds.sample_ids, ds.populations) if pop == args.target_pop]
        if not query_ids_from_q:
            sys.exit(f"Error: No individuals found for target population '{args.target_pop}'")
        
        source_pop = args.source_pop if args.source_pop else args.target_pop
        source_pops = [source_pop]
        q_df = pd.DataFrame({"ID": query_ids_from_q, source_pop: 1.0})

    if args.epmr and (args.ancestry_prop and len(source_pops) > 1):
        sys.exit("Error: -e/--epmr can only be used when -q is omitted or Q-file contains exactly one source population.")

    ploidy_arg = args.ploidy.strip()
    ploidy_map = {}
    if ploidy_arg == "1":
        for p in source_pops: ploidy_map[p] = "pseudohaploid"
    elif ploidy_arg == "2":
        for p in source_pops: ploidy_map[p] = "diploid"
    else:
        parts = [p.strip() for p in ploidy_arg.split(",")]
        if len(parts) != len(source_pops):
            sys.exit(f"Error: Number of ploidy values ({len(parts)}) does not match number of source populations ({len(source_pops)}).")
        for pop, p_val in zip(source_pops, parts):
            if p_val == "1":
                ploidy_map[pop] = "pseudohaploid"
            elif p_val == "2":
                ploidy_map[pop] = "diploid"
            else:
                sys.exit(f"Error: Invalid ploidy value '{p_val}'. Must be 1 or 2.")

    manual_source_sizes = {}
    if args.nind or args.geno_rate:
        if not (args.nind and args.geno_rate):
            sys.exit("Error: Both -n/--nind and -g/--geno-rate must be provided together.")
        
        nind_parts = [p.strip() for p in args.nind.split(",")]
        geno_parts = [p.strip() for p in args.geno_rate.split(",")]
        if len(nind_parts) != len(source_pops):
            sys.exit(f"Error: Number of --nind values ({len(nind_parts)}) does not match source populations ({len(source_pops)}).")
        if len(geno_parts) != len(source_pops):
            sys.exit(f"Error: Number of --geno-rate values ({len(geno_parts)}) does not match source populations ({len(source_pops)}).")
        for pop, n_str, g_str in zip(source_pops, nind_parts, geno_parts):
            manual_source_sizes[pop] = float(n_str) * float(g_str)
    elif external_freq:
        sys.exit("Error: -n/--nind and -g/--geno-rate are required when using an external frequency file.")

    fixed_epmr = None
    if args.epmr:
        if args.epmr.lower() == "median":
            source_idx = [i for i, pop in enumerate(ds.populations) if pop == source_pops[0]]
            if len(source_idx) < 2:
                sys.exit(f"Error: Source population '{source_pops[0]}' has fewer than 2 individuals. Cannot compute median EPMR.")
            from .model import compute_median_pmr
            fixed_epmr = compute_median_pmr(ds.genotypes, source_idx)
        else:
            try:
                fixed_epmr = float(args.epmr)
            except ValueError:
                sys.exit("Error: -e/--epmr must be 'median' or a numeric float value.")

    query_ids_from_q = q_df["ID"].tolist()
    ds_sample_ids_set = set(ds.sample_ids)
    query_ids = [sid for sid in query_ids_from_q if sid in ds_sample_ids_set]
    missing_query_ids = sorted(set(query_ids_from_q) - ds_sample_ids_set)
    if missing_query_ids:
        print(f"Warning: The following individuals from Q-file were not found in the dataset: {', '.join(missing_query_ids)}")

    if not query_ids:
        raise ValueError("No individuals found in dataset to analyze")
    if len(query_ids) < 2:
        raise ValueError("At least two individuals are required to calculate kinship")

    query_idx = [ds.sample_ids.index(x) for x in query_ids]

    if not external_freq:
        missing_pops = [p for p in source_pops if p not in ds.populations]
        if missing_pops:
            sys.exit(f"Error: Source populations not found in the dataset: {', '.join(missing_pops)}")

        source_freqs, source_sizes, source_callrates = estimate_source_freqs(
            genotypes=ds.genotypes,
            sample_ids=ds.sample_ids,
            populations=ds.populations,
            source_pops=source_pops,
        )

        if manual_source_sizes:
            print("Warning: -n/--nind and -g/--geno-rate were provided but will be ignored. Since allele frequencies are estimated internally from the dataset, sizes are automatically estimated as well.")

        src = pd.DataFrame({"ID": ds.snp_ids}).join(source_freqs)
        source_freqs = src[source_pops]

    else:
        source_freq_df = read_source_freq_file(args.source_freq, source_pops)
        source_freqs = pd.DataFrame({"ID": ds.snp_ids}).merge(source_freq_df, on="ID", how="left")
        source_freqs = source_freqs[source_pops]
        source_sizes = manual_source_sizes
        source_callrates = {}
        
    epmr_mode = "overlap" if args.overlap else "all"

    results = run_kinship(
        query_ids=query_ids,
        sample_ids=ds.sample_ids,
        query_genotypes=ds.genotypes,
        query_sample_indices=query_idx,
        query_snp_ids=ds.snp_ids,
        source_freqs=source_freqs,
        source_pops=source_pops,
        q_df=q_df,
        source_sizes=source_sizes,
        source_ploidies=ploidy_map,
        source_callrates=source_callrates,
        epmr_mode=epmr_mode,
        fixed_epmr=fixed_epmr,
    )

    out_path = Path(args.out)
    results.to_csv(out_path, sep="\t", index=False)
    print(f"MixKin finished. Wrote {len(results)} pairs to {out_path}")
    return 0

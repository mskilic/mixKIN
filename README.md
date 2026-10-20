![logo](mixKIN.png)

mixKIN is a Python package for ancestry-aware kinship estimation from low-coverage genomes. It estimates pair-specific expected Pairwise Mismatch Ratio (ePMR) using individual ancestry proportions and source-population allele frequencies. This approach improves kinship estimation in the admixed, heterogeneous, as well as homogeneous populations.

## Table of Contents
- [Installation](#installation)
- [Data Requirements](#data-requirements)
- [Usage Examples](#usage-examples)
- [Parameters](#parameters)
- [Understanding the Output](#understanding-the-output)
- [Citation](#citation)

## Installation

We recommend using a dedicated Conda or Python virtual environment to keep mixKIN's dependencies separate to avoid dependency conflicts.

### Using Conda
```bash
conda create -n mixkin_env python=3.10
conda activate mixkin_env
git clone https://github.com/mskilic/mixKIN.git
cd mixKIN
pip install .
```

### Using Python venv

```bash
git clone https://github.com/mskilic/mixKIN.git
cd mixKIN
python -m venv mixkin_env
source mixkin_env/bin/activate  # On Windows use `mixkin_env\Scripts\activate`
pip install .
```

## Data Requirements

### Genotype Dataset

Supported formats are 

- **PLINK binary** (`.bed`, `.bim`, `.fam`)
- **Unpacked EIGENSTRAT** (`.geno`, `.snp`, `.ind`).

Specify the input file prefix without an extension, mixKIN automatically identifies the files mathching that prefix.

Individual IDs must match between the genotype dataset and ancestry-proportion file. Source-population labels must also match the labels used in the ancestry file and frequency table. 

### Ancestry Proportions

For ancestry-aware analysis, provide a Q-matrix file containing individual admixture proportions. The first column contains individual IDs, and subsequent columns corresponds to source populations.

The following example illustrates a two-source model:

```text
ID        JPT   TSI
sample_1  0.50  0.50
sample_2  0.75  0.25
sample_3  1.00  0.00
```

Ancestry proportions should lie between 0 and 1 ans sum to 1 for each individual. Source columns should have consistent names across input files.

For a single homogeneous population, the Q-file can be omitted by specifying target and source population labels, as shown below. 

### Source reference and ploidy

Source allele frequencies can be calculated internally from reference individuals/populations in the genotype dataset or provided in an external table.

Specify reference ploidy using `-y`: use `1` for pseudohaploid references and `2` for diploid references. Reference sample size, call rate, and ploidy determine the mean observed allele count used in the finite-sample correction.

**The default reference ploidy is 1. Set `-y 2` when using diploid source references.**

### External allele frequencies

The external frequency table contains an `ID` column of SNP IDs followed by one frequency column per source population. For example:

```text
ID	JPT	TSI
rs_example_1	0.20	0.60
rs_example_2	0.70	0.40
```

For each SNP, make sure the external frequency refers to the allele counted in the genotype dataset. Matching SNP IDs is not enough: allele coding and strand orientation must also be consistent between the two inputs.

When external frequencies are supplied, also provide source sample sizes (`-n`), mean genotype call rates (`-g`), and appropriate reference ploidy (`-y`).

## Usage Examples

### Internal frequencies

By default, mixKIN estimates source allele frequencies internally. Use this approach when the input dataset includes target individuals and the required source references.

For pseudohaploid source references:

```bash
mixkin \
  -i /path/to/dataset_prefix \
  -q /path/to/ancestry_q.tsv \
  -y 1 \
  -o kinship_results.tsv
```

For diploid source references, replace `-y 1` with `-y 2`.

### External frequencies

Provide a precomputed source-frequency table to bypass internal frequency estimation. This example assumes frequency columns ordered as JPT and TSI, with 20 diploid reference individuals per source:

```bash
mixkin \
  -i /path/to/dataset_prefix \
  -q /path/to/ancestry_q.tsv \
  -p /path/to/source_freqs.tsv \
  -n 20,20 \
  -g 0.98,0.95 \
  -y 2 \
  -o kinship_results.tsv
```

The sample sizes and call rates are illustrative; replace them with values for the references used to estimate your frequencies. Call rates are proportions between 0 and 1, not percentages.

### Homogeneous population without a Q-file

For a single homogeneous population, specify the target and source population labels:

```bash
mixkin \
  -i /path/to/dataset_prefix \
  -t TargetPop \
  -s SourcePop \
  -y 2 \
  -o kinship_results.tsv
```

Here, `-y 2` assumes diploid source references. If `-s` is omitted, the source label defaults to the target label.

### Pair-specific overlapping SNPs

Use `-v` to calculate ePMR on SNPs called in both individuals in each pair:

```bash
mixkin \
  -i /path/to/dataset_prefix \
  -q /path/to/ancestry_q.tsv \
  -y 2 \
  -v \
  -o kinship_results.tsv
```


## Parameters

### Core arguments

| Argument                | Description                                             |
| ----------------------- | ------------------------------------------------------- |
| `-i`, `--input`         | Genotype file prefix; required.                         |
| `-q`, `--ancestry-prop` | Q-matrix containing individual ancestry proportions.    |
| `-t`, `--target-pop`    | Target population label; required when `-q` is omitted. |
| `-s`, `--source-pop`    | Source population label; defaults to `-t` when omitted. |
| `-p`, `--source-freq`   | External source allele-frequency table.                 |
| `-o`, `--out`           | Output file name; default: `mixkin.out`.                |

### Reference and baseline settings

| Argument            | Description                                                                                              |
| ------------------- | -------------------------------------------------------------------------------------------------------- |
| `-n`, `--nind`      | Comma-separated reference individual counts for source populations; required with `-p`.                  |
| `-g`, `--geno-rate` | Comma-separated mean genotype call rates for source populations; required with `-p`.                     |
| `-y`, `--ploidy`    | Reference ploidy; default: `1`. Accepts a single value or a comma-separated list for source populations. |
| `-v`, `--overlap`   | Calculate ePMR using SNPs called in both members of each pair.                                           |
| `-e`, `--epmr`      | Alternative mismatch baseline; valid only when `-q` is omitted or contains one source population.        |

For `-e`, supported inputs are:

- **`median`:** use the median observed PMR across pairs within the source population.
- **A numeric value**, such as `0.22`: use a fixed expected mismatch baseline.

Median or fixed baselines do not provide the multi-source ancestry weighting used in the default ancestry-aware model.

## Understanding the Output

The output is a tab-separated file containing:

| Column            | Description                                                                               |
| ----------------- | ----------------------------------------------------------------------------------------- |
| `ID1`, `ID2`      | IDs of the two individuals being compared.                                                |
| `Theta`           | Estimated kinship coefficient.                                                            |
| `PMR`             | Observed pairwise mismatch ratio (oPMR).                                                  |
| `EPMR`            | Expected mismatch ratio for an unrelated pair with the corresponding ancestry profiles.   |
| `OverlappingSNPs` | Number of eligible SNPs called in both individuals.                                       |
| `Relationship`    | Inferred category: `Unrelated`, `Third`, `Second`, `First`, or `MZT` (monozygotic twins). |

Kinship is estimated by comparing observed and expected mismatch:

$$
\widehat{\theta}=1-\frac{\mathrm{oPMR}}{\mathrm{ePMR}}.
$$

For non-inbred pedigree relationships, expected kinship coefficients are 0.25 for first-degree, 0.125 for second-degree, and 0.0625 for third-degree relatives. These are expectations, not exact values for every pair.

## Citation

If you use mixKIN in your research, please cite the manuscript:

> Kılıç, M. S., Malaspinas, A.-S., & Koptekin, D. **mixKIN: An ancestry-aware kinship estimation method for low-coverage genomes.**
> https://www.biorxiv.org/content/10.64898/2026.10.13.758828v1

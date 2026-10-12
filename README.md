![logo](mixKIN.png)

mixKIN is a Python package for ancestry-aware kinship estimation from low-coverage genomes. By calculating pair-specific Pairwise Mismatch Ratio (PMR) baselines based on individual ancestry proportions and source allele frequencies, mixKIN improves relatedness inference across admixed, heterogeneous, and homogeneous populations.

## Table of Contents
- [Installation](#installation)
- [Data Requirements](#data-requirements)
- [Usage Examples](#usage-examples)
- [Parameters](#parameters)
- [Understanding the Output](#understanding-the-output)
- [Citation](#citation)

## Installation

We recommend installing mixKIN inside a virtual environment to avoid dependency conflicts.

### Using Conda
```bash
conda create -n mixkin_env python=3.10
conda activate mixkin_env
git clone https://github.com/mskilic/mixKIN.git
cd mixKIN
pip install .
```

### Using standard Virtualenv

```bash
git clone https://github.com/mskilic/mixKIN.git
cd mixKIN
python -m venv mixkin_env
source mixkin_env/bin/activate  # On Windows use `mixkin_env\Scripts\activate`
pip install .
```

Verify the installation:
```bash
mixkin --help
```

## Data Requirements

**Genotype Dataset**: Supported formats are PLINK (`.bed`, `.bim`, `.fam`) and unpacked EIGENSTRAT (`.geno`, `.snp`, `.ind`). The program will auto-detect the format based on the prefix.

**Ancestry Proportions**: A Q-matrix file containing admixture proportions. This file is required for analyzing admixed individuals. The first column must contain individual IDs, and subsequent columns represent source populations.

## Usage Examples

### Internal Frequencies (Default)
By default, mixKIN estimates the source-population allele frequencies internally from the input dataset. Use this approach when your dataset contains both the query individuals and the source populations.

```bash
mixkin \
  -i /path/to/dataset_prefix \
  -q /path/to/ancestry_q.tsv \
  -o kinship_results.tsv
```

### External Frequencies
You can bypass internal estimation by providing a pre-computed external allele frequency table. This must be a tabular file with an `ID` column (SNP IDs) followed by columns for each source population containing their allele frequencies.

```bash
mixkin \
  -i /path/to/dataset_prefix \
  -q /path/to/ancestry_q.tsv \
  -p /path/to/source_freqs.tsv \
  -n 50,45,60 \
  -g 0.98,0.95,0.99 \
  -o kinship_results.tsv
```

### Homogeneous Populations (No Q-file)
If you are analyzing a homogeneous population, you can omit the Q-file and specify the target population directly.

```bash
mixkin \
  -i /path/to/dataset_prefix \
  -t TargetPop \
  -s SourcePop \
  -o kinship_results.tsv
```

## Parameters

### Core Arguments
- `-i`, `--input`: Genotype file prefix (required).
- `-q`, `--ancestry-prop`: Q-matrix file containing admixture proportions.
- `-t`, `--target-pop`: Target population label (required if `-q` is omitted).
- `-s`, `--source-pop`: Source population label (defaults to `-t` if omitted).
- `-p`, `--source-freq`: External source frequency file.
- `-o`, `--out`: Output file name (default: `mixkin.out`).

### Advanced Configurations
- `-n`, `--nind`: Comma-separated list of the number of individuals for each source population. Required when `-p` is used.
- `-g`, `--geno-rate`: Comma-separated list of the average genotyping rates for each source population. Required when `-p` is used.
- `-e`, `--epmr`: EPMR calculation override. Valid only if `-q` is omitted or contains a single source population.
  - `median`: Uses the median PMR of all pairs within the source population.
  - *Float value* (e.g., `0.22`): Bypasses calculations and uses this exact numeric value.
- `-y`, `--ploidy`: Numeric ploidy values (default: 1). Accepts a single value (e.g., `2`) or a comma-separated sequence mapping to each source population.
- `-v`, `--overlap`: Calculate EPMR using strictly the overlapping (non-missing) SNPs between each specific pair of individuals.

## Understanding the Output

The output is a tab-separated (TSV) file with the following columns:

| Column | Description |
|---|---|
| `ID1`, `ID2` | The pair of individuals being compared |
| `Theta` | The estimated kinship coefficient |
| `PMR` | Pairwise Mismatch Rate observed between the two individuals |
| `EPMR` | Expected Pairwise Mismatch Rate given ancestry proportions and frequencies |
| `OverlappingSNPs` | Number of shared non-missing SNPs used for the calculation |
| `Relationship` | Inferred relatedness degree (`Unrelated`, `Third`, `Second`, `First`, `MZT`) |

## Citation

If you use mixKIN in your research, please cite our paper:

> **mixKIN: An ancestry-aware kinship estimation method for low coverage genomes**  
> Kılıç, M. S., Malaspinas, A.-S, Koptekin, D.
> 

# mixKIN

**mixKIN** is a powerful Python package for admixture-aware kinship estimation using PMR (Pairwise Mismatch Rate) normalization. It accurately estimates kinship coefficients and infers relatedness degrees in admixed populations by taking into account ancestry proportions and source population allele frequencies.

---

## 📖 Table of Contents
1. [Installation](#installation)
2. [Data Requirements](#data-requirements)
3. [Tutorial: How to Run mixKIN](#tutorial-how-to-run-mixkin)
   - [Mode 1: `whole` (Internal Frequencies)](#mode-1-whole-internal-frequencies)
   - [Mode 2: `target` (External Frequencies)](#mode-2-target-external-frequencies)
4. [Advanced Configuration](#advanced-configuration)
5. [Understanding the Output](#understanding-the-output)

---

## 🛠️ Installation

We highly recommend installing `mixKIN` inside a virtual environment (e.g., using `conda` or `venv`) to avoid dependency conflicts.

### Option 1: Using Conda (Recommended)

1. Create and activate a new conda environment:
   ```bash
   conda create -n mixkin_env python=3.10
   conda activate mixkin_env
   ```
2. Clone the GitHub repository and install:
   ```bash
   git clone https://github.com/YourUsername/mixKIN.git
   cd mixKIN
   pip install .
   ```

### Option 2: Using standard Pip / Virtualenv

```bash
git clone https://github.com/YourUsername/mixKIN.git
cd mixKIN
python -m venv mixkin_env
source mixkin_env/bin/activate  # On Windows use `mixkin_env\Scripts\activate`
pip install .
```

You can verify the installation by running the command in your terminal. You should see the welcome banner and the help menu:
```bash
mixkin --help
```

---

## 📁 Data Requirements

To run `mixKIN`, you will need:
- **Genotype Dataset**: Supported formats are **PLINK** (`.bed`, `.bim`, `.fam`) and **EIGENSTRAT** (`.geno`, `.snp`, `.ind`).
- **Ancestry Proportions (`--q-file`)**: A file containing the admixture proportions for each target individual. The file must have an `ID` column and subsequent columns for each source population (e.g., values summing to 1.0 per individual).

Optional:
- **Target IDs File (`--target-ids`)**: A simple text file with one individual ID per line, specifying which individuals to calculate kinship for.
- **Source Frequencies (`--source-freq-file`)**: Required only if calculating frequencies externally (see Mode 2).
- **Source Metadata (`--source-meta`)**: A CSV/TSV containing `population`, `size`, and optionally `ploidy` (`diploid` or `pseudohaploid`) for source populations.

---

## 🚀 Tutorial: How to Run mixKIN

`mixKIN` operates in two primary modes: `whole` and `target`.

### Mode 1: `whole` (Internal Frequencies)
Use `whole` mode when your dataset contains **both** your target individuals and the source populations. `mixKIN` will internally estimate the source-population allele frequencies from the dataset.

**Example Command:**
```bash
mixkin \
  --input-mode whole \
  --dataset-type eigenstrat \
  --prefix /path/to/whole_dataset \
  --source-pops PopA,PopB,PopC \
  --q-file /path/to/ancestry_q.tsv \
  --target-pop TARGET_POP \
  --out kinship_results.tsv
```

**Parameters Explained:**
- `--input-mode whole`: Tells mixKIN to calculate allele frequencies internally.
- `--dataset-type eigenstrat`: The format of your genotype data (`eigenstrat` or `plink`).
- `--prefix`: The base path to your dataset (without file extensions).
- `--source-pops`: A comma-separated list of the source population names present in your dataset.
- `--q-file`: The ancestry proportions file.
- `--target-pop`: Calculate kinship for all individuals belonging to this population in the dataset. Alternatively, use `--target-ids target_list.txt` to specify a list of IDs.
- `--out`: The name of the output TSV file.

---

### Mode 2: `target` (External Frequencies)
Use `target` mode when your dataset contains **only** the target individuals, and you want to provide a pre-computed allele frequency table for the source populations.

**Example Command:**
```bash
mixkin \
  --input-mode target \
  --dataset-type plink \
  --prefix /path/to/target_dataset \
  --source-freq-file /path/to/source_freqs.tsv \
  --source-meta /path/to/source_meta.tsv \
  --q-file /path/to/ancestry_q.tsv \
  --target-ids /path/to/target_ids.txt \
  --out kinship_results.tsv
```

**Parameters Explained:**
- `--input-mode target`: Tells mixKIN to use the external frequency table.
- `--source-freq-file`: A tabular file with an `ID` column (SNP IDs) and columns for each source population containing allele frequencies.
- `--source-meta`: **Required** when using an external frequency table. A file containing `population`, `size`, and `ploidy` for the sources.
- `--target-ids`: A text file containing the IDs of individuals to test.

> **Note:** In `target` mode, `--source-pops` can optionally be omitted, as mixKIN can infer the source populations automatically from the header of the `--source-freq-file`.

---

## ⚙️ Advanced Configuration

- **Source Metadata (`--source-meta`)**: 
  You can fine-tune source population sample sizes and ploidy. The file should have columns: `population`, `size`, `ploidy`. If not provided (in `whole` mode), mixKIN uses sizes inferred from the dataset and defaults to the ploidy specified by `--default-source-ploidy`.
  
- **Ploidy Default (`--default-source-ploidy`)**:
  Defaults to `pseudohaploid`. You can set this to `diploid` to treat all source populations as diploid globally without a meta file.

- **EPMR Mode (`--epmr-mode`)**:
  - `overlap` (Default): Calculates Expected PMR (EPMR) using only SNPs that overlap (are non-missing) between the target pair.
  - `all`: Calculates EPMR using all available SNPs from the sources.

---

## 📊 Understanding the Output

mixKIN generates a tab-separated (TSV) output file with the following columns:

| Column | Description |
| :--- | :--- |
| **`ID1`, `ID2`** | The pair of individuals being compared. |
| **`Theta`** | The estimated kinship coefficient. |
| **`PMR`** | Pairwise Mismatch Rate observed between the two individuals. |
| **`EPMR`** | Expected Pairwise Mismatch Rate given their ancestry proportions and source frequencies. |
| **`OverlappingSNPs`** | Number of shared non-missing SNPs used for the calculation. |
| **`Relationship`** | Inferred relatedness degree (`Unrelated`, `Third`, `Second`, `First`, `MZT` (Monozygotic Twins / Identical)). |

---
## 📜 Citation

If you use `mixKIN` in your research, please cite our paper:

> **[Insert Paper Title Here]**
> Muhammed Siddik Kilic et al. (Year)
> [Link to the paper](https://doi.org/10...)

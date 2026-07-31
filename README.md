# CGD — Comprehensive Guide RNA Design

[![CI](https://github.com/vipinmenon1989/CGD/actions/workflows/ci.yml/badge.svg)](https://github.com/vipinmenon1989/CGD/actions/workflows/ci.yml)

**Author:** Vipin Menon, BIG Lab, Hanyang University (HYU)
**Contact:** a.vipin.menon@gmail.com
**Original:** September 2019 | **Modernized:** 2026
**Language:** Python 3.10+
**License:** MIT

---

## Overview

CGD is a Python-based on-target scoring tool for multiple CRISPR systems. It predicts guide RNA activity using **ENLOR (Elastic Net Logistic Regression)** models trained on experimental data.

| Score  | CRISPR System              | PAM       |
|--------|---------------------------|-----------|
| CGDi   | CRISPRi (interference)    | NGG       |
| CGDa   | CRISPRa (activation)      | NGG       |
| CGD9   | CRISPR-Cas9 canonical     | NGG       |
| CGD9NG | CRISPR-Cas9 non-canonical | NGA/NGC/NGT |
| CGD12a | CRISPR-Cas12a (Cpf1)      | TTTV      |

Each model combines: RNA free energy · sequence entropy · GC content · position-dependent nucleotide features.

> **Web interface:** http://big.hanyang.ac.kr:2195/CGD

---

## Installation

### 1. Clone the repository
```bash
git clone https://github.com/vipinmenon1989/CGD.git
cd CGD
```

### 2. Create the CGD runtime environment
```bash
conda env create -f envs/environment.yaml
conda activate cgd
```

This installs Python, NumPy, and the ViennaRNA Python bindings together in one
environment. For a pip-only local installation, `pip install -r requirements.txt`
is also supported when a compatible ViennaRNA wheel is available for the
operating system.

---

## Input Format

A standard **FASTA file** with one or more sequences, each between **100 and 10 000 nt**.

```
>XM_030244935
CGGCGCGGAGTGCGCCGGCGCGTCGTCGGGGACGCCGGGTCCAGGATCTTGCTAGGGAA
CCAGTGTTGTCGCGTCGTCCCGCCCCCTCGGGGCTTTTGCTCCCGTTAACTGTCGGCGG
...
```

---

## Usage

### Comprehensive scoring (all CRISPR systems)
```bash
python CGD.py -a input.fa
```
Output: `CGD.txt` — tab-separated with columns `ID, Start, End, Strand, Sequence, CGDi, CGDa, CGD9, CGDNG, CGD12a`

### CRISPRi only
```bash
python CGD.py -b input.fa
```
Output: `CGDi.txt`

### CRISPRa only
```bash
python CGD.py -c input.fa
```
Output: `CGDa.txt`

### CRISPR-Cas9 (canonical NGG) only
```bash
python CGD.py -d input.fa
```
Output: `CGD9.txt`

### CRISPR-Cas12a only
```bash
python CGD.py -e input.fa
```
Output: `CGD12a.txt`

### CRISPR-Cas9 non-canonical (NGA/NGC/NGT) only
```bash
python CGD.py -f input.fa
```
Output: `CGD9NG.txt`

---

## Output Format

All output files are tab-separated with a header row:

| Column     | Description                              |
|------------|------------------------------------------|
| ID         | Sequence identifier from FASTA header    |
| Start      | Guide start position (0-indexed)         |
| End        | Guide end position                       |
| Strand     | `+` (forward) or `-` (reverse)          |
| Sequence   | 30 or 34 bp guide sequence               |
| Score      | ENLOR-based activity score in [0, 1]     |

Guides with scores **< 0.5** are considered **efficient** by the CGD model.

---

### Custom output path

Every mode accepts `-o/--output` to write results to a specific path instead of
the mode's default filename in the current directory:

```bash
python CGD.py -a input.fa -o results/my_sample_CGD.txt
```

---

## Snakemake Workflow

For batch scoring of multiple FASTA files across multiple CRISPR systems, use
the included Snakemake workflow instead of calling `CGD.py` by hand.

### 1. Install the workflow environment

Install Miniforge, Mamba, or Conda on the workstation or HPC login node, then:

```bash
conda env create -f workflow-env.yaml
conda activate cgd-workflow
```

This environment contains Snakemake 9, Conda for rule-specific software
deployment, and the official SLURM executor plugin. CGD itself runs inside the
isolated environment declared in `envs/environment.yaml`.

### 2. Configure samples and modes
Edit `config/config.yaml`:
```yaml
samples:
  demo: input.fa          # sample_name: path/to/input.fa
  sample_2: /shared/project/sample_2.fa

modes:
  - comprehensive          # CGD.py -a
  - crispri                # CGD.py -b
  - crispra                # CGD.py -c
  - cas9                   # CGD.py -d
  - cas12a                 # CGD.py -e
  - cas9_ng                # CGD.py -f

outdir: results
```

Sample names may contain letters, numbers, dots, underscores, and hyphens.
Input paths are resolved relative to the directory where Snakemake is launched;
absolute input paths are recommended on an HPC system.

### 3. Validate and run locally
```bash
# Validate configuration and preview all jobs
snakemake --cores 1 --software-deployment-method conda --dry-run

# Run up to four jobs locally; dependencies are created automatically
snakemake --cores 4 --software-deployment-method conda
```

Do not run the second command on an HPC login node unless local computation is
explicitly allowed by the cluster administrators.

### 4. Run on a SLURM HPC cluster

The bundled profile submits each CGD scoring job through `sbatch` and requests
one CPU, 2 GB memory, and 30 minutes per job by default:

```bash
# Inspect the planned jobs without submitting anything
snakemake --profile profiles/slurm --dry-run

# Submit up to 50 jobs concurrently
snakemake --profile profiles/slurm
```

Override the default limits for a particular cluster or run without editing the
workflow:

```bash
snakemake --profile profiles/slurm \
  --jobs 100 \
  --set-resources cgd_score:mem_mb=4000 runtime=60
```

Output: `results/<sample>/<SCORE>.txt` (e.g. `results/demo/CGD.txt`,
`results/demo/CGDi.txt`, …), with per-job logs under `results/logs/`.

The submission host must provide `sbatch`, `sacct`, and a shared filesystem that
is visible from the compute nodes. Cluster-specific account and partition values
can be supplied on the command line, for example
`--slurm-account ACCOUNT --slurm-partition PARTITION`, or added to a personal
Snakemake profile outside this repository.

---

## Continuous Integration

Every push and pull request is validated by GitHub Actions
(`.github/workflows/ci.yml`), which:

- Installs dependencies (including ViennaRNA bindings) on Python 3.11 and 3.12
- Lints for syntax errors / undefined names
- Runs the unit test suite (`tests/`)
- Runs `CGD.py` directly on the bundled example FASTA
- Validates the Snakefile and configuration schema with a complete dry-run
- Creates the rule-specific Conda environment
- Runs the Conda-managed Snakemake workflow end to end and checks that every
  expected output file was produced

---

## Project Structure

```
CGD/
├── CGD.py                        # Main scoring tool
├── get_sequence.py               # Reverse complement utility
├── requirements.txt              # Python dependencies
├── workflow-env.yaml              # Snakemake and SLURM executor environment
├── README.md                     # This file
├── .gitignore                    # Git ignore rules
├── Snakefile                     # Snakemake workflow definition
├── config/
│   ├── config.yaml               # Snakemake sample/mode configuration
│   └── config.schema.yaml        # Configuration validation schema
├── envs/
│   └── environment.yaml          # Conda environment for --use-conda
├── profiles/slurm/
│   └── config.yaml               # Portable SLURM execution defaults
├── tests/
│   └── test_cgd.py               # Unit / regression tests
├── .github/workflows/
│   └── ci.yml                    # GitHub Actions CI pipeline
├── Test_dataset/                 # Example input FASTA files
└── Training_data/                # Training data used to build models
```

---

## Citation

Menon AV, Sohn JI, Nam JW. CGD: Comprehensive guide designer for CRISPR-Cas systems. Comput Struct Biotechnol J. 2020 Mar 25;18:814-820. doi: 10.1016/j.csbj.2020.03.020. PMID: 32308928; PMCID: PMC7152703.

---

## License

MIT License — see `LICENSE` file for details.

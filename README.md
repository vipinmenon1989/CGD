# CGD — Comprehensive Guide RNA Design

[![CI](https://github.com/vipinmenon1989/CGD/actions/workflows/ci.yml/badge.svg)](https://github.com/vipinmenon1989/CGD/actions/workflows/ci.yml)

**Author:** Vipin Menon, BIG Lab, Hanyang University (HYU)
**Contact:** a.vipin.menon@gmail.com
**Original:** September 2019 | **Modernized:** 2026
**Language:** Python 3.8+
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

### 2. Create a virtual environment (recommended)
```bash
python3 -m venv venv
source venv/bin/activate       # macOS / Linux
venv\Scripts\activate          # Windows
```

### 3. Install dependencies
```bash
pip install -r requirements.txt
```

> **ViennaRNA** must be installed separately via conda (recommended):
> ```bash
> conda install -c bioconda viennarna
> ```

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

### 1. Install Snakemake (in addition to `requirements.txt`)
```bash
pip install snakemake "pulp<2.8"
```

### 2. Configure samples and modes
Edit `config/config.yaml`:
```yaml
samples:
  demo: input.fa          # sample_name: path/to/input.fa

modes:
  - comprehensive          # CGD.py -a
  - crispri                # CGD.py -b
  - crispra                # CGD.py -c
  - cas9                   # CGD.py -d
  - cas12a                 # CGD.py -e
  - cas9_ng                # CGD.py -f

outdir: results
```

### 3. Run
```bash
# Preview the jobs Snakemake will run
snakemake --cores 1 -n

# Run using dependencies already on PATH
snakemake --cores 4

# Or let Snakemake manage an isolated conda environment (envs/environment.yaml)
snakemake --cores 4 --use-conda
```

Output: `results/<sample>/<SCORE>.txt` (e.g. `results/demo/CGD.txt`,
`results/demo/CGDi.txt`, …), with per-job logs under `results/logs/`.

---

## Continuous Integration

Every push and pull request is validated by GitHub Actions
(`.github/workflows/ci.yml`), which:

- Installs dependencies (including ViennaRNA bindings) on Python 3.10 and 3.11
- Lints for syntax errors / undefined names
- Runs the unit test suite (`tests/`)
- Runs `CGD.py` directly on the bundled example FASTA
- Runs the full Snakemake workflow end to end and checks that every expected
  output file was produced

---

## Project Structure

```
CGD/
├── CGD.py                        # Main scoring tool
├── get_sequence.py               # Reverse complement utility
├── requirements.txt              # Python dependencies
├── README.md                     # This file
├── .gitignore                    # Git ignore rules
├── Snakefile                     # Snakemake workflow definition
├── config/
│   └── config.yaml               # Snakemake sample/mode configuration
├── envs/
│   └── environment.yaml          # Conda environment for --use-conda
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

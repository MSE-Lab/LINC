# LINC

**LINC** (**L**arge-scale **I**ncremental **N**etwork-based **C**lustering) is a modern, high-performance bioinformatics framework for **incremental comparative genomics, Sequence Similarity Network (SSN) construction, ortholog family clustering (MCL)**, and **pan-genome analysis**.

---

## 🌟 Key Features

- **⚡ True Incremental Clustering**: Add new genomes without recomputing all-vs-all alignments across historical datasets. Order-D sorting and dual-cutoff decision trees ensure efficient, scalable joining.
- **🎯 Full Cross-Run Reproducibility**: Rigorously eliminates non-deterministic factors (thread scheduling ties, aligner output row jitter, hash randomization) to guarantee bit-level reproducible clustering results.
- **🔌 Pluggable Aligners & Clusterers**: Supports **MMseqs2** and **DIAMOND** as ultra-fast alignment backends, and system **MCL** for robust orthogroup refinement.
- **🔄 Fault-Tolerant Checkpoint & Resume**: Built-in stage-level state serialization (pickle + gzip) supports seamless recovery (`--resume`) and chainable incremental dataset expansion (`-a / --append`).
- **🧬 Built-in Pan-Genome Analysis**: Includes the `linc pan` toolkit to directly calculate gene presence/absence matrices, core/soft-core/shell/cloud partitions, and per-genome breakdown tables.

---

## 🏗️ Architecture & Workflow

LINC employs a **Hexagonal Architecture (Ports & Adapters)** decoupled into 5 sequential pipeline stages:

```
[FASTA Genomes]
       │
       ▼
┌─────────────────────────┐
│ 1. GenomeLoadStage      │  Load FASTA files, assign G00001 IDs & track original headers
└────────────┬────────────┘
             │
             ▼
┌─────────────────────────┐
│ 2. IntraGenomeStage     │  Self-alignment & k-center greedy coverage to pick local reps
└────────────┬────────────┘
             │
             ▼
┌─────────────────────────┐
│ 3. GlobalJoinStage      │  Order-D sorting & dual-cutoff decision tree against Global Rep Library
└────────────┬────────────┘
             │
             ▼
┌─────────────────────────┐
│ 4. SSNBuilderStage      │  Construct Sequence Similarity Network (none / rep-member / full)
└────────────┬────────────┘
             │
             ▼
┌─────────────────────────┐
│ 5. FamilyRefineStage    │  Connected components & MCL graph clustering -> Orthogroups.txt
└────────────┬────────────┘
             │
             ▼
    [Orthogroups.txt & gene_genome.tsv]
             │
             ▼ (optional)
┌─────────────────────────┐
│       linc pan          │  Pan-genome matrix & core/accessory classification
└─────────────────────────┘
```

---

## 📦 Installation & Deployment

### Recommended: Conda / Mamba Environment

LINC requires Python >= 3.11 and external bioconda binaries (`mmseqs2`, `diamond`, `mcl`). The fastest way to deploy on servers is via `environment.yml`:

```bash
# 1. Clone repository
git clone https://github.com/your-username/LINC.git
cd LINC

# 2. Create conda environment (contains mmseqs2, diamond, mcl, python dependencies)
conda env create -f environment.yml

# 3. Activate environment
conda activate linc

# 4. Install LINC in editable mode
pip install -e .
```

### Manual Installation

If you manage dependencies manually:

```bash
# Install external tools via conda/system package manager
conda install -c bioconda -c conda-forge mmseqs2 diamond mcl

# Install Python package
pip install -e .
```

Verify installation:
```bash
linc --help
```

---

## 🚀 Quick Start

### 1. Run Protein Family Clustering (`linc run`)

Given a directory of protein FASTA files (`.fa`, `.fasta`, `.faa`):

```bash
# Standard run using MMseqs2 backend
linc run -i /path/to/genomes -o output/ --aligner mmseqs --threads 8

# Using DIAMOND backend
linc run -i /path/to/genomes -o output_diamond/ --aligner diamond --threads 8
```

### 2. Resume Interrupted Runs (`--resume`)

If a run was stopped midway, resume from the last completed stage:

```bash
linc run -i /path/to/genomes -o output/ --resume
```

### 3. Incrementally Append New Genomes (`-a / --append`)

Add newly sequenced genomes into an existing clustering result without full re-computation:

```bash
linc run -i /path/to/genomes -a /path/to/new_genomes -o output/
```

### 4. Pan-Genome Analysis (`linc pan`)

Generate presence/absence matrix and core/accessory genome partition summary:

```bash
linc pan -g output/Orthogroups.txt -m output/gene_genome.tsv -o output/pan_genome/
```

---

## ⚙️ CLI Reference

### `linc run` Parameters

| Option | Default | Description |
| :--- | :--- | :--- |
| `-i, --input-dir` | *Required* | Directory containing historical/base FASTA genomes |
| `-o, --output-dir` | *Required* | Working directory for state cache and final output |
| `-a, --append` | `None` | Directory containing new genomes to append incrementally |
| `--resume / --no-resume` | `False` | Resume execution from existing checkpoint state |
| `--force / --no-force` | `False` | Force rerun on parameter mismatch during resume/append |
| `--aligner` | `mmseqs` | Aligner backend: `mmseqs` or `diamond` |
| `--threads` | `0` | Number of CPU threads (0 = auto-detect all cores) |
| `--cutoff-i` | `0.5` | Loose similarity threshold for intra-genome candidate grouping |
| `--cutoff-ii` | `0.8` | Strict high-confidence threshold for representative joining |
| `--cutoff-c` | `0.5` | Coverage cutoff for high-confidence match decisions |
| `--cov-mode` | `2` | Coverage mode: `0` (bidirectional), `1` (target), `2` (query) |
| `--max-per-cluster` | `5` | Maximum number of representatives allowed per local cluster |
| `--edge-completion-mode` | `representative_member` | SSN edge mode: `none`, `representative_member`, or `full` |
| `--inflation` | `1.5` | MCL inflation parameter (-I) for family refinement |
| `--edge-weight-mode` | `identity_x_coverage` | SSN edge weight mode: `identity`, `coverage`, `identity_x_coverage`, `bitscore` |

### `linc pan` Parameters

| Option | Default | Description |
| :--- | :--- | :--- |
| `-g, --orthogroups` | *Required* | Path to `Orthogroups.txt` produced by `linc run` |
| `-m, --gene-genome` | *Required* | Path to `gene_genome.tsv` produced by `linc run` |
| `-o, --output-dir` | *Required* | Destination directory for pan-genome analysis results |

---

## 📊 Output File Specifications

Running `linc run` generates the following files in `--output-dir`:

- **`Orthogroups.txt`**: Standard tab-delimited orthogroup table (Format: `OG0000000: gene_1 gene_2 ...`), directly compatible with OrthoFinder downstream pipelines.
- **`gene_genome.tsv`**: TSV mapping table mapping each original protein ID to its source genome FASTA filename.
- **`load_summary.json`**: Gene and genome counts loaded across datasets.
- **`intra_summary.json`**: Compression statistics and local cluster counts per genome.
- **`join_summary.json`**: Global join statistics and representative library metrics.
- **`ssn_summary.json`**: SSN node and edge completion statistics.
- **`family_summary.json`**: MCL clustering metrics and family size distribution.

Running `linc pan` outputs:

- **`presence_absence_matrix.tsv`**: Gene copy-number matrix across all genomes for each orthogroup.
- **`pan_genome_summary.json`**: Pan-genome classification breakdown:
  - **Core**: Present in $\ge 99\%$ of genomes
  - **Soft-Core**: Present in $[95\%, 99\%)$ of genomes
  - **Shell**: Present in $[15\%, 95\%)$ of genomes
  - **Cloud**: Present in $< 15\%$ of genomes
  - **Singleton**: Single-gene families
- **`per_genome_breakdown.tsv`**: Number of genes per category for each individual genome.

---

## 🧪 Testing & Verification

Run the comprehensive unit and integration test suite:

```bash
pytest tests/ -v
```

---

## 📚 Technical Documentation

For in-depth architectural and algorithmic specifications, refer to:

- [docs/LINC_技术框架设计_v2.md](docs/LINC_技术框架设计_v2.md) — Architectural design, contracts, and stage specifications
- [docs/开发进度.md](docs/开发进度.md) — Detailed development log and changelog
- [docs/归档/](docs/归档/) — Algorithm design notes and historical revisions

---

## 📄 License

This project is licensed under the MIT License.

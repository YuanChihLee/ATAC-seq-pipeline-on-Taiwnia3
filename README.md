# ATAC-seq-pipeline-on-Taiwnia3

This document describes an automated ATAC-seq data analysis pipeline. The pipeline starts with raw FASTQ files and performs a series of steps including sequence preprocessing, alignment, filtering, and finally generates analysis-ready BAM files, visualization-ready BigWig files, and peak files identifying open chromatin regions.

---

## 📜 Table of Contents
- [✨ Features](#-features)
- [⚙️ System Requirements](#️-system-requirements)
- [🚀 Installation and Setup](#-installation-and-setup)
  - [1. Software Dependencies](#1-software-dependencies)
  - [2. Reference Data Setup](#2-reference-data-setup)
  - [3. Script Configuration](#3-script-configuration)
- [▶️ How to Run](#️-how-to-run)
- [📂 Output Directory Structure](#-output-directory-structure)
- [📊 Quality Control](#-quality-control)
- [👁️ How to Visualize Results](#️-how-to-visualize-results)

---

## ✨ Features

* **Automated Workflow**: End-to-end processing from raw FASTQ to final peaks and signal tracks.
* **Parallel Processing**: Utilizes `GNU Parallel` to efficiently process multiple samples, leveraging the power of HPC environments.
* **Robust Filtering**: Implements a multi-step filtering process to remove low-quality reads, PCR duplicates, mitochondrial DNA, and reads mapping to blacklisted regions.
* **Comprehensive QC**: Generates a summary report with key quality metrics for each sample.
* **Reproducibility**: Designed for High-Performance Computing (HPC) environments using SLURM and environment modules.

---

## ⚙️ System Requirements

This pipeline is designed to run on a Linux-based HPC cluster with a SLURM workload manager.

* **Operating System**: Linux
* **Workload Manager**: SLURM
* **Shell**: Bash
* **Core Utilities**: `GNU Parallel`

---

## 🚀 Installation and Setup

### 1. Software Dependencies

The pipeline relies on several bioinformatics tools, which are loaded via environment modules in the `init_env` function. Ensure the following modules are available on your system:

* `Trimmomatic/0.39`
* `SAMTOOLS/1.18`
* `bowtie2/2.4.2`
* `BEDTOOLS/2.31.1`
* `Picard/2.27.4`
* `UCSC_Utilities` (specifically for `bedGraphToBigWig`)
* `MACS3` (must be installed locally, e.g., via `pip install macs3`)

### 2. Reference Data Setup

You need to prepare a dedicated directory (`pkg_dir`) containing the following reference files:

* **Bowtie2 Genome Index**: A pre-built index for your reference genome (e.g., mm39).
* **Chromosome Sizes File**: A two-column text file (`<chromosome_name>\t<size>`).
* **Blacklist Regions File**: A BED file containing regions known to produce artifacts.
* **Adapter Sequences**: FASTA files for adapter trimming, organized by index pairs (e.g., `v2_Ad1_.../v2_Ad2_....fa`).

### 3. Script Configuration

Before running, you must configure the global parameters at the top of the script:

```bash
# === Global Parameters ===
# --- User-defined paths ---
proj_root="/path/to/your/ATAC-seq/project"
pkg_dir="/path/to/your/reference_data"
raw_dir="$proj_root/ATAC"
filelist="$proj_root/filelist.txt"

# --- Reference Genome and Blacklist Files (mm39) ---
GENOME_INDEX="$pkg_dir/index/mouse/mm39"
GENOME_SIZE="mm" # Use 'hs' for human
CHROM_SIZES="$pkg_dir/mm39.chrom.sizes"
BLACKLIST_FILE="$pkg_dir/mm39.excluderanges_2.bed"

# --- Local MACS3 Installation ---
export MACS3_PATH="$HOME/.local/bin/macs3"
```

Create a `filelist.txt` with a tab-separated format: `SampleID\ti5_index\ti7_index`. For example:

```
Sample1	A1	B1
Sample2	A2	B2
```

---

## ▶️ How to Run

Submit the pipeline script to the SLURM scheduler using the `sbatch` command:

```bash
sbatch your_pipeline_script.sh
```

The script will automatically create the necessary output directories and process all samples listed in `filelist.txt` in parallel.

---

## 📂 Output Directory Structure

The pipeline generates a clean and organized directory structure within your `proj_root`:

```
.
├── 00_tmp/             # Temporary files for sorting
├── 01_trimmed/         # Adapter-trimmed and quality-filtered FASTQ files
├── 02_mapping/         # Intermediate SAM and BAM files
├── 03_final_bam/       # Final, analysis-ready BAM files
├── 04_peaks/           # Peak files called by MACS3
├── 05_bigwig/          # CPM-normalized signal tracks for visualization
├── logs/               # Log files for each step and sample
└── QC_reports/         # QC metrics and the final summary report
```

### Key Output Files Explained

* **`03_final_bam/{sample}.final.bam`**: The most important BAM file. It has been filtered to remove duplicates, mitochondrial DNA, and blacklisted regions. Use this for downstream analysis.
* **`04_peaks/{sample}_peaks.narrowPeak`**: Regions of open chromatin identified by MACS3.
* **`05_bigwig/{sample}.bw`**: CPM-normalized signal track. Ideal for visualizing chromatin accessibility in a genome browser like IGV.
* **`QC_reports/QC_summary.tsv`**: A summary table of key quality metrics for all samples.

---

## 📊 Quality Control

The pipeline automatically calculates several important QC metrics, which are compiled into **`QC_reports/QC_summary.tsv`**. This file includes:

* **AlignmentRate**: The percentage of reads successfully aligned to the genome.
* **DuplicationRate**: The fraction of reads marked as PCR duplicates by Picard.
* **MitoContaminationRate**: The percentage of reads mapping to the mitochondrial genome. High values can indicate excessive cell permeabilization.
* **FrIP_Score** (Fraction of Reads in Peaks): The percentage of all filtered reads that fall within the called peaks. A higher score generally indicates better signal-to-noise ratio.

---

## 👁️ How to Visualize Results

The generated BigWig (`.bw`) files are perfect for visualization in a genome browser. We recommend using the **[Integrative Genomics Viewer (IGV)](https://software.broadinstitute.org/software/igv/)**.

1.  **Launch IGV**.
2.  **Load Genome**: Select the appropriate reference genome (e.g., `mm39`) from the top-left dropdown menu.
3.  **Load Tracks**: Drag and drop the `.bw` files from the `05_bigwig/` directory onto the IGV window, or use the menu `File > Load from File...`.

This allows you to visually inspect chromatin accessibility at specific genomic loci and compare signals between different samples.

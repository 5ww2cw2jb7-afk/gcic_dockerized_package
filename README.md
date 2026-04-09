# GCIC Shiny App

This repository provides a Shiny-based web application for GCIC analysis.
The application is distributed as a Docker container to ensure full
reproducibility across operating systems.

---

## Associated Paper

This repository accompanies the following paper:

**"Gene-centered identification of cis-regulatory islands highlights regulatory landscapes complementary to motif-centric approaches"**  
*BMC Genomics, 2026 (accepted)*

---

## Overview

- Interactive Shiny application (R)
- Computational pipeline implemented in Python
- OS-independent execution via Docker
- Designed for reproducible research and publication

---

## Reproducibility (Docker)

### Prerequisites
- Docker Desktop (macOS / Windows)  
  or Docker Engine (Linux)

Docker Desktop must be running before executing the commands below.

### Build and Run

```bash
docker compose up --build -d
```

Open the application in a web browser:

```
http://localhost:3838/gcic/
```

### Outputs
- All result files are written to the host directory:
  ```
  ./results/
  ```
  (mounted into the container at runtime)

---

## Recorded Environment

- Docker: 29.1.3  
- Docker Compose: v2.40.3-desktop.1  
- Base image: rocker/shiny:4.3.2

### R
- R version 4.3.2 (2023-10-31)

### Python
- Python 3.10.12
- Python dependencies are fully pinned in `requirements.txt`

---

## Data Files

This repository does not include large input data files.

Required input FASTA and BED files (examples):

- gene_regions.5prime.fa
- gene_regions.bed
- ATgene_regions.5prime.fa
- ATgene_regions.bed

Please download these files from Zenodo (see below) and place them in the
repository root directory before running the application.

---

## Input data (FASTA/BED)

The FASTA and BED files required to run this application are available on Zenodo:

[Zenodo dataset (FASTA/BED files)](https://doi.org/10.5281/zenodo.18111654)

Please download and place them in the repository root directory before running the application.

---

## License

This software is released under the MIT License.

---

## Citation
If you use this software in your research, please cite it as described in
`CITATION.cff`.

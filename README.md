# GCIC Shiny App

This repository provides a Shiny-based web application for GCIC analysis.
The application is distributed as a Docker container to ensure full
reproducibility across operating systems.

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

### Build and Run

```bash
docker compose build --no-cache
docker compose up -d
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

Required input files (example):

- gene_regions.5prime.fa
- gene_regions.bed
- ATgene_regions.5prime.fa
- ATgene_regions.bed
- resources/
- analysis/

These files must be placed in the repository root directory.

---

## License
Specify your license here (e.g., GPL-3.0-or-later).

---

## Citation
If you use this software in your research, please cite it as described in
`CITATION.cff`.

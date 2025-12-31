\
# GCIC Single-Gene Shiny App - Reproducible container
# Base: R + shiny-server
FROM rocker/shiny:4.3.2

# --- Linux system dependencies ---
RUN apt-get update && apt-get install -y --no-install-recommends \
    python3 python3-pip \
 && rm -rf /var/lib/apt/lists/*

# --- App files ---
# shiny-server serves /srv/shiny-server/*
WORKDIR /srv/shiny-server/gcic
COPY . .

# --- R dependencies (minimal; can be replaced with renv for strict pinning) ---
RUN R -e "install.packages(c('shiny','readr','DT','fs'), repos='https://cloud.r-project.org')"

# --- Python dependencies ---
RUN python3 -m pip install --no-cache-dir -r requirements.txt

# Shiny Server default port
EXPOSE 3838

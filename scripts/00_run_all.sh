#!/bin/bash
# 00 — run the whole pipeline in order.
#
# Scripts sharing a number are ALTERNATIVES at the same pipeline position:
#   04_*  two stability probes (seeds, subsampling)
#   05_*  three parameter sweeps (resolution, KNN k, number of PCs)
#   06_*  four clustering algorithms on one fixed embedding
# 03_baseline_louvain.R must run before 05_* and 06_*, which reuse its
# Harmony embedding and SNN graph.
set -euo pipefail

export DATA_DIR="${DATA_DIR:-data}"
export OUT_DIR="${OUT_DIR:-results/tables}"
export FIG_DIR="${FIG_DIR:-results/figures}"
export COMMON_R="${COMMON_R:-scripts/common.R}"
RSCRIPT="${RSCRIPT:-Rscript}"
PYTHON="${PYTHON:-python3}"

bash   scripts/01_fetch_data.sh
$RSCRIPT scripts/02_build_label_tables.R
$RSCRIPT scripts/03_baseline_louvain.R

$RSCRIPT scripts/04_stability_seeds.R
$RSCRIPT scripts/04_stability_subsample.R

$RSCRIPT scripts/05_sweep_resolution.R
$RSCRIPT scripts/05_sweep_knn_k.R
$RSCRIPT scripts/05_sweep_npcs.R

$RSCRIPT scripts/06_alt_leiden.R
$RSCRIPT scripts/06_alt_kmeans.R
$RSCRIPT scripts/06_alt_hclust_ward.R
$RSCRIPT scripts/06_alt_gmm_mclust.R

$RSCRIPT scripts/07_collect_results.R
$PYTHON  scripts/08_figures.py
echo "pipeline complete"

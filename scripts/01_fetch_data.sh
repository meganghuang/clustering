#!/bin/bash
# 01 — fetch the authors' deposited metadata from Zenodo.
#
# Zheng et al. 2021, Science 374:eabe6474 — code and data: 10.5281/zenodo.5461803
# Only data.metaInfo.tar.gz (60 MB) is required. It contains:
#   CD8.pca.txt.gz          11,972 miniclusters x 15 PCs  (the authors' own PCA)
#   CD8.miniInfo.txt.gz     minicluster -> dataset / metacluster
#   int.CD8.S35.meta.tb.rds 234,550 cells x 21 cols, incl. ClusterID.harmony
#   int.CD8.genes.txt       the 1,500 informative genes used for integration
# data.expression.tar.gz (14 GB) and data.external.tar.gz (5.2 GB) are NOT needed.
set -euo pipefail
DATA_DIR="${DATA_DIR:-data}"
mkdir -p "$DATA_DIR"
cd "$DATA_DIR"
if [ ! -s data.metaInfo.tar.gz ]; then
  curl -L -o data.metaInfo.tar.gz \
    https://zenodo.org/api/records/5461803/files/data.metaInfo.tar.gz/content
fi
tar xzf data.metaInfo.tar.gz --strip-components=1   # -> ./metaInfo/
ls -l metaInfo/CD8.pca.txt.gz metaInfo/int.CD8.S35.meta.tb.rds metaInfo/int.CD8.genes.txt
echo "01 done"

## 02 — build the minicluster and cell label tables from the authors' metadata.
##
## ClusterID.harmony is the authors' RAW global Louvain partition (C00-C35);
## meta.cluster is their hand-curated merge of those 36 into 17 named groups
## (see w.step0.01.ann.CD8.R in the Zenodo code tarball). Both are kept, so the
## reproduction can be scored against the algorithm's output, not only the
## curated relabelling.
source(Sys.getenv("COMMON_R", unset = "scripts/common.R"))

x <- as.data.table(readRDS(file.path(MI, "int.CD8.S35.meta.tb.rds")))
cat("cells:", nrow(x), " miniclusters:", uniqueN(x$miniCluster),
    " datasets:", uniqueN(x$dataset), " cancerTypes:", uniqueN(x$cancerType), "\n")

lab <- unique(x[, .(miniCluster, dataset, dataset.tech, ClusterID.harmony, meta.cluster)])
stopifnot(nrow(lab) == uniqueN(lab$miniCluster))

## dominant cancer type per minicluster (a few span more than one)
ct <- x[, .N, by = .(miniCluster, cancerType)][order(miniCluster, -N)][
        , .SD[1], by = miniCluster][, .(miniCluster, cancerType_major = cancerType)]
lab <- merge(lab, ct, by = "miniCluster")
lab <- merge(lab, x[, .(n_cells = .N), by = miniCluster], by = "miniCluster")

fwrite(lab, file.path(OUT_DIR, "CD8.minicluster_labels.csv"))
fwrite(x[, .(cellID.uniq, miniCluster, dataset, cancerType, ClusterID.harmony, meta.cluster)],
       file.path(OUT_DIR, "CD8.cell_labels.csv.gz"), compress = "gzip")

## sanity: the curated merge must be many-to-one on the raw partition
chk <- unique(x[, .(ClusterID.harmony, meta.cluster)])[
         , .(n_meta = uniqueN(meta.cluster)), by = ClusterID.harmony]
cat("raw clusters mapping to >1 metacluster:", sum(chk$n_meta > 1), "of", nrow(chk), "\n")
cat("02 done -", nrow(lab), "miniclusters,", sum(lab$n_cells), "cells\n")

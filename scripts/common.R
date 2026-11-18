## common.R — shared setup for the CD8 clustering reproduction
## Sourced by every numbered stage script. Defines inputs, the pipeline
## function, and ARI helpers. Parameters match Zheng et al. 2021 exactly.

if (nzchar(Sys.getenv("R_LIBS_USER"))) .libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths()))
suppressMessages({
  library(data.table); library(Seurat); library(harmony)
  library(igraph); library(aricode); library(Matrix)
})

## aricode::ARI coerces character vectors to integer, turning labels like
## "C00".."C35" into NA. Mask it with a factor-safe wrapper.
ARI <- function(a, b) aricode::ARI(as.integer(factor(a)), as.integer(factor(b)))

DATA_DIR <- Sys.getenv("DATA_DIR", unset = "data")
OUT_DIR  <- Sys.getenv("OUT_DIR",  unset = "results/tables")
MI       <- file.path(DATA_DIR, "metaInfo")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

## ---- published parameters (Zenodo 10.5281/zenodo.5461803) ----
PAR <- list(npc = 15, k = 20, res = 2, seed = 1, algorithm = 1)  # Louvain = algorithm 1

load_inputs <- function() {
  pca <- fread(file.path(MI, "CD8.pca.txt.gz"))
  setnames(pca, 1, "miniCluster")
  stopifnot(ncol(pca) == 16, all(paste0("PC_", 1:15) %in% colnames(pca)))
  emb <- as.matrix(pca[, paste0("PC_", 1:15), with = FALSE])
  rownames(emb) <- pca$miniCluster

  lab <- fread(file.path(OUT_DIR, "CD8.minicluster_labels.csv"), header = TRUE)
  setkey(lab, "miniCluster")
  lab <- lab[rownames(emb), ]
  stopifnot(identical(lab$miniCluster, rownames(emb)), !any(is.na(lab$dataset)))

  list(emb   = emb,
       meta  = data.frame(dataset = lab$dataset, row.names = lab$miniCluster,
                          stringsAsFactors = FALSE),
       raw36 = as.character(lab$ClusterID.harmony),   # authors' raw Louvain partition
       meta17= as.character(lab$meta.cluster),        # authors' curated metaclusters
       wts   = as.numeric(lab$n_cells))               # cells per minicluster
}

## Exact cell-level ARI from a weighted contingency table. All cells in a
## minicluster share a label, so this is identical to expanding to 234,550 rows.
ari_tab <- function(tab) {
  ch2 <- function(x) x * (x - 1) / 2
  n <- sum(tab); a <- rowSums(tab); b <- colSums(tab)
  sij <- sum(ch2(tab)); sa <- sum(ch2(a)); sb <- sum(ch2(b)); e <- sa * sb / ch2(n)
  (sij - e) / (0.5 * (sa + sb) - e)
}
ari_cell <- function(cl, ref, w) ari_tab(as.matrix(xtabs(w ~ cl + ref)))

## One pass of the published pipeline, starting from the authors' minicluster PCA.
## The seed enters at RunHarmony: the PCA is inherited from the authors rather
## than recomputed, so Harmony's k-means initialisation is the earliest
## stochastic step. Seeding FindClusters alone on a fixed graph returns ARI ~ 1.
run_pipe <- function(D, npc = PAR$npc, k = PAR$k, res = PAR$res,
                     seed = PAR$seed, idx = NULL, algorithm = PAR$algorithm) {
  e <- D$emb[, 1:npc, drop = FALSE]; m <- D$meta
  if (!is.null(idx)) { e <- e[idx, , drop = FALSE]; m <- m[idx, , drop = FALSE] }
  set.seed(seed)
  hm <- harmony::RunHarmony(e, meta_data = m, vars_use = "dataset", verbose = FALSE)
  nn <- Seurat::FindNeighbors(hm, k.param = k, verbose = FALSE)
  cl <- Seurat::FindClusters(nn$snn, resolution = res, algorithm = algorithm,
                             random.seed = seed, verbose = FALSE)
  list(hm = hm, snn = nn$snn, cl = as.character(cl[[1]]))
}

score <- function(cl, D, arm, setting, idx = NULL) {
  r36 <- if (is.null(idx)) D$raw36  else D$raw36[idx]
  r17 <- if (is.null(idx)) D$meta17 else D$meta17[idx]
  w   <- if (is.null(idx)) D$wts    else D$wts[idx]
  data.table(arm = arm, setting = setting, n_clusters = uniqueN(cl),
             ari_raw36_mini  = ARI(cl, r36),
             ari_meta17_mini = ARI(cl, r17),
             ari_raw36_cell  = ari_cell(cl, r36, w),
             ari_meta17_cell = ari_cell(cl, r17, w))
}

write_rows <- function(dt, file) {
  fwrite(dt, file.path(OUT_DIR, file))
  cat("wrote", file.path(OUT_DIR, file), "-", nrow(dt), "rows\n")
}

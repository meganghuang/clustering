## 05 — PARAMETER SWEEP (b): k in the KNN graph.
## Rebuilds FindNeighbors only; the Harmony embedding is reused from the baseline.
source(Sys.getenv("COMMON_R", unset = "scripts/common.R"))
D <- load_inputs()
hm <- as.matrix(fread(file.path(OUT_DIR, "harmony_embedding.csv.gz"))[, -"miniCluster"])
rownames(hm) <- rownames(D$emb)

rows <- list()
for (k in c(5, 10, 20, 30, 50, 100)) {
  nn <- Seurat::FindNeighbors(hm, k.param = k, verbose = FALSE)
  v  <- as.character(Seurat::FindClusters(nn$snn, resolution = PAR$res,
                      algorithm = PAR$algorithm, random.seed = PAR$seed, verbose = FALSE)[[1]])
  rows[[length(rows)+1]] <- score(v, D, "sweep_k", paste0("k=", k))
  cat(sprintf("  k=%3d clusters=%d ARI_raw36=%.4f\n", k, uniqueN(v), ARI(v, D$raw36)))
}
write_rows(rbindlist(rows), "ari_sweep_k.csv")
cat("05-k done\n")

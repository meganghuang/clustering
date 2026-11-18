## 03 — BASELINE: reproduce the published clustering exactly.
##
## Published settings (Zenodo code):
##   RunPCA npcs = 15                      [inherited: the authors' CD8.pca.txt.gz]
##   RunHarmony(vars_use = "dataset")
##   FindNeighbors(dims = 1:15, k.param = 20)
##   FindClusters(resolution = 2, algorithm = 1)   <- Louvain, the authors' choice
## The authors selected RNA_snn_res.2 for CD8 (res 2.6 for CD4).
source(Sys.getenv("COMMON_R", unset = "scripts/common.R"))

D <- load_inputs()
cat(sprintf("miniclusters=%d PCs=%d datasets=%d raw36=%d meta17=%d cells=%d\n",
            nrow(D$emb), ncol(D$emb), uniqueN(D$meta$dataset),
            uniqueN(D$raw36), uniqueN(D$meta17), sum(D$wts)))

t0 <- Sys.time()
base <- run_pipe(D)
cat("wall:", round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1), "s\n")
cat("clusters:", uniqueN(base$cl), "\n")

row <- score(base$cl, D, "baseline", "npc=15,k=20,res=2,seed=1")
print(row)
write_rows(row, "ari_baseline.csv")

fwrite(data.table(miniCluster = rownames(D$emb), cluster = base$cl,
                  ClusterID.harmony = D$raw36, meta.cluster = D$meta17, n_cells = D$wts),
       file.path(OUT_DIR, "baseline_assignments.csv"))
fwrite(as.data.table(as.matrix(xtabs(D$wts ~ base$cl + D$meta17)), keep.rownames = "cluster"),
       file.path(OUT_DIR, "confusion_baseline_vs_meta17_cells.csv"))
fwrite(as.data.table(as.matrix(xtabs(D$wts ~ base$cl + D$raw36)), keep.rownames = "cluster"),
       file.path(OUT_DIR, "confusion_baseline_vs_raw36_cells.csv"))

hm <- as.data.table(base$hm); hm[, miniCluster := rownames(D$emb)]
fwrite(hm, file.path(OUT_DIR, "harmony_embedding.csv.gz"), compress = "gzip")
saveRDS(base$snn, file.path(OUT_DIR, "baseline_snn.rds"))

if (requireNamespace("uwot", quietly = TRUE)) {
  set.seed(1); um <- uwot::umap(base$hm, n_neighbors = 30, min_dist = 0.3)
  fwrite(data.table(miniCluster = rownames(D$emb), UMAP_1 = um[,1], UMAP_2 = um[,2]),
         file.path(OUT_DIR, "umap_harmony.csv.gz"), compress = "gzip")
}
cat("03 done\n")

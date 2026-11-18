## 06 — ALTERNATIVE CLUSTERING METHOD: k-means (centroid-based)
## All 06_* scripts are alternatives occupying the SAME pipeline position:
## they cluster the identical Harmony embedding produced by 03_baseline_louvain.R,
## so any difference is attributable to the algorithm alone.
source(Sys.getenv("COMMON_R", unset = "scripts/common.R"))
D  <- load_inputs()
hm <- as.matrix(fread(file.path(OUT_DIR, "harmony_embedding.csv.gz"))[, -"miniCluster"])
rownames(hm) <- rownames(D$emb)
K  <- uniqueN(D$raw36)   # 36, matched to the authors' raw partition

set.seed(1)
v <- as.character(stats::kmeans(hm, centers = K, nstart = 25, iter.max = 100)$cluster)
cat("kmeans clusters:", uniqueN(v), "\n")
write_rows(score(v, D, "algorithm", "kmeans"), "ari_alg_kmeans.csv")
fwrite(data.table(miniCluster = rownames(hm), kmeans = v),
       file.path(OUT_DIR, "assign_kmeans.csv.gz"), compress = "gzip")

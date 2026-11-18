## 06 — ALTERNATIVE CLUSTERING METHOD: Ward hierarchical (connectivity-based)
## All 06_* scripts are alternatives occupying the SAME pipeline position:
## they cluster the identical Harmony embedding produced by 03_baseline_louvain.R,
## so any difference is attributable to the algorithm alone.
source(Sys.getenv("COMMON_R", unset = "scripts/common.R"))
D  <- load_inputs()
hm <- as.matrix(fread(file.path(OUT_DIR, "harmony_embedding.csv.gz"))[, -"miniCluster"])
rownames(hm) <- rownames(D$emb)
K  <- uniqueN(D$raw36)   # 36, matched to the authors' raw partition

v <- as.character(stats::cutree(stats::hclust(dist(hm), method = "ward.D2"), k = K))
cat("hclust clusters:", uniqueN(v), "\n")
write_rows(score(v, D, "algorithm", "hclust_ward"), "ari_alg_hclust_ward.csv")
fwrite(data.table(miniCluster = rownames(hm), hclust_ward = v),
       file.path(OUT_DIR, "assign_hclust_ward.csv.gz"), compress = "gzip")

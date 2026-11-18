## 06 — ALTERNATIVE CLUSTERING METHOD: Gaussian mixture model (model-based)
## All 06_* scripts are alternatives occupying the SAME pipeline position:
## they cluster the identical Harmony embedding produced by 03_baseline_louvain.R,
## so any difference is attributable to the algorithm alone.
source(Sys.getenv("COMMON_R", unset = "scripts/common.R"))
D  <- load_inputs()
hm <- as.matrix(fread(file.path(OUT_DIR, "harmony_embedding.csv.gz"))[, -"miniCluster"])
rownames(hm) <- rownames(D$emb)
K  <- uniqueN(D$raw36)   # 36, matched to the authors' raw partition

suppressMessages(library(mclust))
## mclust initialises from a RANDOM subset when n > mclust.options("subset")
## (default 2000; here n = 11,972), so Mclust is seed-dependent. Seed it.
set.seed(1)
## VVI: diagonal, varying volume and shape. VVV is intractable at G = 36 in 15 dims.
g <- mclust::Mclust(hm, G = K, modelNames = "VVI", verbose = FALSE)
v <- as.character(g$classification)
cat("gmm clusters:", uniqueN(v), "\n")
write_rows(score(v, D, "algorithm", "gmm_mclust_VVI"), "ari_alg_gmm.csv")
fwrite(data.table(miniCluster = rownames(hm), gmm_mclust_VVI = v),
       file.path(OUT_DIR, "assign_gmm.csv.gz"), compress = "gzip")

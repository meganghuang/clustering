## 06 — ALTERNATIVE CLUSTERING METHOD: Leiden (graph / modularity)
## All 06_* scripts are alternatives occupying the SAME pipeline position:
## they cluster the identical Harmony embedding produced by 03_baseline_louvain.R,
## so any difference is attributable to the algorithm alone.
source(Sys.getenv("COMMON_R", unset = "scripts/common.R"))
D  <- load_inputs()
hm <- as.matrix(fread(file.path(OUT_DIR, "harmony_embedding.csv.gz"))[, -"miniCluster"])
rownames(hm) <- rownames(D$emb)
K  <- uniqueN(D$raw36)   # 36, matched to the authors' raw partition

## cluster_leiden's refinement is stochastic; seed it so the result is reproducible.
set.seed(1)
snn <- readRDS(file.path(OUT_DIR, "baseline_snn.rds"))
gr  <- igraph::graph_from_adjacency_matrix(snn, mode = "undirected", weighted = TRUE, diag = FALSE)
## igraph >= 2.1.0 renamed resolution_parameter -> resolution; support both.
leiden_args <- names(formals(igraph::cluster_leiden))
le <- if ("resolution" %in% leiden_args) {
  igraph::cluster_leiden(gr, objective_function = "modularity",
                         resolution = PAR$res, n_iterations = 10)
} else {
  igraph::cluster_leiden(gr, objective_function = "modularity",
                         resolution_parameter = PAR$res, n_iterations = 10)
}
v <- as.character(igraph::membership(le))
cat("leiden clusters:", uniqueN(v), "\n")
write_rows(score(v, D, "algorithm", "leiden_igraph"), "ari_alg_leiden.csv")
fwrite(data.table(miniCluster = rownames(hm), leiden_igraph = v),
       file.path(OUT_DIR, "assign_leiden.csv.gz"), compress = "gzip")

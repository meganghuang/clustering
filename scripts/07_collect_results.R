## 07 — collect every stage's ARI rows into one table, plus between-method ARI.
source(Sys.getenv("COMMON_R", unset = "scripts/common.R"))
D <- load_inputs()

parts <- c("ari_baseline.csv", "ari_seeds.csv", "ari_subsample.csv",
           "ari_sweep_resolution.csv", "ari_sweep_k.csv", "ari_sweep_npc.csv",
           "ari_alg_leiden.csv", "ari_alg_kmeans.csv",
           "ari_alg_hclust_ward.csv", "ari_alg_gmm.csv")
have <- parts[file.exists(file.path(OUT_DIR, parts))]
if (length(have) < length(parts))
  warning("missing stage outputs: ", paste(setdiff(parts, have), collapse = ", "))
all <- rbindlist(lapply(file.path(OUT_DIR, have), fread), fill = TRUE)

## the baseline Louvain run is itself one of the compared algorithms
bl <- all[arm == "baseline"][, `:=`(arm = "algorithm", setting = "louvain_baseline")]
all <- rbind(all, bl)
write_rows(all, "ari_results.csv")

## pairwise ARI between the methods that ran
base_as <- fread(file.path(OUT_DIR, "baseline_assignments.csv"))
alg <- list(louvain_baseline = setNames(base_as$cluster, base_as$miniCluster))
alt_files <- c(leiden_igraph  = "assign_leiden.csv.gz",
               kmeans         = "assign_kmeans.csv.gz",
               hclust_ward    = "assign_hclust_ward.csv.gz",
               gmm_mclust_VVI = "assign_gmm.csv.gz")
for (nm_i in names(alt_files)) {
  p <- file.path(OUT_DIR, alt_files[[nm_i]])
  if (file.exists(p)) { t <- fread(p); alg[[nm_i]] <- setNames(t[[2]], t$miniCluster) }
}

nm <- names(alg); ref <- names(alg[[1]])
pa <- CJ(a = nm, b = nm)[a < b]
pa[, ari := mapply(function(x, y) ARI(alg[[x]][ref], alg[[y]][ref]), a, b)]
write_rows(pa, "algorithm_pairwise_ari.csv")

cat("\n== headline ==\n")
print(all[arm %in% c("baseline", "algorithm"),
          .(setting, n_clusters, ari_raw36_mini, ari_meta17_mini, ari_raw36_cell)])
## export the co-assignment matrix ordered by published metacluster, for 08_figures.py
p <- file.path(OUT_DIR, "coassignment_2000.rds")
if (file.exists(p)) {
  L <- readRDS(p); ord <- order(as.character(L$meta17))
  fwrite(data.table(idx = seq_along(ord), meta.cluster = as.character(L$meta17)[ord]),
         file.path(OUT_DIR, "coassign_order_labels.csv"))
  con <- gzfile(file.path(OUT_DIR, "coassign_ordered.csv.gz"), "w")
  write.table(round(L$coassign[ord, ord], 4), con, sep = ",",
              row.names = FALSE, col.names = FALSE)
  close(con)
  cat("co-assignment matrix exported\n")
}

writeLines(capture.output(sessionInfo()), file.path(OUT_DIR, "sessionInfo.txt"))
cat("07 done\n")

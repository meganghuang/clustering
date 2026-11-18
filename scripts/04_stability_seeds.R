## 04 — STABILITY (a): random seeds, data held fixed.
##
## The seed is varied from RunHarmony onward. Reseeding FindClusters alone on a
## fixed graph returns ARI ~ 1 and measures nothing; Harmony's k-means
## initialisation is the earliest stochastic step once the PCA is inherited.
source(Sys.getenv("COMMON_R", unset = "scripts/common.R"))
D <- load_inputs()

seeds <- 1:10
cls <- list(); rows <- list()
for (s in seeds) {
  r <- run_pipe(D, seed = s)
  cls[[as.character(s)]] <- r$cl
  rows[[length(rows)+1]] <- score(r$cl, D, "seed", paste0("seed=", s))
  cat(sprintf("  seed %2d: clusters=%d ARI_raw36=%.4f\n", s, uniqueN(r$cl), ARI(r$cl, D$raw36)))
}
write_rows(rbindlist(rows), "ari_seeds.csv")

pw <- CJ(i = seeds, j = seeds)[i < j]
pw[, ari := mapply(function(a, b) ARI(cls[[as.character(a)]], cls[[as.character(b)]]), i, j)]
write_rows(pw, "seed_pairwise_ari.csv")
cat(sprintf("seed pairwise ARI: min=%.4f median=%.4f max=%.4f\n",
            min(pw$ari), median(pw$ari), max(pw$ari)))
cat("04-seeds done\n")

## 04 — STABILITY (b): subsampling WITHOUT replacement.
##
## Bootstrap resampling is deliberately NOT used: duplicated observations sit at
## distance zero and corrupt the SNN graph, inflating apparent stability.
## Resampling unit = MINICLUSTERS. Cells cannot be the unit here, because the
## minicluster definition is inherited from the authors and cannot be re-derived
## without the 14 GB expression deposit.
## Pairwise ARI is computed only on miniclusters shared between a pair of replicates.
source(Sys.getenv("COMMON_R", unset = "scripts/common.R"))
D <- load_inputs()

nrep <- 10; frac <- 0.8
cls <- list(); rows <- list()
for (r in 1:nrep) {
  set.seed(1000 + r)
  idx <- sort(sample.int(nrow(D$emb), size = floor(frac * nrow(D$emb)), replace = FALSE))
  rr <- run_pipe(D, seed = 1000 + r, idx = idx)
  cls[[r]] <- setNames(rr$cl, rownames(D$emb)[idx])
  rows[[length(rows)+1]] <- score(rr$cl, D, "subsample",
                                  paste0("rep=", r, ",frac=", frac), idx = idx)
  cat(sprintf("  rep %2d: n=%d clusters=%d ARI_raw36=%.4f\n",
              r, length(idx), uniqueN(rr$cl), ARI(rr$cl, D$raw36[idx])))
}
write_rows(rbindlist(rows), "ari_subsample.csv")

pw <- CJ(i = 1:nrep, j = 1:nrep)[i < j]
pw[, ari := mapply(function(a, b) {
      sh <- intersect(names(cls[[a]]), names(cls[[b]])); ARI(cls[[a]][sh], cls[[b]][sh]) }, i, j)]
pw[, n_shared := mapply(function(a, b) length(intersect(names(cls[[a]]), names(cls[[b]]))), i, j)]
write_rows(pw, "subsample_pairwise_ari.csv")
cat(sprintf("subsample pairwise ARI: min=%.4f median=%.4f max=%.4f\n",
            min(pw$ari), median(pw$ari), max(pw$ari)))

## co-assignment on a random 2,000-minicluster subset (all 11,972 would be 143M pairs)
set.seed(7); keep <- sort(sample.int(nrow(D$emb), 2000)); kn <- rownames(D$emb)[keep]
co <- matrix(0, 2000, 2000, dimnames = list(kn, kn)); de <- co
for (r in 1:nrep) {
  v <- cls[[r]][kn]; ok <- !is.na(v)
  de[ok, ok] <- de[ok, ok] + 1
  co[ok, ok] <- co[ok, ok] + outer(v[ok], v[ok], "==")
}
saveRDS(list(coassign = ifelse(de > 0, co / de, NA_real_),
             meta17 = D$meta17[keep], raw36 = D$raw36[keep], names = kn),
        file.path(OUT_DIR, "coassignment_2000.rds"))
cat("04-subsample done\n")

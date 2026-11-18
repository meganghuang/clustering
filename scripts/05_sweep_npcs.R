## 05 — PARAMETER SWEEP (c): number of principal components.
## The only sweep that forces a full redo: Harmony, the graph and the clustering
## are all recomputed at each setting. Capped at 15, the number of PCs deposited.
source(Sys.getenv("COMMON_R", unset = "scripts/common.R"))
D <- load_inputs()

rows <- list()
for (np in c(5, 8, 10, 12, 15)) {
  r <- run_pipe(D, npc = np)
  rows[[length(rows)+1]] <- score(r$cl, D, "sweep_npc", paste0("npc=", np))
  cat(sprintf("  npc=%2d clusters=%d ARI_raw36=%.4f\n", np, uniqueN(r$cl), ARI(r$cl, D$raw36)))
}
write_rows(rbindlist(rows), "ari_sweep_npc.csv")
cat("05-npc done\n")

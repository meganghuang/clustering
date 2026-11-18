## 05 — PARAMETER SWEEP (a): clustering resolution.
## Nearly free: FindClusters accepts a vector of resolutions on one fixed graph,
## so the Harmony run and the SNN graph are reused from the baseline.
source(Sys.getenv("COMMON_R", unset = "scripts/common.R"))
D <- load_inputs()
snn <- readRDS(file.path(OUT_DIR, "baseline_snn.rds"))

resv <- seq(0.2, 5, 0.2)
cl <- Seurat::FindClusters(snn, resolution = resv, algorithm = PAR$algorithm,
                           random.seed = PAR$seed, verbose = FALSE)
rows <- lapply(seq_along(resv), function(i)
          score(as.character(cl[[i]]), D, "sweep_resolution", paste0("res=", resv[i])))
write_rows(rbindlist(rows), "ari_sweep_resolution.csv")
cat("05-resolution done -", length(resv), "settings\n")

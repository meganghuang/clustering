# Clustering methods report

Reproduction and stress-test of the CD8 T cell clustering in
**Zheng et al. 2021, *Science* 374:eabe6474**, "Pan-cancer single cell landscape of
tumor-infiltrating T cells".

---

## 1. Scope

| | |
|---|---|
| Compartment | CD8 T cells |
| Cells | **234,550** |
| Clustering units (miniclusters) | **11,972** |
| Datasets | 35 |
| Cancer types | 22 |
| Published raw clusters (`ClusterID.harmony`) | **36** (C00–C35) |
| Published curated metaclusters (`meta.cluster`) | **17** (`CD8.c01` – `CD8.c17`) |

All clustering operates on **miniclusters, not cells**. This is the authors' own design:
cells are sub-clustered within each dataset at very high resolution (`res.hi = 50`) and each
minicluster is reduced to its mean expression vector. A cell inherits its minicluster's label.
Minicluster size here is 2–130 cells (median 16).

### Two reference levels, and why both are reported

`meta.cluster` is **not** a clustering output. In the authors' `w.step0.01.ann.CD8.R` the 36
Louvain clusters are merged **by hand** into 17 named groups:

```r
sce$meta.cluster[sce$ClusterID %in% c("C16","C04")]                   <- "CD8.c01.Tn.MAL"
sce$meta.cluster[sce$ClusterID %in% c("C15","C10","C33","C06","C28")] <- "CD8.c02.Tm.IL7R"
```

GEO (GSE156728) deposits only the merged `meta.cluster`. The **raw** 36-cluster partition is
deposited on Zenodo, as the `ClusterID.harmony` column of `int.CD8.S35.meta.tb.rds`. We verified
that the raw→curated map is many-to-one for **36 of 36** clusters, i.e. a pure merge with no
cluster split across metaclusters.

Because the raw partition is available, **the headline number is ARI against the 36-cluster
algorithmic output**. No majority-vote merge to 17 groups is used anywhere; such a merge would be
defined using the very labels it is then scored against.

---

## 2. Starting point and published parameters

We start from the authors' **own deposited minicluster PCA** (`CD8.pca.txt.gz`,
11,972 × 15 PCs) rather than re-deriving it from counts. This inherits their normalization,
minicluster construction and gene selection exactly, so differences measure the **clustering**
step rather than preprocessing. The 14 GB expression deposit is therefore not required.

Parameters, taken from the authors' deposited code (Zenodo 10.5281/zenodo.5461803) and their
`scPip` package, not from the paper's prose:

| Stage | Setting |
|---|---|
| Per-dataset Seurat | 15 PCs, Louvain resolution 2 |
| Minicluster construction | within-dataset sub-clustering at `res.hi = 50`, mean expression |
| Informative genes | 1,500; occupancy ≥ 0.85; blacklist `exclude.gene.misc.human.v4` + `^RP[LS]` + `MALAT1` |
| Covariates regressed | `S.Score`, `G2M.Score`, `DIG.Score1` |
| PCA | `npcs = 15` |
| Batch correction | `RunHarmony(vars_use = "dataset")` |
| Graph | `FindNeighbors(dims = 1:15, k.param = 20)` |
| Clustering | `FindClusters(resolution = 2, algorithm = 1)` → Louvain |

The authors selected `RNA_snn_res.2` for CD8 (and `res.2.6` for CD4).

---

## 3. Clustering methods compared

Five methods, all run on the **same Harmony embedding**, so every difference is attributable to
the algorithm alone. Non-graph methods are fixed at 36 clusters to match the published count.

| Method | Family | Script | Key settings | Clusters | ARI vs published 36 | ARI vs curated 17 |
|---|---|---|---|---|---|---|
| **Louvain** (the authors' method) | graph / modularity | `03_baseline_louvain.R` | SNN, k=20, res=2 | 36 | **0.496** | 0.362 |
| **Leiden** | graph / modularity | `06_alt_leiden.R` | same SNN graph, resolution 2, 10 iterations, seeded | 34 | **0.527** | 0.377 |
| **k-means** | centroid | `06_alt_kmeans.R` | k=36, `nstart=25` | 36 | 0.344 | 0.247 |
| **Ward hierarchical** | connectivity | `06_alt_hclust_ward.R` | Euclidean, `ward.D2`, cut at 36 | 36 | 0.318 | 0.240 |
| **Gaussian mixture** | model-based | `06_alt_gmm_mclust.R` | `mclust`, `G=36`, `modelNames="VVI"`, seeded | 36 | 0.298 | 0.245 |

**The method family dominates every parameter choice.** The two graph methods agree with each
other at ARI 0.936, but with k-means, Ward and GMM at only 0.26–0.44. Leiden actually
recovers the authors' Louvain partition **better than our own Louvain run does**
(0.527 vs 0.496) — the two optimise the same modularity
objective, and Leiden's refinement step makes it the less erratic of the pair.

GMM is the weakest (0.298). A diagonal-covariance mixture in 15 Harmony
dimensions assumes roughly ellipsoidal components of comparable scale, which this embedding does
not provide — the published metaclusters span three orders of magnitude in size
(1,427 to 38,415 cells).

---

## 4. Reproduction result

Re-running the published pipeline at the published settings gives **36 clusters,
exactly matching the authors' 36**, but only partial agreement on membership:

| Comparison | Minicluster level | Cell level |
|---|---|---|
| vs published raw 36 clusters | **0.4956** | **0.5126** |
| vs published curated 17 metaclusters | 0.3623 | 0.3768 |

An ARI near 0.5 means the partitions share broad structure but disagree substantially on
boundaries. Section 5 shows this gap is **not** evidence of an error: it is the same size as the
run-to-run variation of the method itself.

---

## 5. Stability

**Seeds (data fixed).** 10 seeds, varied from `RunHarmony` onward — Harmony's k-means
initialisation is the earliest stochastic step once the PCA is inherited. Reseeding
`FindClusters` alone on a fixed graph returns ARI ≈ 1 and measures nothing.

- cluster count ranged **32–37**
- ARI vs published: 0.4607–0.5134
- **pairwise ARI between two of our own runs: 0.4268–0.5204, median 0.4777**

**Subsampling without replacement.** 10 replicates at 80% of miniclusters (n = 9,577).
Bootstrap resampling is deliberately excluded: duplicated observations sit at distance zero and
corrupt the SNN graph, inflating apparent stability. The resampling unit is the **minicluster** —
cells cannot be the unit, because the minicluster definition is inherited from the authors.
Pairwise ARI uses only miniclusters shared between a pair of replicates (7,633–7,713).

- cluster count ranged **34–40**
- pairwise ARI: 0.3724–0.4882, median 0.4331

**The key comparison.** Our reproduction agrees with the published partition at
0.496. Two of our own runs, differing only in random seed, agree with each other at a
median of 0.478. **The reproduction gap is no larger than the method's own
irreproducibility**, so the published partition is one draw from a distribution of partitions the
pipeline produces, not a uniquely determined result.

---

## 6. Parameter sensitivity

| Parameter | Published | Range swept | ARI vs raw 36 | Best setting |
|---|---|---|---|---|
| Clustering resolution | 2 | 0.2 – 5.0 (25 settings) | 0.233 – 0.518 | res=2.6 (40 clusters) |
| k in the KNN graph | 20 | 5 – 100 | 0.449 – 0.497 | k=30 (31 clusters) |
| Number of PCs | 15 | 5 – 15 | 0.314 – 0.496 | npc=15 |

**Number of PCs is the most damaging parameter.** ARI falls monotonically from
0.496 at 15 PCs to 0.314 at 5 — dropping representation costs more than
any clustering setting. The sweep is capped at 15 because that is what the authors deposited.

**Resolution and k are flat-topped.** ARI vs the raw partition stays within ~0.03 of its maximum
across resolution 1.2–3.4 and k 20–30. The published settings are not finely tuned; they sit on a
plateau.

**The two reference levels are optimised by opposite settings.** Against the raw 36 clusters,
the best resolution is 2.6 (40 clusters). Against the
curated 17 metaclusters, it is **0.6** — giving 18 clusters,
ARI 0.526, far better than resolution 2 achieves (0.362).
The same holds for k: ARI vs the 17 metaclusters is maximised at k = 100
(15 clusters, 0.517).

This is the most interesting methodological result here. **The biological granularity the authors
report is roughly a quarter of the algorithmic resolution they ran.** They clustered at resolution
2 to get 36 fine clusters and then merged by hand down to 17; a single run at resolution ~0.6
would have landed on that granularity directly. The manual merge is doing real work that the
clustering parameters alone were never set to do — which is a point in favour of their curation,
not against it.

---

## 7. Platform reproducibility

The same scripts, the same seed and **identical package versions** (R 4.4.3, Seurat 5.5.1,
harmony 2.0.5, igraph 2.3.4, aricode 1.1.0) were run on two machines:

| Platform | Clusters | ARI vs published raw 36 |
|---|---|---|
| Linux (UCLA Hoffman2, `pi_shared`) | 36 | 0.4956 |
| Windows (local, conda) | 34 | 0.4773 |

Only the operating system and BLAS differ. This independently confirms the authors' own warning,
written in `w.step0.00.run.CD8.sh`, that rerunning their wrapper
*"will yield result very similar with original result, but not exact the same because of
differences of the computer enviroment, such as R version, Seurat version etc."*

The canonical results in `results/` are the **Linux** run.

---

## 8. Limitations

1. **Preprocessing is inherited, not reproduced.** Starting from the authors' deposited PCA means
   normalization, minicluster construction, gene selection and covariate regression are taken on
   trust. A from-counts reproduction would test those too, but would mostly measure preprocessing
   and could not use all 35 datasets, since GEO carries only the 10 newly sequenced 10X cancer types.
2. **Resampling unit is the minicluster, not the cell.** The minicluster definition cannot be
   re-derived without the 14 GB expression deposit, so cell-level resampling was not possible.
3. **The PC sweep is capped at 15**, the number deposited.
4. **GMM used a diagonal (`VVI`) covariance.** Full covariance (`VVV`) is intractable at G = 36 in
   15 dimensions; the GMM result is therefore a floor, not that family's best achievable score.
5. **Cluster-count matching for non-graph methods is a choice.** k-means, Ward and GMM were fixed
   at 36 to match the published count; letting each select its own k by BIC or silhouette would
   give different, probably lower, agreement.
6. **Assumption diagnostics were not run** (no silhouette, gap statistic or modularity profiling);
   the comparison is confined to ARI against the authors' labels, as the assignment specifies.

---

## 9. Reproducibility of this repository

Every number above comes from the numbered scripts in `scripts/`, run end-to-end on Linux
(UCLA Hoffman2, `pi_shared`). The split pipeline was checked against the original single-script
run and reproduces it exactly — baseline ARI 0.4956281 in both.

Two stochastic steps were found unseeded during that check and are now seeded:

- `igraph::cluster_leiden` — its refinement is stochastic (`set.seed(1)` in `06_alt_leiden.R`).
- `mclust::Mclust` — initialises from a **random subset** when n exceeds
  `mclust.options("subset")` (default 2,000; here n = 11,972), so it is seed-dependent
  (`set.seed(1)` in `06_alt_gmm_mclust.R`).

This is the same class of problem the report documents in the source paper: a published number
that moves when you re-run it. Both are now pinned.

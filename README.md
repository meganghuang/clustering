# clustering

Reproduction and stress-test of the CD8 T cell clustering in
**Zheng et al. 2021, *Science* 374:eabe6474** — *Pan-cancer single cell landscape of
tumor-infiltrating T cells*.

The published pipeline is re-run at the authors' exact parameters, then perturbed along four
axes — random seed, data subsampling, three key parameters, and the clustering algorithm itself —
with every comparison made by Adjusted Rand Index against the authors' own labels.

**→ [`report/clustering_methods_report.md`](report/clustering_methods_report.md)** is the write-up.

---

## Headline result

Re-running the published pipeline at the published settings reproduces the **cluster count
exactly** (36) but the **membership only in part**:

| Comparison | minicluster level | cell level |
|---|---|---|
| vs the authors' raw 36-cluster partition | **0.4956** | **0.5126** |
| vs the authors' curated 17 metaclusters | 0.3623 | 0.3768 |

That gap looks alarming until it is put next to the method's own noise: two of *our* runs,
differing only in random seed, agree with each other at a median ARI of **0.478**.
The reproduction gap is no larger than the pipeline's run-to-run variation.

---

## Data

Everything comes from the authors' deposits. Nothing is simulated.

| Source | What | Size |
|---|---|---|
| Zenodo [10.5281/zenodo.5461803](https://doi.org/10.5281/zenodo.5461803) — `data.metaInfo.tar.gz` | minicluster PCA, label tables, informative genes | 60 MB |
| GEO [GSE156728](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE156728) | processed expression (not needed by this pipeline) | — |

`01_fetch_data.sh` downloads the 60 MB tarball. The 14 GB `data.expression.tar.gz` and the
5.2 GB `data.external.tar.gz` are **not** required: the clustering operates on miniclusters, and
the authors deposited their own minicluster PCA (11,972 × 15).

Scope: **234,550 CD8 T cells / 11,972 miniclusters / 35 datasets / 22 cancer types**.

### Two sets of published labels

- `ClusterID.harmony` — the authors' **raw** Louvain partition, 36 clusters (C00–C35). On Zenodo only.
- `meta.cluster` — their **hand-curated** merge of those 36 into 17 named groups. This is the only
  one deposited in GEO.

The merge is many-to-one for 36 of 36 clusters. Scoring is primarily against the raw partition,
because that is the algorithm's actual output.

---

## Repository layout

```
scripts/     numbered pipeline stages (see below)
results/
  tables/    ARI tables, assignments, confusion matrices, co-assignment matrix
  figures/   fig1–fig4 PNGs
report/      clustering_methods_report.md
```

## Scripts

Run in numeric order. **Scripts sharing a number are alternatives at the same pipeline position.**

| Script | Purpose |
|---|---|
| `00_run_all.sh` | driver — runs everything in order |
| `common.R` | shared inputs, the pipeline function, ARI helpers, published parameters |
| `01_fetch_data.sh` | download + extract the Zenodo metadata tarball |
| `02_build_label_tables.R` | minicluster and cell label tables from `int.CD8.S35.meta.tb.rds` |
| `03_baseline_louvain.R` | **baseline** — Harmony + SNN + Louvain at the published settings |
| `04_stability_seeds.R` | 10 random seeds, data fixed |
| `04_stability_subsample.R` | 10 subsamples at 80%, without replacement, + co-assignment |
| `05_sweep_resolution.R` | resolution 0.2 – 5.0 |
| `05_sweep_knn_k.R` | k = 5 – 100 in the KNN graph |
| `05_sweep_npcs.R` | 5 – 15 principal components |
| `06_alt_leiden.R` | Leiden (graph / modularity) |
| `06_alt_kmeans.R` | k-means (centroid) |
| `06_alt_hclust_ward.R` | Ward hierarchical (connectivity) |
| `06_alt_gmm_mclust.R` | Gaussian mixture (model-based) |
| `07_collect_results.R` | merge stage outputs into `ari_results.csv` + between-method ARI |
| `08_figures.py` | the four figures |

`03` must run before `05_*` and `06_*`, which reuse its Harmony embedding and SNN graph.
This is deliberate: resolution is nearly free on a fixed graph, changing k rebuilds only the
neighbour graph, and only the PC count forces a full Harmony redo.

### Running

```bash
export DATA_DIR=data OUT_DIR=results/tables FIG_DIR=results/figures COMMON_R=scripts/common.R
bash scripts/00_run_all.sh
```

Requires R ≥ 4.4 with `Seurat`, `harmony`, `igraph`, `aricode`, `mclust`, `data.table`, `uwot`,
and Python with `pandas`/`matplotlib`. Total runtime is a few minutes on 8 cores — the clustering
units are 11,972 miniclusters, not 234,550 cells.

---

## Published parameters

Taken from the authors' deposited code and their `scPip` package, not from the paper's prose:

```
PCA            npcs = 15                       (inherited from CD8.pca.txt.gz)
Harmony        RunHarmony(vars_use = "dataset")
Graph          FindNeighbors(dims = 1:15, k.param = 20)
Clustering     FindClusters(resolution = 2, algorithm = 1)      # Louvain
Genes          1,500; occupancy >= 0.85; blacklist exclude.gene.misc.human.v4 + ^RP[LS] + MALAT1
Covariates     S.Score, G2M.Score, DIG.Score1
```

The authors selected `RNA_snn_res.2` for CD8 (`res.2.6` for CD4).

---

## Three findings worth the read

1. **Leiden beats Louvain at reproducing Louvain.** Leiden recovers the authors' partition at
   ARI 0.527, better than our own Louvain run (0.496).
2. **The published granularity is ~4× coarser than the resolution they ran.** ARI against the 17
   curated metaclusters peaks at resolution **0.6**
   (18 clusters, 0.526), not at the
   resolution 2 the authors used. Their manual 36→17 merge is doing work the parameters never did.
3. **Same code and versions, different OS, different answer.** Linux gave 36 clusters / ARI
   0.4956; Windows gave 34 / 0.4773 with identical package versions and seed — confirming the
   authors' own note that a rerun is *"very similar ... but not exact the same"*.

Canonical results in `results/` are from the Linux run (UCLA Hoffman2).

#!/usr/bin/env python3
"""08 - figures for the CD8 clustering reproduction.

Reads results/tables/*.csv written by stages 03-07 and writes
results/figures/fig1..fig4 PNGs. Pure matplotlib; no Seurat dependency.
"""
import os
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt

TAB = os.environ.get("OUT_DIR", "results/tables")
FIG = os.environ.get("FIG_DIR", "results/figures")
os.makedirs(FIG, exist_ok=True)
mpl.rcParams.update({"font.size": 8, "axes.titlesize": 8, "axes.labelsize": 8,
                     "xtick.labelsize": 6, "ytick.labelsize": 6,
                     "axes.spines.top": False, "axes.spines.right": False,
                     "figure.dpi": 150, "savefig.dpi": 300})

ari  = pd.read_csv(f"{TAB}/ari_results.csv")
base = ari[ari.arm == "baseline"].iloc[0]
sp   = pd.read_csv(f"{TAB}/seed_pairwise_ari.csv")
ss   = pd.read_csv(f"{TAB}/subsample_pairwise_ari.csv")
algp = pd.read_csv(f"{TAB}/algorithm_pairwise_ari.csv")

# ---------------- Figure 1: baseline vs published, on the authors' embedding
um = pd.read_csv(f"{TAB}/umap_harmony.csv.gz")
ba = pd.read_csv(f"{TAB}/baseline_assignments.csv")
d  = um.merge(ba, on="miniCluster")
dd = d.iloc[np.random.default_rng(0).permutation(len(d))]
sz = np.clip(dd.n_cells.values / 6, 1.0, 9)

fig, axes = plt.subplots(1, 3, figsize=(11.5, 4.0))
for ax, col, ttl in [
        (axes[0], "ClusterID.harmony", "Authors' published partition\n(36 Louvain clusters)"),
        (axes[1], "cluster",           "Our reproduction\n(36 Louvain clusters)"),
        (axes[2], "meta.cluster",      "Authors' curated metaclusters\n(17, hand-merged from the 36)")]:
    codes = pd.Categorical(dd[col]).codes
    cmap = mpl.colormaps["tab20"] if dd[col].nunique() <= 20 else mpl.colormaps["gist_ncar"]
    ax.scatter(dd.UMAP_1, dd.UMAP_2, c=codes, cmap=cmap, s=sz,
               linewidths=0, alpha=0.85, rasterized=True)
    ax.set_title(ttl, fontsize=8); ax.set_xticks([]); ax.set_yticks([])
    for s_ in ax.spines.values():
        s_.set_visible(False)
    ax.margins(0.04)
for xy, xytext in [((0.12, 0.04), (0.02, 0.04)), ((0.02, 0.14), (0.02, 0.04))]:
    axes[0].annotate("", xy=xy, xytext=xytext, xycoords="axes fraction",
                     arrowprops=dict(arrowstyle="->", lw=0.8, color="0.35"))
axes[0].text(0.135, 0.035, "UMAP 1", transform=axes[0].transAxes, fontsize=6, color="0.35", va="center")
axes[0].text(0.015, 0.155, "UMAP 2", transform=axes[0].transAxes, fontsize=6, color="0.35")
fig.suptitle("Rerunning the published pipeline recovers the cluster count exactly, "
             "the partition only in part", fontsize=9, y=0.99)
fig.text(0.5, 0.015,
         "11,972 miniclusters (234,550 CD8 T cells, 35 datasets) on the authors' own 15-PC embedding. "
         f"ARI vs the published 36-cluster partition = {base.ari_raw36_mini:.3f} (minicluster) / "
         f"{base.ari_raw36_cell:.3f} (cell).", ha="center", fontsize=6.5, color="0.3")
fig.tight_layout(rect=[0, 0.045, 1, 0.96])
fig.savefig(f"{FIG}/fig1_baseline_vs_published.png"); plt.close(fig)

# ---------------- Figure 2: stability
co   = pd.read_csv(f"{TAB}/coassign_ordered.csv.gz", header=None, na_values=["NA"]).to_numpy(float)
lab2 = pd.read_csv(f"{TAB}/coassign_order_labels.csv")["meta.cluster"].values
names = ["Seed only\n(same data)", "Subsample 80%\n(no replacement)"]
cols  = ["#1f6fb4", "#d1642f"]

fig = plt.figure(figsize=(11.5, 4.2))
gs = fig.add_gridspec(1, 3, width_ratios=[1.0, 1.0, 1.25], wspace=0.34)
ax = fig.add_subplot(gs[0, 0])
for i, (v, c) in enumerate(zip([sp.ari.values, ss.ari.values], cols)):
    ax.scatter(i + np.random.default_rng(i).uniform(-0.1, 0.1, len(v)), v,
               s=14, color=c, alpha=0.75, linewidths=0, zorder=3)
    ax.plot([i - 0.22, i + 0.22], [np.median(v)] * 2, color="0.15", lw=1.8, zorder=4)
ax.axhline(base.ari_raw36_mini, ls="--", lw=0.9, color="0.45", zorder=2)
ax.text(1.52, base.ari_raw36_mini + 0.006, "our run vs published", fontsize=6, color="0.4", va="bottom")
ax.set_xticks([0, 1]); ax.set_xticklabels(names); ax.set_xlim(-0.5, 2.05)
ax.set_ylabel("ARI between two runs"); ax.set_ylim(0.30, 0.75)
ax.set_title("Re-running alone moves the partition\nas much as the reproduction gap", fontsize=8)

ax = fig.add_subplot(gs[0, 1])
for i, (v, c) in enumerate(zip([ari[ari.arm == "seed"].n_clusters.values,
                                ari[ari.arm == "subsample"].n_clusters.values], cols)):
    ax.scatter(i + np.random.default_rng(i + 5).uniform(-0.1, 0.1, len(v)), v,
               s=14, color=c, alpha=0.75, linewidths=0, zorder=3)
    ax.plot([i - 0.22, i + 0.22], [np.median(v)] * 2, color="0.15", lw=1.8, zorder=4)
ax.axhline(36, ls="--", lw=0.9, color="0.45")
ax.text(1.52, 36.18, "published: 36", fontsize=6, color="0.4", va="bottom")
ax.set_xticks([0, 1]); ax.set_xticklabels(names); ax.set_xlim(-0.5, 2.05)
ax.set_ylabel("Number of clusters"); ax.set_ylim(31.3, 41.2)
ax.set_title("Cluster count is not stable either", fontsize=8)

ax = fig.add_subplot(gs[0, 2])
im = ax.imshow(co, cmap="magma", vmin=0, vmax=1, interpolation="nearest")
for b in np.where(lab2[1:] != lab2[:-1])[0] + 1:
    ax.axhline(b, color="white", lw=0.3); ax.axvline(b, color="white", lw=0.3)
ax.set_xticks([]); ax.set_yticks([])
ax.set_title("Co-assignment over 10 subsamples\n(2,000 miniclusters, blocks = published metaclusters)", fontsize=8)
cb = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.02)
cb.set_label("Replicate fraction\nin the same cluster", fontsize=6.5); cb.ax.tick_params(labelsize=6)
for s_ in ax.spines.values():
    s_.set_visible(False)
fig.suptitle("The published partition sits inside a cloud of equally valid runs", fontsize=9.5, y=0.995)
fig.text(0.5, 0.012, "Resampling unit is the minicluster; pairwise ARI uses only miniclusters shared "
         f"between a pair of replicates ({ss.n_shared.min():,}-{ss.n_shared.max():,} of 9,577).",
         ha="center", fontsize=6.3, color="0.35")
fig.savefig(f"{FIG}/fig2_stability.png", bbox_inches="tight"); plt.close(fig)

# ---------------- Figure 3: parameter sweeps
sw_r = ari[ari.arm == "sweep_resolution"].copy(); sw_r["x"] = sw_r.setting.str.replace("res=", "").astype(float)
sw_k = ari[ari.arm == "sweep_k"].copy();          sw_k["x"] = sw_k.setting.str.replace("k=", "").astype(int)
sw_p = ari[ari.arm == "sweep_npc"].copy();        sw_p["x"] = sw_p.setting.str.replace("npc=", "").astype(int)
C36, C17 = "#1f6fb4", "#c43d6a"
fig, axes = plt.subplots(1, 3, figsize=(11.5, 3.8), sharey=True)
for ax, df, xl, pub, logx in [(axes[0], sw_r, "Clustering resolution", 2.0, False),
                              (axes[1], sw_k, "k in the KNN graph", 20, True),
                              (axes[2], sw_p, "Number of PCs", 15, False)]:
    df = df.sort_values("x")
    ax.plot(df.x, df.ari_raw36_mini, "-o", color=C36, ms=3.2, lw=1.4, zorder=3)
    ax.plot(df.x, df.ari_meta17_mini, "-o", color=C17, ms=3.2, lw=1.4, zorder=3)
    ax.axvline(pub, ls="--", lw=0.9, color="0.5", zorder=1)
    if logx:
        ax.set_xscale("log"); ax.set_xticks([5, 10, 20, 30, 50, 100])
        ax.set_xticklabels([5, 10, 20, 30, 50, 100])
    ax.set_xlabel(xl); ax.margins(0.06)
    b36 = df.loc[df.ari_raw36_mini.idxmax()]; b17 = df.loc[df.ari_meta17_mini.idxmax()]
    ax.plot(b36.x, b36.ari_raw36_mini, "o", mfc="none", mec=C36, ms=9, mew=1.4, zorder=4)
    ax.plot(b17.x, b17.ari_meta17_mini, "o", mfc="none", mec=C17, ms=9, mew=1.4, zorder=4)
axes[0].set_ylabel("ARI vs the authors' labels"); axes[0].set_ylim(0.18, 0.60)
axes[0].text(2.1, 0.195, "published\nsetting", fontsize=6, color="0.45")
axes[1].text(21.5, 0.195, "published\nsetting", fontsize=6, color="0.45")
axes[2].text(13.1, 0.195, "published\nsetting", fontsize=6, color="0.45", ha="right")
axes[0].text(2.75, sw_r.ari_raw36_mini.max() + 0.012, "vs the 36 raw clusters", color=C36, fontsize=7)
axes[0].text(0.35, sw_r.ari_meta17_mini.max() + 0.012, "vs the 17 metaclusters", color=C17, fontsize=7)
fig.suptitle("The two reference levels are optimised by opposite settings; circles mark each curve's best",
             fontsize=9, y=1.0)
fig.text(0.5, -0.03, "Open circles mark the maximum of each curve. The authors ran resolution 2, k = 20, 15 PCs (dashed).",
         ha="center", fontsize=6.3, color="0.35")
fig.savefig(f"{FIG}/fig3_parameter_sweeps.png", bbox_inches="tight"); plt.close(fig)

# ---------------- Figure 4: algorithms
alg = ari[ari.arm == "algorithm"].copy()
pretty = {"louvain_baseline": "Louvain\n(the authors' method)", "leiden_igraph": "Leiden\n(graph, modularity)",
          "kmeans": "k-means\n(centroid)", "hclust_ward": "Ward hierarchical\n(connectivity)",
          "gmm_mclust_VVI": "Gaussian mixture\n(model-based)"}
fam = {"louvain_baseline": "#1f6fb4", "leiden_igraph": "#1f6fb4", "kmeans": "#8a8f98",
       "hclust_ward": "#8a8f98", "gmm_mclust_VVI": "#8a8f98"}
alg["label"] = alg.setting.map(pretty)
alg = alg.sort_values("ari_raw36_mini")
fig, axes = plt.subplots(1, 2, figsize=(11.0, 3.9), gridspec_kw={"width_ratios": [1.25, 1.0], "wspace": 0.30})
ax = axes[0]
for yi, (_, r) in enumerate(alg.iterrows()):
    c = fam[r.setting]
    ax.plot([0, r.ari_raw36_mini], [yi, yi], color=c, lw=1.0, alpha=0.5, zorder=2)
    ax.scatter(r.ari_raw36_mini, yi, s=52, color=c, zorder=3)
    ax.text(r.ari_raw36_mini + 0.012, yi, f"{r.ari_raw36_mini:.3f}", va="center", fontsize=7, color="0.2")
ax.set_yticks(range(len(alg))); ax.set_yticklabels(alg.label)
ax.set_xlabel("ARI vs the authors' 36-cluster partition"); ax.set_xlim(0, 0.62)
ax.set_title("Only graph-based methods come close to the published partition", fontsize=8)
ax.text(0.545, 0.22, "higher = better", fontsize=6, color="0.45", transform=ax.transAxes)

order = [o for o in ["louvain_baseline", "leiden_igraph", "kmeans", "hclust_ward", "gmm_mclust_VVI"]
         if o in set(algp.a) | set(algp.b)]
M = np.full((len(order), len(order)), np.nan)
for _, r in algp.iterrows():
    if r.a in order and r.b in order:
        i, j = order.index(r.a), order.index(r.b); M[i, j] = M[j, i] = r.ari
np.fill_diagonal(M, 1.0)
ax = axes[1]
im = ax.imshow(M, cmap="viridis", vmin=0.2, vmax=1.0)
short = {"louvain_baseline": "Louvain", "leiden_igraph": "Leiden", "kmeans": "k-means",
         "hclust_ward": "Ward", "gmm_mclust_VVI": "GMM"}
lbl = [short[o] for o in order]
ax.set_xticks(range(len(order))); ax.set_xticklabels(lbl, rotation=45, ha="right")
ax.set_yticks(range(len(order))); ax.set_yticklabels(lbl)
for i in range(len(order)):
    for j in range(len(order)):
        ax.text(j, i, f"{M[i, j]:.2f}", ha="center", va="center", fontsize=6.5,
                color="white" if M[i, j] < 0.72 else "0.1")
ax.set_title("Louvain and Leiden agree with each other\nfar more than with anything else", fontsize=8)
cb = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.03)
cb.set_label("ARI between methods", fontsize=6.5); cb.ax.tick_params(labelsize=6)
for s_ in ax.spines.values():
    s_.set_visible(False)
fig.suptitle("Method family matters more than any parameter choice", fontsize=9.5, y=1.0)
fig.text(0.5, -0.06, "All methods run on the same Harmony embedding; k-means, Ward and GMM fixed at 36 "
         "clusters to match the published count.", ha="center", fontsize=6.3, color="0.35")
fig.savefig(f"{FIG}/fig4_algorithms.png", bbox_inches="tight"); plt.close(fig)
print("08 done - 4 figures written to", FIG)

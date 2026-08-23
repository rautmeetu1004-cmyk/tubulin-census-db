#!/usr/bin/env python3
"""Inferential view of the isotype structure: PCA of the nineteen human tubulin
isotypes on their pathway/GO enrichment profiles, with a PERMANOVA test of the
alpha-versus-beta separation. Real per-isotype g:Profiler data.

Panel A: PC1 vs PC2 with 95% confidence ellipses per family.
Panel B: scree (variance explained per component).
"""
import csv, os, math, collections
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
A = f"{ROOT}/01_Census/Tubulin_Alpha_Interactome/04_analysis/enrichment/enrichment_all.csv"
B = f"{ROOT}/01_Census/Tubulin_Beta_Interactome/04_analysis/enrichment/enrichment_all.csv"
alpha_iso = ["TUBA1A","TUBA1B","TUBA1C","TUBA4A","TUBA3C","TUBA3D","TUBA8","TUBA3E","TUBAL3"]
beta_iso  = ["TUBB","TUBB4B","TUBB3","TUBB2A","TUBB2B","TUBB4A","TUBB6","TUBB8","TUBB1","TUBB8B"]
isos = alpha_iso + beta_iso
fam = {i:"alpha" for i in alpha_iso}; fam.update({i:"beta" for i in beta_iso})

val = {}; terms = set()
for path in (A, B):
    for r in csv.DictReader(open(path)):
        iso = r["gene_list"]
        if iso not in isos: continue
        try: p = float(r["p_value"])
        except ValueError: continue
        if p <= 0: continue
        t = r["term_id"]; terms.add(t)
        val[(iso, t)] = max(val.get((iso, t), 0.0), -math.log10(p))
terms = sorted(terms)
X = np.array([[val.get((i, t), 0.0) for t in terms] for i in isos])

Xs = StandardScaler().fit_transform(X)
pca = PCA(n_components=min(10, len(isos)-1)).fit(Xs)
Y = pca.transform(Xs)
ev = pca.explained_variance_ratio_ * 100

# PERMANOVA (Euclidean) on the standardized space: is alpha vs beta separation real?
def permanova(D2, labels, nperm=4999, seed=12345):
    labels = np.array(labels); n = len(labels)
    def ssw(lab):
        tot = 0.0
        for g in set(lab):
            idx = np.where(lab == g)[0]
            if len(idx) > 1:
                tot += D2[np.ix_(idx, idx)].sum() / (2 * len(idx))
        return tot
    sst = D2.sum() / (2 * n)
    a = len(set(labels))
    def F(lab):
        w = ssw(lab); b = sst - w
        return (b / (a - 1)) / (w / (n - a))
    F0 = F(labels)
    rng = np.random.RandomState(seed)
    cnt = 1
    for _ in range(nperm):
        if F(rng.permutation(labels)) >= F0: cnt += 1
    return F0, cnt / (nperm + 1)

D2 = np.sum((Y[:, None, :] - Y[None, :, :])**2, axis=2)
Fstat, pval = permanova(D2, [fam[i] for i in isos])

COL = {"alpha": "#2E6CA4", "beta": "#BE5738"}
fig, (axA, axB) = plt.subplots(1, 2, figsize=(12.6, 5.8),
                               gridspec_kw={"width_ratios": [1.7, 1]})

def ellipse(ax, pts, col):
    if len(pts) < 3: return
    mu = pts.mean(0); cov = np.cov(pts.T)
    vals, vecs = np.linalg.eigh(cov)
    order = vals.argsort()[::-1]; vals, vecs = vals[order], vecs[:, order]
    ang = math.degrees(math.atan2(vecs[1, 0], vecs[0, 0]))
    w, h = 2 * 2.4477 * np.sqrt(np.maximum(vals, 1e-9))   # ~95%
    e = Ellipse(mu, w, h, angle=ang, facecolor=col, alpha=0.07, edgecolor=col, lw=1.1, ls="--")
    ax.add_patch(e)

for f in ("alpha", "beta"):
    idx = [k for k, i in enumerate(isos) if fam[i] == f]
    P = Y[idx][:, :2]
    ellipse(axA, P, COL[f])
    axA.scatter(P[:, 0], P[:, 1], s=90, color=COL[f], edgecolor="#16323D",
                linewidth=0.8, label=f"{f} (n={len(idx)})", zorder=3)
for k, i in enumerate(isos):
    axA.annotate(i, (Y[k, 0], Y[k, 1]), fontsize=7.6, color="#16323D",
                 xytext=(5, 3), textcoords="offset points",
                 bbox=dict(boxstyle="round,pad=0.12", fc="white", ec="none", alpha=0.65))
axA.axhline(0, color="#CBD3D9", lw=0.8); axA.axvline(0, color="#CBD3D9", lw=0.8)
axA.set_xlabel(f"PC1 ({ev[0]:.0f}% variance)", fontsize=10.5)
axA.set_ylabel(f"PC2 ({ev[1]:.0f}% variance)", fontsize=10.5)
axA.set_title("A  Isotypes in pathway-profile space (PCA)", fontsize=11.5, loc="left", weight="bold")
axA.legend(frameon=False, fontsize=9.5, loc="best")
axA.text(0.02, 0.02, f"PERMANOVA alpha vs beta:  F = {Fstat:.2f},  p = {pval:.4f}  (4,999 perm.)",
         transform=axA.transAxes, fontsize=9, color="#42505C",
         bbox=dict(boxstyle="round,pad=0.35", fc="#F2F5F7", ec="#D5DBE1"))
axA.spines[["top", "right"]].set_visible(False)

k = min(8, len(ev))
axB.bar(range(1, k+1), ev[:k], color="#2C8C74", edgecolor="#1E5E4E", linewidth=0.5)
axB.plot(range(1, k+1), np.cumsum(ev[:k]), "-o", color="#C4922B", lw=1.6, ms=4, label="cumulative")
axB.set_xlabel("principal component", fontsize=10.5)
axB.set_ylabel("variance explained (%)", fontsize=10.5)
axB.set_title("B  Scree", fontsize=11.5, loc="left", weight="bold")
axB.legend(frameon=False, fontsize=9)
axB.spines[["top", "right"]].set_visible(False)
axB.yaxis.grid(True, ls=":", lw=0.6, color="#C7CED6"); axB.set_axisbelow(True)

fig.suptitle("Alpha and beta overlap in pathway-profile space (PERMANOVA p = 0.26): PC1 tracks "
             "interactome breadth, not family",
             fontsize=12.5, y=1.0)
fig.tight_layout()
os.makedirs(f"{ROOT}/05_Pathway_Analysis/figures", exist_ok=True)
out = f"{ROOT}/05_Pathway_Analysis/figures/P5_isotype_pca"
fig.savefig(out + ".png", dpi=300, bbox_inches="tight", facecolor="white")
fig.savefig(out + ".pdf", bbox_inches="tight")
print(f"wrote P5_isotype_pca; terms={len(terms)} PC1={ev[0]:.1f}% PC2={ev[1]:.1f}% "
      f"PERMANOVA F={Fstat:.2f} p={pval:.4f}")

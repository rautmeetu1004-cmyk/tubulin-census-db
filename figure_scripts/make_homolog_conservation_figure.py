#!/usr/bin/env python3
"""Cross-species conservation of the tubulin interactome (the homolog study).

Left  : per species, the fraction of the 73 curated interactome partners that
        retains a detectable ortholog (OrthoDB presence/absence), ordered by
        distance from human. The set thins with distance to an ancient core.
Right : mean species-breadth (out of 9 non-human species) by functional module,
        showing the folding machinery is universal while the division and
        ciliary partners are progressively younger.
Real data only; no synthetic background.
"""
import csv, os, collections
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ARM = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
rows = list(csv.DictReader(open(f"{ARM}/data/taxonomy_conservation.csv")))
partners = [r for r in rows if r["group"] != "core_tubulin_isoform"]
N = len(partners)

species = ["M.musculus", "G.gallus", "X.tropicalis", "D.rerio", "D.melanogaster",
           "C.elegans", "Dictyostelium", "A.thaliana", "S.cerevisiae"]
slabels = ["M. musculus", "G. gallus", "X. tropicalis", "D. rerio", "D. melanogaster",
           "C. elegans", "D. discoideum", "A. thaliana", "S. cerevisiae"]
def frac(sp): return sum(int(r[sp]) for r in partners) / N * 100
sv = [frac(sp) for sp in species]

# module breadth
mod = collections.defaultdict(list)
for r in partners:
    mod[r["family_module"]].append(sum(int(r[s]) for s in species))
mod_order = ["folding;prefoldin", "CCT_chaperonin;folding", "folding",
             "spindle", "centrosome", "cilium"]
mlabels = {"folding;prefoldin": "Prefoldin", "CCT_chaperonin;folding": "CCT / TriC",
           "folding": "Other folding", "spindle": "Spindle",
           "centrosome": "Centrosome", "cilium": "Cilium"}
mod_order = [m for m in mod_order if m in mod]
mv = [np.mean(mod[m]) for m in mod_order]
mn = [len(mod[m]) for m in mod_order]

fig, (axL, axR) = plt.subplots(1, 2, figsize=(12.6, 5.4),
                               gridspec_kw={"width_ratios": [1.25, 1]})

# ---- left: per-species partner conservation ----
y = np.arange(len(species))[::-1]
cols = plt.cm.YlOrRd(np.linspace(0.75, 0.35, len(species)))
bars = axL.barh(y, sv, 0.66, color=cols, edgecolor="#7A3B12", linewidth=0.5)
for bar, v in zip(bars, sv):
    axL.text(v + 1.3, bar.get_y() + bar.get_height()/2, f"{v:.0f}%",
             va="center", ha="left", fontsize=9, color="#333")
axL.set_yticks(y); axL.set_yticklabels(slabels, fontstyle="italic", fontsize=10.5)
axL.set_xlim(0, 100); axL.set_xlabel("Interactome partners with an ortholog (%)", fontsize=10.5)
axL.set_title(f"A  Conservation across species (n = {N} partners)", fontsize=11, loc="left", weight="bold")
_univ = round(sum(1 for r in partners if all(r[s] == "1" for s in species)) / N * 100)
axL.axvline(_univ, ls="--", lw=1, color="#555")
axL.text(_univ, len(species)-0.35, f"  {_univ}% universal\n  (all 9 species)", fontsize=8.5, color="#555", va="top")
axL.spines[["top", "right"]].set_visible(False)
axL.xaxis.grid(True, ls=":", lw=0.6, color="#C7CED6"); axL.set_axisbelow(True)

# ---- right: module breadth as violins (distribution, not just the mean) ----
ym = np.arange(len(mod_order))[::-1]
mcol = ["#1E7A34", "#2FA34B", "#6DBE45", "#2E6DB4", "#4C8FC0", "#C0392B"][:len(mod_order)]
data = [mod[m] for m in mod_order]
vp = axR.violinplot(data, positions=ym, vert=False, widths=0.85,
                    showmeans=False, showextrema=False)
for body, col in zip(vp["bodies"], mcol):
    body.set_facecolor(col); body.set_edgecolor("#243B4A"); body.set_alpha(0.55); body.set_linewidth(0.6)
rng = np.random.RandomState(0)
for yi, vals, col in zip(ym, data, mcol):
    jit = yi + (rng.rand(len(vals)) - 0.5) * 0.28
    axR.scatter(vals, jit, s=14, color=col, edgecolor="#243B4A", linewidth=0.3, alpha=0.9, zorder=3)
    axR.scatter([np.mean(vals)], [yi], marker="|", s=420, color="#16323D", zorder=4, linewidth=2)
    axR.text(9.5, yi, f"mean {np.mean(vals):.1f}", va="center", ha="right", fontsize=8.5, color="#555")
axR.set_yticks(ym)
axR.set_yticklabels([f"{mlabels[m]} (n={n})" for m, n in zip(mod_order, mn)], fontsize=10)
axR.set_xlim(0, 9.7); axR.set_ylim(-0.6, len(mod_order)-0.4)
axR.set_xlabel("Species-breadth per protein (of 9 species)", fontsize=10.5)
axR.set_title("B  Conservation depth by functional module", fontsize=11, loc="left", weight="bold")
axR.spines[["top", "right"]].set_visible(False)
axR.xaxis.grid(True, ls=":", lw=0.6, color="#C7CED6"); axR.set_axisbelow(True)

fig.suptitle("The tubulin interactome contracts with evolutionary distance to an ancient folding core",
             fontsize=12.5, y=1.005)
fig.tight_layout()
os.makedirs(f"{ARM}/figures", exist_ok=True)
fig.savefig(f"{ARM}/figures/H1_homolog_conservation.png", dpi=300, bbox_inches="tight", facecolor="white")
fig.savefig(f"{ARM}/figures/H1_homolog_conservation.pdf", bbox_inches="tight")
print("wrote H1_homolog_conservation")
print("species%:", [f"{v:.0f}" for v in sv])
print("modules:", list(zip([mlabels[m] for m in mod_order], [round(v,1) for v in mv])))

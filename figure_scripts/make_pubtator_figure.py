#!/usr/bin/env python3
"""Literature text-mining arm (PubTator3), better represented.

Specificity of a partner = tubulin co-mentions / total mentions (entrez-keyed),
exactly the ratio the census ranks by. Panel A: the most tubulin-specific
partners (>= 40 total mentions, tubulin seeds excluded) - the microtubule
machinery rises to the top. Panel B: every partner placed by how well studied
it is (total mentions) against how tubulin-specific it is; the specific core
separates from the well-studied but promiscuous bulk.
Real PubTator counts from partner_external_evidence.csv (both families).
"""
import csv, os, math, re
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
FILES = {
    "alpha": f"{ROOT}/01_Census/Tubulin_Alpha_Interactome/04_analysis/external_evidence/partner_external_evidence.csv",
    "beta":  f"{ROOT}/01_Census/Tubulin_Beta_Interactome/04_analysis/external_evidence/partner_external_evidence.csv",
}
SEED = re.compile(r"^TUB[AB]")   # exclude the alpha/beta seeds themselves

rec = {}
for fam, path in FILES.items():
    for r in csv.DictReader(open(path)):
        g = r["gene"]
        try:
            co = float(r["lit_comention_tubulin_entrez"]); tot = float(r["lit_total_mentions_entrez"])
        except (ValueError, KeyError):
            continue
        if tot <= 0:
            continue
        if g not in rec or tot > rec[g][2]:
            rec[g] = (g, fam, tot, co, co/tot)
rows = list(rec.values())

COL = {"alpha": "#2E6DB4", "beta": "#C0392B"}

# ---- Panel A: top specificity (machinery), seeds excluded, >=40 mentions ----
cand = [r for r in rows if r[2] >= 40 and not SEED.match(r[0])]
topA = sorted(cand, key=lambda r: -r[4])[:16]

fig, (axA, axB) = plt.subplots(1, 2, figsize=(13.2, 6.4),
                               gridspec_kw={"width_ratios": [1.0, 1.15]})
y = np.arange(len(topA))[::-1]
axA.barh(y, [r[4]*100 for r in topA], 0.7,
         color=[COL[r[1]] for r in topA], edgecolor="#222", linewidth=0.4)
for yi, r in zip(y, topA):
    axA.text(r[4]*100 + 1, yi, f"{r[4]*100:.0f}%  ({int(r[3])}/{int(r[2])})",
             va="center", ha="left", fontsize=8.3, color="#333")
axA.set_yticks(y); axA.set_yticklabels([r[0] for r in topA], fontsize=9.5)
axA.set_xlim(0, 100); axA.set_xlabel("Tubulin-specificity of the literature (co-mentions / total, %)", fontsize=10)
axA.set_title("A  Most tubulin-specific partners (seeds excluded)", fontsize=11, loc="left", weight="bold")
axA.spines[["top", "right"]].set_visible(False)
axA.xaxis.grid(True, ls=":", lw=0.6, color="#C7CED6"); axA.set_axisbelow(True)
from matplotlib.patches import Patch
axA.legend(handles=[Patch(color=COL["alpha"], label="alpha partner"),
                    Patch(color=COL["beta"], label="beta partner")],
           loc="lower right", frameon=False, fontsize=8.5)

# ---- Panel B: mentions vs specificity ----
xs = np.array([r[2] for r in rows]); ys = np.array([r[4]*100 for r in rows])
cs = [COL[r[1]] for r in rows]
axB.scatter(xs, ys, s=14, c=cs, alpha=0.45, edgecolor="none")
axB.set_xscale("log")
label_set = {"TUBGCP4", "TUBGCP3", "TUBGCP2", "TUBG1", "NEDD1", "CEP97", "IFT56",
             "SVBP", "CCT2", "CCT5", "SQSTM1", "MAPT", "STMN1", "CLASP1"}
for r in rows:
    if r[0] in label_set:
        axB.annotate(r[0], (r[2], r[4]*100), fontsize=7.6, color="#16323D",
                     xytext=(3, 2), textcoords="offset points")
axB.axhline(50, ls="--", lw=1, color="#888")
axB.text(axB.get_xlim()[1]*0.55, 52, "tubulin-specific", fontsize=8.5, color="#666")
axB.set_xlabel("How well studied (total literature mentions, log)", fontsize=10)
axB.set_ylabel("Tubulin-specificity (%)", fontsize=10)
axB.set_title("B  Specific core vs the well-studied but promiscuous bulk", fontsize=11, loc="left", weight="bold")
axB.spines[["top", "right"]].set_visible(False)
axB.grid(True, ls=":", lw=0.5, color="#D5DBE1"); axB.set_axisbelow(True)

fig.suptitle(f"Literature corroboration (PubTator3): partner specificity across {len(rows):,} tubulin partners",
             fontsize=12.5, y=1.01)
fig.tight_layout()
OUT = f"{os.path.dirname(os.path.dirname(os.path.abspath(__file__)))}/figures"
os.makedirs(OUT, exist_ok=True)
fig.savefig(f"{OUT}/L1_pubtator_specificity.png", dpi=300, bbox_inches="tight", facecolor="white")
fig.savefig(f"{OUT}/L1_pubtator_specificity.pdf", bbox_inches="tight")
print("wrote L1_pubtator_specificity; top:",
      ", ".join(f"{r[0]} {r[4]*100:.0f}%" for r in topA[:6]))

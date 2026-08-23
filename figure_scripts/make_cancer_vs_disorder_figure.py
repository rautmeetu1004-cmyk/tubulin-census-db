#!/usr/bin/env python3
"""Inference from the census: cancer versus constitutional (non-cancer) disorders.

The disease-linked tubulin partners are largely pleiotropic, implicated in both
neoplastic and constitutional disease, but the two disease classes reach them by
different mutational routes: cancer associations carry somatic-mutation evidence,
while non-cancer disorders are dominated by germline (genetic) evidence.

Data: Open Targets association export (disease_associations.csv). Real values;
percentages are association-level, evidence flagged as somatic_score>0 or
genetic_score>0.
"""
import csv, os, collections
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ARM = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
rows = list(csv.DictReader(open(f"{ARM}/analysis/disease_associations.csv")))
CANCER = "cancer or benign tumor"
def f(x):
    try: return float(x)
    except: return 0.0

genes = set(r["gene"] for r in rows)
gc = collections.defaultdict(bool); go = collections.defaultdict(bool)
for r in rows:
    if CANCER in r["therapeutic_areas"]: gc[r["gene"]] = True
    else: go[r["gene"]] = True
both = sum(1 for g in genes if gc[g] and go[g])
conly = sum(1 for g in genes if gc[g] and not go[g])
oonly = sum(1 for g in genes if go[g] and not gc[g])

cancer = [r for r in rows if CANCER in r["therapeutic_areas"]]
other = [r for r in rows if CANCER not in r["therapeutic_areas"]]
def pos(l, k): return sum(1 for r in l if f(r[k]) > 0) / len(l) * 100
c_som, c_gen = pos(cancer, "somatic_score"), pos(cancer, "genetic_score")
o_som, o_gen = pos(other, "somatic_score"), pos(other, "genetic_score")

fig, (axA, axB) = plt.subplots(1, 2, figsize=(12.4, 5.0),
                               gridspec_kw={"width_ratios": [1.15, 1]})

# ---- Panel A: overlap of disease-linked partners ----
segs = [("Cancer only", conly, "#C0392B"),
        ("Both cancer &\nconstitutional", both, "#8E44AD"),
        ("Constitutional\ndisorder only", oonly, "#2E6DB4")]
left = 0
for lab, val, col in segs:
    axA.barh(0, val, left=left, height=0.55, color=col, edgecolor="white")
    axA.text(left + val/2, 0, f"{val:,}", va="center", ha="center",
             color="white", fontsize=13, weight="bold")
    left += val
axA.set_xlim(0, len(genes)); axA.set_ylim(-0.7, 0.9)
axA.set_yticks([])
axA.set_xlabel(f"Disease-linked tubulin partners (n = {len(genes):,})", fontsize=10.5)
axA.set_title("A  The same partners span cancer and constitutional disease",
              fontsize=11, loc="left", weight="bold")
# legend
handles = [plt.Rectangle((0, 0), 1, 1, color=c) for _, _, c in segs]
axA.legend(handles, [s[0].replace("\n", " ") for s in segs],
           loc="upper center", bbox_to_anchor=(0.5, -0.16), ncol=3,
           frameon=False, fontsize=9)
axA.spines[["top", "right", "left"]].set_visible(False)
pct = both / len(genes) * 100
axA.text(len(genes)/2, 0.72, f"{pct:.0f}% of disease-linked partners are implicated in both",
         ha="center", fontsize=9.5, color="#444", style="italic")

# ---- Panel B: mutational route contrast ----
x = np.arange(2); w = 0.36
som = [c_som, o_som]; gen = [c_gen, o_gen]
b1 = axB.bar(x - w/2, som, w, color="#C0392B", label="Somatic-mutation evidence")
b2 = axB.bar(x + w/2, gen, w, color="#2E6DB4", label="Germline (genetic) evidence")
for bars, vals in ((b1, som), (b2, gen)):
    for bar, v in zip(bars, vals):
        axB.text(bar.get_x() + bar.get_width()/2, v + 0.9, f"{v:.1f}%",
                 ha="center", fontsize=9.5, color="#333")
axB.set_xticks(x)
axB.set_xticklabels([f"Cancer\nassociations\n(n={len(cancer):,})",
                     f"Constitutional\ndisorders\n(n={len(other):,})"], fontsize=9.5)
axB.set_ylabel("Associations with evidence of type (%)", fontsize=10.5)
axB.set_ylim(0, 60)
axB.set_title("B  Cancer is somatic-driven, constitutional disease germline-driven",
              fontsize=11, loc="left", weight="bold")
axB.legend(frameon=False, fontsize=9, loc="upper center", bbox_to_anchor=(0.5, 1.0))
axB.spines[["top", "right"]].set_visible(False)
axB.yaxis.grid(True, ls=":", lw=0.6, color="#C7CED6"); axB.set_axisbelow(True)

fig.suptitle("Inference from the census: tubulin partners bridge cancer and inherited disease "
             "through distinct mutational routes", fontsize=12.5, y=1.02)
fig.tight_layout()
os.makedirs(f"{ARM}/figures", exist_ok=True)
fig.savefig(f"{ARM}/figures/D2_cancer_vs_disorder.png", dpi=300, bbox_inches="tight", facecolor="white")
fig.savefig(f"{ARM}/figures/D2_cancer_vs_disorder.pdf", bbox_inches="tight")
print(f"wrote D2_cancer_vs_disorder; both={both} conly={conly} oonly={oonly} "
      f"cancer som/gen={c_som:.1f}/{c_gen:.1f} other som/gen={o_som:.1f}/{o_gen:.1f}")

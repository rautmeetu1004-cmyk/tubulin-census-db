#!/usr/bin/env python3
"""P1: drug-pocket conservation vs isotype-variability bubble scatter.
Recomputed from analysis/pocket_conservation.csv (current, nine-alpha isotypes)."""
import csv, os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ARM = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
rows = list(csv.DictReader(open(f"{ARM}/analysis/pocket_conservation.csv")))
os.makedirs(f"{ARM}/figures", exist_ok=True)
COL = {"Pironetin": "#E9A23B", "Vinca": "#4C72B0", "Taxane": "#DD8452",
       "Colchicine": "#5E7A86", "Laulimalide": "#2A9D8F"}

fig, ax = plt.subplots(figsize=(7.6, 5.2))
ax.axvline(0.15, ls="--", lw=1, color="#B0BEC5")
for r in rows:
    x = float(r["mean_isotype_entropy"]); y = float(r["ET_percentile"]); n = int(r["n_contact"])
    ax.scatter(x, y, s=40 + 22 * n, color=COL.get(r["pocket"], "#888"),
               alpha=0.75, edgecolor="#16323D", linewidth=0.9, zorder=3)
    ax.annotate(f'{r["pocket"]}\n({r["location"]})', (x, y), fontsize=8.5, color="#16323D",
                xytext=(6, 6), textcoords="offset points", weight="bold")
ax.set_xlabel("human-isotype variability (mean Shannon entropy across 9 alpha / 10 beta)  →  resistance-prone", fontsize=9)
ax.set_ylabel("evolutionary conservation (ET percentile)  ↑", fontsize=9)
ax.text(0.06, ax.get_ylim()[0] + 0.4, "more durable", color="#2A9D8F", fontsize=8.5)
ax.text(0.27, ax.get_ylim()[0] + 0.4, "more vulnerable", color="#C44E52", fontsize=8.5)
ax.set_title("Tubulin drug pockets: conservation vs isotype-variability\n"
             "(bubble = contact residues; left = durable pan-isotype target, right = isotype-vulnerable)",
             fontsize=10)
ax.spines[["top", "right"]].set_visible(False)
fig.savefig(f"{ARM}/figures/P1_pocket_conservation_variability.png", dpi=200, bbox_inches="tight", facecolor="white")
fig.savefig(f"{ARM}/figures/P1_pocket_conservation_variability.pdf", bbox_inches="tight")
print("wrote P1_pocket_conservation_variability;",
      ", ".join(f'{r["pocket"]} ET{r["ET_percentile"]}/ent{r["mean_isotype_entropy"]}' for r in rows))

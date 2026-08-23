#!/usr/bin/env python3
"""Pathway / GO enrichment across all nineteen human tubulin isotypes
(nine alpha + ten beta). Value = -log10(p) from the per-isotype g:Profiler
over-representation runs; blank where a term is not enriched in that isotype.
Real data pooled from both family enrichment_all.csv tables.
"""
import csv, os, math, collections
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
A = f"{ROOT}/01_Census/Tubulin_Alpha_Interactome/04_analysis/enrichment/enrichment_all.csv"
B = f"{ROOT}/01_Census/Tubulin_Beta_Interactome/04_analysis/enrichment/enrichment_all.csv"

alpha_iso = ["TUBA1A", "TUBA1B", "TUBA1C", "TUBA4A", "TUBA3C", "TUBA3D", "TUBA8", "TUBA3E", "TUBAL3"]
beta_iso  = ["TUBB", "TUBB4B", "TUBB3", "TUBB2A", "TUBB2B", "TUBB4A", "TUBB6", "TUBB8", "TUBB1", "TUBB8B"]
isos = alpha_iso + beta_iso

# (isotype, term) -> -log10 p ; term -> display name
val = {}
termname = {}
seen_iso_term = collections.defaultdict(set)
for path in (A, B):
    for r in csv.DictReader(open(path)):
        iso = r["gene_list"]
        if iso not in isos:
            continue
        try:
            p = float(r["p_value"])
        except ValueError:
            continue
        if p <= 0:
            continue
        term = r["term_id"]
        termname[term] = r["term_name"]
        s = -math.log10(p)
        val[(iso, term)] = max(val.get((iso, term), 0), s)
        seen_iso_term[term].add(iso)

# choose rows: terms enriched in the most isotypes, then by peak score; keep a
# readable set and drop near-duplicate names.
def peak(term): return max(val[(i, term)] for i in isos if (i, term) in val)
cands = sorted(seen_iso_term, key=lambda t: (len(seen_iso_term[t]), peak(t)), reverse=True)
rows, used_names = [], set()
for t in cands:
    nm = termname[t].strip()
    key = nm.lower()[:34]
    if key in used_names:
        continue
    used_names.add(key)
    rows.append(t)
    if len(rows) == 18:
        break
# order rows by peak score descending for a clean gradient
rows.sort(key=peak, reverse=True)

M = np.full((len(rows), len(isos)), np.nan)
for ri, t in enumerate(rows):
    for ci, iso in enumerate(isos):
        if (iso, t) in val:
            M[ri, ci] = val[(iso, t)]

fig, ax = plt.subplots(figsize=(12.2, 7.6))
cmap = plt.cm.magma.copy()
cmap.set_bad("#EDEFF2")
vmax = np.nanpercentile(M, 98)
im = ax.imshow(M, aspect="auto", cmap=cmap, vmin=0, vmax=vmax)

ax.set_xticks(range(len(isos)))
ax.set_xticklabels(isos, rotation=55, ha="right", fontsize=9)
for lab, iso in zip(ax.get_xticklabels(), isos):
    lab.set_color("#B0281A" if iso in beta_iso else "#1F4E96")
def short(s):
    s = s.replace("Formation of tubulin folding intermediates by CCT/TriC", "Tubulin folding intermediates (CCT/TriC)")
    return s if len(s) <= 46 else s[:44] + "…"
ax.set_yticks(range(len(rows)))
ax.set_yticklabels([short(termname[t]) for t in rows], fontsize=9)

# divider between alpha and beta blocks
ax.axvline(len(alpha_iso) - 0.5, color="white", lw=3)
ax.axvline(len(alpha_iso) - 0.5, color="#444", lw=1, ls="--")
ax.text((len(alpha_iso)-1)/2, -1.15, "alpha (9)", ha="center", color="#1F4E96", fontsize=11, weight="bold")
ax.text(len(alpha_iso) + (len(beta_iso)-1)/2, -1.15, "beta (10)", ha="center", color="#B0281A", fontsize=11, weight="bold")

cb = fig.colorbar(im, ax=ax, fraction=0.025, pad=0.015)
cb.set_label(r"$-\log_{10}\,p$  (grey = not enriched)", fontsize=10)
ax.set_title("Pathway and GO enrichment across all nineteen human tubulin isotypes",
             fontsize=13, pad=26)
fig.tight_layout()
os.makedirs(f"{os.path.dirname(os.path.dirname(os.path.abspath(__file__)))}/figures", exist_ok=True)
out = f"{os.path.dirname(os.path.dirname(os.path.abspath(__file__)))}/figures/P4_isotype_enrichment_heatmap"
fig.savefig(out + ".png", dpi=300, bbox_inches="tight", facecolor="white")
fig.savefig(out + ".pdf", bbox_inches="tight")
print("wrote P4_isotype_enrichment_heatmap; rows:", len(rows), "isotypes:", len(isos))
print("terms:", [short(termname[t]) for t in rows])

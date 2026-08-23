#!/usr/bin/env python3
"""Sequence logos for the headline MEME motifs, selected data-drivenly from the
current per-family runs (motifs_by_set.csv assigns each motif a family, so the
picks follow the data and never go stale on a re-run). Renders bits logos with
logomaker. Real MEME output only.
"""
import os, re, csv
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import logomaker

ARM = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ROOT = os.path.dirname(ARM)
AA = list("ACDEFGHIKLMNPQRSTVWY")
MEME_PATH = {
    "α partners": f"{ARM}/motifs_by_family/alpha/meme.txt",
    "β partners": f"{ARM}/motifs_by_family/beta/meme.txt",
    "all TAPs (pooled)": f"{ARM}/motifs/all_classic/meme.txt",
    "C. elegans partners": f"{ROOT}/14_Celegans_MEME/motifs/worm_partners_noseeds/meme.txt",
}

def parse_meme(path):
    txt = open(path).read().splitlines()
    out, i = [], 0
    while i < len(txt):
        m = re.search(r"letter-probability matrix: alength= 20 w= (\d+) nsites= (\d+) E= (\S+)", txt[i])
        if m:
            w, ns, E = int(m.group(1)), int(m.group(2)), m.group(3)
            rows = [[float(x) for x in txt[j].split()] for j in range(i+1, i+1+w)]
            out.append({"w": w, "nsites": ns, "E": E, "prob": pd.DataFrame(rows, columns=AA)})
            i += w
        i += 1
    return out

# map (set, motif-index) from the assignment table
rows = list(csv.DictReader(open(f"{ARM}/analysis/motifs_by_set.csv")))
def find(set_name, family, nth=0):
    hits = [r for r in rows if r["set"] == set_name and r["assigned_family"] == family]
    if not hits:
        return None
    r = hits[nth if nth < len(hits) else 0]
    return set_name, int(r["motif"].split("-")[1]) - 1, r["assigned_family"], r["sites"], r["E_value"]

# headline picks that exist in the current run (kinase / P-loop machinery)
wanted = [
    find("α partners", "kinase catalytic loop"),
    find("α partners", "kinase activation seg."),
    find("α partners", "kinase Gly-rich loop"),
    find("α partners", "Walker A P-loop"),
    find("β partners", "kinase catalytic loop"),
    find("C. elegans partners", "kinesin motor") or find("C. elegans partners", "Walker A P-loop"),
]
picks = [w for w in wanted if w]

cache = {}
def matrix(set_name, idx):
    if set_name not in cache:
        cache[set_name] = parse_meme(MEME_PATH[set_name])
    return cache[set_name][idx]

LABEL = {"kinase catalytic loop": "Kinase catalytic loop (HRD...N)",
         "kinase activation seg.": "Kinase activation segment (DFG...APE)",
         "kinase Gly-rich loop": "Kinase glycine-rich loop (GxGxxG)",
         "Walker A P-loop": "Walker A P-loop (GxxxxGKT)",
         "kinesin motor": "Kinesin motor (P-loop NTPase)"}
SETLAB = {"α partners": "alpha partners", "β partners": "beta partners",
          "C. elegans partners": "C. elegans partners", "all TAPs (pooled)": "all TAPs"}

n = len(picks)
fig, axes = plt.subplots(n, 1, figsize=(12.5, 2.05*n))
if n == 1:
    axes = [axes]
for ax, (set_name, idx, fam, sites, E) in zip(axes, picks):
    mat = matrix(set_name, idx)
    info = logomaker.transform_matrix(mat["prob"], from_type="probability", to_type="information")
    logomaker.Logo(info, ax=ax, color_scheme="chemistry", stack_order="big_on_top")
    ax.set_ylim(0, np.log2(20)); ax.set_ylabel("bits", fontsize=9); ax.set_xticks([])
    ax.spines[["top", "right"]].set_visible(False)
    ax.set_title(f"{LABEL.get(fam, fam)}    —    {SETLAB.get(set_name, set_name)}, "
                 f"{mat['nsites']} sites, E = {E}", fontsize=10.5, loc="left",
                 weight="bold", color="#16323D")

fig.suptitle("MEME motif logos: the recurring signal in the tubulin partner proteomes is the "
             "kinase / P-loop nucleotide-binding fold", fontsize=12.5, y=1.005)
fig.tight_layout()
os.makedirs(f"{ARM}/figures", exist_ok=True)
fig.savefig(f"{ARM}/figures/M8_motif_logos.png", dpi=300, bbox_inches="tight", facecolor="white")
fig.savefig(f"{ARM}/figures/M8_motif_logos.pdf", bbox_inches="tight")
print("wrote M8_motif_logos with", n, "logos:",
      ", ".join(f"{LABEL.get(f,f)}({s})" for _, _, f, _, s in [(a,b,c,d,e) for (a,b,c,d,e) in picks]))

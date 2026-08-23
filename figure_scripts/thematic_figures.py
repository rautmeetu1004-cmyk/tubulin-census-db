#!/usr/bin/env python3
"""
Thematic manuscript figures — data panels + conceptual schematics
==================================================================
  T1  Structural role of tubulin, proven by the data (3 data panels):
        A conservation profile of the 32-tubulin alignment with the invariant
          GTP-binding motifs marked; B per-isoform ortholog constraint;
        C folding/chaperone partners are the highest-degree hub.
  T2  alpha/beta functional convergence + cross-species conservation (3 panels):
        A partner overlap (6.09x over chance); B functional-annotation Jaccard;
        C conserved worm->human core with disease genes highlighted.
  S1  Conceptual schematic: the alpha-beta heterodimer, N-/E-site GTP, invariant
        motifs, C-terminal tubulin-code tails.
  S2  Conceptual schematic: why C. elegans — human maps -> DIOPT orthology ->
        conserved disease core -> genetically tractable testbed.
All data panels read live CSVs; schematics are drawn, not data (labelled as such).
"""
from pathlib import Path
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Rectangle, Circle

def _project_root():
    import os
    from pathlib import Path as _P
    env = os.environ.get("TUBULIN_ROOT")
    if env:
        return _P(env)
    here = _P(__file__).resolve()
    for parent in here.parents:
        if (parent / "00_Publication_Package").exists() or (parent / ".git").exists():
            return parent
    return here.parent


ROOT = _project_root()
MSA = ROOT / "03_Structure_Evolution/Tubulin_MSA_Conservation"
CORR = ROOT / "04_Networks/Tubulin_Correlation_Analysis/data"
CMP = ROOT / "01_Census/Tubulin_AlphaBeta_Comparison/analysis"
CROSS = ROOT / "01_Census/Tubulin_AlphaBeta_Comparison/crossspecies/conserved_core.csv"
FIG = Path(__file__).resolve().parents[1] / "figures"
plt.rcParams.update({"font.size": 9, "svg.fonttype": "none",
                     "axes.spines.top": False, "axes.spines.right": False})
CA, CB = "#0072B2", "#D55E00"
OI = {"folding": "#0072B2", "motor": "#D55E00", "spindle": "#009E73",
      "centrosome": "#CC79A7", "cilium": "#56B4E9", "PTM": "#E69F00"}
# grounded disease/function labels for conserved-core genes (well-established only)
DISEASE = {"PRKN": "Parkinson's (Parkin)", "TBCE": "tubulinopathy (HRD)",
           "DYNC1H1": "dynein neuropathy", "TUBA4A": "ALS / motor neuron",
           "PLK1": "mitosis / cancer target", "TUBB4B": "tubulinopathy"}


# ------------------------------------------------------------------ T1
def fig_t1():
    fig = plt.figure(figsize=(12, 4.2))
    gs = fig.add_gridspec(1, 3, width_ratios=[1.5, 1, 1], wspace=0.32)

    # A: conservation profile + GTP motifs
    axA = fig.add_subplot(gs[0])
    prof = pd.read_csv(MSA / "analysis/all_tubulins_columns.csv")
    axA.fill_between(prof["column"], prof["identity"] * 100, color="#BFD6EC", lw=0)
    smooth = prof["identity"].rolling(7, center=True, min_periods=1).mean() * 100
    axA.plot(prof["column"], smooth, color="#185FA5", lw=1.1)   # eye-guide
    for (a, b, name), yl in zip([(79, 83, "VDLEP"), (152, 158, "GGGTGSG"),
                                 (191, 197, "VVEPYN")], [107, 107, 100]):
        axA.axvspan(a, b, color="#009E73", alpha=0.3)
        axA.annotate(name, xy=((a + b) / 2, 100), xytext=((a + b) / 2, yl),
                     ha="center", fontsize=7, color="#0F6E56",
                     arrowprops=dict(arrowstyle="-", color="#0F6E56", lw=0.5))
    axA.set_xlabel("alignment column"); axA.set_ylabel("% identity (34 tubulins)")
    axA.set_ylim(0, 112)
    axA.set_title("A  Invariant structural core\n(green = GTP-binding motifs)", fontsize=9.5)

    # B: per-isoform ortholog constraint
    axB = fig.add_subplot(gs[1])
    con = pd.read_csv(MSA / "ortholog_arm/analysis/per_isoform_constraint.csv")
    con = con.sort_values("mean_pairwise_identity")
    axB.barh(range(len(con)), con["mean_pairwise_identity"],
             color=[CA if f == "alpha" else CB for f in con["family"]])
    axB.set_yticks(range(len(con))); axB.set_yticklabels(con["isoform"], fontsize=6.5)
    _lo = float(con["mean_pairwise_identity"].min()); _hi = float(con["mean_pairwise_identity"].max())
    axB.set_xlim(max(0, _lo - 4), _hi + 2); axB.set_xlabel("mean ortholog % identity")
    for _i, _v in enumerate(con["mean_pairwise_identity"]):
        axB.text(_v + 0.2, _i, f"{_v:.0f}", va="center", fontsize=6)
    axB.set_title("B  Per-isoform\nevolutionary constraint", fontsize=9.5)

    # C: folding partners = highest degree
    axC = fig.add_subplot(gs[2])
    med = []
    for fam in ["alpha", "beta"]:
        f = pd.read_csv(CORR / f"partner_feature_table_{fam}.csv")
        f = f[(f["reactome_category"] != "other") & (f["degree"] > 0)]
        m = f.groupby("reactome_category")["degree"].median()
        med.append(m.rename(fam))
    md = pd.concat(med, axis=1).reindex([c for c in OI if c in pd.concat(med, axis=1).index])
    x = np.arange(len(md))
    axC.bar(x - 0.2, md["alpha"], 0.4, color=CA, label="α")
    axC.bar(x + 0.2, md["beta"], 0.4, color=CB, label="β")
    axC.set_xticks(x); axC.set_xticklabels(md.index, rotation=45, ha="right", fontsize=7)
    axC.set_ylabel("median degree"); axC.legend(frameon=False, fontsize=7)
    axC.set_title("C  Folding = pan-tubulin\nstructural hub", fontsize=9.5)

    fig.suptitle("The structural role of tubulin, proven by the data", y=1.04, fontsize=12)
    for e in ("pdf", "png"):
        fig.savefig(FIG / f"T1_structural_role.{e}", bbox_inches="tight", dpi=200)
    plt.close(fig); print("T1 done")


# ------------------------------------------------------------------ T2
def fig_t2():
    stats = json.load(open(CMP / "comparison_stats.json"))
    fig = plt.figure(figsize=(12, 4.4))
    gs = fig.add_gridspec(1, 3, width_ratios=[1, 1, 1.4], wspace=0.3)

    # A: partner overlap
    axA = fig.add_subplot(gs[0])
    po = stats["partner_overlap"]
    axA.bar([0, 1, 2], [po["alpha_only"], po["shared"], po["beta_only"]],
            color=[CA, "#7B4FA0", CB])
    axA.set_xticks([0, 1, 2]); axA.set_xticklabels(["α-only", "shared", "β-only"])
    axA.set_ylabel("partners")
    for i, v in enumerate([po["alpha_only"], po["shared"], po["beta_only"]]):
        axA.text(i, v, str(v), ha="center", va="bottom", fontsize=8)
    axA.set_title(f"A  Partner convergence\n{po['fold_enrichment']:.1f}× over chance, p≈0", fontsize=9.5)

    # B: functional-annotation Jaccard
    axB = fig.add_subplot(gs[1])
    fs = stats["functional_similarity"]
    cats = list(fs.keys()); jac = [fs[c]["jaccard"] for c in cats]
    axB.bar(range(len(cats)), jac, color="#7B4FA0")
    axB.set_xticks(range(len(cats))); axB.set_xticklabels(cats, rotation=45, ha="right", fontsize=7)
    axB.set_ylabel("α/β Jaccard"); axB.set_ylim(0, 1)
    axB.set_title("B  Shared functional\nannotations", fontsize=9.5)

    # C: conserved worm->human core, disease highlighted
    # always include disease genes (e.g. PRKN, 7 methods) even if below the top cutoff
    axC = fig.add_subplot(gs[2])
    ccf = pd.read_csv(CROSS)
    top = ccf.sort_values("diopt_methods", ascending=False).head(18)
    dis = ccf[ccf["human_ortholog"].isin(DISEASE)]
    cc = pd.concat([top, dis]).drop_duplicates("human_ortholog")
    cc = cc.sort_values("diopt_methods")
    colors = ["#B03060" if g in DISEASE else "#888" for g in cc["human_ortholog"]]
    axC.barh(range(len(cc)), cc["diopt_methods"], color=colors)
    axC.set_yticks(range(len(cc)))
    labels = [f"{w}→{h}" + (" ★" if h in DISEASE else "") for w, h in
              zip(cc["worm"], cc["human_ortholog"])]
    axC.set_yticklabels(labels, fontsize=6)
    axC.set_xlabel("DIOPT method support (max 9)")
    axC.set_title("C  Conserved worm→human core\n(★ = disease gene)", fontsize=9.5)

    fig.suptitle("α and β converge functionally; the core is conserved to C. elegans",
                 y=1.04, fontsize=12)
    for e in ("pdf", "png"):
        fig.savefig(FIG / f"T2_convergence_conservation.{e}", bbox_inches="tight", dpi=200)
    plt.close(fig); print("T2 done")


# ------------------------------------------------------------------ S1 schematic
def fig_s1():
    fig, ax = plt.subplots(figsize=(8, 6))
    ax.set_xlim(0, 10); ax.set_ylim(0, 10); ax.axis("off")
    # alpha (top) and beta (bottom) monomers
    ax.add_patch(FancyBboxPatch((3, 5.4), 4, 2.2, boxstyle="round,pad=0.1",
                                fc="#CFE3F5", ec=CA, lw=2))
    ax.add_patch(FancyBboxPatch((3, 2.9), 4, 2.2, boxstyle="round,pad=0.1",
                                fc="#F6DCCB", ec=CB, lw=2))
    ax.text(5, 6.5, "α-tubulin", ha="center", fontsize=13, color=CA, weight="bold")
    ax.text(5, 4.0, "β-tubulin", ha="center", fontsize=13, color=CB, weight="bold")
    # GTP sites
    ax.add_patch(Circle((5, 5.4), 0.28, fc="#009E73", ec="k"))
    ax.text(5.5, 5.4, "N-site: GTP locked (structural)", fontsize=8, va="center")
    ax.add_patch(Circle((5, 2.9), 0.28, fc="#E69F00", ec="k"))
    ax.text(5.5, 2.9, "E-site: GTP↔GDP switch (dynamics)", fontsize=8, va="center")
    # invariant motifs annotation
    ax.text(1.4, 4.15, "invariant\nmotifs\nVDLEP\nGGGTGSG\nVVEPYN", ha="center",
            fontsize=7.5, color="#0F6E56",
            bbox=dict(boxstyle="round", fc="#E0F2EA", ec="#009E73"))
    ax.annotate("", xy=(3, 4.4), xytext=(2.3, 4.4),
                arrowprops=dict(arrowstyle="->", color="#009E73"))
    # C-terminal tails (tubulin code)
    for y, c in [(6.4, CA), (3.9, CB)]:
        ax.plot([7, 8.6], [y, y + 0.3], color=c, lw=1.5)
        ax.plot([7, 8.6], [y, y - 0.2], color=c, lw=1.5)
    ax.text(8.8, 5.2, "C-terminal tails\n= 'tubulin code'\n(PTMs, partner\nrecruitment)",
            fontsize=7.5, va="center")
    # heterodimer bracket
    ax.annotate("", xy=(2.6, 3.0), xytext=(2.6, 7.5),
                arrowprops=dict(arrowstyle="-", color="k", lw=1))
    ax.text(2.2, 5.25, "αβ-heterodimer\n(microtubule\nbuilding block)",
            rotation=90, ha="center", va="center", fontsize=8)
    ax.set_title("S1  Structure of the αβ-tubulin heterodimer  (schematic)",
                 fontsize=11)
    for e in ("pdf", "png"):
        fig.savefig(FIG / f"S1_heterodimer_schematic.{e}", bbox_inches="tight", dpi=200)
    plt.close(fig); print("S1 done")


# ------------------------------------------------------------------ S2 schematic
def fig_s2():
    fig, ax = plt.subplots(figsize=(10, 5.5))
    ax.set_xlim(0, 12); ax.set_ylim(0, 8); ax.axis("off")

    def box(x, y, w, h, text, fc, ec, fs=9):
        ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.08",
                                    fc=fc, ec=ec, lw=1.8))
        ax.text(x + w / 2, y + h / 2, text, ha="center", va="center", fontsize=fs)

    box(0.3, 5.3, 2.6, 1.8, "Human α-map\n1,801 partners", "#CFE3F5", CA)
    box(0.3, 1.0, 2.6, 1.8, "Human β-map\n2,100 partners", "#F6DCCB", CB)
    box(4.2, 3.1, 2.6, 1.8, "DIOPT\northology\n(9 methods)", "#EDE7F6", "#7B4FA0")
    box(8.0, 5.0, 3.4, 1.9, "Conserved core\n33 worm→human\northologs", "#E0F2EA", "#009E73")
    box(8.0, 1.0, 3.4, 1.9, "Disease testbed\nPRKN (Parkinson),\nTBCE, TUBA4A…", "#FBE3E8", "#B03060")
    for (x0, y0), (x1, y1) in [((2.9, 6.2), (4.2, 4.4)), ((2.9, 1.9), (4.2, 3.6)),
                               ((6.8, 4.4), (8.0, 5.6)), ((9.7, 5.0), (9.7, 2.9))]:
        ax.add_patch(FancyArrowPatch((x0, y0), (x1, y1), arrowstyle="->",
                                     mutation_scale=14, color="#555", lw=1.4))
    ax.text(6.0, 0.2, "Why C. elegans: the machinery is conserved AND the animal is "
            "genetically tractable — knock out pdr-1/PRKN, ebp-1/EB1, plk-1 and test "
            "the human hypotheses in vivo.", ha="center", fontsize=8.5, style="italic")
    ax.set_title("S2  Why C. elegans is the right model  (schematic)", fontsize=11)
    for e in ("pdf", "png"):
        fig.savefig(FIG / f"S2_why_celegans_schematic.{e}", bbox_inches="tight", dpi=200)
    plt.close(fig); print("S2 done")


if __name__ == "__main__":
    # T1/T2 depend on the MSA-conservation arm, which was archived on 2026-07-17
    # (_archive/removed_2026-07-17/); both figures were pruned from the publication
    # package, so a missing MSA input is expected, not an error. S1/S2 are the kept
    # schematics and have no such dependency.
    for name, fn in (("T1", fig_t1), ("T2", fig_t2)):
        try:
            fn()
        except FileNotFoundError as e:
            print(f"  skipped {name}: needs the archived MSA arm ({e.filename})")
    fig_s1(); fig_s2()
    print(f"Thematic figures -> {FIG}")

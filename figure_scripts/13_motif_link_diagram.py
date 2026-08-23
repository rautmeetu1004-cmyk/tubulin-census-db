#!/usr/bin/env python3
"""The alpha / beta / C. elegans motif diagram: every motif, and every match between them.

Script 12 produces the counts. This produces the picture they are counted from -- the three
sets side by side with a line drawn for every Tomtom match at q < 0.05, so "four shared
between alpha and beta" is something you can see rather than something you are told.

WHY NOT A VENN

A Venn needs clean set membership, and the data do not have it. alpha's Gly-rich loop
matches beta's Gly-rich loop at q = 3.4e-31 AND C. elegans' Walker A P-loop at q = 0.024 --
both are glycine-rich nucleotide loops, so the "sets" overlap partially, at different
strengths. Drawing edges keeps that visible; a Venn would have to pick a compartment and
throw the q-values away.

JOIN KEY

Tomtom's Query_ID / Target_ID are the MEME motif IDs verbatim, so they join to the
`MOTIF <id> MEME-n` lines exactly. Do NOT join on Query_consensus: Tomtom writes resolved
residues (…NILLDEDGHVKI…) where meme.txt writes ambiguity codes (…NILLBEDGHVKJ…), and a
consensus join silently matches nothing.

Run:  python3 13_TAP_AllProtein_MEME/scripts/13_motif_link_diagram.py
Needs 12_cross_set_motifs.py to have run (it writes tomtom_cross/).
"""
import csv
import os
import re
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import FancyBboxPatch, PathPatch
from matplotlib.path import Path as MPath


def _project_root():
    env = os.environ.get("TUBULIN_ROOT")
    if env:
        return Path(env)
    here = Path(__file__).resolve()
    for parent in here.parents:
        if (parent / "00_Publication_Package").exists() or (parent / ".git").exists():
            return parent
    return here.parent


ROOT = _project_root()
ARM = ROOT / "13_TAP_AllProtein_MEME"
AN, FIG, TT = ARM / "analysis", ARM / "figures", ARM / "tomtom_cross"

SETS = [
    ("alpha", "α partners", ARM / "motifs_by_family/alpha/meme.txt", 1867),
    ("beta", "β partners", ARM / "motifs_by_family/beta/meme.txt", 2429),
    ("worm", "C. elegans partners",
     ROOT / "14_Celegans_MEME/motifs/worm_partners_noseeds/meme.txt", 123),
]
PAIRS = [("alpha", "beta"), ("alpha", "worm"), ("beta", "worm")]

COL = {"kinase catalytic loop": "#B03030", "kinesin motor": "#E06666",
       "kinase Gly-rich loop": "#D9743B", "kinase activation seg.": "#D9A441",
       "Walker A P-loop": "#8E5BA6", "poly-proline": "#6A8D3A", "acidic run": "#4F9DA6",
       "basic / histone-like": "#7B6A94", "Cys-loop channel": "#2A9D8F",
       "other": "#B8C2C7"}
INK, MUTE = "#16323D", "#5E7A86"

plt.rcParams.update({"font.family": "DejaVu Sans", "axes.edgecolor": "#4A4A4A",
                     "figure.dpi": 150})


def family(c):
    """Same classifier as 12_cross_set_motifs.py -- kept identical on purpose."""
    u = c.upper()
    if "PYRDSK" in u or "HRDSK" in u:
        return "kinesin motor"
    if "HRDLK" in u:
        return "kinase catalytic loop"
    if re.search(r"GKG[SA]FG", u):
        return "kinase Gly-rich loop"
    if "DVWS" in u:
        return "kinase activation seg."
    if re.search(r"G[PQ].G[TS]GKT|GPPG.GKT", u):
        return "Walker A P-loop"
    if u.count("P") >= len(u) * 0.5:
        return "poly-proline"
    if u.count("E") >= len(u) * 0.45:
        return "acidic run"
    if u.count("K") + u.count("R") >= len(u) * 0.4:
        return "basic / histone-like"
    if "SWVSFW" in u or "CFIDLASY" in u:
        return "Cys-loop channel"
    return "other"


def parse(path):
    out = []
    for line in open(path):
        m = re.match(r"MOTIF\s+(\S+)\s+MEME-(\d+).*?sites =\s*(\d+).*?E-value = (\S+)", line)
        if m:
            out.append({"id": m.group(1), "rank": int(m.group(2)),
                        "sites": int(m.group(3)), "evalue": m.group(4),
                        "family": family(m.group(1))})
    return out


def main():
    missing = [str(p) for _, _, p, _ in SETS if not Path(p).exists()]
    if missing:
        sys.exit(f"ABORT: missing motif files {missing}. Run 11 and the C. elegans arm.")
    if not TT.exists():
        sys.exit(f"ABORT: {TT} not found. Run 12_cross_set_motifs.py first -- this "
                 f"script draws its Tomtom output, it does not recompute it.")
    FIG.mkdir(parents=True, exist_ok=True)

    motifs = {key: parse(path) for key, _, path, _ in SETS}
    for key, ms in motifs.items():
        if not ms:
            sys.exit(f"ABORT: no motifs parsed for {key}")

    # ---- edges, joined on motif ID
    edges = []
    for a, b in PAIRS:
        tsv = TT / f"{a}_vs_{b}" / "tomtom.tsv"
        if not tsv.exists():
            sys.exit(f"ABORT: missing {tsv}; re-run 12_cross_set_motifs.py")
        ida = {m["id"] for m in motifs[a]}
        idb = {m["id"] for m in motifs[b]}
        for r in csv.DictReader(open(tsv), delimiter="\t"):
            if not r.get("Query_ID") or r["Query_ID"].startswith("#"):
                continue
            if float(r["q-value"]) >= 0.05:
                continue
            if r["Query_ID"] not in ida or r["Target_ID"] not in idb:
                sys.exit(f"ABORT: Tomtom id {r['Query_ID']}/{r['Target_ID']} is absent "
                         f"from the parsed {a}/{b} motif lists -- the join key is wrong "
                         f"or the motif files were regenerated after Tomtom ran.")
            edges.append({"a": a, "b": b, "ida": r["Query_ID"], "idb": r["Target_ID"],
                          "q": float(r["q-value"]), "overlap": int(r["Overlap"])})
    if not edges:
        sys.exit("ABORT: zero cross-set matches parsed. The kinase loop is present in "
                 "both human sets, so this is a parsing failure, not a result.")

    # ---- layout: three columns, motifs ordered by MEME rank
    XC = {"alpha": 0.0, "beta": 1.0, "worm": 2.0}
    ORDER = {"alpha": 0, "beta": 1, "worm": 2}
    BW, BH, GAP = 0.62, 0.062, 0.028
    GUT = -0.70                                 # left gutter, keeps clear of alpha column
    pos = {}
    for key, _, _, _ in SETS:
        ms = sorted(motifs[key], key=lambda m: m["rank"])
        top = 0.94
        for i, m in enumerate(ms):
            y = top - i * (BH + GAP)
            pos[(key, m["id"])] = (XC[key], y)

    fig, ax = plt.subplots(figsize=(17.8, 7.4))
    ax.set_xlim(-1.05, 2.45)
    ax.set_ylim(0.02, 1.10)
    ax.axis("off")

    # edges first, so boxes sit on top of them
    used_slots = []
    for e in sorted(edges, key=lambda e: e["q"]):
        x0, y0 = pos[(e["a"], e["ida"])]
        x1, y1 = pos[(e["b"], e["idb"])]
        x0 += BW / 2
        x1 -= BW / 2
        col = COL.get(family(e["ida"]), "#B8C2C7")
        strong = e["q"] < 1e-3
        lw = 3.4 if e["q"] == 0 else (2.4 if strong else 1.2)
        mid = (x0 + x1) / 2
        verts = [(x0, y0), (mid, y0), (mid, y1), (x1, y1)]
        ax.add_patch(PathPatch(MPath(verts, [MPath.MOVETO, MPath.CURVE4,
                                             MPath.CURVE4, MPath.CURVE4]),
                               fc="none", ec=col, lw=lw,
                               alpha=0.85 if strong else 0.45, zorder=2))
        # Label placement: adjacent columns have empty space at the midpoint, but an
        # alpha->C.elegans span midpoint lands ON the beta column, so those labels are
        # pinned just before the target box instead.
        adjacent = abs(ORDER[e["b"]] - ORDER[e["a"]]) == 1
        if adjacent:
            lx, ly, ha = mid, (y0 + y1) / 2 + 0.012, "center"
        else:
            lx, ly, ha = x1 - 0.05, y1 + 0.014, "right"
        while any(abs(lx - ux) < 0.30 and abs(ly - uy) < 0.030 for ux, uy in used_slots):
            ly -= 0.034
        used_slots.append((lx, ly))
        qtxt = "q = 0" if e["q"] == 0 else f"q = {e['q']:.0e}".replace("e-0", "e−")
        ax.text(lx, ly, f"{qtxt} · {e['overlap']} aa", ha=ha, va="bottom", fontsize=7.4,
                color=col if strong else MUTE, zorder=5,
                bbox=dict(fc="white", ec="none", pad=1.0, alpha=0.9))

    # boxes
    for key, label, _, n in SETS:
        ms = sorted(motifs[key], key=lambda m: m["rank"])
        xc = XC[key]
        ax.text(xc, 1.035, label, ha="center", fontsize=13.5, color=INK, weight="bold")
        ax.text(xc, 1.005, f"{n:,} sequences · {len(ms)} motifs",
                ha="center", fontsize=9.5, color=MUTE)
        for m in ms:
            _, y = pos[(key, m["id"])]
            col = COL.get(m["family"], "#B8C2C7")
            named = m["family"] != "other"
            ax.add_patch(FancyBboxPatch((xc - BW / 2, y - BH / 2), BW, BH,
                                        boxstyle="round,pad=0.004,rounding_size=0.012",
                                        fc=col if named else "#F2F5F6",
                                        ec=col, lw=1.1,
                                        alpha=1.0 if named else 0.9, zorder=3))
            ax.text(xc, y + 0.011, m["id"][:29], ha="center", va="center",
                    fontsize=6.9, family="monospace", zorder=4,
                    color="white" if named else MUTE)
            ax.text(xc, y - 0.014, f"MEME-{m['rank']} · {m['sites']} sites · "
                                   f"{100 * m['sites'] / n:.1f}% · E {m['evalue']}",
                    ha="center", va="center", fontsize=6.0, zorder=4,
                    color="#F0F4F5" if named else MUTE)

    # the pooled run, in the left gutter as context
    gl = GUT - BW / 2
    ax.text(gl, 1.035, "pooled run", fontsize=13.5, color=INK, weight="bold")
    ax.text(gl, 1.005, "2,734 TAPs · 1 motif", fontsize=9.5, color=MUTE)
    ax.add_patch(FancyBboxPatch((gl, 0.94 - BH / 2), BW, BH,
                                boxstyle="round,pad=0.004,rounding_size=0.012",
                                fc=COL["kinase catalytic loop"],
                                ec=COL["kinase catalytic loop"], lw=1.1, zorder=3))
    ax.text(GUT, 0.951, "IIHRDLKPENILLDS", ha="center", va="center",
            fontsize=6.9, family="monospace", color="white", zorder=4)
    ax.text(GUT, 0.926, "MEME-1 · 129 sites · 4.7% · E 1.3e-133",
            ha="center", va="center", fontsize=6.0, color="#F0F4F5", zorder=4)
    ax.text(gl, 0.885,
            "One motif only, because 125 kinases\n"
            "outvote everyone else. Split by family,\n"
            "α and β each return ten — both still led\n"
            "by that same catalytic loop.",
            fontsize=8.4, color=MUTE, ha="left", va="top")

    ax.text(gl, 0.70,
            "Lines are Tomtom matches at q < 0.05;\n"
            "thick = q < 1e−3. Grey boxes resolve to no\n"
            "catalogued domain. Consensus strings differ\n"
            "between sets, so an exact-string test finds\n"
            "nothing — the motifs are the same all the same.",
            fontsize=8.4, color=MUTE, ha="left", va="top")

    ax.text(gl, 0.52,
            f"α–β  {sum(1 for e in edges if e['b'] == 'beta')} matches\n"
            f"α–C. elegans  {sum(1 for e in edges if e['a'] == 'alpha' and e['b'] == 'worm')}\n"
            f"β–C. elegans  {sum(1 for e in edges if e['a'] == 'beta')}\n\n"
            "The kinase / P-loop nucleotide-binding fold\n"
            "is the only signal present in all three.",
            fontsize=8.6, color=INK, ha="left", va="top")

    used = [k for k in COL if any(m["family"] == k for ms in motifs.values() for m in ms)]
    fig.legend(handles=[plt.Line2D([], [], marker="s", ls="", color=COL[k],
                                   markersize=10, label=k) for k in used],
               loc="lower center", bbox_to_anchor=(0.5, -0.035), ncol=5,
               frameon=False, fontsize=9)

    fig.suptitle("Motif discovery run separately on human α partners, human β partners and "
                 "C. elegans partners — and every match between them",
                 fontsize=13.5, x=0.012, ha="left", y=0.995)
    fig.tight_layout()
    out = FIG / "M7_alpha_beta_worm_links.png"
    fig.savefig(out, dpi=200, bbox_inches="tight", facecolor="white")
    plt.close(fig)

    with open(AN / "motif_links_alpha_beta_worm.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["a", "b", "ida", "idb", "q", "overlap"])
        w.writeheader(); w.writerows(edges)

    n_ab = sum(1 for e in edges if e["a"] == "alpha" and e["b"] == "beta")
    n_aw = sum(1 for e in edges if e["b"] == "worm" and e["a"] == "alpha")
    n_bw = sum(1 for e in edges if e["b"] == "worm" and e["a"] == "beta")
    print(f"  α–β {n_ab} · α–C.elegans {n_aw} · β–C.elegans {n_bw}  (q < 0.05)")
    print(f"  -> {out}")
    print(f"  -> {AN / 'motif_links_alpha_beta_worm.csv'}")


if __name__ == "__main__":
    main()

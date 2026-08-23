#!/usr/bin/env python3
"""Publication-quality PPI network figures: human alpha, human beta, C. elegans.

Human networks come from the replication-filtered induced subgraph (an edge needs two
independent publications or two databases), restricted to high-confidence partners, so
what is drawn is the network the study's claims actually rest on.

The worm network is built here from the worm BioGRID and HINT bulk files, because the
worm arm never had an induced network of its own.

Two design decisions worth stating, both driven by what the data actually looks like:

1. Nodes are coloured by NETWORK COMMUNITY (Newman-Girvan greedy modularity), not by
   Reactome module. Reactome covers only 5% of the drawn nodes, so module colouring
   left 95% grey and said nothing. Communities are derived from the network itself,
   cover every node, and make the modularity result (Q about 0.32) visible.

2. Layout runs on the largest connected component. Disconnected singletons (TUBA3D in
   alpha) otherwise push spring_layout apart and crush the real network into a sliver.
   Such nodes are reported in the caption rather than silently dropped.
"""
import csv
from collections import Counter, defaultdict
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import networkx as nx
from networkx.algorithms.community import greedy_modularity_communities, modularity
from adjustText import adjust_text
import matplotlib.patheffects as pe

import sys as _sys
from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parents[2] / "_scripts"))
from evidence_tiers import (ALL_TIERS, ASSAY_TIERS,  # noqa: E402,F401
                            load_partner_union, load_partner_union_df)

# ASSAY_TIERS: this arm makes a claim about biology, and a degree, denominator or
# network edge computed over a text-mined STRING edge overstates its confidence.
TIERS_USED = ASSAY_TIERS


def _project_root():
    """Walk up from this file to the repo root. Honours TUBULIN_ROOT if set."""
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


PROJECT_ROOT = _project_root()


S = PROJECT_ROOT
CEN = S / "01_Census"; NET = S / "04_Networks"; ANN = S / "02_Annotation"; OUT = S / "07_Figures"
csv.field_size_limit(10_000_000)

# Okabe-Ito colour-blind-safe palette for communities
PAL = ["#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00", "#56B4E9", "#8C564B", "#4C72B0"]
GREY = "#C9D5DB"; SEEDC = "#B5179E"; INK = "#16323D"
plt.rcParams.update({"font.family": "DejaVu Sans", "figure.dpi": 150})

SEEDS = {
    "alpha": {"TUBA1A", "TUBA1B", "TUBA1C", "TUBA3C", "TUBA3D", "TUBA4A", "TUBA8", "TUBA3E", "TUBAL3"},
    "beta": {"TUBB", "TUBB1", "TUBB2A", "TUBB2B", "TUBB3", "TUBB4A", "TUBB4B",
             "TUBB6", "TUBB8", "TUBB8B"},
}
WORM_SEEDS = {"tba-1", "tba-2", "tba-4", "tba-5", "tba-6", "tba-7", "tba-8", "tba-9",
              "mec-12", "tbb-1", "tbb-2", "tbb-4", "tbb-6", "mec-7", "ben-1"}


def loc_class():
    """broad localisation per gene, used to name each community"""
    out = {}
    KEY = [("cytoskeleton / MT", ("microtubule", "cytoskelet", "spindle", "centrosome", "cilium", "axoneme")),
           ("nucleus", ("nucleus", "nucleopl", "nucleol", "chromat")),
           ("membrane", ("plasma membrane", "membrane")),
           ("mitochondrion", ("mitochond",)),
           ("secreted / vesicle", ("exosome", "extracellular", "vesicle", "golgi", "endoplasmic"))]
    per = defaultdict(Counter)
    p = ANN / "Tubulin_Partner_Annotation/analysis/localization.csv"
    if p.exists():
        for r in csv.DictReader(open(p)):
            if r["source"] != "GO_CC":
                continue
            t = r["term"].lower()
            for lab, keys in KEY:
                if any(k in t for k in keys):
                    per[r["gene"]][lab] += 1
    for g, c in per.items():
        out[g] = c.most_common(1)[0][0]
    return out


def human_graph(fam):
    hc = {r["partner_canonical"] for r in csv.DictReader(
        open(CEN / f"Tubulin_{fam.capitalize()}_Interactome/03_parsed_data/partner_highconfidence.csv"))}
    keep = hc | SEEDS[fam]
    G = nx.Graph()
    for r in csv.DictReader(open(NET / f"Tubulin_Network_Topology/{fam}/data/edges.csv")):
        if r["confidence"] != "high":
            continue
        a, b = r["source_node"], r["target_node"]
        if a in keep and b in keep:
            G.add_edge(a, b)
    G.add_nodes_from(SEEDS[fam] & keep)
    return G


def worm_graph():
    RAW = CEN / "Tubulin_Celegans_Interactome/02_raw_data"
    partners = {r["partner"] for r in load_partner_union(
        CEN / "Tubulin_Celegans_Interactome/03_parsed_data/partner_union.csv",
        TIERS_USED, label="worm")[0]}
    keep = partners | WORM_SEEDS
    G = nx.Graph()
    f = next(RAW.glob("biogrid/*Caenorhabditis*.tab3.txt"), None)
    if f:
        with open(f, encoding="utf-8", errors="ignore") as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                if r.get("Experimental System Type") != "physical":
                    continue
                a = r.get("Official Symbol Interactor A"); b = r.get("Official Symbol Interactor B")
                if a in keep and b in keep and a != b:
                    G.add_edge(a, b)
    h = RAW / "hint/CaenorhabditisElegans_binary_hq.txt"
    if h.exists():
        with open(h) as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                a = r.get("Gene_A") or r.get("Alias_A"); b = r.get("Gene_B") or r.get("Alias_B")
                if a in keep and b in keep and a != b:
                    G.add_edge(a, b)
    G.add_nodes_from(WORM_SEEDS & keep)
    # Both bulk files are raw downloads and are gitignored, so a clone -- or a cleanup
    # that removed 02_raw_data/biogrid and 02_raw_data/hint -- leaves this graph with
    # nodes and no edges. modularity() then divides by a zero degree sum and the run
    # dies inside networkx with a bare ZeroDivisionError, which says nothing about the
    # cause. Fail here instead, naming the files that are missing.
    if G.number_of_edges() == 0:
        have_b = f is not None
        have_h = h.exists()
        _sys.exit(
            "ABORT: the C. elegans graph has 0 edges, so PPI_celegans.png cannot be "
            "drawn.\n"
            f"  BioGRID  {RAW / 'biogrid'}/*Caenorhabditis*.tab3.txt : "
            f"{'found' if have_b else 'MISSING'}\n"
            f"  HINT     {h} : {'found' if have_h else 'MISSING'}\n"
            "  Both are bulk downloads, gitignored in the parent repo. Restore them with\n"
            "    python3 01_Census/Tubulin_Celegans_Interactome/01_scripts/pipeline/"
            "celegans_fetch_bulk_ppi.py\n"
            "  which pins BioGRID 5.0.259 and HINT 2024-06. Do NOT ship the stale PNG as "
            "if it were current.")
    return G


def _layout(LCC, comms, seed=11):
    """ForceAtlas2, with spring_layout as a fallback on older networkx.

    Three layouts were rendered and compared on the alpha graph before choosing:

        spring_layout       one hairball plus two satellite fans, ~55% of the canvas empty
        kamada_kawai        fills the canvas best but mixes every community into one blob,
                            so the Q the caption reports is invisible
        forceatlas2         communities separate cleanly -- the CCT/TriC chaperonin block
                            reads as a block -- with far less dead space

    FA2 was picked on the second criterion, not the first. All three are layouts of the
    SAME graph: nothing here changes a node, an edge, a community or Q. Only which of many
    equally valid force-directed minima gets drawn.
    """
    if hasattr(nx, "forceatlas2_layout"):
        pos = nx.forceatlas2_layout(LCC, max_iter=500, scaling_ratio=2.0,
                                    gravity=1.0, seed=seed)
    else:                                    # networkx < 3.4
        pos = nx.spring_layout(LCC, seed=seed, iterations=600,
                               k=1.9 / np.sqrt(LCC.number_of_nodes()))
    # Rescale into [0, 1]. The raw output occupies a box far smaller than its own range,
    # which is where the white margins came from. One scale for both axes: never distort.
    xs = np.array([p[0] for p in pos.values()])
    ys = np.array([p[1] for p in pos.values()])
    span = max(float(np.ptp(xs)) or 1.0, float(np.ptp(ys)) or 1.0)
    return {node: np.array([(p[0] - xs.min()) / span, (p[1] - ys.min()) / span])
            for node, p in pos.items()}


def draw(G, seeds, title, out, locs=None, n_comm=6, label_hubs=10):
    locs = locs or {}
    comps = sorted(nx.connected_components(G), key=len, reverse=True)
    dropped = [n for c in comps[1:] for n in c]
    LCC = G.subgraph(comps[0]).copy()

    comms = sorted(greedy_modularity_communities(LCC), key=len, reverse=True)
    Q = modularity(LCC, comms)
    cmap = {}
    for i, c in enumerate(comms):
        for n in c:
            cmap[n] = i if i < n_comm else -1

    deg = dict(LCC.degree())
    pos = _layout(LCC, comms)
    col = {n: (PAL[cmap[n]] if cmap[n] >= 0 else GREY) for n in LCC}

    fig, ax = plt.subplots(figsize=(15, 12.4))

    # Edges split by community. Drawn as two passes so the modularity is legible instead
    # of being a uniform grey wash: within-community edges take the module colour, edges
    # that cross modules stay grey and sit underneath.
    intra = [(u, v) for u, v in LCC.edges() if cmap[u] == cmap[v] and cmap[u] >= 0]
    inter = [(u, v) for u, v in LCC.edges() if (u, v) not in set(intra)]
    nx.draw_networkx_edges(LCC, pos, edgelist=inter, ax=ax, edge_color="#C8D4DA",
                           width=0.35, alpha=0.55)
    nx.draw_networkx_edges(LCC, pos, edgelist=intra, ax=ax,
                           edge_color=[col[u] for u, _ in intra],
                           width=0.55, alpha=0.34)

    ns = [n for n in LCC if n not in seeds]
    nx.draw_networkx_nodes(LCC, pos, nodelist=ns, ax=ax,
                           node_color=[col[n] for n in ns],
                           node_size=[26 + 6.5 * deg[n] for n in ns],
                           linewidths=0.4, edgecolors="white", alpha=0.95)
    sd = [n for n in LCC if n in seeds]
    if sd:
        nx.draw_networkx_nodes(LCC, pos, nodelist=sd, ax=ax, node_color=SEEDC,
                               node_shape="*", node_size=[300 + 6 * deg[n] for n in sd],
                               linewidths=1.1, edgecolors=INK)

    # Labels. Every seed plus the top hubs, each with a white halo so it stays readable
    # over a dense patch, then adjustText with every node passed as an obstacle -- without
    # x/y the labels only avoid each other and still land on top of the stars.
    want = list(sd) + [n for n, _ in sorted(((n, deg[n]) for n in ns),
                                            key=lambda x: -x[1])[:label_hubs]]
    halo = [pe.withStroke(linewidth=2.6, foreground="white")]
    texts = [ax.text(pos[n][0], pos[n][1], n, fontsize=8.2, fontweight="bold",
                     color=INK if n not in seeds else SEEDC, ha="center", va="center",
                     path_effects=halo, zorder=6) for n in want]
    px = np.array([pos[n][0] for n in LCC])
    py = np.array([pos[n][1] for n in LCC])
    adjust_text(texts, x=px, y=py, ax=ax, expand=(1.35, 1.55),
                force_text=(0.4, 0.6), force_static=(0.25, 0.4),
                arrowprops=dict(arrowstyle="-", color="#8FA3AC", lw=0.55), time_lim=25)

    # Name each community by its highest-degree member. "Most common localisation"
    # was tried first and was useless: GO:CC "membrane" is so common that three of
    # five communities came out "mostly membrane". A lead protein is concrete and
    # lets a reader recognise the module (e.g. the CCT chaperonin).
    names = []
    for i, c in enumerate(comms[:n_comm]):
        inner = [n for n in c if n not in seeds] or list(c)
        lead = max(inner, key=lambda n: deg.get(n, 0))
        ll = Counter(locs.get(n) for n in c if locs.get(n))
        tag = ll.most_common(1)[0][0] if ll else ""
        names.append(f"C{i+1}  n={len(c)}, led by {lead}" + (f" ({tag})" if tag else ""))

    ax.set_title(title, fontsize=17, loc="left", color=INK, fontweight="bold", pad=28)
    sub = (f"{LCC.number_of_nodes()} proteins, {LCC.number_of_edges()} replicated interactions. "
           f"Colour = network community ({len(comms)} found, Q = {Q:.2f}); size = degree; star = tubulin.")
    if dropped:
        sub += f"  {len(dropped)} node(s) with no internal edge not shown: {', '.join(sorted(dropped)[:4])}."
    ax.text(0, 1.012, sub, transform=ax.transAxes, fontsize=10, color="#5E7A86", va="bottom")
    handles = [Line2D([], [], marker="o", ls="", ms=9, mfc=PAL[i], mec="white", label=names[i])
               for i in range(min(n_comm, len(comms)))]
    if len(comms) > n_comm:
        handles.append(Line2D([], [], marker="o", ls="", ms=9, mfc=GREY, mec="white",
                              label=f"{len(comms)-n_comm} smaller communities"))
    handles.append(Line2D([], [], marker="*", ls="", ms=16, mfc=SEEDC, mec=INK, label="tubulin"))
    ax.legend(handles=handles, fontsize=9.6, frameon=False, loc="upper left",
              bbox_to_anchor=(0, -0.005), ncol=3, handletextpad=0.4, columnspacing=1.6)
    ax.axis("off")
    ax.set_xlim(-0.04, max(p[0] for p in pos.values()) + 0.04)
    ax.set_ylim(-0.04, max(p[1] for p in pos.values()) + 0.04)
    ax.set_aspect("equal")
    fig.tight_layout()
    fig.savefig(out, dpi=260, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"  {Path(out).name}: {LCC.number_of_nodes()} nodes, {LCC.number_of_edges()} edges, "
          f"{len(comms)} communities Q={Q:.3f}, {len(dropped)} isolated")


if __name__ == "__main__":
    locs = loc_class()
    for fam, nice in (("alpha", "α"), ("beta", "β")):
        draw(human_graph(fam), SEEDS[fam],
             f"Human {nice}-tubulin high-confidence interactome",
             OUT / f"PPI_human_{fam}.png", locs=locs)
    draw(worm_graph(), WORM_SEEDS, "C. elegans tubulin interactome",
         OUT / "PPI_celegans.png", locs={}, n_comm=5, label_hubs=12)
    print("done ->", OUT)

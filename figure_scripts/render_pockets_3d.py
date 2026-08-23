#!/usr/bin/env python3
"""PyMOL renders of the tubulin drug pockets on the alpha/beta dimer.

One ray-traced panel per pocket, from the PDB that actually contains the ligand:

    Colchicine    1SA0  CN2   alpha/beta intradimer interface
    Vinca         1Z2B  VLB   interdimer interface
    Laulimalide   4O4H  LLM   beta, external site
    Pironetin     5FNV  X3H   alpha

Alpha and beta chains are coloured distinctly in every panel, so the reader can see
which chain each drug actually binds -- the point of the slide. Chain->family follows
the same T2R / T2R-TTL convention the analysis uses (A,C = alpha; B,D = beta).

    Taxane        1JFF  TA1   beta, luminal site

All five clinical pockets are rendered. 1JFF is the docetaxel/taxol-bound structure the
analysis already uses for the taxane contacts, downloaded here so the render comes from
the same structure the numbers do.

Run:  python3 15_Drug_Pocket_Conservation/scripts/render_pockets_3d.py
Needs pymol on PATH.
"""
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path


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
ARM = ROOT / "15_Drug_Pocket_Conservation"
STRUCT = ARM / "structures"
FIG = ARM / "figures"
RENDER = ARM / "renders"

ALPHA_COL, BETA_COL = "0x6E9BD6", "0xE0A06A"
LIG_COL, CONTACT_COL = "0xC0392B", "0xF2C744"

POCKETS = [
    ("Taxane",      "1jff", "TA1", "beta, luminal site"),
    ("Colchicine",  "1sa0", "CN2", "alpha/beta intradimer interface"),
    ("Vinca",       "1z2b", "VLB", "interdimer interface"),
    ("Laulimalide", "4o4h", "LLM", "beta, external site"),
    ("Pironetin",   "5fnv", "X3H", "alpha"),
]


def pml(name, pdb, lig):
    """PyMOL script for one pocket. Kept as text so the render is auditable."""
    out = RENDER / f"{name.lower()}.png"
    return f"""
load {STRUCT / (pdb + '.pdb')}, m
hide everything
remove solvent
remove resn HOH+GOL+EDO+SO4+PO4+ACT+MG+CA+ZN+GTP+GDP
bg_color white
set ray_opaque_background, 1
set cartoon_transparency, 0.55
set ray_shadows, 0
set antialias, 2
set specular, 0.2

# Tubulin only. T2R/T2R-TTL constructs also contain stathmin-like RB3 and TTL, which
# are crystallisation scaffold, not the drug target -- drawn faint grey, never coloured
# as if they were tubulin.
show cartoon, polymer
color grey80, polymer
color {ALPHA_COL}, polymer and chain A+C
color {BETA_COL},  polymer and chain B+D

# ONE ligand copy. These structures carry two symmetry copies of the drug; orienting on
# both zooms out to the whole filament and the pocket becomes invisible.
select ligall, resn {lig}
python
from pymol import cmd
ch = sorted({{a.chain for a in cmd.get_model("ligall").atom}})
cmd.select("lig", f"resn {lig} and chain {{ch[0]}}" if ch else "ligall")
python end

show spheres, lig
color {LIG_COL}, lig
set sphere_scale, 0.40, lig
select pocket, byres (polymer within 4.5 of lig)
show sticks, pocket and sidechain
color {CONTACT_COL}, pocket and elem C
set stick_radius, 0.14
orient lig
zoom lig, 6
turn x, -10
set field_of_view, 20
ray 1400, 1400
png {out}, dpi=200
"""


def main():
    if not shutil.which("pymol"):
        sys.exit("ABORT: pymol not on PATH. brew install pymol")
    missing = [p for _, p, _, _ in POCKETS if not (STRUCT / f"{p}.pdb").exists()]
    if missing:
        sys.exit(f"ABORT: missing structures: {missing}")
    RENDER.mkdir(parents=True, exist_ok=True)
    FIG.mkdir(parents=True, exist_ok=True)

    made = []
    for name, pdb, lig, loc in POCKETS:
        with tempfile.NamedTemporaryFile("w", suffix=".pml", delete=False) as fh:
            fh.write(pml(name, pdb, lig))
            script = fh.name
        r = subprocess.run(["pymol", "-cq", script], capture_output=True, text=True)
        os.unlink(script)
        out = RENDER / f"{name.lower()}.png"
        if not out.exists() or out.stat().st_size == 0:
            print(f"  FAILED {name}: {(r.stderr or r.stdout)[-300:]}")
            continue
        made.append((name, loc, out))
        print(f"  rendered {name:12s} {pdb.upper()} / {lig}  -> {out.name}")

    if not made:
        sys.exit("ABORT: no panel rendered; refusing to write an empty figure.")

    # assemble
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.image as mpimg
    import numpy as np

    # Grid, not a single row. Five panels in one row is 4:1 on the page, which leaves
    # each pocket about 2 cm wide on a slide -- small enough that readers have counted
    # three panels and asked where the other two went. A 3-over-2 grid doubles the size
    # of every panel for the same slide width.
    n = len(made)
    ncol = 3 if n > 3 else n
    nrow = -(-n // ncol)
    fig, axes = plt.subplots(nrow, ncol, figsize=(4.6 * ncol, 5.3 * nrow))
    axes = np.atleast_1d(axes).ravel()
    for ax, (name, loc, path) in zip(axes, made):
        ax.imshow(mpimg.imread(path))
        ax.set_title(f"{name}\n{loc}", fontsize=15, color="#16323D")
        ax.axis("off")
    for ax in axes[n:]:
        ax.axis("off")
    handles = [
        plt.Line2D([], [], marker="s", ls="", color="#6E9BD6", markersize=11, label="α-tubulin"),
        plt.Line2D([], [], marker="s", ls="", color="#E0A06A", markersize=11, label="β-tubulin"),
        plt.Line2D([], [], marker="o", ls="", color="#C0392B", markersize=10, label="ligand"),
        plt.Line2D([], [], marker="s", ls="", color="#F2C744", markersize=11, label="contact residues (4.5 Å)"),
    ]
    fig.legend(handles=handles, loc="lower center", ncol=4, frameon=False, fontsize=13,
               bbox_to_anchor=(0.5, -0.005))
    fig.tight_layout(rect=[0, 0.04, 1, 1])
    out = FIG / "P3_pockets_3d.png"
    fig.savefig(out, dpi=200, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"  {len(made)}/{len(POCKETS)} panels -> {out}")


if __name__ == "__main__":
    main()

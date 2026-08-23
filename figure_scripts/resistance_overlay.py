#!/usr/bin/env python3
"""Overlay literature-curated tubulin drug-resistance mutations onto the pockets, and test
whether resistance positions coincide with the isotype-variable pocket residues.
Resistance positions are LITERATURE-CURATED (Giannakakou 1997/2000; Kavallaris; Berrieman 2004;
Kanakkanthara 2011; Prota 2013 pironetin) — COSMIC is login-gated, so we do not use it here."""
import csv, math, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
import os as _os
from pathlib import Path as _P


def _project_root():
    """Walk up to the repo root. Honours TUBULIN_ROOT. Replaces a hardcoded absolute
    path, which broke the moment the project was copied or shared."""
    env = _os.environ.get("TUBULIN_ROOT")
    if env:
        return str(_P(env))
    here = _P(__file__).resolve()
    for parent in here.parents:
        if (parent / "00_Publication_Package").exists() or (parent / ".git").exists():
            return str(parent)
    return str(here.parent)


ROOT = _project_root()
SE = f"{ROOT}/03_Structure_Evolution/Tubulin_Structural_Evolution"
ARM = f"{ROOT}/15_Drug_Pocket_Conservation"

# per-residue human-isotype entropy (same as pocket study), keyed (family, human_resno)
def read_fasta(p):
    s={},; s={}; n=None
    for ln in open(p):
        ln=ln.rstrip("\n")
        if ln.startswith(">"): n=ln[1:].split()[0]; s[n]=""
        elif n: s[n]+=ln
    return s
def ent(col):
    a=[c for c in col if c not in "-."]
    if len(a)<2: return 0.0
    h=0.0
    for x in set(a):
        pf=a.count(x)/len(a); h-=pf*math.log2(pf)
    return h
def var_by_resno(aln, ref_key):
    s=read_fasta(aln); ref=[k for k in s if ref_key in k]; ref=ref[0] if ref else list(s)[0]
    mat=[s[k] for k in s]; E=[ent([m[c] for m in mat]) for c in range(len(s[ref]))]
    out={}; rn=0
    for c,ch in enumerate(s[ref]):
        if ch not in "-.": rn+=1; out[rn]=E[c]
    return out
VA=var_by_resno(f"{SE}/alignments/alpha_human.aln.fasta","TUBA1B")
VB=var_by_resno(f"{SE}/alignments/beta_human.aln.fasta","TUBB")

# map pdb_resno -> human_resno per family (from the ET map)
p2h={}
for r in csv.DictReader(open(f"{SE}/analysis/structure_et_map.csv")):
    if r["human_resno"].isdigit(): p2h[(r["family"],int(r["pdb_resno"]))]=int(r["human_resno"])

# pocket contact residues (from pocket_conservation run — re-derive taxane; others via CSV note)
# reuse the same detection by importing contact residues saved implicitly: recompute quickly here
import numpy as np
STR=f"{ARM}/structures"; CHAINMAP={"A":"alpha","B":"beta","C":"alpha","D":"beta"}
def contacts(pdb,lig,cut=4.5):
    L=[];
    for ln in open(f"{STR}/{pdb}.pdb"):
        if ln.startswith("HETATM") and ln[17:20].strip()==lig and ln[76:78].strip()!="H":
            L.append((float(ln[30:38]),float(ln[38:46]),float(ln[46:54])))
    L=np.array(L); prot={}
    for ln in open(f"{STR}/{pdb}.pdb"):
        if ln.startswith("ATOM") and ln[76:78].strip()!="H":
            fam=CHAINMAP.get(ln[21]);
            if not fam: continue
            try: rn=int(ln[22:26])
            except: continue
            prot.setdefault((fam,rn),[]).append((float(ln[30:38]),float(ln[38:46]),float(ln[46:54])))
    return sorted(k for k,c in prot.items() if np.sqrt(((np.array(c)[:,None]-L[None])**2).sum(-1)).min()<=cut)
TAX=[("beta",r) for r in [22,23,26,27,217,226,229,230,233,236,272,274,275,276,277,278,360,369,370]]
POCK={"Taxane":TAX,"Colchicine":contacts("1sa0","CN2"),"Vinca":contacts("1z2b","VLB"),
      "Laulimalide":contacts("4o4h","LLM"),"Pironetin":contacts("5fnv","X3H")}

# LITERATURE-CURATED resistance mutation positions (human numbering)
RES={"Taxane":{("beta",p) for p in [26,185,215,217,228,270,274,282,292,364]},
     "Colchicine":{("beta",p) for p in [236,238,239,248,350]}|{("alpha",316)},
     "Vinca":{("beta",p) for p in [390]},
     "Laulimalide":{("beta",p) for p in [306,312,313]},
     "Pironetin":{("alpha",316)}}

rows=[]
for name,res in POCK.items():
    # convert pocket residues to human numbering
    hres=set()
    for fam,rn in res:
        h=p2h.get((fam,rn))
        if h: hres.add((fam,h))
    rmark=RES.get(name,set())
    hit=hres & rmark
    var_hit=[ (VA if f=="alpha" else VB).get(h,0) for f,h in hit]
    var_oth=[ (VA if f=="alpha" else VB).get(h,0) for f,h in (hres-rmark)]
    rows.append(dict(pocket=name, n_pocket=len(hres), n_resistance=len(hit),
        pct=round(100*len(hit)/max(len(hres),1)),
        var_resistance=round(np.mean(var_hit),3) if var_hit else None,
        var_other=round(np.mean(var_oth),3) if var_oth else None))
with open(f"{ARM}/analysis/resistance_overlay.csv","w",newline="") as f:
    w=csv.DictWriter(f,fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)

# ---- figure ----
fig,ax=plt.subplots(1,2,figsize=(13,5))
names=[r["pocket"] for r in rows]; col={'Taxane':'#DD8452','Colchicine':'#6D6875','Vinca':'#4C72B0','Laulimalide':'#2A9D8F','Pironetin':'#E9A23B'}
ax[0].bar(names,[r["n_resistance"] for r in rows],color=[col[n] for n in names],edgecolor="k")
for i,r in enumerate(rows): ax[0].text(i,r["n_resistance"]+0.05,f"{r['n_resistance']}/{r['n_pocket']}",ha="center",fontsize=9)
ax[0].set_ylabel("pocket contact residues with a known resistance mutation")
ax[0].set_title("A  Documented resistance mutations per pocket\n(literature-curated)",fontsize=10,loc="left")
ax[0].tick_params(axis='x',rotation=20)
# B: taxane per-residue variability with resistance flagged
tax_h=sorted({p2h[("beta",rn)] for (_,rn) in TAX if ("beta",rn) in p2h})
vals=[VB.get(h,0) for h in tax_h]; isres=[("beta",h) in RES["Taxane"] for h in tax_h]
bars=ax[1].bar(range(len(tax_h)),vals,color=["#B0303A" if x else "#c9d3d8" for x in isres],edgecolor="k",lw=0.3)
ax[1].set_xticks(range(len(tax_h))); ax[1].set_xticklabels([f"β{h}" for h in tax_h],rotation=90,fontsize=7)
ax[1].set_ylabel("human-isotype variability (entropy)")
ax[1].set_title("B  Taxane pocket, per residue — red = known resistance position\n(resistance clusters on the isotype-variable residues)",fontsize=10,loc="left")
from matplotlib.patches import Patch
ax[1].legend(handles=[Patch(color="#B0303A",label="known resistance mutation"),Patch(color="#c9d3d8",label="other contact")],fontsize=8,frameon=False)
fig.tight_layout(); out=f"{ARM}/figures/P2_resistance_overlay.png"; fig.savefig(out,dpi=150,facecolor="white"); print("wrote",out)
print(f"{'pocket':12}{'n_res/n':>10}{'var_res':>10}{'var_other':>11}")
for r in rows: print(f"{r['pocket']:12}{str(r['n_resistance'])+'/'+str(r['n_pocket']):>10}{str(r['var_resistance']):>10}{str(r['var_other']):>11}")

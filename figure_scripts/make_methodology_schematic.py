#!/usr/bin/env python3
"""Methodology workflow schematic for the tubulin interactome census.
Card-based design: inputs -> resource -> analysis. Every number is taken from
the verified census. Renders 07_Figures/methodology_schematic.svg and .png .
"""
import os, cairosvg

OUT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
W, H = 860, 940

INK, SUB, BG = "#17262F", "#5A6B76", "#F4F7F9"
CARD, BORDER = "#FFFFFF", "#E3E9ED"
ACC1, TINT1 = "#2E6CA4", "#EAF2F8"   # inputs  – steel blue
ACC2, TINT2 = "#BE5738", "#F9EEE9"   # resources – terracotta
ACC3, TINT3 = "#2C8C74", "#E6F2EE"   # outputs – emerald
GOLD = "#C4922B"

def esc(s): return s.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")

def txt(x, y, s, size, col, weight="normal", anchor="middle", spacing=0, font="Helvetica,Arial,sans-serif"):
    sp = f' letter-spacing="{spacing}"' if spacing else ""
    return (f'<text x="{x}" y="{y}" font-family="{font}" font-size="{size}" fill="{col}" '
            f'font-weight="{weight}" text-anchor="{anchor}"{sp}>{esc(s)}</text>')

def rrect(x, y, w, h, fill, rx=14, stroke="none", sw=0, shadow=False):
    st = f' stroke="{stroke}" stroke-width="{sw}"' if stroke != "none" else ""
    fl = ' filter="url(#sh)"' if shadow else ""
    return f'<rect x="{x}" y="{y}" width="{w}" height="{h}" rx="{rx}" fill="{fill}"{st}{fl}/>'

def arrow(cx, y0, y1, col):
    w, hw, hh = 12, 30, 16
    p = [(cx-w/2, y0), (cx+w/2, y0), (cx+w/2, y1-hh), (cx+hw/2, y1-hh),
         (cx, y1), (cx-hw/2, y1-hh), (cx-w/2, y1-hh)]
    return f'<polygon points="{" ".join(f"{a:.1f},{b:.1f}" for a,b in p)}" fill="{col}" opacity="0.9"/>'

def chip(x, y, s, size=11.5, fill="#EEF2F5", col=INK, pad=11):
    w = len(s) * size * 0.56 + pad * 2
    return (rrect(x, y, w, size + 12, fill, rx=(size+12)/2) +
            txt(x + w/2, y + size + 1, s, size, col, "600"), w)

el = ['<?xml version="1.0" encoding="UTF-8"?>',
      f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H}" viewBox="0 0 {W} {H}">',
      '<defs><filter id="sh" x="-20%" y="-20%" width="140%" height="140%">'
      '<feDropShadow dx="0" dy="2" stdDeviation="5" flood-color="#1B2A36" flood-opacity="0.10"/></filter>'
      f'<linearGradient id="badge" x1="0" y1="0" x2="1" y2="1">'
      f'<stop offset="0" stop-color="{ACC1}"/><stop offset="1" stop-color="{ACC3}"/></linearGradient></defs>',
      rrect(0, 0, W, H, BG, rx=0)]

# ---------- title ----------
el += [txt(W/2, 40, "HOW THE CENSUS IS BUILT", 13, ACC1, "bold", spacing=3),
       txt(W/2, 62, "inputs feed an evidence-tiered resource, queried through three analysis arms", 12.5, SUB)]

def card_header(x, y, w, acc, tint, num, title):
    e = [rrect(x, y, w, 300, CARD, rx=16, stroke=BORDER, sw=1, shadow=True),
         rrect(x, y, 5, 300, acc, rx=2.5),                    # left accent stripe
         f'<circle cx="{x+30}" cy="{y+30}" r="14" fill="{tint}"/>',
         txt(x+30, y+35, num, 15, acc, "bold"),
         txt(x+52, y+35, title, 15.5, INK, "bold", anchor="start")]
    return e

def stat_tile(x, y, w, h, tint, acc, num, label):
    e = [rrect(x, y, w, h, tint, rx=11), txt(x+16, y+38, num, 27, acc, "bold", anchor="start")]
    words = label.split(); lines, cur = [], ""
    for wd in words:
        if len(cur+" "+wd) > 20: lines.append(cur); cur = wd
        else: cur = (cur+" "+wd).strip()
    lines.append(cur)
    for j, ln in enumerate(lines):
        e.append(txt(x+16, y+58+j*13, ln, 10.5, SUB, anchor="start"))
    return e

def cat_tile(x, y, w, h, tint, acc, head, src):
    e = [rrect(x, y, w, h, tint, rx=11),
         f'<rect x="{x}" y="{y+14}" width="4" height="{h-28}" rx="2" fill="{acc}"/>',
         txt(x+16, y+30, head, 14, acc, "bold", anchor="start")]
    words = src.split(", "); lines, cur = [], ""
    for wd in words:
        if len(cur+", "+wd) > 26: lines.append(cur); cur = wd
        else: cur = (cur+", "+wd).strip(", ").strip()
    lines.append(cur)
    for j, ln in enumerate(lines):
        e.append(txt(x+16, y+50+j*14, ln, 10.5, SUB, anchor="start"))
    return e

# ---------- TIER 1 : inputs ----------
lx, lw = 30, 375
el += card_header(lx, 82, lw, ACC1, TINT1, "1", "Tubulin seeds & interaction screens")
tw, th = 165, 92
tx = [lx+22, lx+22+tw+12]; ty = [130, 130+th+12]
blue = [("34", "Seeds: 9 alpha, 10 beta, 15 worm"), ("19", "Human isotypes"),
        ("7", "Interaction databases + HIPPIE"), ("27,622", "Interaction records")]
for i, (n, l) in enumerate(blue):
    el += stat_tile(tx[i%2], ty[i//2], tw, th, TINT1, ACC1, n, l)

rx0, rw = 455, 375
el += card_header(rx0, 82, rw, ACC2, TINT2, "2", "Orthogonal data resources")
rtx = [rx0+22, rx0+22+tw+12]
red = [("Annotation", "UniProt, InterPro, GO"), ("Pathways", "KEGG, Reactome"),
       ("Disease & drug", "OMIM, Open Targets, ClinVar, GDSC"), ("Structure & context", "PDB 1JFF, DepMap, GTEx")]
for i, (hd, sr) in enumerate(red):
    el += cat_tile(rtx[i%2], ty[i//2], tw, th, TINT2, ACC2, hd, sr)

el += [arrow(lx+lw/2, 388, 430, ACC1), arrow(rx0+rw/2, 388, 430, ACC2)]

# ---------- TIER 2 : resource ----------
cy = 436
el += [rrect(30, cy, W-60, 150, CARD, rx=18, stroke=BORDER, sw=1, shadow=True),
       rrect(30, cy, W-60, 6, GOLD, rx=3)]
# monogram badge
el += [rrect(58, cy+34, 84, 84, "url(#badge)", rx=18),
       txt(100, cy+74, "TIC", 34, "#FFFFFF", "bold"),
       txt(100, cy+96, "census", 12, "#EAF3F0", "600")]
el += [txt(180, cy+52, "Tubulin Interactome Census", 25, INK, "bold", anchor="start"),
       txt(180, cy+76, "3,187 human partners + 144 C. elegans, each evidence-tiered", 13, SUB, anchor="start")]
# stat pills
px = 180; pillrow = cy+96
for s, fill, col in [("physical + pfam2go + pathway routes", TINT1, ACC1),
                     ("MongoDB resource", TINT3, ACC3),
                     ("109-claim verifier", TINT2, ACC2),
                     ("mapped onto structure", "#F1ECDD", GOLD)]:
    piece, wpill = chip(px, pillrow, s, 11, fill, col)
    el.append(piece); px += wpill + 10

el += [arrow(180, 590, 632, ACC3), arrow(W/2, 590, 632, ACC3), arrow(W-180, 590, 632, ACC3)]

# ---------- TIER 3 : analysis arms ----------
bw = (W-60-2*20)/3
bx = [30, 30+bw+20, 30+2*(bw+20)]
titles = ["Census & evidence tiers", "Domain, pathway & annotation", "Structure, drug & validation"]
bodies = [
    ["Physical interaction route", "6 primary databases + HIPPIE",
     "Tiers: experimental / compiled / predicted", "High-confidence (>=2 DBs): 423 a, 533 b",
     "Per-isotype census (all 19)"],
    ["pfam2go: 34 Pfam families", "Reactome ORA: 895 a, 942 b pathways",
     "GO function & co-annotation", "Cellular localisation & PTM",
     "Cross-species conservation & homologs"],
    ["Evolutionary Trace on PDB 1JFF", "5 drug pockets on structure",
     "Complex Portal benchmark (109)", "MCODE, hu.MAP, PubTator",
     "DepMap co-dependency arm", "Queryable MongoDB database"],
]
bh = 292
for i in range(3):
    x = bx[i]
    el += [rrect(x, 640, bw, bh, CARD, rx=16, stroke=BORDER, sw=1, shadow=True),
           rrect(x, 640, bw, 46, ACC3, rx=16),
           rrect(x, 664, bw, 22, ACC3),
           f'<circle cx="{x+26}" cy="{663}" r="12" fill="#FFFFFF" opacity="0.22"/>',
           txt(x+26, 668, str(i+1), 13, "#FFFFFF", "bold"),
           txt(x+46, 668, titles[i], 12.5, "#FFFFFF", "bold", anchor="start")]
    yy = 712
    for ln in bodies[i]:
        el.append(f'<circle cx="{x+20}" cy="{yy-4}" r="2.6" fill="{ACC3}"/>')
        # wrap long lines
        words = ln.split(); wl, cur = [], ""
        for wd in words:
            if len(cur+" "+wd) > 30: wl.append(cur); cur = wd
            else: cur = (cur+" "+wd).strip()
        wl.append(cur)
        for k, seg in enumerate(wl):
            el.append(txt(x+32, yy + k*15, seg, 11, INK, anchor="start"))
        yy += 15*len(wl) + 9

el.append("</svg>")
svg = "\n".join(el)
open(os.path.join(OUT, "methodology_schematic.svg"), "w").write(svg)
cairosvg.svg2png(bytestring=svg.encode(),
                 write_to=os.path.join(OUT, "methodology_schematic.png"), output_width=1720)
print("wrote methodology_schematic .svg/.png")

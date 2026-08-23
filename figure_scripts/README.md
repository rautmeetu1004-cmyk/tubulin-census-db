# Figure generator scripts

The Python scripts that render the manuscript figures. Each reads from the
project data files (or the `raw_data/` mirror) and writes PNG/PDF. Run with
Python 3 + matplotlib (some also need logomaker, cairosvg, numpy, scipy,
scikit-learn, pandas, networkx; render_pockets_3d.py needs PyMOL).

| Script | Manuscript figure |
|---|---|
| `make_methodology_schematic.py` | Fig 1 — methodology workflow schematic |
| `make_isotype_enrichment_heatmap.py` | Fig 6 — pathway/GO enrichment across 19 isotypes |
| `make_cancer_vs_disorder_figure.py` | Fig 8 — cancer vs constitutional disorders |
| `make_homolog_conservation_figure.py` | Fig 10 — cross-species conservation (+ violin) |
| `thematic_figures.py` (`fig_t1`) | Fig 11 — evolutionary trace / structural role |
| `make_pubtator_figure.py` | Fig 15 — PubTator literature specificity (+ scatter) |
| `make_ppi_networks.py` | Figs 16–18 — alpha / beta / worm high-confidence maps |
| `13_motif_link_diagram.py` | Fig 23 — cross-set motif links (M7) |
| `make_motif_logos.py` | Fig 24 — MEME sequence logos (M8) |
| `render_pockets_3d.py` | Fig 27 — five drug pockets on structure |
| `make_pocket_figure.py` | Fig 28 — drug-pocket durability |
| `resistance_overlay.py` | Fig 29 — resistance mutations per pocket |
| `make_isotype_pca.py` | isotype PCA + PERMANOVA (supplementary/inferential) |

The remaining manuscript figures (three-route Venn, Pfam treemap, domain-PPI
recovery, pathway dumbbell, disease enrichment, cross-species core, correlation
heatmap, MCODE, Complex Portal benchmark, annotation coverage, GO chord,
localisation, PTM, DepMap, GTEx, graphical abstract) are produced by their
own arm scripts in the full study pipeline.

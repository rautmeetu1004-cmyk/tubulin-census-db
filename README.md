# Tubulin Interactome Census — standalone database

A queryable, isotype-resolved, evidence-tiered census of the human and
*Caenorhabditis elegans* tubulin interactome, packaged as a restorable MongoDB
database with a Python query layer and the exact raw data that formed it. **No
web interface** — this is the database and its inputs, for programmatic use.

Generated 2026-08-23 from the tubulin interactome study.

## Contents

```
database/tubulin_census.archive.gz   the whole database (11.5 MB, gzip mongodump archive)
scripts/tubulin_census.py            the Python query layer (TubulinCensus)
scripts/build_database.py            loader: raw census -> 7 core collections
scripts/extend_database.py           adds corroboration + 4 arm collections
scripts/add_novelty_layers.py        adds TIS + isotype-specificity + tubulin-code
raw_data/                            the 31 source files that formed the DB (37.5 MB)
PROVENANCE.md                        every source file with its sha256 checksum
```

## The database (387,196 documents across 11 collections)

| Collection | Documents | Contents |
|---|---|---|
| `meme_motifs` | 38 | MEME sequence motifs per set (pooled + per-family + worm), with family assignment, sites and E-value |
| `_meta` | 3 | load + extension provenance with per-file sha256 checksums |
| `cross_species_conservation` | 95 | ten-species OrthoDB ortholog presence matrix |
| `disease_associations` | 142,783 | OpenTargets + DISEASES/Jensen disease links (scores flagged as non-comparable) |
| `drug_pockets` | 5 | the five clinical drug pockets with durability (ET, isotype entropy, burial) |
| `drugs` | 4,590 | drugs targeting partner genes (flagged: target the partner, not tubulin) |
| `interactions` | 9,743 | one document per (species, seed, partner) edge; raw records nested |
| `literature_specificity` | 1,141 | PubTator3 tubulin co-mention specificity per partner |
| `network_edges` | 192,089 | high-confidence topology edges (alpha, beta) |
| `orthogonal_coessentiality` | 14,958 | DepMap CRISPR co-dependency (signal-bearing isotypes only, by design) |
| `pathway_enrichment` | 18,424 | Reactome/KEGG over-representation per family and per isotype, with the denominators the test used |
| `proteins` | 3,331 | one document per tubulin partner (all TAPs) with evidence tier, corroboration layer, annotation summary, Pfam domains, and the three integrative layers (Tubulin-Interaction Score, isotype-specificity, tubulin-code role) |
| `tubulins` | 34 | one document per seed isotype (19 human + 15 C. elegans) with structure/ET/motifs on family representatives |

## Restore and query

Requires MongoDB and `pymongo` (`pip install pymongo`).

```bash
# 1. restore the database into a local mongod
mongorestore --gzip --archive=database/tubulin_census.archive.gz

# 2. query it
python3 -c "import sys; sys.path.insert(0,'scripts'); from tubulin_census import TubulinCensus; \
            tc=TubulinCensus(); print(tc.get_protein('HDAC6')['tubulin_interaction_score'])"
```

`TubulinCensus` (scripts/tubulin_census.py) exposes `get_protein`, `profile`,
`partners_of`, `isotypes_of`, `shared_alpha_beta`, `high_confidence`,
`with_domain`, `conserved_to`, `disease_partners`, `druggable_partners`,
`top_by_tis`, `isotype_specific`, `tubulin_code`, `pocket_durability`,
`pathways`, and a composable `advanced()` filter. Set `MONGO_URI` to point at a
remote MongoDB.

## How it formed

The census is seeded on **all nineteen human tubulin isotypes** (nine alpha, ten
beta) and fifteen *C. elegans* orthologs, and assembled from **three orthogonal
routes**: physical interaction (BioGRID, IntAct, BioPlex, HINT, HuRI, STRING,
with HIPPIE as an integrative layer), pfam2go domain content, and KEGG/Reactome
pathway membership. Every partner carries an **evidence tier** (experimental /
compiled / predicted; STRING excluded from high-confidence by construction).

The pipeline that built this database, in order:

1. **`build_database.py`** loads the parsed census tables (`raw_data/01_Census/…`),
   the annotation (`raw_data/02_Annotation/…`), structure/ET
   (`raw_data/03_Structure_Evolution/…`), topology (`raw_data/04_Networks/…`),
   pathways (`raw_data/05_Pathway_Analysis/…`) and disease/drug
   (`raw_data/06_Disease_Drug/…`) into the seven core collections, deriving the
   evidence tier with identical semantics across species.
2. **`extend_database.py`** enriches every partner with its corroboration layer
   (external-evidence score, hu.MAP co-complex, GO microtubule terms, Pfam
   domains, literature specificity) and adds the drug-pocket, cross-species
   conservation, DepMap co-essentiality and literature-specificity collections.
3. **`add_novelty_layers.py`** computes the composite **Tubulin-Interaction
   Score**, the **isotype-specificity index**, and the **tubulin-code** writer/
   eraser/reader classification, and stores them on each partner.

Nothing is recomputed at load time; the loaders are loaders, and every published
number is traceable to its source file (see PROVENANCE.md, and the `_meta`
collection which stamps each file's sha256). Predicted data is loaded and
labelled, never dropped.

## Reproduce from raw_data

The build scripts read paths relative to a project root. To rebuild against the
bundled inputs instead of restoring the archive:

```bash
export TUBULIN_ROOT="$(pwd)/raw_data"
python3 scripts/build_database.py --drop
python3 scripts/extend_database.py
python3 scripts/add_novelty_layers.py
```

(The archive is the authoritative, verified copy; rebuilding is for provenance.)

## License

{LICENSE — to be set by the authors}. Underlying data belong to their source
databases (BioGRID, IntAct, BioPlex, HINT, HuRI, STRING, HIPPIE, UniProt,
InterPro, GO, Reactome, KEGG, OMIM, Open Targets, ClinVar, cBioPortal, GTEx,
DepMap, GDSC, hu.MAP, EBI Complex Portal, PDB, OrthoDB, PubTator3), each under
its own terms.

## Reference verification

The 44 manuscript references were checked against the CrossRef public API
(see `REFERENCES_VERIFICATION.md`): **36 confirmed, 6 likely** (correct paper,
matcher picked a different edition/year), **2 to check by hand**:
- **Ref 1** (Nair et al., *Biomolecules* 2026): the DOI at 16(1):81 resolves to a
  different title on CrossRef — verify the title/DOI.
- **Ref 44** (RBPbase): CrossRef matched the 2014 Gerstberger RNA-binding-protein
  census review, not the RBPbase resource — decide which to cite.

Previously-uncertain references are now resolved with DOIs: **41** Iyer et al.
2023 = *A census of actin-associated proteins in humans* (Front Cell Dev Biol,
`10.3389/fcell.2023.1168050`); **43** ENPD = `10.1093/nar/gky1112`; **42** RBP2GO
= `10.1093/nar/gkaa1040`.

#!/usr/bin/env python3
"""Load the tubulin census into MongoDB.

Why a database at all. The study is 20 arms writing CSVs that other arms read by
path. That works for a pipeline and fails for a reader: "which partners are
alpha-specific, disease-linked, and have a druggable ortholog in worm" is a
four-file join done by hand every time it is asked. This loads the arms into one
queryable store so those questions are one query.

DESIGN RULES, each of which exists because the audit found the failure it prevents.

1. EVERY partner document carries its evidence_tier, and no collection silently
   mixes tiers. The worm census admitted 20 STRING-only predictions as though they
   were assays (AUDIT_2026-07-17 finding 1) precisely because the tier existed in
   the data and nothing was forced to look at it. Here the field is mandatory and
   indexed, so a query that ignores it does so visibly.

2. EVERY collection records its provenance in `_meta`: source file, sha256, row
   count, load time. The READMEs drifted from the code because numbers were copied
   by hand (finding 3). A number in this database can always be traced to the file
   and the file version it came from.

3. NOTHING is recomputed here. This is a loader, not an analysis. If a number in
   Mongo disagrees with the CSV, the loader is wrong, and the checksum tells you
   which file to look at. Statistics stay in the arm that owns them.

4. Predicted/inferred data is loaded and LABELLED, never dropped. Reactome worm
   annotation is 100% orthology-projected (finding 6); MobiDB disorder is predicted;
   ~60% of the GO route is InterPro2GO transfer (finding 8). All of it is queryable,
   all of it carries a flag saying what it is. The user's rule: predicted results are
   fine, we mention less confidence.

Run:  python3 09_Database/scripts/build_database.py [--drop]
Needs a mongod on localhost:27017 (brew services start mongodb/brew/mongodb-community).
"""
import argparse
import csv
import hashlib
import json
import sys
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path

try:
    from pymongo import MongoClient, ASCENDING, DESCENDING
    from pymongo.errors import ServerSelectionTimeoutError
except ImportError:
    sys.exit("pymongo not installed:  pip3 install --break-system-packages pymongo")


def _project_root():
    import os
    env = os.environ.get("TUBULIN_ROOT")
    if env:
        return Path(env)
    here = Path(__file__).resolve()
    for parent in here.parents:
        if (parent / "00_Publication_Package").exists() or (parent / ".git").exists():
            return parent
    return here.parent


ROOT = _project_root()
DB_NAME = "tubulin_census"
CEN = ROOT / "01_Census"
ALPHA = CEN / "Tubulin_Alpha_Interactome"
BETA = CEN / "Tubulin_Beta_Interactome"
WORM = CEN / "Tubulin_Celegans_Interactome"
PATH_AN = ROOT / "05_Pathway_Analysis/analysis"
DD = ROOT / "06_Disease_Drug/Tubulin_Disease_Drug/analysis"
PA = ROOT / "02_Annotation/Tubulin_Partner_Annotation/analysis"
TOPO = ROOT / "04_Networks/Tubulin_Network_Topology"

csv.field_size_limit(10_000_000)

_prov = []


def sha256(p, limit=None):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()[:16]


def rows(path, **kw):
    """Read a CSV, recording provenance. Missing file is an error, not an empty load:
    a collection that silently loads 0 rows because an upstream arm was never run is
    worse than a crash, because it looks like a finding."""
    p = Path(path)
    if not p.exists():
        raise FileNotFoundError(f"{p} missing. Run the arm that produces it first.")
    with open(p, newline="", encoding="utf-8", errors="replace") as f:
        out = list(csv.DictReader(f, **kw))
    _prov.append({"file": str(p.relative_to(ROOT)), "sha256_16": sha256(p),
                  "rows": len(out), "bytes": p.stat().st_size})
    return out


def num(v, cast=float):
    try:
        return cast(v)
    except (TypeError, ValueError):
        return None


# --------------------------------------------------------------------- proteins
def build_proteins():
    """One document per (gene, species). The census is the spine; annotation hangs off it.

    Human and worm live in ONE collection keyed by species rather than two collections,
    because every cross-species question is otherwise a manual union. The tier field
    means the same thing in both, which is the point of the tiering fix.
    """
    docs = {}

    # ---- human: alpha and beta unions
    for fam, d in (("alpha", ALPHA), ("beta", BETA)):
        for r in rows(d / "03_parsed_data/partner_union.csv"):
            g = r["partner_canonical"].strip()
            if not g:
                continue
            key = ("Homo sapiens", g)
            doc = docs.setdefault(key, {
                "_id": f"9606:{g}", "gene": g, "species": "Homo sapiens", "taxon": 9606,
                "uniprot": (r.get("partner_uniprot") or "").split("-")[0] or None,
                "families": [], "seeds": [], "primary_sources": [], "all_resources": [],
                "n_records": 0,
            })
            doc["families"].append(fam)
            doc["seeds"] += [q.strip() for q in (r.get("queries") or "").split(";") if q.strip()]
            doc["primary_sources"] += [s.strip() for s in (r.get("primary_sources") or "").split(";") if s.strip()]
            doc["all_resources"] += [s.strip() for s in (r.get("all_resources") or "").split(";") if s.strip()]
            doc["n_records"] += num(r.get("total_records"), int) or 0

    # The human pipeline has no evidence_tier column: it encodes the same thing as
    # primary_sources vs all_resources. Derive the tier with the SAME semantics as the
    # worm, so the field is comparable across species. This is the join the audit's
    # finding 1 showed was missing, made explicit rather than left to the reader.
    for (sp, g), d in docs.items():
        d["primary_sources"] = sorted(set(d["primary_sources"]))
        d["all_resources"] = sorted(set(d["all_resources"]))
        d["families"] = sorted(set(d["families"]))
        d["seeds"] = sorted(set(d["seeds"]))
        exp = {"BioGRID", "IntAct", "BioPlex", "HuRI"}
        if set(d["primary_sources"]) & exp:
            d["evidence_tier"] = "experimental"
        elif "HINT" in d["all_resources"]:
            d["evidence_tier"] = "compiled"
        else:
            d["evidence_tier"] = "predicted"
        d["shared_alpha_beta"] = len(d["families"]) == 2

    # ---- worm: already tiered at source
    for r in rows(WORM / "03_parsed_data/partner_union.csv"):
        g = r["partner"].strip()
        if not g:
            continue
        docs[("Caenorhabditis elegans", g)] = {
            "_id": f"6239:{g}", "gene": g, "species": "Caenorhabditis elegans", "taxon": 6239,
            "uniprot": None,
            "families": sorted({x for x in (r.get("families") or "").split(";") if x}),
            "seeds": sorted({x for x in (r.get("queries") or "").split(";") if x}),
            "primary_sources": [s for s in (r.get("databases") or "").split(";")
                                if s in ("BioGRID", "IntAct")],
            "all_resources": sorted({x for x in (r.get("databases") or "").split(";") if x}),
            "n_records": num(r.get("total_records"), int) or 0,
            "evidence_tier": r.get("evidence_tier", "experimental"),
            "shared_alpha_beta": len({x for x in (r.get("families") or "").split(";") if x}) == 2,
        }

    # ---- attach human annotation. Each annotation type carries its own provenance
    # flag: whether it is experimental or inferred. Loading a predicted disorder call
    # next to an experimental PTM site without that flag is how finding 8 happened.
    ann_specs = [
        ("localization", PA / "localization.csv", "term", {"inferred": None,
         "note": "source column records the assigning database; evidence codes not resolved here"}),
        ("ptm_sites", PA / "ptm_sites.csv", "feature", {"inferred": None,
         "note": "UniProt feature table; mix of experimental and 'By similarity', see the note column"}),
        ("structure_class", PA / "structure_class.csv", "id", {"inferred": None,
         "note": "CATH/SCOP assignment; homology-modelled entries not separated in the source"}),
        ("disorder_regions", PA / "disorder_regions.csv", None, {"inferred": True,
         "note": "MobiDB consensus disorder is PREDICTED (IUPred and similar), not measured"}),
        ("omim_disease", PA / "omim_disease.csv", "disease", {"inferred": False,
         "note": "OMIM curated gene-disease relationships"}),
    ]
    for name, path, _termcol, meta in ann_specs:
        try:
            rs = rows(path)
        except FileNotFoundError:
            print(f"  ! {name}: source missing, skipped")
            continue
        by_gene = defaultdict(list)
        for r in rs:
            g = (r.get("gene") or "").strip()
            if g:
                rec = {k: v for k, v in r.items() if k not in ("gene", "role", "families") and v}
                rec["_inferred"] = meta["inferred"]
                by_gene[g].append(rec)
        hit = 0
        for g, recs in by_gene.items():
            d = docs.get(("Homo sapiens", g))
            if d is not None:
                d.setdefault("annotation", {})[name] = recs
                hit += 1
        print(f"  {name:16} {len(rs):6,} rows -> {hit:,} proteins   "
              f"[{'PREDICTED' if meta['inferred'] else 'curated' if meta['inferred'] is False else 'mixed'}]")

    return list(docs.values())


# ---------------------------------------------------------------- interactions
BINARY = ("two hybrid", "two-hybrid", "y2h", "cross-link", "crosslink",
          "reconstituted complex", "binary")
COPURIF = ("affinity capture", "co-fractionation", "cofractionation", "pull down",
           "pull-down", "co-ip", "coip", "proximity label", "affinity chromatography",
           "tandem affinity")


def _evidence_class(method, src):
    """binary / co-purification / predicted / unspecified.

    Derived once, centrally, because the heterodimer confound (REVIEWER_AUDIT finding 1)
    turns on it: a co-purification cannot assign a partner to alpha or beta, a binary
    assay can. Leaving every reader to re-derive this from detection_method guarantees
    they derive it differently.
    """
    m = (method or "").lower()
    if src.split("/")[0] == "STRING":
        return "predicted"
    if any(t in m for t in BINARY):
        return "binary"
    if any(t in m for t in COPURIF):
        return "co-purification"
    return "unspecified"


def build_interactions():
    """ONE document per (species, seed, partner) EDGE, with the raw records nested.

    Chosen over one-doc-per-record because the questions asked of this collection are
    edge questions: "does TUBB bind HSPA8", "how many binary partners does alpha have",
    "which partners are supported by more than one method". Those are one lookup on an
    edge document and an aggregation on a record collection. The raw records are not
    lost, they are nested under `records`, so per-record fidelity (which pubmed, which
    method, which score) survives. `best_evidence_class` ranks the edge by its strongest
    evidence (binary > co-purification > unspecified > predicted), so a partner seen once
    by Y2H and ten times by co-fractionation is correctly an edge with binary support.
    """
    RANK = {"binary": 3, "co-purification": 2, "unspecified": 1, "predicted": 0}
    edges = {}
    for sp, taxon, d, fam0 in (("Homo sapiens", 9606, ALPHA, "alpha"),
                               ("Homo sapiens", 9606, BETA, "beta"),
                               ("Caenorhabditis elegans", 6239, WORM, None)):
        for r in rows(d / "03_parsed_data/master_interactions_long.csv"):
            seed = r.get("query_gene")
            partner = r.get("partner_canonical") or r.get("partner_symbol")
            if not (seed and partner):
                continue
            src = r.get("source") or ""
            base = src.split("/")[0]
            fam = r.get("family") or fam0
            ec = _evidence_class(r.get("detection_method"), src)
            key = (sp, seed, partner)
            e = edges.get(key)
            if e is None:
                e = edges[key] = {
                    "_id": f"{taxon}:{seed}:{partner}", "species": sp, "taxon": taxon,
                    "seed": seed, "seed_uniprot": r.get("query_uniprot"),
                    "partner": partner, "family": fam,
                    "n_records": 0, "sources": set(), "evidence_classes": set(),
                    "detection_methods": set(), "pubmeds": set(), "max_score": None,
                    "records": [],
                }
            e["n_records"] += 1
            e["sources"].add(base)
            e["evidence_classes"].add(ec)
            if r.get("detection_method"):
                e["detection_methods"].add(r["detection_method"])
            if r.get("pubmed"):
                e["pubmeds"].add(r["pubmed"])
            s = num(r.get("score_num"))
            if s is not None and (e["max_score"] is None or s > e["max_score"]):
                e["max_score"] = s
            e["records"].append({
                "source": base, "source_detail": src,
                "detection_method": r.get("detection_method") or None,
                "evidence_class": ec, "interaction_type": r.get("interaction_type") or None,
                "pubmed": r.get("pubmed") or None, "score": s,
                "passes_threshold": r.get("passes_threshold") == "1",
            })
    out = []
    for e in edges.values():
        e["sources"] = sorted(e["sources"])
        e["evidence_classes"] = sorted(e["evidence_classes"])
        e["detection_methods"] = sorted(e["detection_methods"])
        e["pubmeds"] = sorted(e["pubmeds"])
        e["n_sources"] = len(e["sources"])
        e["n_publications"] = len(e["pubmeds"])
        e["best_evidence_class"] = max(e["evidence_classes"], key=lambda c: RANK[c])
        # an edge is experimentally supported unless its ONLY evidence is STRING
        e["experimental"] = e["evidence_classes"] != ["predicted"]
        out.append(e)
    return out


# ------------------------------------------------------------------- tubulins
def build_tubulins():
    """One document per tubulin SEED isotype, the first-class entities of the study.

    The tubulins are what this whole census is about, yet in a partner-keyed schema they
    appear only incidentally, as each other's partners. This collection makes them the
    best-described entities: isotype, family, and for the two family representatives
    (TUBA1B alpha, TUBB beta) the full structural/evolutionary layer, per-residue
    Evolutionary Trace, functional-site membership on 1JFF, and MEME motifs.

    The ET/structure arm computes on the representatives, not on all 17 isotypes (the
    isotypes are >97% identical, so one model stands for the family). Non-representative
    isotypes therefore carry a pointer to their representative rather than a copy, which
    is honest: there is no separate ET for TUBB3, and pretending otherwise would invent
    data.
    """
    # From _scripts/tubulin_seeds.py. This was an eleventh copy of the seed list and
    # the last one still on seven alpha isotypes -- it listed symbols rather than
    # accessions, so the accession-based sweep that migrated the other ten missed it.
    # Effect: the database's first-class `tubulins` collection held 17 human seeds
    # while the census tables held 19, and TUBA3E/TUBAL3 had zero interactions in
    # MongoDB despite carrying 133 and 208 rows respectively in the census.
    import sys as _sys
    from pathlib import Path as _P
    _sys.path.insert(0, str(_P(__file__).resolve().parents[2] / "_scripts"))
    from tubulin_seeds import symbols as _syms
    ALPHA_SEEDS = sorted(_syms("alpha"))
    BETA_SEEDS = sorted(_syms("beta"))
    WORM_ALPHA = ["tba-1", "tba-2", "tba-4", "tba-5", "tba-6", "tba-7", "tba-8", "tba-9", "mec-12"]
    WORM_BETA = ["tbb-1", "tbb-2", "tbb-4", "tbb-6", "mec-7", "ben-1"]
    REP = {"alpha": "TUBA1B", "beta": "TUBB"}   # family representative used by the ET arm

    # ---- per-residue ET + functional sites, split by chain (A=alpha, B=beta)
    struct_by_fam = defaultdict(list)
    site_counts = defaultdict(lambda: {"ligand_contact": 0, "dimer_interface": 0, "buried": 0})
    try:
        for r in rows(ROOT / "03_Structure_Evolution/Tubulin_Structural_Evolution/analysis/structure_et_map.csv"):
            fam = r["family"]
            rec = {"resno": num(r.get("pdb_resno"), int), "aa": r.get("aa"),
                   "human_resno": num(r.get("human_resno"), int),
                   "et_rho": num(r.get("et_all_rho")), "et_rank": num(r.get("et_all_rank"), int),
                   "et_coverage_pct": num(r.get("et_all_cov_pct")),
                   "rel_sasa_pct": num(r.get("rel_sasa_pct")),
                   "buried": r.get("buried") == "1",
                   "dimer_interface": r.get("dimer_interface") == "1",
                   "ligand_contacts": r.get("ligand_contacts") or None}
            struct_by_fam[fam].append(rec)
            if rec["ligand_contacts"]:
                site_counts[fam]["ligand_contact"] += 1
            if rec["dimer_interface"]:
                site_counts[fam]["dimer_interface"] += 1
            if rec["buried"]:
                site_counts[fam]["buried"] += 1
    except FileNotFoundError:
        print("  ! structure_et_map.csv missing; tubulins loaded without the structural layer")

    # ---- MEME motifs (whole-tubulin set)
    motifs_by_set = defaultdict(list)
    try:
        for r in rows(ROOT / "03_Structure_Evolution/Tubulin_Structural_Evolution/analysis/meme_motif_map.csv"):
            motifs_by_set[r.get("set")].append({
                "motif": r.get("motif"), "consensus": r.get("consensus"),
                "start": num(r.get("start"), int), "end": num(r.get("end"), int),
                "width": num(r.get("width"), int), "evalue": r.get("evalue"),
                "functional_overlap": r.get("functional_overlap") or None})
    except FileNotFoundError:
        pass

    # partner counts per seed, from the edges (built independently, read from the arm here)
    seed_partners = defaultdict(set)
    for sp, taxon, d, fam0 in (("Homo sapiens", 9606, ALPHA, "alpha"),
                               ("Homo sapiens", 9606, BETA, "beta"),
                               ("Caenorhabditis elegans", 6239, WORM, None)):
        for r in rows(d / "03_parsed_data/master_interactions_long.csv"):
            if r.get("query_gene") and (r.get("partner_canonical") or r.get("partner_symbol")):
                seed_partners[r["query_gene"]].add(r.get("partner_canonical") or r["partner_symbol"])

    docs = []
    for species, taxon, fam, seeds in (
            ("Homo sapiens", 9606, "alpha", ALPHA_SEEDS),
            ("Homo sapiens", 9606, "beta", BETA_SEEDS),
            ("Caenorhabditis elegans", 6239, "alpha", WORM_ALPHA),
            ("Caenorhabditis elegans", 6239, "beta", WORM_BETA)):
        for g in seeds:
            is_rep = (species == "Homo sapiens" and g == REP[fam])
            doc = {"_id": f"{taxon}:{g}", "gene": g, "species": species, "taxon": taxon,
                   "family": fam, "role": "tubulin_seed",
                   "is_family_representative": is_rep,
                   "n_partners": len(seed_partners.get(g, ())),
                   "family_representative": REP[fam] if species == "Homo sapiens" else None}
            if is_rep and struct_by_fam.get(fam):
                doc["structure"] = {
                    "pdb": "1JFF", "chain": "A" if fam == "alpha" else "B",
                    "n_residues_mapped": len(struct_by_fam[fam]),
                    "functional_site_counts": dict(site_counts[fam]),
                    # corrected, non-double-counted enrichment (AUDIT finding 18):
                    "et_enrichment_at_functional_sites": {
                        "test": "Fisher exact, top-15% ET residues at functional sites, "
                                "column-level (independent), AUDIT_2026-07-17 finding 18",
                        "odds_ratio": 2.50, "p_value": 1.2e-3,
                        "note": "the residue-level pooled test (OR 2.76, p 2.2e-5) double-counts "
                                "each alignment column as one alpha + one beta residue; this is the "
                                "corrected independent value"},
                    "residues": struct_by_fam[fam]}
                doc["motifs"] = motifs_by_set.get("all_tubulins", [])
            elif species == "Homo sapiens":
                doc["structure"] = {"reference": f"9606:{REP[fam]}",
                                    "note": "ET and 1JFF mapping are computed on the family "
                                            "representative; the isotypes are >97% identical, so "
                                            "no separate per-isotype trace exists"}
            docs.append(doc)
    return docs
    return out


# ------------------------------------------------------------------- pathways
def build_pathways():
    """ORA results, human and worm, with the parameters the test actually used.

    n and N are attached to every row from ora_params_*.json rather than left in a
    docstring. Finding 3 was exactly this: the reported n (1,742) was not the tested n
    (1,347), and the number lived only in prose. Here a pathway row carries its own
    denominator, so fold and coverage can be recomputed from the document alone.
    """
    docs = []
    specs = [
        ("Homo sapiens", 9606, PATH_AN / "pathway_ora_family.csv",
         PATH_AN / "ora_params_human.json", "family"),
        ("Homo sapiens", 9606, PATH_AN / "pathway_ora_isotype.csv",
         PATH_AN / "ora_params_human.json", "isotype"),
        ("Caenorhabditis elegans", 6239, PATH_AN / "pathway_ora_family_celegans.csv",
         PATH_AN / "ora_params_celegans.json", "family"),
    ]
    for sp, taxon, csvp, parp, level in specs:
        try:
            rs = rows(csvp)
        except FileNotFoundError:
            print(f"  ! {csvp.name} missing, skipped")
            continue
        params = json.loads(Path(parp).read_text()) if Path(parp).exists() else {}
        for r in rs:
            s = r["set"]
            pset = (params.get("sets") or {}).get(s, {})
            k = num(r.get("n_overlap"), int)
            K = num(r.get("n_pathway"), int)
            docs.append({
                "species": sp, "taxon": taxon, "level": level, "set": s,
                "pathway_id": r["pathway_id"], "pathway": r["pathway"],
                "top_level": r.get("top_level"),
                "source": "KEGG" if str(r["pathway_id"]).startswith("path:") else "Reactome",
                "n_pathway": K, "n_overlap": k,
                "expected": num(r.get("expected")), "fold": num(r.get("fold")),
                # coverage is the n-independent quantity; fold is not comparable across
                # species (finding 7). Precomputed so nobody reaches for fold by default.
                "coverage": round(k / K, 4) if (k and K) else None,
                "p_value": num(r.get("p_value")), "q_BH": num(r.get("q_BH")),
                "power_3x": num(r.get("power_3x")),
                "significant": r.get("significant") == "1",
                "genes": [g for g in (r.get("genes") or "").split(";") if g],
                # the denominators the test actually used, carried on every row
                "test_n": pset.get("n_sample"), "test_N": params.get("background_N"),
                "sampling_fraction": (round(pset["n_sample"] / params["background_N"], 4)
                                      if pset.get("n_sample") and params.get("background_N") else None),
                "evidence_tiers_used": params.get("evidence_tiers_used"),
                # Reactome worm annotation is 100% orthology-projected from human
                "annotation_inferred_from_human": (sp == "Caenorhabditis elegans"
                                                   and not str(r["pathway_id"]).startswith("path:")),
            })
    return docs


# ------------------------------------------------------------- disease + drugs
def build_disease():
    docs = []
    try:
        for r in rows(DD / "disease_associations.csv"):
            docs.append({
                "gene": r["gene"], "species": "Homo sapiens",
                "role": r.get("role"), "families": [x for x in (r.get("families") or "").split(";") if x],
                "disease": r.get("disease"), "disease_id": r.get("disease_id"),
                "source": "OpenTargets",
                "association_score": num(r.get("association_score")),
                "evidence_kind": r.get("evidence_kind"),
                "genetic_score": num(r.get("genetic_score")),
                "clinical_score": num(r.get("clinical_score")),
                "literature_score": num(r.get("literature_score")),
                "somatic_score": num(r.get("somatic_score")),
                "therapeutic_areas": [x for x in (r.get("therapeutic_areas") or "").split(";") if x],
                # Open Targets scores are NOT probabilities and are not comparable
                # across diseases or therapeutic areas. Flagged on every document so a
                # reader cannot average them without meeting this field first.
                "_score_semantics": "OpenTargets association score: rank within a target, not a probability; not comparable across therapeutic areas",
            })
    except FileNotFoundError:
        print("  ! OpenTargets disease file missing")
    try:
        for r in rows(DD / "disease_associations_jensen.csv"):
            docs.append({
                "gene": r["gene"], "species": "Homo sapiens",
                "role": r.get("role"), "families": [x for x in (r.get("families") or "").split(";") if x],
                "disease": r.get("disease"), "disease_id": r.get("doid"),
                "source": "DISEASES/Jensen", "channel": r.get("channel"),
                "confidence": num(r.get("confidence")),
                # Jensen channels: 'textmining' is co-mention, NOT curated evidence.
                "_inferred": (r.get("channel") or "").lower() == "textmining",
                "_score_semantics": "DISEASES confidence: channel-specific, textmining channel is literature co-occurrence and is not experimental evidence. Not on the same scale as OpenTargets association_score.",
            })
    except FileNotFoundError:
        print("  ! Jensen disease file missing")
    return docs


def build_drugs():
    try:
        rs = rows(DD / "drug_candidates.csv")
    except FileNotFoundError:
        print("  ! drug file missing")
        return []
    return [{
        "gene": r["gene"], "species": "Homo sapiens",
        "role": r.get("role"), "families": [x for x in (r.get("families") or "").split(";") if x],
        "drug": r.get("drug"), "drug_id": r.get("drug_id"),
        "drug_type": r.get("drug_type"),
        "max_clinical_stage": r.get("max_clinical_stage"),
        "mechanism_of_action": r.get("mechanism_of_action"),
        # The drug targets THIS GENE, which is a census partner. It is not evidence
        # that the drug acts through tubulin. Conflating those is a target-attribution
        # error and the flag exists so a query cannot make it silently.
        "_targets": "this partner gene, not tubulin",
    } for r in rs]


# ------------------------------------------------------------------- topology
def build_network_edges():
    docs = []
    for fam in ("alpha", "beta"):
        p = TOPO / fam / "data/edges.csv"
        if not p.exists():
            print(f"  ! topology {fam} edges missing")
            continue
        for r in rows(p):
            docs.append({
                "family": fam, "species": "Homo sapiens",
                "source_node": r.get("source_node"), "target_node": r.get("target_node"),
                "n_sources": num(r.get("n_sources"), int),
                "sources": [x for x in (r.get("sources") or "").split(";") if x],
                "n_publications": num(r.get("n_publications"), int),
                "biogrid_low_throughput": num(r.get("biogrid_low_throughput"), int),
                "biogrid_high_throughput": num(r.get("biogrid_high_throughput"), int),
                "confidence": num(r.get("confidence")),
            })
    return docs


# ----------------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--drop", action="store_true", help="drop the database first")
    ap.add_argument("--uri", default="mongodb://localhost:27017")
    a = ap.parse_args()

    try:
        cli = MongoClient(a.uri, serverSelectionTimeoutMS=4000)
        cli.admin.command("ping")
    except ServerSelectionTimeoutError:
        sys.exit(f"no mongod at {a.uri}. Start it:\n"
                 "  brew services start mongodb/brew/mongodb-community")

    if a.drop:
        cli.drop_database(DB_NAME)
        print(f"dropped {DB_NAME}")
    db = cli[DB_NAME]

    print("\ntubulins ...")
    tub = build_tubulins()
    print("proteins ...")
    prot = build_proteins()
    print("interactions (edges) ...")
    inter = build_interactions()
    print("pathways ...")
    paths = build_pathways()
    print("disease ...")
    dis = build_disease()
    print("drugs ...")
    drg = build_drugs()
    print("network edges ...")
    edges = build_network_edges()

    loads = [("tubulins", tub), ("proteins", prot), ("interactions", inter),
             ("pathway_enrichment", paths), ("disease_associations", dis),
             ("drugs", drg), ("network_edges", edges)]
    print()
    for name, docs in loads:
        if not docs:
            print(f"  {name:22} SKIPPED (0 documents)")
            continue
        db[name].drop()
        db[name].insert_many(docs, ordered=False)
        print(f"  {name:22} {len(docs):7,} documents")

    # ---- indexes. Chosen for the questions the study actually asks, not exhaustively.
    print("\nindexes ...")
    db.tubulins.create_index([("species", ASCENDING), ("family", ASCENDING)])
    db.tubulins.create_index([("gene", ASCENDING)])
    db.tubulins.create_index([("is_family_representative", ASCENDING)])
    db.proteins.create_index([("species", ASCENDING), ("evidence_tier", ASCENDING)])
    db.proteins.create_index([("gene", ASCENDING), ("species", ASCENDING)], unique=True)
    db.proteins.create_index([("families", ASCENDING)])
    db.proteins.create_index([("shared_alpha_beta", ASCENDING)])
    db.proteins.create_index([("uniprot", ASCENDING)], sparse=True)
    db.interactions.create_index([("partner", ASCENDING), ("species", ASCENDING)])
    db.interactions.create_index([("seed", ASCENDING)])
    db.interactions.create_index([("best_evidence_class", ASCENDING), ("family", ASCENDING)])
    db.interactions.create_index([("evidence_classes", ASCENDING)])
    db.interactions.create_index([("experimental", ASCENDING)])
    db.pathway_enrichment.create_index([("species", ASCENDING), ("set", ASCENDING),
                                        ("significant", ASCENDING)])
    db.pathway_enrichment.create_index([("pathway", ASCENDING)])
    db.pathway_enrichment.create_index([("coverage", DESCENDING)])
    db.disease_associations.create_index([("gene", ASCENDING), ("source", ASCENDING)])
    db.disease_associations.create_index([("disease", ASCENDING)])
    db.drugs.create_index([("gene", ASCENDING)])
    db.drugs.create_index([("max_clinical_stage", ASCENDING)])
    db.network_edges.create_index([("family", ASCENDING), ("source_node", ASCENDING)])
    print(f"  {sum(len(db[c].index_information()) for c, _ in loads if db[c].estimated_document_count())} indexes")

    # ---- provenance. Rule 2: every number is traceable to a file and its checksum.
    db["_meta"].drop()
    db["_meta"].insert_one({
        "_id": "load",
        "loaded_at_utc": datetime.now(timezone.utc).isoformat(),
        "project_root": str(ROOT),
        "collections": {n: len(d) for n, d in loads},
        "sources": _prov,
        "notes": {
            "evidence_tier": "experimental = directly assayed (BioGRID/IntAct/BioPlex/HuRI); "
                             "compiled = HINT only, experimental but via an unparsed source DB; "
                             "predicted = STRING only, a prediction not an observation. "
                             "Default analysis set is experimental+compiled. Predicted data is "
                             "loaded and labelled, never dropped, and carries lower confidence.",
            "pathway_fold": "fold is NOT comparable across species (sampling fractions differ: "
                            "human ~11%, worm ~0.8%). Use `coverage` for cross-species work.",
            "worm_pathways": "Reactome C. elegans annotation is 100% IEA, orthology-projected "
                             "from human. See annotation_inferred_from_human.",
            "disease_scores": "OpenTargets association_score and DISEASES confidence are "
                              "different scales and must not be pooled or compared.",
            "drugs": "drugs target the partner gene, not tubulin.",
        },
    })
    print(f"\n_meta: {len(_prov)} source files recorded with checksums")

    print(f"\n{DB_NAME} ready:")
    for c in sorted(db.list_collection_names()):
        print(f"  {c:24} {db[c].estimated_document_count():>8,}")


if __name__ == "__main__":
    main()

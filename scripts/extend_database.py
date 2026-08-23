#!/usr/bin/env python3
"""Extend the tubulin_census database so it matches the current manuscript.

1. Enrich every TAP in `proteins` with its full corroboration + annotation layer
   (external-evidence score, hu.MAP co-complex, PubTator specificity, GO-MT terms,
   Pfam domains, annotation counts). This is the part the base loader left on the
   census tables.
2. Add collections for the newer arms: drug_pockets, cross_species_conservation,
   orthogonal_coessentiality, literature_specificity. (meme_motifs is loaded by a
   separate step once the arm-13 re-run finishes.)

Idempotent: re-running overwrites the same fields / collections.
Needs mongod on localhost:27017.
"""
import csv, os, sys, hashlib
from datetime import datetime, timezone
from pathlib import Path
from pymongo import MongoClient, ASCENDING, DESCENDING

ROOT = Path(__file__).resolve().parents[2]
DB = MongoClient("mongodb://localhost:27017", serverSelectionTimeoutMS=4000)["tubulin_census"]
csv.field_size_limit(10_000_000)
prov = []

def rows(rel):
    p = ROOT / rel
    if not p.exists():
        raise FileNotFoundError(p)
    with open(p, newline="", encoding="utf-8", errors="replace") as f:
        out = list(csv.DictReader(f))
    h = hashlib.sha256(p.read_bytes()).hexdigest()[:16]
    prov.append({"file": rel, "sha256_16": h, "rows": len(out)})
    return out

def fnum(v):
    try: return float(v)
    except (TypeError, ValueError): return None
def inum(v):
    try: return int(float(v))
    except (TypeError, ValueError): return None

# ---------------------------------------------------------------- 1. enrich TAPs
def enrich_proteins():
    # augmented union tables carry the corroboration fields for EVERY census partner
    # (the external_evidence table only covers the ~960 high-value ones).
    EE = {"alpha": "01_Census/Tubulin_Alpha_Interactome/03_parsed_data/partner_union_augmented.csv",
          "beta":  "01_Census/Tubulin_Beta_Interactome/03_parsed_data/partner_union_augmented.csv"}
    corr = {}
    for fam, rel in EE.items():
        for r in rows(rel):
            g = r["gene"].strip()
            co, tot = fnum(r.get("lit_comention_tubulin")), fnum(r.get("lit_total_mentions"))
            spec = round(co / tot, 4) if (co and tot) else None
            rec = {
                "n_primary_sources": inum(r.get("n_primary_sources")),
                "external_evidence_score": fnum(r.get("external_evidence_score")),
                "humap2_cocomplex_with_tubulin": r.get("humap2_cocomplex_with_tubulin") in ("1", "True", "true"),
                "humap2_complex_ids": [x for x in (r.get("humap2_complex_ids") or "").split(";") if x],
                "humap2_max_confidence": fnum(r.get("humap2_max_confidence")),
                "n_go_mt_terms": inum(r.get("n_go_mt_terms")),
                "go_tubulin_binding": r.get("go_tubulin_binding") in ("1", "True", "true"),
                "go_microtubule_binding": r.get("go_microtubule_binding") in ("1", "True", "true"),
                "lit_comention_tubulin": inum(r.get("lit_comention_tubulin")),
                "lit_total_mentions": inum(r.get("lit_total_mentions")),
                "lit_specificity": spec,
            }
            # keep the row with more literature if a gene appears in both families
            if g not in corr or (rec["lit_total_mentions"] or 0) > (corr[g]["lit_total_mentions"] or 0):
                corr[g] = rec

    # annotation-count summary
    annsum = {}
    for r in rows("02_Annotation/Tubulin_Partner_Annotation/analysis/annotation_master.csv"):
        annsum[r["gene"].strip()] = {
            "length": inum(r.get("length")), "n_omim": inum(r.get("n_omim")),
            "has_disease": r.get("has_disease") in ("1", "True", "true"),
            "n_go_cc": inum(r.get("n_go_cc")), "has_subcell_cc": r.get("has_subcell_cc") in ("1", "True", "true"),
            "n_ptm": inum(r.get("n_ptm")), "n_lipid": inum(r.get("n_lipid")),
            "n_glyco": inum(r.get("n_glyco")), "n_cath": inum(r.get("n_cath")),
            "n_scop": inum(r.get("n_scop")), "n_pfam": inum(r.get("n_pfam"))}

    # pfam domains per protein
    pfam = {}
    for rel in ("02_Annotation/Pfam_GO_Analysis/Alpha/data/pfam_associated_proteins.csv",
                "02_Annotation/Pfam_GO_Analysis/Beta/data/pfam_associated_proteins.csv"):
        try: rs = rows(rel)
        except FileNotFoundError: continue
        for r in rs:
            g = r["gene"].strip()
            doms = [d for d in (r.get("domains") or "").split(";") if d]
            pfam.setdefault(g, set()).update(doms)

    n = 0
    genes = set(corr) | set(annsum) | set(pfam)
    for g in genes:
        upd = {}
        if g in corr: upd["corroboration"] = corr[g]
        if g in annsum: upd["annotation_summary"] = annsum[g]
        if g in pfam: upd["pfam_domains"] = sorted(pfam[g])
        if upd:
            res = DB.proteins.update_one({"_id": f"9606:{g}"}, {"$set": upd})
            n += res.modified_count
    print(f"  enriched {n} protein documents (corroboration + annotation_summary + pfam_domains)")

# ------------------------------------------------------------- 2. new collections
def load_drug_pockets():
    docs = []
    for r in rows("15_Drug_Pocket_Conservation/analysis/pocket_conservation.csv"):
        docs.append({"_id": r["pocket"], "pocket": r["pocket"], "location": r.get("location"),
                     "n_contact_residues": inum(r.get("n_contact")), "n_mapped": inum(r.get("n_mapped")),
                     "mean_ET": fnum(r.get("mean_ET")), "ET_percentile": fnum(r.get("ET_percentile")),
                     "mean_isotype_entropy": fnum(r.get("mean_isotype_entropy")),
                     "frac_buried": fnum(r.get("frac_buried")),
                     "durability_score": fnum(r.get("durability_score")),
                     "_note": "durability combines ET percentile, isotype Shannon entropy across the "
                              "19 isotypes and burial; higher = harder for an isotype switch to remodel"})
    return docs

def load_cross_species():
    sp = ["S.cerevisiae","A.thaliana","Dictyostelium","C.elegans","D.melanogaster",
          "D.rerio","X.tropicalis","G.gallus","M.musculus","H.sapiens"]
    docs = []
    for r in rows("02_Annotation/Taxonomy_Conservation_Analysis/data/taxonomy_conservation.csv"):
        present = {s: (r.get(s) == "1") for s in sp}
        docs.append({"gene": r["gene"], "uniprot": r.get("uniprot_acc"),
                     "group": r.get("group"), "family_module": r.get("family_module"),
                     "n_orthologs_orthodb": inum(r.get("n_orthologs")),
                     "species_present": [s for s in sp if present[s]],
                     "n_species": sum(present.values()),
                     "deepest_clade": r.get("deepest_clade"),
                     "conservation_depth_rank": inum(r.get("conservation_depth_rank")),
                     "_note": "OrthoDB 1:1 ortholog presence/absence across ten species"})
    return docs

def load_coessentiality():
    docs = []
    for r in rows("11_Orthogonal_Evidence/DepMap_Coessentiality/analysis/tap_tubulin_coessentiality.csv"):
        docs.append({"gene": r["gene"], "isotype": r.get("isotype"),
                     "pearson_r": fnum(r.get("r")), "n_lines": inum(r.get("n_lines")),
                     "p_value": fnum(r.get("p_value")), "q_value": fnum(r.get("q_value")),
                     "is_tap": r.get("is_tap") in ("1", "True", "true"),
                     "evidence_tier": r.get("evidence_tier"), "role": r.get("role"),
                     "_note": "DepMap CRISPR co-dependency after >=200 jointly-screened-line floor "
                              "and principal-component correction"})
    return docs

def load_literature():
    docs = []
    for fam, rel in (("alpha","01_Census/Tubulin_Alpha_Interactome/04_analysis/external_evidence/partner_external_evidence.csv"),
                     ("beta","01_Census/Tubulin_Beta_Interactome/04_analysis/external_evidence/partner_external_evidence.csv")):
        for r in rows(rel):
            co, tot = inum(r.get("lit_comention_tubulin_entrez")), inum(r.get("lit_total_mentions_entrez"))
            if not tot: continue
            docs.append({"gene": r["gene"], "family": fam,
                         "co_mentions_tubulin": co, "total_mentions": tot,
                         "specificity": round(co / tot, 4) if co else 0.0,
                         "_source": "PubTator3 (entrez-keyed)",
                         "_note": "specificity = tubulin co-mentions / total mentions"})
    return docs

def load(name, docs, indexes):
    DB[name].drop()
    if docs:
        DB[name].insert_many(docs, ordered=False)
    for idx in indexes:
        DB[name].create_index(idx)
    print(f"  {name:28} {len(docs):6,} documents")

def main():
    print("enriching proteins ...")
    enrich_proteins()
    DB.proteins.create_index([("corroboration.external_evidence_score", DESCENDING)])
    DB.proteins.create_index([("corroboration.humap2_cocomplex_with_tubulin", ASCENDING)])
    DB.proteins.create_index([("corroboration.lit_specificity", DESCENDING)])
    DB.proteins.create_index([("pfam_domains", ASCENDING)])

    print("\nnew collections ...")
    load("drug_pockets", load_drug_pockets(),
         [[("durability_score", DESCENDING)]])
    load("cross_species_conservation", load_cross_species(),
         [[("gene", ASCENDING)], [("n_species", DESCENDING)], [("family_module", ASCENDING)]])
    load("orthogonal_coessentiality", load_coessentiality(),
         [[("gene", ASCENDING)], [("isotype", ASCENDING)], [("pearson_r", DESCENDING)]])
    load("literature_specificity", load_literature(),
         [[("gene", ASCENDING)], [("specificity", DESCENDING)], [("family", ASCENDING)]])

    DB["_meta"].update_one({"_id": "extension"}, {"$set": {
        "_id": "extension",
        "extended_at_utc": datetime.now(timezone.utc).isoformat(),
        "added_collections": ["drug_pockets", "cross_species_conservation",
                              "orthogonal_coessentiality", "literature_specificity"],
        "proteins_enriched_with": ["corroboration", "annotation_summary", "pfam_domains"],
        "meme_motifs": "pending arm-13 re-run",
        "sources": prov,
    }}, upsert=True)
    print("\ncollections now:")
    for c in sorted(DB.list_collection_names()):
        print(f"  {c:28} {DB[c].estimated_document_count():>9,}")

if __name__ == "__main__":
    main()

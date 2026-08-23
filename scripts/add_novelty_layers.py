#!/usr/bin/env python3
"""Three novel layers on the census, stored on every TAP in `proteins`.

1. tubulin_interaction_score (TIS)  -- an RBP2GO-style composite confidence score
   (0-100) that sums independent, orthogonal support for a partner. Missing
   evidence contributes zero, so more independent arms = higher score. The
   per-component breakdown is stored so the number is never a black box.

2. isotype_specificity -- which isotype(s) a partner prefers, from the per-isotype
   interaction records. Turns the isotype-agnostic partner list into a recruitment
   map: preferred isotype, #alpha/#beta isotypes bound, and a 0-1 specificity index
   (1 = binds one isotype only, 0 = binds them all evenly).

3. tubulin_code_role -- writer / eraser / reader of the tubulin post-translational
   code, from curated enzyme lists (writers, erasers) and mark-responsive domains
   (readers). A functional layer specific to tubulin biology.

Idempotent. Needs mongod. Run after extend_database.py.
"""
import math
from collections import defaultdict
from pymongo import MongoClient, DESCENDING

DB = MongoClient("mongodb://localhost:27017", serverSelectionTimeoutMS=4000)["tubulin_census"]

ALPHA = ["TUBA1A","TUBA1B","TUBA1C","TUBA3C","TUBA3D","TUBA3E","TUBA4A","TUBA8","TUBAL3"]
BETA = ["TUBB","TUBB1","TUBB2A","TUBB2B","TUBB3","TUBB4A","TUBB4B","TUBB6","TUBB8","TUBB8B"]
NALL = len(ALPHA) + len(BETA)

# ---- curated tubulin-code enzymes (established literature) --------------------
WRITERS = {
    # polyglutamylases / polyglycylases (TTLL family)
    "TTLL1","TTLL2","TTLL4","TTLL5","TTLL6","TTLL7","TTLL9","TTLL11","TTLL13",
    "TTLL3","TTLL8","TTLL10",
    "TTL",       # tyrosination
    "ATAT1",     # alpha-tubulin acetyltransferase (MEC-17)
    "SETD2",     # tubulin methylation
}
ERASERS = {
    # deglutamylases (cytosolic carboxypeptidases)
    "AGTPBP1","AGBL1","AGBL2","AGBL3","AGBL4","AGBL5","CCP110",
    "VASH1","VASH2","SVBP",      # detyrosination (VASH-SVBP)
    "HDAC6","SIRT2",             # deacetylation
}
# note: CCP110 is centriolar CP110, not a deglutamylase; keep the true CCPs only
ERASERS.discard("CCP110")
READER_DOMAINS = {"CAP-Gly domain", "Doublecortin (DCX)", "EB1-like C-terminal motif",
                  "Kinesin motor domain", "XMAP215_CLASP_TOG", "HOOK", "CAMSAP_CKK"}


def pathway_supported():
    genes = set()
    for p in DB.pathway_enrichment.find({"significant": True}, {"genes": 1}):
        genes.update(p.get("genes", []))
    return genes


def coessential_genes():
    return {d["gene"] for d in DB.orthogonal_coessentiality.find(
        {"q_value": {"$lte": 0.05}}, {"gene": 1})}


def conservation_depth():
    return {d["gene"]: d.get("n_species", 0)
            for d in DB.cross_species_conservation.find({}, {"gene": 1, "n_species": 1})}


def isotype_records():
    """partner -> {seed: n_records} for human edges."""
    rec = defaultdict(dict)
    for e in DB.interactions.find({"taxon": 9606}, {"partner": 1, "seed": 1, "n_records": 1}):
        rec[e["partner"]][e["seed"]] = e.get("n_records", 0)
    return rec


def tis_components(d, pw, coess, cons):
    c = d.get("corroboration", {}) or {}
    tier = d.get("evidence_tier")
    comp = {}
    comp["tier"] = {"experimental": 30, "compiled": 18, "predicted": 6}.get(tier, 0)
    comp["databases"] = min(c.get("n_primary_sources") or 0, 5) * 4          # <=20
    comp["records"] = round(min(math.log1p(d.get("n_records") or 0) / math.log1p(100), 1) * 10, 1)
    comp["humap"] = 10 if c.get("humap2_cocomplex_with_tubulin") else 0
    comp["domain_go"] = 10 if (c.get("go_tubulin_binding") or c.get("go_microtubule_binding")
                               or d.get("pfam_domains")) else 0
    comp["pathway"] = 6 if d["gene"] in pw else 0
    comp["conservation"] = round((cons.get(d["gene"], 0) / 9) * 8, 1)         # <=8
    comp["coessential"] = 8 if d["gene"] in coess else 0
    comp["literature"] = round((c.get("lit_specificity") or 0) * 8, 1)        # <=8
    total = round(min(sum(comp.values()), 100), 1)
    return total, comp


def isotype_specificity(seedrec):
    if not seedrec:
        return None
    tot = sum(seedrec.values()) or 1
    a = [s for s in seedrec if s in ALPHA]
    b = [s for s in seedrec if s in BETA]
    pref = max(seedrec, key=seedrec.get)
    # normalised entropy over isotypes actually bound
    ps = [v / tot for v in seedrec.values() if v > 0]
    H = -sum(p * math.log(p) for p in ps) if len(ps) > 1 else 0.0
    Hmax = math.log(len([v for v in seedrec.values() if v > 0])) or 1
    spec = round(1 - (H / Hmax if Hmax else 0), 3)
    return {"n_alpha_isotypes": len(a), "n_beta_isotypes": len(b),
            "n_isotypes_total": len([v for v in seedrec.values() if v > 0]),
            "preferred_isotype": pref, "preferred_family": "alpha" if pref in ALPHA else "beta",
            "specificity_index": spec,
            "records_by_isotype": dict(sorted(seedrec.items(), key=lambda x: -x[1]))}


def code_role(d):
    g = d["gene"]
    if g in WRITERS:
        return "writer", "adds a tubulin-code mark (glutamylation/glycylation/tyrosination/acetylation/methylation)"
    if g in ERASERS:
        return "eraser", "removes a tubulin-code mark (deglutamylation/detyrosination/deacetylation)"
    doms = set(d.get("pfam_domains") or [])
    if doms & READER_DOMAINS:
        return "reader", f"carries a mark-responsive microtubule-binding domain ({', '.join(sorted(doms & READER_DOMAINS))})"
    return None, None


def main():
    pw, coess, cons, irec = pathway_supported(), coessential_genes(), conservation_depth(), isotype_records()
    n_tis = n_iso = n_role = 0
    roles = defaultdict(int)
    for d in DB.proteins.find({"species": "Homo sapiens"}):
        upd = {}
        tis, comp = tis_components(d, pw, coess, cons)
        upd["tubulin_interaction_score"] = tis
        upd["tis_components"] = comp
        n_tis += 1
        iso = isotype_specificity(irec.get(d["gene"]))
        if iso:
            upd["isotype_specificity"] = iso; n_iso += 1
        role, detail = code_role(d)
        if role:
            upd["tubulin_code_role"] = role
            upd["tubulin_code_detail"] = detail
            roles[role] += 1; n_role += 1
        DB.proteins.update_one({"_id": d["_id"]}, {"$set": upd})

    DB.proteins.create_index([("tubulin_interaction_score", DESCENDING)])
    DB.proteins.create_index([("isotype_specificity.preferred_isotype", 1)])
    DB.proteins.create_index([("isotype_specificity.specificity_index", DESCENDING)])
    DB.proteins.create_index([("tubulin_code_role", 1)])

    print(f"TIS set on {n_tis} proteins; isotype_specificity on {n_iso}; tubulin_code_role on {n_role}")
    print("code roles:", dict(roles))
    print("\nTop 10 partners by Tubulin-Interaction Score:")
    for d in DB.proteins.find({"species": "Homo sapiens"}).sort("tubulin_interaction_score", DESCENDING).limit(10):
        print(f"  {d['gene']:10} TIS={d['tubulin_interaction_score']:5}  tier={d['evidence_tier']:12} "
              f"code={d.get('tubulin_code_role') or '-'}")
    print("\nMost isotype-specific partners (>=3 records, bind one isotype):")
    q = DB.proteins.find({"isotype_specificity.specificity_index": 1.0,
                          "n_records": {"$gte": 3}}).limit(8)
    for d in q:
        s = d["isotype_specificity"]
        print(f"  {d['gene']:10} only {s['preferred_isotype']} ({s['records_by_isotype'][s['preferred_isotype']]} records)")
    print("\nTubulin-code writers/erasers found in the census:")
    for r in ("writer", "eraser"):
        gs = sorted(x["gene"] for x in DB.proteins.find({"tubulin_code_role": r}, {"gene": 1}))
        print(f"  {r}s ({len(gs)}): {', '.join(gs)}")


if __name__ == "__main__":
    main()

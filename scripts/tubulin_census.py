#!/usr/bin/env python3
"""Query API for the tubulin interactome census MongoDB database.

A thin, documented layer over the `tubulin_census` collections so the canonical
questions of the study are one call instead of a hand-written multi-collection
join. Import it, or run it directly to print the worked examples:

    python3 09_Database/tubulin_census.py                # runs example_queries()
    >>> from tubulin_census import TubulinCensus
    >>> tc = TubulinCensus()
    >>> tc.get_protein("HDAC6")["corroboration"]["lit_specificity"]
    >>> tc.partners_of("TUBA1A", tier="experimental")[:5]

Collections: tubulins, proteins (all TAPs + corroboration + annotation_summary +
pfam_domains), interactions (edges), pathway_enrichment, disease_associations,
drugs, network_edges, drug_pockets, cross_species_conservation,
orthogonal_coessentiality, literature_specificity.  meme_motifs is added once the
arm-13 re-run completes.

Every score keeps the semantics the loader stamped: evidence_tier is mandatory,
OpenTargets and DISEASES scores are not comparable, drugs target the partner gene
(not tubulin), pathway `fold` is not comparable across species (use `coverage`).
"""
from pymongo import MongoClient, DESCENDING

HUMAN, WORM = 9606, 6239


class TubulinCensus:
    def __init__(self, uri=None, db="tubulin_census"):
        import os
        uri = uri or os.environ.get("MONGO_URI", "mongodb://localhost:27017")
        self.cli = MongoClient(uri, serverSelectionTimeoutMS=4000)
        self.cli.admin.command("ping")
        self.db = self.cli[db]

    # ---- single entities -------------------------------------------------
    def get_protein(self, gene, species="Homo sapiens"):
        """Full partner record: census tier, corroboration, annotation, pfam."""
        taxon = HUMAN if species.startswith("Homo") else WORM
        return self.db.proteins.find_one({"_id": f"{taxon}:{gene}"})

    def get_tubulin(self, isotype, species="Homo sapiens"):
        taxon = HUMAN if species.startswith("Homo") else WORM
        return self.db.tubulins.find_one({"_id": f"{taxon}:{isotype}"})

    def profile(self, gene):
        """Everything the database knows about one partner, across collections."""
        p = self.get_protein(gene) or {}
        return {
            "protein": p,
            "interactions": list(self.db.interactions.find({"partner": gene})),
            "disease": list(self.db.disease_associations.find({"gene": gene})),
            "drugs": list(self.db.drugs.find({"gene": gene})),
            "conservation": self.db.cross_species_conservation.find_one({"gene": gene}),
            "literature": self.db.literature_specificity.find_one({"gene": gene}),
        }

    # ---- reverse queries -------------------------------------------------
    def partners_of(self, isotype, tier=None, best_evidence=None, species="Homo sapiens"):
        """Every partner of a seed isotype (edge view), optionally filtered by tier."""
        taxon = HUMAN if species.startswith("Homo") else WORM
        q = {"seed": isotype, "taxon": taxon}
        if best_evidence:
            q["best_evidence_class"] = best_evidence
        edges = list(self.db.interactions.find(q))
        if tier:
            genes = {e["partner"] for e in edges}
            keep = {d["gene"] for d in self.db.proteins.find(
                {"gene": {"$in": list(genes)}, "taxon": taxon, "evidence_tier": tier})}
            edges = [e for e in edges if e["partner"] in keep]
        return edges

    def isotypes_of(self, partner, species="Homo sapiens"):
        """Which seed isotypes bind a chosen partner."""
        taxon = HUMAN if species.startswith("Homo") else WORM
        return sorted({e["seed"] for e in self.db.interactions.find(
            {"partner": partner, "taxon": taxon})})

    # ---- set queries -----------------------------------------------------
    def shared_alpha_beta(self, tier=None):
        q = {"shared_alpha_beta": True, "species": "Homo sapiens"}
        if tier:
            q["evidence_tier"] = tier
        return list(self.db.proteins.find(q))

    def by_tier(self, tier, species="Homo sapiens"):
        return list(self.db.proteins.find({"evidence_tier": tier, "species": species}))

    def high_confidence(self, species="Homo sapiens"):
        """Partners reported by >=2 primary databases (STRING excluded by definition)."""
        return list(self.db.proteins.find(
            {"species": species, "corroboration.n_primary_sources": {"$gte": 2}}))

    def with_domain(self, pfam):
        return list(self.db.proteins.find({"pfam_domains": pfam}))

    def conserved_to(self, species_label, group=None):
        """Partners with a detectable ortholog in a given species (OrthoDB)."""
        q = {"species_present": species_label}
        if group:
            q["group"] = group
        return list(self.db.cross_species_conservation.find(q))

    def disease_partners(self, disease=None, source="OpenTargets", min_score=None):
        q = {"source": source}
        if disease:
            q["disease"] = {"$regex": disease, "$options": "i"}
        if min_score is not None:
            q["association_score"] = {"$gte": min_score}
        return list(self.db.disease_associations.find(q))

    def druggable_partners(self, max_stage=None):
        q = {}
        if max_stage:
            q["max_clinical_stage"] = max_stage
        return list(self.db.drugs.find(q))

    def top_literature_specific(self, n=20, min_mentions=40, exclude_seeds=True):
        cur = self.db.literature_specificity.find(
            {"total_mentions": {"$gte": min_mentions}}).sort("specificity", DESCENDING)
        out, seen = [], set()
        for d in cur:
            if d["gene"] in seen:
                continue
            if exclude_seeds and (d["gene"].startswith("TUBA") or d["gene"].startswith("TUBB")):
                continue
            seen.add(d["gene"])
            out.append(d)
            if len(out) >= n:
                break
        return out

    def coessential_with(self, isotype, min_r=0.1):
        return list(self.db.orthogonal_coessentiality.find(
            {"isotype": isotype, "pearson_r": {"$gte": min_r}}).sort("pearson_r", DESCENDING))

    def pocket_durability(self):
        return list(self.db.drug_pockets.find().sort("durability_score", DESCENDING))

    def pathways(self, sample="alpha", significant=True, species="Homo sapiens"):
        q = {"set": sample, "species": species}
        if significant:
            q["significant"] = True
        return list(self.db.pathway_enrichment.find(q).sort("coverage", DESCENDING))

    # ---- the flexible "advanced search" ---------------------------------
    def advanced(self, family=None, tier=None, shared=None, has_disease=None,
                 druggable=None, conserved_in=None, min_ext_evidence=None,
                 domain=None):
        """Compose the filters the paper's advanced query needs. Returns gene list."""
        q = {"species": "Homo sapiens"}
        if family:
            q["families"] = family
        if tier:
            q["evidence_tier"] = tier
        if shared is not None:
            q["shared_alpha_beta"] = shared
        if min_ext_evidence is not None:
            q["corroboration.external_evidence_score"] = {"$gte": min_ext_evidence}
        if domain:
            q["pfam_domains"] = domain
        genes = {d["gene"] for d in self.db.proteins.find(q, {"gene": 1})}
        if has_disease:
            dz = {d["gene"] for d in self.db.disease_associations.find({}, {"gene": 1})}
            genes &= dz
        if druggable:
            dr = {d["gene"] for d in self.db.drugs.find({}, {"gene": 1})}
            genes &= dr
        if conserved_in:
            cs = {d["gene"] for d in self.db.cross_species_conservation.find(
                {"species_present": conserved_in}, {"gene": 1})}
            genes &= cs
        return sorted(genes)

    # ---- novelty layers --------------------------------------------------
    def top_by_tis(self, n=20, exclude_seeds=True, tier=None):
        """Partners ranked by the composite Tubulin-Interaction Score (0-100)."""
        q = {"species": "Homo sapiens"}
        if tier:
            q["evidence_tier"] = tier
        cur = self.db.proteins.find(q).sort("tubulin_interaction_score", DESCENDING)
        out = []
        for d in cur:
            if exclude_seeds and (d["gene"].startswith("TUBA") or d["gene"].startswith("TUBB")):
                continue
            out.append(d)
            if len(out) >= n:
                break
        return out

    def tis_of(self, gene):
        d = self.get_protein(gene) or {}
        return {"gene": gene, "tubulin_interaction_score": d.get("tubulin_interaction_score"),
                "components": d.get("tis_components")}

    def isotype_specific(self, isotype=None, min_index=1.0, min_records=3):
        """Partners that concentrate on one isotype (specificity_index high)."""
        q = {"isotype_specificity.specificity_index": {"$gte": min_index},
             "n_records": {"$gte": min_records}}
        if isotype:
            q["isotype_specificity.preferred_isotype"] = isotype
        return list(self.db.proteins.find(q))

    def tubulin_code(self, role=None):
        """Writers / erasers / readers of the tubulin post-translational code."""
        q = {"tubulin_code_role": role} if role else {"tubulin_code_role": {"$exists": True}}
        return list(self.db.proteins.find(q))

    # ---- summary ---------------------------------------------------------
    def stats(self):
        return {c: self.db[c].estimated_document_count()
                for c in sorted(self.db.list_collection_names())}


def example_queries():
    tc = TubulinCensus()
    print("=== database contents ===")
    for c, n in tc.stats().items():
        print(f"  {c:28} {n:>9,}")

    print("\n=== Q1. HDAC6: the shared disease hub, full corroboration ===")
    d = tc.get_protein("HDAC6") or {}
    c = d.get("corroboration", {})
    print(f"  tier={d.get('evidence_tier')} families={d.get('families')} "
          f"ext_evidence={c.get('external_evidence_score')} "
          f"lit_specificity={c.get('lit_specificity')} pfam={d.get('pfam_domains')}")

    print("\n=== Q2. experimental partners of TUBA1A (top 5 by records) ===")
    e = sorted(tc.partners_of("TUBA1A", tier="experimental"),
               key=lambda x: -x.get("n_records", 0))[:5]
    for x in e:
        print(f"  {x['partner']:12} records={x['n_records']} best={x['best_evidence_class']}")

    print("\n=== Q3. which isotypes bind CCT2 ===")
    print("  ", tc.isotypes_of("CCT2"))

    print("\n=== Q4. most tubulin-specific partners in the literature (seeds removed) ===")
    for x in tc.top_literature_specific(n=6):
        print(f"  {x['gene']:10} specificity={x['specificity']:.2f} "
              f"({x['co_mentions_tubulin']}/{x['total_mentions']})")

    print("\n=== Q5. drug-pocket durability ranking ===")
    for x in tc.pocket_durability():
        print(f"  {x['pocket']:12} durability={x['durability_score']} "
              f"ET%={x['ET_percentile']} isotype_entropy={x['mean_isotype_entropy']}")

    print("\n=== Q6. advanced: alpha, experimental, disease-linked, druggable, "
          "conserved to C. elegans ===")
    g = tc.advanced(family="alpha", tier="experimental", has_disease=True,
                    druggable=True, conserved_in="C.elegans")
    print(f"  {len(g)} partners; e.g. {g[:12]}")

    print("\n=== Q7. shared alpha/beta high-confidence partners ===")
    hc = [d for d in tc.shared_alpha_beta() if d.get("corroboration", {}).get("n_primary_sources", 0) >= 2]
    print(f"  {len(hc)} shared partners with >=2 primary sources")

    print("\n=== Q8. top partners by Tubulin-Interaction Score (seeds removed) ===")
    for d in tc.top_by_tis(n=6):
        print(f"  {d['gene']:10} TIS={d['tubulin_interaction_score']:5}  "
              f"components={d.get('tis_components')}")

    print("\n=== Q9. partners specific to a tissue-restricted isotype (TUBA8) ===")
    for d in tc.isotype_specific(isotype="TUBA8")[:6]:
        print(f"  {d['gene']:10} {d['isotype_specificity']['records_by_isotype']}")

    print("\n=== Q10. tubulin-code enzymes in the census ===")
    for role in ("writer", "eraser", "reader"):
        gs = sorted(d["gene"] for d in tc.tubulin_code(role))
        print(f"  {role+'s':8} ({len(gs)}): {', '.join(gs[:12])}{' ...' if len(gs) > 12 else ''}")


if __name__ == "__main__":
    example_queries()

#!/usr/bin/env python3
"""Single source of truth for the human tubulin seed sets.

Why this exists. Before 2026-08-10 the alpha seed map was copy-pasted into
thirteen scripts, each in a slightly different shape -- `{uid: symbol}`,
`{symbol: uid}`, `{uid: (symbol, description)}`, a bare `set()`. When the census
moved from seven alpha isotypes to nine (adding TUBA3E and TUBAL3), only
build_master.py and fetch_string.py were updated. The other ten kept the
seven-isotype map and silently under-reported: the census tables were correct
because build_master drives them, but enrichment, methodology, annotation,
diagrams and the Complex Portal benchmark all still described a seven-seed study.
Nothing failed; the numbers were just quietly wrong in a subset of outputs.

Import from here rather than re-declaring. A seed set that lives in one place can
be wrong once; a seed set that lives in thirteen places is wrong in a way no
single review can catch.

Membership follows Nair, Saroj & Kunwar, Biomolecules 2026;16:81, Table 1:
"There are nine alpha- and ten beta-tubulin isotypes in humans."

TUBA4B (Q9H853) is deliberately excluded. It is absent from that table, and
UniProt annotates it "Tubulin-like protein alpha-4B" at 241 aa -- roughly half
the ~450 aa of a genuine alpha-tubulin chain.
"""

# uniprot -> (gene symbol, UniProt protein name)
ALPHA = {
    "Q71U36": ("TUBA1A", "Tubulin alpha-1A chain"),
    "P68363": ("TUBA1B", "Tubulin alpha-1B chain"),
    "Q9BQE3": ("TUBA1C", "Tubulin alpha-1C chain"),
    "P0DPH7": ("TUBA3C", "Tubulin alpha-3C chain"),
    "P0DPH8": ("TUBA3D", "Tubulin alpha-3D chain"),
    "Q6PEY2": ("TUBA3E", "Tubulin alpha-3E chain"),
    "P68366": ("TUBA4A", "Tubulin alpha-4A chain"),
    "Q9NY65": ("TUBA8", "Tubulin alpha-8 chain"),
    "A6NHL2": ("TUBAL3", "Tubulin alpha chain-like 3"),
}

BETA = {
    "P07437": ("TUBB", "Tubulin beta chain"),
    "Q9H4B7": ("TUBB1", "Tubulin beta-1 chain"),
    "Q13885": ("TUBB2A", "Tubulin beta-2A chain"),
    "Q9BVA1": ("TUBB2B", "Tubulin beta-2B chain"),
    "Q13509": ("TUBB3", "Tubulin beta-3 chain"),
    "P04350": ("TUBB4A", "Tubulin beta-4A chain"),
    "P68371": ("TUBB4B", "Tubulin beta-4B chain"),
    "Q9BUF5": ("TUBB6", "Tubulin beta-6 chain"),
    "Q3ZCM7": ("TUBB8", "Tubulin beta-8 chain"),
    "A6NNZ2": ("TUBB8B", "Tubulin beta-8B chain"),
}

# Ensembl gene id -> uniprot, for the HuRI parse (HuRI is Ensembl-keyed).
# Resolved from the Ensembl REST xrefs/symbol endpoint, not hand-written.
# TUBA3C and TUBA3D are indistinguishable at the Ensembl level and both map to
# P0DPH8; TUBA3E carries two current gene IDs.
ALPHA_ENSEMBL = {
    "ENSG00000075886": "P0DPH8", "ENSG00000198033": "P0DPH8",
    "ENSG00000123416": "P68363", "ENSG00000127824": "P68366",
    "ENSG00000167552": "Q71U36", "ENSG00000167553": "Q9BQE3",
    "ENSG00000183785": "Q9NY65",
    "ENSG00000152086": "Q6PEY2", "ENSG00000292120": "Q6PEY2",
    "ENSG00000178462": "A6NHL2",
}

# TUBAL3 is a real seed but a sequence outlier: 70.7% identity to its closest
# C. elegans alpha (mec-12) against 86.7-93.1% for the other eight. Any
# conservation-weighted claim that pools all nine should say so explicitly.
# See 01_Census/Tubulin_AlphaBeta_Comparison/crossspecies/seed_conservation_celegans.py
DIVERGENT_ALPHA = {"TUBAL3"}


def uid_to_symbol(family="alpha"):
    """{uniprot: symbol}"""
    src = ALPHA if family == "alpha" else BETA
    return {u: s for u, (s, _) in src.items()}


def symbol_to_uid(family="alpha"):
    """{symbol: uniprot}"""
    src = ALPHA if family == "alpha" else BETA
    return {s: u for u, (s, _) in src.items()}


def symbols(family="alpha"):
    """set of gene symbols"""
    src = ALPHA if family == "alpha" else BETA
    return {s for s, _ in src.values()}


def uid_to_name(family="alpha"):
    """{uniprot: (symbol, protein name)}"""
    return dict(ALPHA if family == "alpha" else BETA)


def all_symbols():
    """Every human tubulin seed symbol, both families (19 total)."""
    return symbols("alpha") | symbols("beta")


if __name__ == "__main__":
    print(f"alpha seeds ({len(ALPHA)}): {', '.join(sorted(symbols('alpha')))}")
    print(f"beta seeds  ({len(BETA)}): {', '.join(sorted(symbols('beta')))}")
    print(f"total: {len(all_symbols())}")
    print(f"divergent (flag in conservation claims): {', '.join(sorted(DIVERGENT_ALPHA))}")

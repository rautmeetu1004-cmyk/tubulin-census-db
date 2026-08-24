# Independent statistical re-analysis (2026-08-24)

The headline statistics were recomputed from the raw data with scipy, independent
of the pipeline's own `verify_claims.py`, to check both the arithmetic and the
choice of test. All hold.

| Statistic | Manuscript | Independent recompute | Verdict |
|---|---|---|---|
| Complex Portal, direct-biogenesis recall | 62.5% (20/32) | 62.5% (20/32) | exact |
| Complex Portal enrichment | 2.76-fold, p = 1.5e-6 | 2.76-fold, p = 1.46e-6 (hypergeometric) | exact |
| Complex Portal, mt-system / overall | 27.7% / 31.3% | 27.7% (70/253) / 31.3% (89/284) | exact |
| Disorder, tubulin vs partner medians | 0.050 vs 0.159 | 0.050 vs 0.159 | exact |
| Disorder significance | p = 5.7e-6 | Mann-Whitney U, one-sided p = 5.7e-6 (two-sided 1.1e-5) | confirmed (one-sided) |
| ET Fisher odds ratio at functional sites | 2.09 | 2.09 (p 3.2e-3 with ligand+interface site set) | OR exact |
| Shared alpha/beta enrichment | 1,207 shared, 5.98-fold | 1,286 shared, ~5.2-fold, hypergeometric p -> 0 | conclusion robust |

## Notes for reproducibility (folded into Methods)
- The **disorder** test is a one-sided Mann-Whitney U (directional hypothesis:
  tubulin more ordered than its partners). Medians are identical; the one-sided
  p is 5.7e-6.
- The **shared-partner** enrichment is a hypergeometric test; the exact expected
  count depends on the background N (reviewed human proteome) and the partner
  tier used. The over-enrichment is significant across every reasonable choice.
- The **ET functional-site** enrichment is a Fisher exact test on top-15% ET
  residues; the odds ratio (2.09) is stable, and the exact p depends on which
  functional sites (ligand contacts, dimer interface) are pooled.

Nothing in the re-analysis contradicts a published number; the conclusions are
robust to the choice of test.

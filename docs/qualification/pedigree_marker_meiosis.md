# Marker-level pedigree meiosis qualification

## Scope

This is the bounded qualification for the byte-matrix Mendelian inheritance
oracle frozen in `docs/design/pedigree_marker_meiosis.md`.  It uses no external
data, runtime, package dependency, broad benchmark, or large pedigree.

Exact fixtures cover founder-ID alignment, later founders, a grandchild,
paternal and maternal phase meanings, one-marker and zero-length chromosomes,
fixed zero/one/multiple crossovers, crossover/marker equality, multiple
chromosomes, canonical traversal, invalid inputs, output options, operational
batches, external pedigree permutations, chromosome subsets, unrelated-animal
addition, and R global-RNG isolation.  Mendelian checking compares every child
allele independently with the corresponding two parental alleles; the observed
inconsistency count is zero.

## Focused relationship experiment

The relationship fixture contains four unrelated founders and four families of
eight offspring: `(S1,D1)`, `(S1,D2)`, `(S2,D1)`, and `(S2,D2)`.  Thus it has
parent-offspring, full-sib, paternal-half-sib, and maternal-half-sib groups.
There are 256 independent, single-marker chromosomes.  Founder alleles are
fixed Bernoulli-0.5 fixture values.  With `z = genotype - 1`, relationships are
calculated independently of meiosis as `2 * tcrossprod(z) / 256`.

Markers are divided into sixteen fixed blocks of sixteen.  For each relationship
group, the Monte Carlo standard error is the standard deviation of its sixteen
block means divided by four.  The predeclared acceptance interval is the
realized mean plus or minus four standard errors and must contain the pedigree
expectation.

| Group | Mean | Monte Carlo SE | Four-SE interval | Pair SD | Expected |
|---|---:|---:|---:|---:|---:|
| Parent-offspring | 0.481323 | 0.029466 | [0.363461, 0.599186] | 0.053108 | 0.50 |
| Full siblings | 0.479422 | 0.029319 | [0.362145, 0.596700] | 0.051567 | 0.50 |
| Paternal half siblings | 0.236389 | 0.030953 | [0.112578, 0.360201] | 0.046949 | 0.25 |
| Maternal half siblings | 0.253479 | 0.029711 | [0.134636, 0.372322] | 0.053130 | 0.25 |
| Unrelated founders | 0.001302 | 0.028430 | [-0.112418, 0.115022] | 0.065442 | 0.00 |

All expectations lie inside their uncertainty-based intervals.  Nonzero
full-sib and half-sib pair standard deviations demonstrate realized segregation
variation rather than fixed pedigree coefficients.

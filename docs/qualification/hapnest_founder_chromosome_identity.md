# Founder chromosome-identity qualification

## Contract correction

The original native founder core keyed SplitMix64 streams by chromosome block
position.  This qualification covers its intentional replacement with exact
UTF-8 chromosome-label identity.  Backward compatibility with the defective
seeded outputs is not retained.  HAPNEST donor, coalescent-age, segment-length,
endpoint, phase, and mutation-filtering semantics remain unchanged.

The decisive fixture has three distinguishable chromosomes whose supplied
order differs from lexical order, phase-asymmetric donor panels, informative
marker patterns, and explicit variant identifiers.  Comparisons realign by
variant ID and require exact byte equality.  It exercises each chromosome
alone, forward and reverse orders, surrounding chromosomes, removal and
addition of unrelated chromosomes, operational individual batches, changed
labels, fixed and different seeds, optional genotype and segment-audit output,
and R global-RNG isolation.

## Acceptance

Every retained chromosome comparison must report zero H1, H2, genotype, and
segment-audit mismatches after coordinate realignment.  A changed chromosome
label and a changed global seed must each produce a nonzero difference on the
informative fixture.  Distinct labels must expose distinct native FNV-1a keys,
and the fixed known FNV-1a byte fixtures must match exactly.  The phase-specific
H1/H2 donor fixture and all existing scientific primitive tests remain exact.

## Results

All values below were measured from an isolated native installation of the
corrected package.  Each comparison reports `(H1, H2, genotype, audit)`
mismatches after alignment by explicit chromosome and variant identity.

| Comparison | `chrZ` | `01` | `CHR1` |
|---|---:|---:|---:|
| chromosome alone versus forward collection | (0, 0, 0, 0) | (0, 0, 0, 0) | (0, 0, 0, 0) |
| reversed versus forward collection | (0, 0, 0, 0) | (0, 0, 0, 0) | (0, 0, 0, 0) |
| surrounded versus forward collection | (0, 0, 0, 0) | (0, 0, 0, 0) | (0, 0, 0, 0) |

Removing chromosome `01` produced `(0, 0, 0, 0)` mismatches on the retained
chromosomes.  Both operational divisions, 5+7 and 3+4+5 individuals, produced
`(0, 0, 0, 0)` mismatches.  Relabelling `chrZ` changed its output, and changing
the seed changed both phases on the informative fixture.

The focused results were 140/140 founder expectations, 80/80 pedigree-meiosis
expectations, 51/51 existing `gsim()` expectations, and 102/102 existing
pedigree/record expectations.  The isolated native package installation and
`git diff --check` also completed successfully.  The known-answer founder
fixture was intentionally replaced because the former values encoded the
defective block-position stream contract.

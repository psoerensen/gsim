# Packed genotypes to phenotypes

This report retains the 0.14.0 dense-path measurements. From 0.15.0 the
default BED Glist path uses bounded native accumulation; qgg supplies
construction but is no longer needed for accumulation. See the
[current phenotype contract](../design/marker_specific_phenotype.md#bounded-bed-accumulation)
and [focused qualification](glist_bounded_accumulation.md).

## Workflow at the measured revision

No genotype/phenotype adapter is required. The shortest supported route is:

```text
phased VCF or existing HAP/BIM/FAM reference
  -> gsim_reference()
  -> gsim_simulate_founders()                 [reusable phased HAP]
  -> gsim_simulate_pedigree(format = "bed")   [BED/BIM/FAM]
  -> qgg::gprep()                              [canonical Glist]
  -> gsim(Glist = ...)                         [phenotypes]
```

`qgg` remains an optional suggested package: it is needed only to construct and
read the Glist in this route. A gsim BED is SNP-major PLINK 1 data with dosage
equal to H1 + H2 and BIM A1 equal to ALT. `qgg::gprep()` preserves BIM marker
order in `Glist$rsids`, FAM sample order in `Glist$ids`, and computes MAF into
`Glist$maf`. `gsim()` explicitly realigns returned genotype rows and columns to
those identities.

The initial Glist preparation summarizes the BED resource. During phenotype
simulation itself, gsim selects components from metadata first and calls
`qgg::getG()` only for causal marker columns unless optional full summary
statistics are requested. Those selected columns are decoded into an
individual-by-causal-marker R matrix. This is not a chromosome-local packed
phenotype engine, and no such claim is made.

`gprep()` does not automatically produce LD scores. Realistic LD scores should
normally be prepared once with the genotype resource through the established
qgg workflow. A complete, uniquely named, positive `ld_score` vector can instead
be supplied directly to `gsim()`; examples using synthetic scores label them as
illustrative.

## Three separate marker-level quantities

For marker j, `causal_probability = q_j` sets

```text
Pr(Z_j != 1) = q_j.
```

Given that the marker is active, `pi_k / sum(pi_active)` chooses its non-null
mixture component. Conditional on active component k,

```text
Cov(beta_j | Z_j = k) = v_k * w_j * R_g,
w_j = [p_j(1-p_j)]^a * r_j^b * s_j^c.
```

Thus `q_j` controls causal membership, active `pi_k` controls conditional
component membership, and `w_j` controls conditional effect magnitude. Effects
scale by `sqrt(w_j)`. The scalar `annotation_score = s_j` is distinct from the
SBayesRC matrix `A` and coefficients `alpha`, which continue to alter component
probabilities. The general gsim model can express a HAPNEST-like variance
parameterization; it is not a separate HAPNEST architecture and the example
exponents are not universal recommendations.

## Alignment and memory

FAM order is the phenotype sample order and BIM order is the phenotype marker
order. Pedigree founder and descendant IDs remain exact FAM IIDs. Every external
non-scalar marker vector requires complete, unique names and is reordered to
canonical marker IDs. MAF from Glist is aligned through `Glist$rsids`; Glist LD
scores, when available, are aligned through `Glist$rsidsLD`.

The packed founder and pedigree stages remain chromosome-local. The subsequent
Glist phenotype stage materializes only its causal genotype columns (or bounded
chunks of all markers when `compute_sumstats = TRUE`). The bounded qualification
fixture independently decodes the complete small BED only to establish exact
H1 + H2 parity; that diagnostic decode is not part of production simulation.

See `inst/examples/end_to_end_phenotype.R` for the complete local-data example.

## Bounded qualification result

The committed local fixture uses 8 phased reference donors, 32 markers on two
chromosomes, a 9-founder reusable base, and a 24-animal three-generation
pedigree containing full siblings, both half-sibling types, a later founder, and
grandchildren. The focused qualification produced:

- zero HAP H1+H2 versus BED decoding mismatches;
- zero BED versus `qgg::getG()` dosage mismatches;
- zero sample-order, marker-order, or pedigree-ID mismatches;
- zero Mendelian inconsistencies;
- exact default versus explicit-zero-exponent components, effects, genetic
  values, residuals, phenotypes, and final R RNG state;
- exact supplied `q_j` and conditional BayesR active proportions;
- weights equal to the hand calculation for MAF/LD and MAF/LD/annotation
  models, with effects differing from unit weights by exactly `sqrt(w_j)`; and
- simultaneous `q_j` and `w_j` operation with both values correctly recorded
  in the causal-marker table.

The two-phase packed pedigree payload is analytically 512 bytes
(`2 * 32 * 8 * ceiling(24/64)`) apart from HAP headers. The default phenotype
run requested 8 causal columns: its dense numeric biological payload is about
1,536 bytes (`24 * 8 * 8`), compared with 6,144 bytes for the deliberately full
`24 * 32` diagnostic decode. R matrix headers and dimnames add bounded overhead.
The full decode exists only in this small independent parity test; ordinary
`gsim(Glist=...)` did not request it.

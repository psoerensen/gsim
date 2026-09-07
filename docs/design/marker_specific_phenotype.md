# Marker-specific phenotype architecture

## Scope and ownership

`gsim()` owns a general marker-specific truth model for phenotype simulation.
It supports separate marker-level control of causal membership and active-effect
variance without changing genotype storage or requiring qgg beyond optional
Glist access. HAPNEST is provenance for one possible parameterization of
marker-specific variance; callers supply the general probabilities and weights
required by their study.

The model does not include dominance, epistasis, GxE, multiallelic effects, or a new
genotype representation.

## Call flow

`gsim()` first establishes canonical marker and sample IDs. In Glist mode this
uses marker metadata only. It builds probabilities and draws components, computes
selected causal genotype statistics in bounded blocks, draws effects, and
accumulates genetic values natively. Residuals, phenotypes, and optional bounded
summary-statistic scans follow. qgg supplies Glist construction through `gprep`;
native BED accumulation from an existing supported object needs only gsim.

The `marker_multipliers` path represents marker-specific
active-effect variance. The `causal_probability` path is deliberately
independent of it.

## Probability contract

Let marker `j` have component `Z_j`, where component 1 is null and components
`2,...,K` are active. Let `pi` be the architecture's normalized baseline
component probabilities and let `q_j` be the supplied causal probability. With
`causal_probability` supplied,

$$
\Pr(Z_j = 1) = 1 - q_j, \qquad
\Pr(Z_j = k) = q_j \frac{\pi_k}{\sum_{\ell=2}^{K}\pi_\ell},
\quad k = 2,\ldots,K.
$$

Thus `q_j` is a genuine Bernoulli non-null probability, not a relative sampling
weight and not an exact causal count. Each marker retains the active-component
mixture implied by `pi`. Values may include zero and one; they must otherwise be
finite and in `[0,1]`, and at least one marker must have positive mass. A scalar
is repeated. A marker-specific vector requires unique, nonempty names whose set
exactly matches the canonical marker IDs, and is reordered once to that order.

Direct probabilities cannot be combined with `n_causal`: conditioning on an
exact count would destroy their marginal-probability meaning. They are also
rejected for fixed effects, clustered active-odds enrichment, and SBayesRC
annotation-driven component probabilities. These models remain available when
`causal_probability = NULL`.

If a Bernoulli draw happens to contain no causal marker, simulation fails with
the existing clear no-causal-marker error; it does not silently condition,
resample, or change the requested probabilities.

## Variance contract

Let `v_k` be the relative variance of active component `k`, `w_j` be the
positive finite `marker_multipliers` value for marker `j`, and `R_g` be the
trait-effect correlation matrix. Before optional realized-variance calibration,

$$
\operatorname{Cov}(\beta_j \mid Z_j = k) = v_k w_j R_g, \quad k > 1,
\qquad \beta_j = 0 \text{ when } k = 1.
$$

The `w_j` values are not clipped or normalized. They affect neither causal
membership nor the conditional active-component probabilities. Under the
MAF-dependent architecture, the existing normalized
`[2p_j(1-p_j)]^maf_exponent` factor additionally multiplies active-effect
variance. When `scale_effects = TRUE`, the established trait-wise post-draw
calibration rescales effects to realized `vg`; it does not change component
membership.

`gsim()` can either accept `w_j` directly through `marker_multipliers`, or derive
it from marker metadata:

```text
w_j = [p_j(1-p_j)]^a * r_j^b * s_j^c
```

Here `p_j` is MAF, `r_j` is LD score, `s_j` is one scalar annotation score, and
`a`, `b`, and `c` are finite scalar exponents. Computation uses the fixed order

```text
log(w_j) = a*log(p_j(1-p_j)) + b*log(r_j) + c*log(s_j)
w_j = exp(log(w_j)).
```

Only terms with nonzero exponents are resolved or validated. Their bases must
be finite and strictly positive; MAF must additionally be strictly below one.
Underflow, overflow, and any non-positive result are errors and identify the
affected marker IDs. There is no flooring, truncation, winsorization, or
normalization. In particular, LD scores are not floored at 0.0001.

The default exponents are zero, so the derived model is not requested and unit
weights remain the default. Supplying `marker_multipliers` directly takes the
place of the derived model; combining direct multipliers with any nonzero
exponent is rejected rather than multiplied twice. A scalar direct multiplier
is repeated, while every non-scalar external marker vector must have complete,
unique, nonempty names whose set exactly matches the canonical marker IDs.

In Glist mode, `p_j` is obtained from `Glist$maf` aligned through
`Glist$rsids`, and `r_j` from `Glist$ldscores` aligned through
`Glist$rsidsLD`, unless an explicit vector is supplied. Metadata are resolved
before causal selection without reading genotypes. Once components have been
drawn, the complete weight vector is subset to causal markers by the existing
effect path, and only causal genotype columns are requested from `getG()`.
For in-memory or internally simulated genotypes, MAF comes from the existing
genotype-based calculation unless supplied; LD score is never calculated by
gsim and must be supplied when `b` is nonzero.

`annotation_score` is exactly one positive scalar score per marker and affects
only conditional effect variance. It is not inferred from, combined with, or
reinterpreted as the SBayesRC marker-by-annotation matrix `A`; `A` and `alpha`
retain their existing role in component probabilities. The two mechanisms may
be used together because they operate on separate parts of the model.

The formula includes a parameterization used by HAPNEST, but it is implemented
here as a general marker-specific variance model integrated with Glist. It is
not a separate architecture and the illustrated exponents are not universal
defaults. The pre-existing `maf_dependent` architecture retains its additional
normalized `[2p_j(1-p_j)]^maf_exponent` factor, including when combined with
direct or derived marker multipliers.

## RNG and returned truth

`gsim()` continues to use R's RNG. A supplied seed calls `set.seed(seed)` once.
The direct probability surface is deterministic and consumes no random draws;
component drawing retains one uniform draw per marker, followed by the existing
effect and residual draws. The default `causal_probability = NULL` branch is
unchanged and introduces no additional draws, so established default seeded
components, effects, and phenotypes remain unchanged.

The result contains the canonical effective `causal_probability` and
`marker_multipliers` vectors. The causal-marker table records each selected
marker's probability and final variance weight. Compact settings record the
variance-model type, exponents, formula, metadata sources, alignment, weight
range and geometric mean; full metadata vectors are not duplicated there.
Causal settings retain the policy, alignment, range, mean, and expected causal
count `sum(q_j)`. The existing full component-probability surface remains
optional through `return_marker_probabilities`.

## Bounded BED accumulation

For BED-backed Glist input, `gsim()` uses `.gsim_bed_plan()` to align canonical
sample and marker IDs to each physical FAM and BIM, `.gsim_glist_statistics()`
to freeze transformations, then `.gsim_glist_accumulate()` calls the registered
private `C_gsim_bed_accumulate`. Its standalone `BedReader` adapts checked record
I/O and the scalar `xmat_scalar` loop from gbits (see
[provenance](../../inst/COPYRIGHTS)). It reads each selected marker for all traits
in canonical causal order, retaining one packed physical BED record. The existing
packed BED diagnostic reader reuses this checked reader. No sibling shared
library, external ABI discovery, adapter layer, or configured path is required.

For marker j, let x be BIM A1 dosage (00 ? 2, 01 ? missing, 10 ? 1,
11 ? 0). Replace missing values with the mean of observed dosages among the
**selected samples**. With `standardize_W = TRUE`, subtract the imputed column
mean and divide by its sample SD (denominator n ? 1 after imputation). With
FALSE, use imputed raw dosages. All-missing causal columns fail; invariant
columns fail when standardizing, and can contribute unstandardized provided
aggregate genetic variance is positive. Glist allele frequencies are not used
for this transformation. Observed causal MAF is mean raw dosage / 2, folded
about 0.5, ignoring missing values. Existing supplied MAF and Glist MAF used by
the derived variance model retain their precedence; observed values only fill
missing causal entries. Conditional probabilities and effect-variance weights
are unchanged.

Raw effect draws keep the same R RNG calls and ordering. Each trait's genetic
value is the sum of transformed dosage ? raw effect over causal markers.
Calibration uses sqrt(target variance / observed genetic variance), multiplying
both effects and the accumulated genetic values. Returned `B` and `B_causal`
are the calibrated effects; `settings$effect_scale` identifies that rescaling.
This replaces a second dense matrix product and allows ordinary floating-point
sum differences, not changes to the scientific model. Residual covariance is
formed from realized genetic variance, target h2 and re; the existing normal
residual draw and Y = G + E remain unchanged. `rg` describes effect covariance;
it is not an assertion that finite-sample genetic correlations equal the target.

### Supported storage and identity

The default route requires SNP-major BED, six-column BIM/FAM, unique sample IDs
and globally unique marker IDs, and `Glist$bedfiles`, `ids`, and per-file `rsids`.
`bimfiles`/`famfiles` can be supplied or inferred from matching `.bed` stems.
Optional `n`, `mchr`, `a1` and `a2` must agree with the declared catalog and
physical allele orientation. Selected rows and columns may be noncontiguous
and reordered, across multiple files with independently ordered FAM rows.
`rsidsLD` may define eligible markers; physical offsets always come from BIM,
not positions in filtered `rsids` or `rsidsLD`. This repairs the old qgg access
assumption for filtered metadata rather than reproducing its wrong offsets.
Each BED header, exact extent and selected offset is checked before access.

Other automatic storage modes are unsupported; the former default qgg reader
also used `bedfiles`. An explicit `getG_fun` retains a bounded custom route with
the established callback signature. It must return deterministic raw dosages,
honor requested IDs/order and bounded requests, and tolerate a statistics pass
plus an accumulation pass (and optional summary pass). gsim cannot bound memory
allocated inside a caller's callback. No dense fallback is used. HAP phenotype
input and broad gsuite Glist compatibility are not established.

### Complete memory bound

Each statistics or summary request is capped at B = min(chunk_size, 64),
independent of causal count C. Several transient raw, imputed, centered and
summary-work matrices can coexist, each at most n ? B; this is O(nB), not a
single-buffer peak claim. The native multiplication itself holds one packed
record of ceil(physical sample count / 4) bytes and an n ? traits result.
Calibration and residual generation create further n ? traits copies.

Existing returned effects and truth include full m ? traits `B`, C ? traits
`B_causal`, m ? mixture-component probability arrays, optional m ? annotation
arrays, and optional m ? traits summary tables (including transient frame
assembly). Metadata includes per-file sample mappings and BIM/FAM catalogs.
Thus total working storage is O(n ? traits + m ? traits + C ? traits +
m ? components + m ? annotations + metadata + nB), not constant memory or
strictly O(C ? traits) while the existing full-marker truth contract is retained.
`settings$genotype_stream` records the backend, decode capacity and packed
record bytes, not measured peak process memory. Glist `return_genotypes = TRUE`
is explicitly rejected; W and internally simulated paths retain dense returns.

See the [focused qualification](../qualification/glist_bounded_accumulation.md)
for baseline, allocation and resource evidence. Historical phenotype reports
retain their original dense-path measurements.

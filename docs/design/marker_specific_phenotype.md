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
uses marker metadata only. It then builds component probabilities, draws marker
components, and asks `getG()` only for the selected causal columns. Effects,
genetic values, residuals, phenotypes, and optional summary statistics follow.
The full Glist genotype panel is read only in chunks when marginal summary
statistics are explicitly requested. Selected causal columns form a dense R
matrix; this is not a chromosome-local packed phenotype engine. qgg is required
for this Glist route and optional for other gsim workflows.

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

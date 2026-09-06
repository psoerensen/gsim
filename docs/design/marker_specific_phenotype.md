# Marker-specific phenotype architecture

## Scope and ownership

`gsim()` owns a general marker-specific truth model for phenotype simulation.
It supports separate marker-level control of causal membership and active-effect
variance without changing genotype storage or requiring qgg beyond optional
Glist access. HAPNEST is provenance for one possible parameterization of
marker-specific variance; callers supply the general probabilities and weights
required by their study.

This milestone adds no dominance, epistasis, GxE, multiallelic effects, or new
genotype representation. The established architectures and their defaults are
unchanged.

## Existing call flow

`gsim()` first establishes canonical marker and sample IDs. In Glist mode this
uses marker metadata only. It then builds component probabilities, draws marker
components, and asks `getG()` only for the selected causal columns. Effects,
genetic values, residuals, phenotypes, and optional summary statistics follow.
The full Glist genotype panel is read only in chunks when marginal summary
statistics are explicitly requested.

The existing `marker_multipliers` path already represents marker-specific
active-effect variance. The new `causal_probability` path is deliberately
independent of it.

## Probability contract

Let marker `j` have component `Z_j`, where component 1 is null and components
`2,...,K` are active. Let `pi` be the architecture's normalized baseline
component probabilities and let `q_j` be the supplied causal probability. With
`causal_probability` supplied,

```text
Pr(Z_j = 1) = 1 - q_j
Pr(Z_j = k) = q_j * pi_k / sum(pi_2,...,pi_K),  k = 2,...,K.
```

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

```text
Cov(beta_j | Z_j = k) = v_k * w_j * R_g,  k > 1,
beta_j = 0,                                      k = 1.
```

The `w_j` values are not clipped or normalized. They affect neither causal
membership nor the conditional active-component probabilities. Under the
MAF-dependent architecture, the existing normalized
`[2p_j(1-p_j)]^maf_exponent` factor additionally multiplies active-effect
variance. When `scale_effects = TRUE`, the established trait-wise post-draw
calibration rescales effects to realized `vg`; it does not change component
membership.

No `a`, `b`, or `c` convenience parameterization is added. A study can calculate
`w_j` from MAF, annotations, or another fixed model and pass the resulting named
vector directly. This keeps the simulation truth explicit and avoids a second
partly overlapping variance interface.

## RNG and returned truth

`gsim()` continues to use R's RNG. A supplied seed calls `set.seed(seed)` once.
The direct probability surface is deterministic and consumes no random draws;
component drawing retains one uniform draw per marker, followed by the existing
effect and residual draws. The default `causal_probability = NULL` branch is
unchanged and introduces no additional draws, so established default seeded
components, effects, and phenotypes remain unchanged.

The result contains the canonical effective `causal_probability` and
`marker_multipliers` vectors. The causal-marker table records each selected
marker's probability and variance weight. Compact settings record the policy,
alignment, range, mean, and expected causal count `sum(q_j)`. The existing full
component-probability surface remains optional through
`return_marker_probabilities`.

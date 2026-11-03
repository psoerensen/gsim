# Experimental healthcare trajectory simulation

Development helper: [tools/healthcare_history_simulation.R](../../tools/healthcare_history_simulation.R).
This is a source-level benchmark helper, not an exported installed R API or a
clinical simulator calibrated to Danish registries. Existing genomic simulation
functions and RNG behavior are unchanged. No new packages are installed.

## Process

`hc_simulate_histories()` generates adults from age 30 to 85 or death. Every
patient has an independent stream derived from the supplied seed and patient
index. Appending patients preserves earlier trajectories. RNG kind and the
caller's RNG state are restored. Splits, sex, landmark age, observation start,
patient/module frailties and independent dropout clocks are drawn before the
event process, without conditioning the split on outcomes.

The default vocabulary has 128 first-diagnosis codes in eight observable
ontology groups; the last code is the synthetic target cancer. Baseline rates
follow a declining frequency spectrum, with explicit target-rate override and
heterogeneous age/sex effects. Global and group-specific Normal frailties induce
correlated disease histories. A diagnosis increases its group's recent burden,
which decays with a three-year time scale and affects subsequent diagnoses.
Diagnoses occur at most once per code; death is terminal. The process continues
after target diagnosis until death or age 85, retaining a full disease trajectory.

At each event or fixed age-grid boundary, rates are refreshed. Between refreshes
intensities are constant. The simulator draws an exponential waiting time with
rate equal to the total and a cause proportional to the individual rates. A
wait crossing the next age boundary is censored there and redrawn with refreshed
rates. Thus recorded event ages are continuous; rates are piecewise constant
by definition. The one-year grid is configurable, not a claim of exact smooth
continuous-time baseline hazards.

For code j, the base log rate includes

    log(base_j) + age_effect_j*(age-50) + sex_effect_j*sex
    + .4*global + .6*group_factor[module_j] - .26
    + .18*recent_module_burden[module_j].

The -.26 term centers the lognormal frailty multiplier marginally. The target
also depends additively on selected previous diagnoses with an eight-year decay.
The three scenario profiles share that base and exogenous streams:

| Profile | Additional mechanism |
| --- | --- |
| additive | Additive observed-history log-rate effects; correlated latent factors remain |
| interaction | Products of recent diagnosis effects for the target and selected other diseases |
| order | Interaction profile plus distinct A-before-B versus B-before-A effects, including a five-year gap condition |

These are diverse conditional-intensity mechanisms, not guarantees that any
fitted model is correctly specified after marginalizing latent states/future
diagnoses. The generating rate function is explicit in `.hc_rates`; tests check
interaction rate ratios and fixed-rate exponential/categorical sampling.

## Observation and outcomes

Landmarks are ages 50--70. Observation begins 10--25 years earlier, bounded below
by age 30. Default dropout is independent with rate .02/year after the landmark;
optional diagnosis-report delays have a separately seeded observation stream.
Death reporting is immediate. With nonzero delay the endpoint means first
AVAILABLE target diagnosis, not biological onset. Primary generation uses zero
report delay. Known target diagnosis or death at/before the landmark excludes
that patient using pre-landmark registry information, including older records
outside the predictor lookback. No later follow-up criterion defines eligibility.

The history table contains only observed diagnosis records whose occurrence
and availability are both at/before the landmark. Full follow-up is limited to
five years, dropout or first target/death. Unknown horizon labels stay NA.
They are never coded as negative. The target endpoint is incidence; this does
not simulate cancer-specific death or a causal screening effect.

`hc_complete_followup()` provides realized five-year outcomes with dropout
removed for simulation evaluation only. They are not conditional-history risk
probabilities. Use held-out rows for final complete-ascertainment scoring; never
fill censored training/validation labels from this truth file. A real-data
censored score would require separately validated censoring adjustment.

## Matched feature preparation

`hc_prepare_views()` freezes an external code vocabulary and fixed windows
(0--1,1--3,3--5,5+ years), plus a cubic basis with interior knots 1,3,5,10,20 and
boundary 0,40 years. The same historical records feed baseline, window, spline
and window-plus-interaction views. Age, sex and observation-span covariates are
shared. Sparse matrices retain zero structure and row/column bindings.

Each history code's average squared basis energy is normalized to one using
eligible TRAINING histories only. This matches average prior scale across code
representations, not the whole induced function prior. Additional interaction
priors must budget their total contribution rather than blindly reuse a
per-coefficient variance across different feature counts.

Interaction candidates are the top training-supported code pairs, with a
prespecified minimum co-occurrence (20) and maximum catalog size (256). Neither
outcomes, validation/test histories nor the true generating interaction graph
select the catalog. Each pair has co-occurrence and two direction indicators.
An exact tie has neither direction. True pairs/coefficients and latent states
are kept in truth, never in X. The public ontology grouping itself is observable
metadata and must be equally available to competing methods.

## Export and checks

`hc_export_benchmark()` refuses an existing directory and separates:

- model_inputs: eligible baseline patients, as-of histories, sparse matrices,
  feature dictionaries, training-derived preprocessing/pairs and observed
  training/validation trajectories and follow-up;
- evaluation: held-out observed futures/follow-up and cohort audit;
- truth: latent generating state and complete-ascertainment realized outcomes.

Matrix Market row order binds to model_inputs/patients.csv; feature dictionaries
bind every matrix column. Training trajectories are truncated at landmark+5.
Generated matrix/CSV/RDS files and hashes belong in a developer cache, not Git.
Checks cover RNG preservation, append stability, increasing times, death as the
last event, no repeat diagnosis, availability/eligibility, censoring labels,
training-only pair selection/scaling, sparse values, interaction identities and
fixed-rate sampling references. A separate complete-export check validates
labels, matrices, pair catalogs, split scope and hashes.

The gsuite benchmark guide owns the model-comparison protocol and commands.
No model is fitted by this helper. Delphi tokenization/training, native BLR
precision qualification, repeated seeds, clinical calibration, and package/release
qualification remain separate. Surrogate validation stays paused.

## Signal qualification and additive positive control

The separate development helper `tools/healthcare_signal_simulation.R` adds
`hc_simulate_strong_additive()` and `hc_latent_forward_risk()`. Source the existing
history helper first. No installed public API is added and the original seeded
trajectory generator is unchanged.

The additive control uses fully observed code-presence effects, annual
multinomial cancer/death probabilities and exact cumulative risks. Four codes
have log-odds effects 1.6, 1.2, -1.4 and 0.9; other codes are null. It intentionally
omits hidden outcome frailty and intermediate future diagnoses. It is a learning
positive control, not a clinically calibrated disease simulator. Exact risks
and coefficients are truth-only, excluded from model inputs.

The forward integrator conditions on the original simulator's full latent and
pre-landmark disease state, preserves rate-update boundaries and marginalizes
new future disease paths. These are hidden-state risks, not observable-history
oracle risks or a precise attainable AUC ceiling. Monte Carlo precision is
reported. Empty histories and an independent piecewise-exponential competing
risk reference are checked.

From gsuite, run `tools/validation/test_healthcare_signal.R` and the development
workflow `tools/examples/gbayes/healthcare/qualify_simulation_signal.R`.
The scientific contract and measured comparison are recorded in
`../gsuite/benchmarks/healthcare_simulation/signal-qualification.md`.
These optional comparison links do not create a build or installed dependency.
Surrogate validation remains paused.

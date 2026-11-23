# Standalone gsim simulation workflow

`gsim simulate simulation.yml --output data/smoke` uses the same native workflow
as `gsim_simulate(settings, output, base_directory)` in R. Standalone requires
Python 3.10+, PyYAML and the installed gsim_native executable, but no R runtime.
The marker/prepared-Gaussian R wrapper has no Python dependency. Full Gaussian model preparation calls the matching public gsuite command. Existing gsim(), gsim_pedigree(),
gsim_records() and their caller-RNG behavior are unchanged. Defaults and final
validation are in the shared native core; neither adapter implements simulation.

The first version intentionally supports a frozen supplied parent-first FAM
pedigree, one-chromosome BIM map, founder phase rows `ID bitstring1 bitstring2`,
and ordered marker-effect rows `marker effect_trait1 ...`. FAM/BIM use gsim's
existing canonical tab/LF format. Identities use ASCII letters/digits/_.-.
Every nonfounder needs two earlier parents; every founder needs two phases.
This is not a replacement for the richer legacy R pedigree/founder generators.

YAML schema `gsim.simulation/v1` has required inputs (pedigree, variants, founders,
effects) and ordered traits. Optional controls: seed (default 1), genotyped and
phenotyped counts (default all), records_per_individual (1), missing_rate (0), residual_covariance
(identity), maximum_memory_bytes (4294967296, 4 GiB). Covariance rows/columns follow
traits. Paths resolve beside simulation.yml; R paths resolve beside the explicit
base_directory. Unknown/duplicate keys, aliases/tags, nonfinite values, invalid
maps/parentage/dimensions/covariance and insufficient numeric workspace fail.
The workspace allowance is an estimated RAM admission limit, not an output disk-size
limit or a process RSS cap; it does not reserve 4 GiB. Set it explicitly to suit
the job allocation. `visits` remains a compatibility alias for `records_per_individual`;
supplying both is an error, even when equal. Records are recording opportunities
before missingness, not a time/correlation or longitudinal-model specification.
Effective settings report records_per_individual and retain visits for old readers. R exact integers are limited
to 2^53-1; the core seed reserves offsets below that limit.

The workflow reuses the existing identity-keyed SplitMix64 streams, meiosis event
sampling, packed inheritance and PLINK writers. Supplied marker effects and
founder dosage centers define genetic truth. Random genotyping membership uses
the existing selection key; output PLINK rows are parent-first. Gaussian residual
streams use record-unit and trait identities; missingness uses a separate domain.
This is marker-derived truth, not an exact Gaussian H draw or covariance-recovery
calibration. Existing R RNG streams are not silently substituted or changed.

CSV/JSON outputs are directly readable from R/Python, with true genetic values
separated from observed phenotype records. effective_settings.json records
resolved controls, simulation.json records dimensions and labelled FNV-1a content
identities (not cryptographic signatures). Evaluation bundles additionally bind
files and sources with SHA-256. Files are published only from a fresh staging
directory; failed stages remain for diagnosis and are never silently overwritten.

Build/install without R:

```sh
cmake -S . -B build/standalone -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=/desired/prefix
cmake --build build/standalone --parallel 2
cmake --install build/standalone
```

Install the R interface using the existing isolated R CMD INSTALL workflow.
Generated help and the namespace export come from R/simulation-workflow.R via
roxygen2. R/standalone outputs and legacy RNG/meiosis/record results were verified
on a bounded Windows fixture; Linux/HPC qualification is a separate check.


## Covariance-matched Gaussian model fixtures

Schema `gsim.gaussian_model/v1` accepts `model`, named `design_inputs`, optional
`seed` (default 1), `maximum_dimension` (1..1024) and `maximum_memory_bytes`
(default 4 GiB). The gsuite frontend validates the fitting YAML and uses greml's
native prior/residual export; gsim performs the shared Gaussian draws. Fixed
means are zero. Pedigree, marker/ancestry data, covariates and missingness are
fixed inputs. REML draws use the declared starting covariance, not fitted values.
Explicit APY and other model approximations remain part of the generated model.

The sampler is deliberately bounded to 1,024 latent/scalar dimensions and 4,096
reported coefficient components. It is for small model-matched checks, not a
large dense simulation method. Invalid dimensions, lengths, nonfinite arrays,
non-SPD covariances, changed inputs and insufficient workspace fail without
publishing a completed dataset. No jitter or eigenvalue clipping is applied.

Both standalone and `gsim_simulate()` accept an already prepared recipe with
`inputs: {specification: path}` instead of model/design_inputs. Preparation uses
only installed public gsuite/greml/gsolve interfaces. No sibling native headers
or libraries are linked into gsim. Existing marker-derived and caller-RNG
simulation behavior is unchanged.

CSV outputs include coefficient truth, labelled scalar observation truth,
standard-normal innovations and phenotypes. Full model mode additionally renders
observations and labelled prediction/component truth through the existing
model output contract. No genotype files are generated in this conditional mode.
A full catalogue and per-pair fitting qualification ship with the evaluation
bundle. R/standalone parity and a reconstructed Gaussian draw are numerical
checks; they are not reliability calibration, variance recovery or biological
realism qualification.

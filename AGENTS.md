# AGENTS.md

## Purpose

`gsim` provides simulation functionality for the broader gsuite ecosystem.

## Scope

Own deterministic simulation primitives, simulation-facing R APIs, and the
private packed genotype storage and file I/O needed by those workflows. Do not
turn gsim into a general-purpose matrix library or place numerical solvers,
mixed-model estimation, general matrix kernels, scoring, or correlation methods
here. qgg is optional generally and required for the supported Glist workflow;
no sibling native project dependency is required.

## Build and test

Install with `R CMD INSTALL --preclean .`. For a focused check, create an empty
isolated library outside the repository and substitute its absolute path for
`<isolated-library>` in both commands below. Run from the repository root with
testthat already available in the R library search path:

```text
R CMD INSTALL --preclean --library="<isolated-library>" .
Rscript --vanilla -e "lib <- normalizePath('<isolated-library>', mustWork = TRUE); .libPaths(c(lib, .libPaths())); library(gsim, lib.loc = lib); testthat::test_file('tests/testthat/test-gsim.R')"
```

This explicitly loads the newly installed package from the isolated library.
See [docs/README.md](docs/README.md) for documentation and reproduction entry
points. Generated `man/*.Rd` files name their roxygen source in `R/`; update
those sources and regenerate with roxygen2 when such changes are authorized,
never edit generated Rd by hand.

## Rules

Preserve unrelated worktree changes. Keep seeded behavior reproducible and update generated R documentation only through the repository's documented workflow. Benchmarks and broad qualification runs require authorization. Never commit generated build, installation, cache, log, or check artifacts.

## Website

Follow [website/README.md](website/README.md). Build with
`python website/build.py`; render only into ignored `website/_site/`, never
into `docs/`. Scientific Markdown and public Rd remain authoritative. Website
builds must not load gsim, execute examples or benchmarks, install dependencies,
or require private repositories. Deployment stays manual-only until authorized.

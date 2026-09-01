# AGENTS.md

## Purpose

`gsim` provides simulation functionality for the broader gsuite ecosystem.

## Scope

Own deterministic simulation primitives and simulation-facing R APIs. Do not place numerical solvers, mixed-model estimation, genotype storage, matrix kernels, scoring, or correlation methods here. Optional integration may use `qgg`; no native project dependency is required.

## Build and test

Install with `R CMD INSTALL --preclean .`. Run the focused package check with `Rscript --vanilla -e "testthat::test_file('tests/testthat/test-gsim.R')"` in an isolated test library.

## Rules

Preserve unrelated worktree changes. Keep seeded behavior reproducible and update generated R documentation only through the repository's documented workflow. Benchmarks and broad qualification runs require authorization. Never commit generated build, installation, cache, log, or check artifacts.

# Genetic extensions: focused qualification

Observed 4 November 2026, Windows, R 4.4.1, GCC 13.2.0/C++17.

The modified source package compiled and installed successfully into a new
isolated temporary R library. Validation explicitly loaded that installation.
The final R-only refresh reused the unchanged verified native objects.

| Test file | Test cases | Passed assertions |
| --- | ---: | ---: |
| test-genetic-extensions.R | 5 | 51 |
| test-gsim.R | 15 | 142 |
| test-glist-streaming.R | 5 | 394 |
| test-breeding-blocks.R | 4 | 60 |
| test-packed-simulation.R | 9 | 72 |
| test-boa-simulation.R | 2 | 55 |
| test-separated-simulation.R | 1 | 26 |
| Total | 41 | 800 |

No test failures or skipped tests occurred in these files. The local
[genetic extension example](../../inst/examples/genetic_extensions.R) also ran
successfully from the isolated installation. An additional smoke check confirmed
that additive-only simulation leaves unrequested dominance components NULL.
R help was regenerated from roxygen source with local links resolved, and the
README was rendered with its example chunks disabled by its existing setup.

The new checks compare supplied dominance and environmental slopes with direct
matrix formulas, common total-genetic variance calibration, pure dominance,
allele-reversal invariance, BED/dense parity with missing calls, fixed-base
persistent values, and record-specific exposures with unchanged residual draws.
Mutation checks cover exact zero/one endpoints, unchanged founders/crossovers,
replay, incremental/full-pedigree parity, a bounded binomial count check,
multigeneration transmission and unchanged breed-origin labels.

## Reproduction

Create an empty isolated library outside the repository, then run from the
package root (with testthat available):

```text
R CMD INSTALL --preclean --library="<isolated-library>" .
```

```r
lib <- normalizePath("<isolated-library>", mustWork=TRUE)
.libPaths(c(lib, .libPaths()))
library(gsim, lib.loc=lib)
source("tests/testthat/helper-glist-bed.R")
files <- c("test-genetic-extensions.R", "test-gsim.R",
           "test-glist-streaming.R", "test-breeding-blocks.R",
           "test-packed-simulation.R", "test-boa-simulation.R",
           "test-separated-simulation.R")
for (f in files)
  testthat::test_file(file.path("tests/testthat", f),
                      reporter="summary", stop_on_failure=TRUE)
source("inst/examples/genetic_extensions.R")
```

## Limits

The initial run above is focused numerical and regression validation. The full
source-package follow-up below adds package-level checks, but neither run is a
public-release qualification, performance benchmark or biological calibration study.
Startup locale warnings and a testthat build-version warning were present;
they did not produce test warnings or failures. The scientific model and units
are specified in [the genetic extension contract](../design/genetic_extensions.md).

## Full source-package check follow-up

Observed 4 November 2026 on the same Windows R 4.4.1/GCC 13.2.0 toolchain.
A source archive was built with `R CMD build` and checked with unabridged
`R CMD check gsim_0.17.0.tar.gz`. Installation, examples, the complete test suite,
compiled-code inspection, documentation/usage checks and the PDF manual all
passed. The test result was **1,643 passed assertions, zero failures, zero test
warnings and zero skipped tests**.

Final package-check status: **0 errors, 0 warnings, 1 NOTE**. The note reports
that the repository's existing LICENSE file is not mentioned in DESCRIPTION;
DESCRIPTION retains its existing GPL-3 declaration. No licensing change was
made as part of these numerical features.

The first attempt stopped during DESCRIPTION inspection because inherited
`C.UTF-8` locale settings are unsupported by this Windows R installation.
The complete successful retry set `LC_ALL=C`, `LANG=C` and `LC_CTYPE=C` for that
process only. TinyTeX's Windows binary directory and Rtools were placed on the
process PATH so the manual was checked rather than omitted. Suggested packages
were available; no install, test or manual stage was disabled. The installed
testthat version emitted its existing R-build-version startup message, while
the test suite itself reported WARN 0.

The checked archive and raw logs remain outside the repository under the local
temporary `gsim-full-check-c66403d4625a4671aa164940a01f97fe` directory. They are
not source dependencies and are not committed. Code and test sources were
unchanged during the successful check; only this validation record and its
index links were updated afterwards.

# Building the gsim website

The site uses Quarto, following sblr's Flatly theme, MathJax, GitHub syntax
highlighting, static examples, and Pages artifact deployment. Unlike sblr's
in-place docs output, gsim renders into `website/_site/`. Quarto never receives
`docs/` as its project or output directory.

## Local commands

Prerequisites already on PATH: Python 3.12, Rscript (base R 4.4.1), and Quarto
1.7.32, matching the initial CI configuration. No gsim installation, R package
installation, qgg, knitr, pkgdown, or biological data is required. The script
reports missing tools rather than installing them.

From the repository root, build once:

```text
python website/build.py
```

The helper generates Quarto inputs from the public sources, converts Rd through
`tools::parse_Rd()` and `tools::Rd2HTML()`, then invokes `quarto render website`.
It supplies the version from DESCRIPTION explicitly, so an older installed gsim
cannot leak into the reference pages. It never loads the package or executes
Rd examples. All code fences are display-only; Quarto also sets `eval: false`.

To reopen the retained preview without another build:

```text
python -m http.server 8000 --bind 127.0.0.1 --directory website/_site
```

Open <http://127.0.0.1:8000/> and stop the server with Ctrl+C. The published URL
is configured as <https://psoerensen.github.io/gsim/>; local asset links remain
relative. MathJax uses Quarto's normal CDN asset, so mathematical typesetting
requires browser access to that CDN. Site generation itself downloads no data.

For an identified page defect only:

```text
python website/build.py --page contracts/marker_specific_phenotype.qmd
```

`--prepare-only` regenerates page sources without rendering. Generated QMD files,
reference HTML fragments, downloaded-script copies, Quarto state/cache, and HTML
output are ignored. The build never recursively cleans source directories.
On Windows only the render subprocess's LOCALAPPDATA is redirected to
`website/.cache/` for Quarto's Sass cache; on Linux it uses XDG_CACHE_HOME there.
No user settings or global environment variables are changed.

## Navigation and source authority

| Site destination | Authoritative source |
| --- | --- |
| Home | README introduction, installation, and first small deterministic example |
| Getting started | README packed-reference, phenotype handoff, and pedigree-record sections |
| Function reference | All NAMESPACE exports plus the print S3 method, from public man/*.Rd |
| Examples | website/examples.md for running instructions; unchanged inst/examples scripts copied byte-for-byte as downloads |
| Scientific contracts | docs/README.md index and all docs/design Markdown contracts |
| Validation and performance | docs/README.md evidence index and all docs/qualification reports |
| Attribution | Unmodified inst/COPYRIGHTS, also available as a text download |

Do not edit generated QMD or HTML. Edit the public Markdown authority or the
roxygen source and regenerate Rd through the documented package workflow when
R documentation changes are authorized. The website only converts existing Rd.
Links between public documents are mapped to site pages; links to source code,
tests, and development instructions point to the public repository. Historical
workstation checkout paths in reports become portable placeholders on the site;
measurements, revisions, hardware, hashes, uncertainty, and limitations remain.
No sibling repository or private archive is read by the build.

`website/` and `.github/` are excluded by .Rbuildignore. Website tooling is not
listed in DESCRIPTION and is not a runtime dependency. Only `_site/` is uploaded
to Pages; build helpers, caches, and source directories are not uploaded.

## Review and eventual publication

The prepared `.github/workflows/pages.yml` runs only by workflow_dispatch, only
for psoerensen/gsim on main. Ordinary pushes and pull requests do not publish.
The build job has contents: read; only the deployment job has pages: write and
id-token: write, with the github-pages environment and serialized deployment.
CI provisions the pinned website tools but installs no gsim package/dependencies,
runs no simulations or test suites, and downloads no biological inputs.

After Peter reviews the local preview:

1. Review and commit the intended gsim changes, then push main when authorized.
2. In repository Settings > Pages, select GitHub Actions as the build source.
3. Review the github-pages environment protection rules and restrict deployment
   to main as appropriate.
4. Run the manual Publish documentation workflow on main, approve any environment
   gate, and inspect its artifact/deployment result and the configured site URL.

No step above was performed in this local milestone. To enable automatic
publication later, deliberately add a push trigger restricted to main and paths
README.md, DESCRIPTION, NAMESPACE, man/**, docs/**, inst/examples/**,
inst/COPYRIGHTS, website/**, and .github/workflows/pages.yml. Retain manual
workflow_dispatch and the main/repository deployment guards. Do this only after
review and publication authorization.

## Local verification record

The initial site milestone completed one full 27-page render with Quarto 1.7.32
and a focused homepage refresh after adding the build-guide link. An initial
cache-opening failure was resolved by the subprocess cache setting above.
Local Edge inspection covered Home, Getting started, a founder API page, the
phenotype mathematics contract, Examples, and Validation and performance.
MathJax rendered both display equations; code and navigation were readable
without page-level horizontal overflow at 1440 pixels. The six review screenshots
are retained locally in ignored `website/.preview/`; they are not site assets.

All 852 local links/anchors/assets and 27 sitemap URLs passed static checks.
Example downloads and COPYRIGHTS match their sources byte-for-byte. Public R/C++,
Rd, tests, executable examples/benchmarks, qualification reports, licenses, and
version were unchanged. No package build, package test, simulation, benchmark,
or R CMD check ran. The preview server and isolated browser were stopped, and
the disposable browser profile and inspection helper were removed. Deployment
has not been triggered or validated against GitHub, and the site is not claimed
live. The first manual deployment remains subject to the review steps above.

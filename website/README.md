# RStudio development and package website

Open `gsim.Rproj` in RStudio. The Build pane uses devtools for package
installation, checking and roxygen documentation, as in qgg/gact.
Edit roxygen comments in `R/` and use **Build > Document** to regenerate help.
Edit `README.Rmd`, then use its Knit button to update `README.md`.
Install devtools and roxygen2 as developer tools if they are missing.
Package checks retain the scope and prerequisites documented in AGENTS.md.

## Build the website

Website development uses pkgdown, knitr and rmarkdown. Install these developer
tools once if missing; RStudio supplies Pandoc. They are not runtime package
dependencies. From the project root, run in the RStudio R console:

```r
pkgdown::build_site(examples = FALSE, install = FALSE)
```

Examples remain display-only. No package installation, simulation, benchmark or
biological download is performed by this command. Output goes to ignored
`website/_site/`; existing scientific sources in `docs/` are preserved.
The pkgdown configuration is `_pkgdown.yml`. The homepage is `README.md`,
function reference comes from `man/`, and website-only articles live under
`vignettes/articles/`. Their small include chunks read the authoritative
Markdown and adapt links; scientific code is never executed. Edit the original
Markdown, not copies or generated HTML. Static downloads in `pkgdown/assets/files/`
retain the example and notice files; refresh those copies when their authoritative
`inst/` files change. Article wrappers are excluded from R
source packages, so documentation tools are not installation requirements.

For command-line builds, use Rscript with the same expression. RStudio users
do not need a custom build.py or build.R, Python or Quarto.

## GitHub Pages

Every push to main automatically builds and publishes the website through
`.github/workflows/pages.yml`. Edit in RStudio, Knit README.Rmd when it changes,
then commit and push. Local website previews are optional. The workflow prepares
R website tools and Pandoc, renders without installing the package or running
examples, and uploads only `website/_site/`.

In repository Settings > Pages, select GitHub Actions. The **Publish documentation**
workflow also retains **Run workflow** for a manual rebuild on main. Environment
protection rules, if configured in GitHub, still apply to deployment.
The configured URL is <https://psoerensen.github.io/gsim/>. A successful
Actions deployment establishes that the published site is available.

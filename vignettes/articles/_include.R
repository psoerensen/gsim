# Include authoritative Markdown in a website-only pkgdown article.
# The caller's working directory is vignettes/articles (knitr's default).
include_document <- function(source) {
  root <- normalizePath("../..", winslash = "/", mustWork = TRUE)
  package <- read.dcf(file.path(root, "DESCRIPTION"))[1, "Package"]
  repository <- paste0("https://github.com/psoerensen/", package, "/blob/main/")
  documents <- list.files(file.path(root, "docs"), "[.]md$", recursive = TRUE)
  sources <- c(paste0("docs/", documents), "website/README.md", "website/examples.md", "inst/COPYRIGHTS")
  article <- function(path) sub("[.]md$", "", gsub("[/_]", "-", path))
  mapping <- as.list(setNames(paste0(article(sources), ".html"), sources))
  mapping[["README.md"]] <- "../index.html"
  help <- list.files(file.path(root, "man"), "[.]Rd$")
  for (name in help) mapping[[paste0("man/", name)]] <- paste0("../reference/", sub("[.]Rd$", ".html", name))
  text <- paste(readLines(file.path(root, source), warn = FALSE, encoding = "UTF-8"), collapse = "\n")
  # The article metadata already supplies the document heading.
  if (grepl("^# ", text)) text <- sub("^[^\n]*\n", "", text)
  if (source == "inst/COPYRIGHTS") text <- paste0("```text\n", text, "\n```")
  # Display code, including any executable fences, without evaluating it.
  text <- gsub("```\\{r[^}]*\\}", "```r", text)
  text <- gsub("C:/Users/[^`\\s]+/synthetic_data", "<hapnest-source>", text, perl = TRUE)
  pattern <- "\\[([^\\]]*)\\]\\(([^\\s)]+)\\)"
  hits <- gregexpr(pattern, text, perl = TRUE)
  matches <- regmatches(text, hits)[[1]]
  replacements <- vapply(matches, function(link) {
    parts <- regmatches(link, regexec(pattern, link, perl = TRUE))[[1]]
    destination <- parts[3]
    if (grepl("^[A-Za-z][A-Za-z0-9+.-]*:|^//|^#", destination)) return(link)
    anchor <- if (grepl("#", destination, fixed = TRUE)) sub("^[^#]*", "", destination) else ""
    path <- sub("#.*$", "", destination)
    resolved <- normalizePath(file.path(root, dirname(source), path), winslash = "/", mustWork = FALSE)
    relative <- substring(resolved, nchar(root) + 2L)
    if (!startsWith(resolved, paste0(root, "/"))) stop("Link escapes repository: ", destination)
    mapped <- mapping[[relative]]
    if (is.null(mapped)) mapped <- paste0(repository, relative)
    paste0("[", parts[2], "](", mapped, anchor, ")")
  }, character(1))
  regmatches(text, hits) <- list(replacements)
  if (startsWith(source, "docs/qualification/")) {
    text <- paste0("Historical evidence: these measurements retain their recorded revisions, fixtures and limitations; website rendering does not rerun scientific qualification.\n\n", text)
  }
  paste0(text, "\n\n[Authoritative source](", repository, source, ")\n")
}

# Convert installed-help sources without installing/loading gsim or running R code.
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 2L)
root <- normalizePath(args[[1L]], winslash = "/", mustWork = TRUE)
destination <- normalizePath(args[[2L]], winslash = "/", mustWork = TRUE)
rd_files <- list.files(file.path(root, "man"), "[.]Rd$", full.names = TRUE)
topics <- sub("[.]Rd$", "", basename(rd_files))
links <- stats::setNames(paste0(topics, ".html"), topics)
links["set.seed"] <- "https://stat.ethz.ch/R-manual/R-devel/library/base/html/Random.html"
version <- read.dcf(file.path(root, "DESCRIPTION"), fields = "Version")[[1L]]
for (i in seq_along(rd_files)) {
  rd <- tools::parse_Rd(rd_files[[i]])
  output <- file.path(destination, paste0(topics[[i]], ".qmd"))
  html <- capture.output(tools::Rd2HTML(
    rd, package = c("gsim", version), fragment = FALSE, Links = links, Links2 = character(),
    texmath = "mathjax", outputEncoding = "UTF-8"
  ))
  start <- grep("<body>", html, fixed = TRUE)
  end <- grep("</body>", html, fixed = TRUE)
  stopifnot(length(start) == 1L, length(end) == 1L)
  html <- html[seq.int(start + 1L, end - 1L)]
  html <- gsub("</?main[^>]*>", "", html)
  html <- gsub('href="00Index.html"', 'href="index.html"', html, fixed = TRUE)
  # Rd2HTML supplies headings, usage, arguments, details, and literal examples.
  # Quarto provides the surrounding navigation and stylesheet.
  writeLines(c(
    "---", paste0('title: "', topics[[i]], '()"'), "---", "",
    paste0("[Authoritative Rd source](https://github.com/psoerensen/gsim/blob/main/man/",
           basename(rd_files[[i]]), ")"), "",
    "Examples are displayed only; the website build does not execute them.", "",
    "```{=html}", html, "```"
  ), output, useBytes = TRUE)
}

"""Prepare public sources and render once. Requires Python 3, base R, Quarto.

Only generated files below website/ are written. No package loading, installation,
example execution, network request, or recursive cleanup occurs here.
"""
from pathlib import Path
import argparse
import json
import os
import re
import shutil
import subprocess
from urllib.parse import urlsplit
import posixpath

SITE = Path(__file__).resolve().parent
ROOT = SITE.parent
REPO = "https://github.com/psoerensen/gsim/blob/main/"
DESIGNS = sorted((ROOT / "docs/design").glob("*.md"))
REPORTS = sorted((ROOT / "docs/qualification").glob("*.md"))
MAPPING = {p.relative_to(ROOT).as_posix(): f"contracts/{p.stem}.qmd" for p in DESIGNS}
MAPPING.update({p.relative_to(ROOT).as_posix(): f"evidence/{p.stem}.qmd" for p in REPORTS})
MAPPING.update({f"man/{p.name}": f"reference/{p.stem}.qmd" for p in (ROOT / "man").glob("*.Rd")})
MAPPING.update({"README.md": "index.qmd", "docs/README.md": "contracts.qmd",
                "inst/COPYRIGHTS": "attribution.qmd"})


def read(path):
    return (ROOT / path).read_text(encoding="utf-8")


def section(text, heading):
    marker = "## " + heading + "\n"
    if marker not in text:
        raise ValueError(f"Authoritative heading missing: {heading}")
    return marker + text.split(marker, 1)[1].split("\n## ", 1)[0].rstrip() + "\n\n"


def adapt(text, source, output):
    def link(match):
        dest = match.group(2)
        if urlsplit(dest).scheme or dest.startswith("#"):
            return match.group(0)
        path, sep, anchor = dest.partition("#")
        original = posixpath.normpath(posixpath.join(posixpath.dirname(source), path))
        if original == "docs/README.md" and anchor == "qualification-and-reproducibility":
            target = "validation.qmd"
        else:
            target = MAPPING.get(original)
        if target:
            dest = posixpath.relpath(target, posixpath.dirname(output) or ".")
        else:
            if not (ROOT / original).is_file():
                raise ValueError(f"Missing linked source: {source}: {original}")
            dest = REPO + original
        return f"[{match.group(1)}]({dest}{sep}{anchor})"

    text = re.sub(r"\[([^\]]*)\]\(([^\s)]+)\)", link, text)
    # Historical reports retain their full measurements. Replace only the
    # historical workstation checkout location, which is not a prerequisite.
    text = re.sub(r"C:/Users/[^`\s]+/synthetic_data", "<hapnest-source>", text)
    if re.search(r"(?i)gfactory|[A-Z]:[/\\](?!/)|/Users/|/home/", text):
        raise ValueError(f"Private or machine-specific path in {source}")
    return text


def page(output, title, body, source=None):
    if source:
        body = adapt(body, source, output)
        body += f"\n[Authoritative source]({REPO}{source})\n"
    path = SITE / output
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(f"---\ntitle: {json.dumps(title)}\n---\n\n" + body.rstrip() + "\n",
                    encoding="utf-8")


def prepare():
    for tool in ("Rscript", "quarto"):
        if not shutil.which(tool):
            raise SystemExit(f"Missing prerequisite: {tool} on PATH. No installation attempted.")
    readme = read("README.md")
    intro = readme.split("\n", 1)[1].split("## Installation", 1)[0]
    example = section(readme, "Examples").split("Simulate from a caller-provided", 1)[0]
    page("index.qmd", "gsim", intro + section(readme, "Installation") + example,
         "README.md")
    start = "Follow the biological stages in order. The sketches below require prepared inputs; "
    start += "use the [complete examples](examples.qmd) for runnable scripts.\n\n"
    # Resolve source-relative links before adding the site-local introduction.
    body = "".join(section(readme, h) for h in
                   ("Packed reference workflow", "Genotypes to phenotypes", "Pedigree and record workloads"))
    page("getting-started.qmd", "Getting started", start + adapt(body, "README.md", "getting-started.qmd"))
    index = read("docs/README.md")
    page("contracts.qmd", "Scientific contracts", section(index, "Authoritative contracts"), "docs/README.md")
    page("validation.qmd", "Validation and performance",
         section(index, "Qualification and reproducibility"), "docs/README.md")
    for source in DESIGNS + REPORTS:
        relative = source.relative_to(ROOT).as_posix()
        heading, body = read(relative).split("\n", 1)
        if source in REPORTS:
            body = ("Historical evidence: the measurements and recommendations below belong "
                    "to their recorded revisions and fixtures. They are not a new qualification "
                    "of this website or the current checkout. Historical local checkout paths "
                    "are displayed as portable placeholders.\n\n" + body)
        page(MAPPING[relative], heading.lstrip("# "), body, relative)
    (SITE / "reference").mkdir(exist_ok=True)
    subprocess.run(["Rscript", "--vanilla", str(SITE / "reference.R"), str(ROOT),
                    str(SITE / "reference")], check=True)
    exports = re.findall(r"^export\(([^)]+)\)", read("NAMESPACE"), re.M)
    topics = exports + ["print.gsim"]
    for topic in topics:
        if not (SITE / "reference" / f"{topic}.qmd").is_file():
            raise ValueError(f"Missing public reference: {topic}")
    page("reference/index.qmd", "Function reference",
         "Generated from public Rd without loading gsim. All exported functions "
         "and the registered print method are listed.\n\n" +
         "\n".join(f"- [`{name}()`]({name}.qmd)" for name in topics))
    files = SITE / "files"
    files.mkdir(exist_ok=True)
    for name in ("1000G_chr22.R", "end_to_end_phenotype.R", "gene_sets.R"):
        shutil.copyfile(ROOT / "inst/examples" / name, files / name)
    shutil.copyfile(ROOT / "inst/COPYRIGHTS", files / "COPYRIGHTS.txt")
    page("attribution.qmd", "Attribution and license notices",
         "[Download the unmodified notices](files/COPYRIGHTS.txt).\n\n```text\n" +
         read("inst/COPYRIGHTS").rstrip() + "\n```\n")
    page("examples.qmd", "Examples", (SITE / "examples.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--prepare-only", action="store_true")
    parser.add_argument("--page", help="Render one generated QMD after a concrete page fix")
    args = parser.parse_args()
    prepare()
    if not args.prepare_only:
        target = SITE
        if args.page:
            target = (SITE / args.page).resolve()
            if not target.is_relative_to(SITE) or target.suffix != ".qmd":
                raise SystemExit("--page must name a generated QMD below website/")
        # Quarto 1.7 uses LOCALAPPDATA on Windows and XDG_CACHE_HOME on Linux.
        # Scope cache overrides to this child process, not the user's settings.
        render_env = os.environ.copy()
        cache = SITE / ".cache"
        cache.mkdir(exist_ok=True)
        if os.name == "nt":
            render_env["LOCALAPPDATA"] = str(cache)
        else:
            render_env["XDG_CACHE_HOME"] = str(cache)
        subprocess.run(["quarto", "render", str(target)], check=True, cwd=ROOT,
                       env=render_env)

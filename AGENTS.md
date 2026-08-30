# AGENTS.md

## Project

This repository contains `gsim`, a small base-R-first package for genomic
simulation used in validation and methodological studies.

## Development rules

- Preserve the public `gsim()` API, defaults, numerical behavior, RNG behavior,
  validation, and return structure unless a task explicitly changes them.
- Keep internal helpers unexported and register S3 methods through roxygen2.
- Keep `qgg` optional; it is required only for real `Glist` operation when the
  caller does not provide `getG_fun`.
- Do not add compiled code or broad dependencies without explicit approval.
- Generate `NAMESPACE` and `man/*.Rd` with roxygen2 rather than editing them by
  hand.
- Keep tests focused and deterministic; do not add large simulations or
  external-data campaigns.
- Check `git status --short` before editing and do not overwrite unrelated user
  changes.
- Do not commit or push unless explicitly requested.

## Windows command discipline

This repository is developed on Windows.

- Use simple commands compatible with the active Windows shell.
- When commands run through `cmd.exe`, do not use POSIX single quotes to quote command bodies or arguments.
- Do not wrap commands as `cmd.exe /c 'command'`.
- For `rg` searches containing alternatives, prefer multiple `-e` arguments instead of a pipe-delimited regular expression.
- Resolve filenames with `rg --files` before referring to a specific source path.
- Do not assume that a guessed filename exists.
- Treat `rg` exit code 1 with no output as “no matches,” not as an execution failure.
- Prefer several short, explicit searches over deeply nested or heavily quoted shell commands.

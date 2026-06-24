# AGENTS.md

## Project overview

This repository contains a Snakemake-based bioinformatics workflow.

Before making changes, inspect the repository structure first. In particular, check:

- `README.md`
- `Snakefile` or `workflow/Snakefile`
- `workflow/rules/`
- `scripts/`
- `config/` or example config files
- `envs/`, `env.yaml`, or `pyproject.toml`
- `tests/`, if present

Do not assume the workflow layout before reading the actual files.

## General working rules

- Make minimal, targeted changes for the requested task.
- Do not rewrite the workflow structure unless explicitly asked.
- Do not change scientific assumptions, input formats, output formats, or default parameters unless explicitly asked.
- Keep code comments and docstrings in English.
- Use complete copyright and license headers in code and workflow files where headers are used; do not add copyright headers to YAML files.
- Preserve existing style, naming conventions, and directory layout.
- Do not add new dependencies unless necessary.
- Never commit generated outputs, temporary files, benchmark files, logs, or large data files.

## Snakemake conventions

- Keep Snakemake rules declarative where possible.
- Put non-trivial Python logic in `scripts/` or package modules instead of large `run:` blocks.
- Prefer explicit `input`, `output`, `params`, `log`, `threads`, and `resources` fields.
- Use `log:` for command logs when possible.
- Avoid hard-coded absolute paths.
- Use config values for user-adjustable paths and parameters.
- Do not put runtime-specific settings such as threads, resources, or output directories into config files unless the existing workflow already does so.
- Do not silently change wildcard behavior.
- When editing rules, check that target expansion and wildcards still resolve correctly.
- Do not use inline `lambda` functions in Snakemake rules. [Snakemake best practices](https://snakemake.readthedocs.io/en/stable/snakefiles/best_practices.html#care-about-code-readability) **explicitly** recommend avoiding `lambda` expressions inside rules and moving helper logic into named functions or scripts. Inline `lambda` logic makes rules harder to read, test, debug, and modify safely.
- Do not use Snakemake checkpoints unless explicitly requested. Checkpoints create data-dependent dynamic DAGs and make dry-runs, CI, wildcard resolution, and reruns harder to reason about. Use regular explicit rules whenever the workflow structure can be known before execution.
- When using `expand()` with patterns that intentionally keep some wildcards unresolved, use `allow_missing=True`.
- Do not add `conda: ../envs/env.yaml` or another shared main environment to rules that only use dependencies from the main workflow environment.
- Add `conda:` only for rules that require a dedicated environment.

Suggested rule field order:

1. `input`
2. `output`
3. `params`
4. `log`
5. `threads`
6. `resources`
7. `conda` / `container`
8. `shell` / `script` / `run`

## Python scripts executed by Snakemake

- Python scripts executed by Snakemake must assume the injected `snakemake` object exists.
- Do not use `globals()` checks such as `if "snakemake" in globals()`.
- Do not add standalone fallback/default parameters for non-Snakemake execution.
- Read values directly from `snakemake.input`, `snakemake.output`, `snakemake.params`, `snakemake.threads`, `snakemake.resources`, and `snakemake.log`.
- Keep reusable logic in functions where appropriate, but pass arguments from the `snakemake` object explicitly.
- Do not rely on global variables to pass workflow state.
- All Python functions must use type hints.
- All Python functions must have NumPy-style docstrings.
- Python docstrings must put the opening triple double quotes on their own first line, then start the summary on the next line.
- If a function returns `None`, do not include a `Returns` section in the docstring.

## Testing and validation

After modifying workflow logic, run at least:

- `snakemake -n -p`

If available, also run:

- `snakemake --lint`
- `snakefmt --check .`

If there is an example config or test dataset, run the smallest available real workflow test, for example:

- `snakemake -n -p --configfile <example_config>`
- `snakemake --cores 1 -p --configfile <example_config>`

For Snakemake commands, use a longer timeout because DAG construction or test runs may take time.

If Python code is changed, run relevant tests if available:

- `pytest`
- or the smallest relevant test file/function

If formatting tools exist, run the project’s configured formatter/linter, for example:

- `ruff check .`
- `ruff format .`
- `pre-commit run --all-files`

Only run commands that match the tools actually present in the repository.

## Minimal test data and CI

- Preserve existing minimal test data if present.
- Do not delete or replace test data unless explicitly asked.
- If CI configuration exists, keep workflow changes compatible with it.
- Do not add new CI jobs unless explicitly asked.

## Files not to commit

Never commit:

- `.snakemake/`
- `.snakemake/log/`
- `*.snakemake.log`
- `slurm-*.out`
- `*.tmp`
- `*.bak`
- generated workflow outputs
- benchmark files
- large input datasets
- downloaded reference data
- local conda/pixi/virtualenv directories

## Git rules

- Show `git diff` before finalizing changes.
- Do not commit or push unless explicitly asked.
- If asked to commit, run `git status` first and include only relevant files.
- Use Conventional Commits for commit messages: https://www.conventionalcommits.org
- Keep commit messages concise and focused on the actual change.
- Do not include AI-authorship lines in commit messages.
- Include a concise description of the changes in the PR in English and Link relevant related PRs or issues.

## Stop conditions

Stop and ask for direction before proceeding if:

- Scientific logic would change.
- Workflow semantics would change.
- Output filenames or formats would change.
- A new dependency is required.
- Runtime behavior depends on cluster-specific assumptions.
- Multiple reasonable designs are possible and a human decision is needed.

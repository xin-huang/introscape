# Introscape v1 SSTAR2 Workflow Plan

## Summary

Reorganize the top-level `introscape` repository into a standard Snakemake workflow project, using `selscape` only as a structural and style reference. Do not modify `selscape/`, `sai-analysis/`, or `sstar-qr-analysis/`.

Version 1 should implement only the `sstar2` workflow. `sai` remains a future method module. The goal is to connect a `selscape`-style upstream data/config layer to an `sstar2` method layer and a lightweight downstream layer that reports candidate regions and candidate genes.

The existing top-level GITA/LR/UNet workflow and old simulation-oriented configs should be removed during implementation and replaced with the new `introscape` structure.

## Key Changes

- Build a top-level Snakemake workflow:
  - `workflow/Snakefile`
  - `workflow/rules/common.smk`
  - `workflow/rules/preprocess.smk`
  - `workflow/rules/sstar2.smk`
  - `workflow/rules/annotation.smk`
  - `workflow/rules/report.smk`
  - `workflow/scripts/`
  - `workflow/envs/`
  - `workflow/schemas/`

- Rebuild configuration around the `selscape` pattern:
  - `config/main.yaml`
  - `config/datasets/example.yaml`
  - `config/methods/sstar2.yaml`
  - explicit analysis units in `main.yaml`

- Keep reference projects untouched:
  - `selscape/` is only the workflow layout and upstream/downstream template.
  - `sstar-qr-analysis/` is only the reference for `sstar2 train`, config rendering, and `sstar2 infer`.
  - `sai-analysis/` is not part of v1 implementation.

## Workflow Design

- Dataset configs follow the useful parts of `selscape`:
  - dataset ID, species, reference genome, chromosomes
  - VCF path pattern
  - sample metadata
  - optional ancestral allele BED files
  - annotation and gene mapping files

- `main.yaml` defines analysis units explicitly:
  - analysis ID
  - dataset ID
  - ref, tgt, and src populations
  - enabled method list, with `sstar2` enabled for v1
  - demes model path and method-specific overrides

- Preprocessing should:
  - extract biallelic SNPs per chromosome
  - generate ref/tgt/src sample lists from metadata
  - optionally polarize variants if ancestral allele data is configured
  - produce stable `results/processed_data/...` inputs for `sstar2`

- `sstar2` should:
  - render one config per analysis unit
  - train an ONNX model from the configured demes model
  - run inference on the processed analysis VCF
  - convert inferred tract BED output into a unified candidate-region table

- Downstream v1 should:
  - use candidate regions as the common output object
  - annotate candidate regions to genes
  - write candidate region TSV, candidate gene TSV, and candidate gene HTML
  - not force Manhattan plots in v1
  - keep enrichment as a later extension unless the minimal fixture supports it cleanly

## Test Data And CI

- Use upstream `xin-huang/selscape` GitHub test fixtures as the source for test data, not the local vendored `selscape/` directory.
- Copy or otherwise vendor the upstream `.tests/ci/data` fixture files into the top-level `introscape` test area during implementation:
  - `chr21.hg19.test.vcf.gz`
  - `test.21.bed.gz`
  - `test.gene2go.gz`
  - `test.gtf.gz`
  - `test_metadata.txt`

- Create a top-level `.tests/ci/config` modeled after upstream `selscape` CI config, adapted for `introscape` and `sstar2`.
- Add the extra minimal files needed by `sstar2`, especially a tiny demes model and method config.
- CI/test commands should start with dry-run validation and then run the smallest feasible real target up to candidate-region conversion.

## Test Plan

- Run:
  - `snakemake -n -p --configfile .tests/ci/config/main.yaml`
  - `snakemake --lint`, if available
  - `snakefmt --check .`, if available

- Add focused tests for adapter scripts using tiny BED/TSV inputs.
- If the `sstar2` environment resolves locally, run the smallest real workflow target through candidate-region conversion.
- Verify no files under `selscape/`, `sai-analysis/`, or `sstar-qr-analysis/` are modified.
- Show `git status` and `git diff` before finalizing implementation.

## Assumptions

- v1 is `sstar2` only.
- `sai` is deferred.
- `selscape` is a template, not an implementation target.
- Upstream `xin-huang/selscape` is the source of test fixture files.
- The user will handle commit and push to GitHub after reviewing local changes.

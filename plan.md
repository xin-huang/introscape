# Introscape v1 Workflow Reorganization Plan

## Summary

Build a new top-level `introscape` Snakemake workflow using `selscape` only as the structural template. Do not modify `selscape/`, `sai-analysis/`, or `sstar-qr-analysis/`. Remove the existing top-level GITA/LR-oriented `workflow/` and `config/` content, then replace it with a standard Snakemake project that runs `sai` and `sstar2` as first-version introgression analysis methods.

The v1 workflow targets real dataset analysis, not simulation benchmarking. Inputs are `selscape`-style dataset configs plus explicit `ref/tgt/src` analysis definitions. Outputs are candidate regions, annotated candidate genes, enrichment tables, enrichment plots, and simple HTML tables. Do not force Manhattan plots in v1.

## Key Changes

- Create a standard top-level workflow layout:
  - `workflow/Snakefile`
  - `workflow/rules/common.smk`
  - `workflow/rules/preprocess.smk`
  - `workflow/rules/sai.smk`
  - `workflow/rules/sstar2.smk`
  - `workflow/rules/annotation.smk`
  - `workflow/rules/enrichment.smk`
  - `workflow/rules/report.smk`
  - `workflow/scripts/`
  - `workflow/envs/`
  - `workflow/schemas/`

- Replace old top-level files:
  - Delete old GITA/LR and empty UNet rules.
  - Delete old top-level GITA simulation scripts and plotting scripts.
  - Delete old `config/scenarios`, `config/features`, and `config/demog`.
  - Recreate `config/` around the new `main + datasets + methods` layout.

- Keep external reference projects untouched:
  - `selscape/` remains read-only reference for workflow organization, preprocessing, annotation, enrichment, and report conventions.
  - `sai-analysis/` remains read-only reference for `sai score`, `sai outlier`, and 1-source config shape.
  - `sstar-qr-analysis/` remains read-only reference for `sstar2 train`, config rendering, and `sstar2 infer`.

## Workflow Design

- Main config:
  - `config/main.yaml` lists dataset config files, method config files, and explicit analysis units.
  - Each analysis unit has an ID, dataset ID, ref population, tgt population, src population, enabled methods, and method-specific overrides where needed.

- Dataset config:
  - Follows the useful parts of `selscape` dataset schema: dataset, species, ref genome, VCF folder/prefix/suffix, metadata, chromosomes, genome annotation, gene2go, repeats, ancestral alleles, ploidy.
  - Metadata format is `Sample` + `Population`; all ref/tgt/src populations must exist in the same multisample VCF for v1.

- Preprocessing:
  - Extract biallelic SNPs per chromosome.
  - Generate per-analysis sample lists for `ref`, `tgt`, and `src`.
  - Optionally polarize variants only when ancestral allele config is present.
  - Produce cleaned VCFs and sample lists in stable `results/processed_data/...` paths.

- SAI method:
  - Implement v1 as generic 1-source analysis.
  - Generate a `sai` config per analysis unit from ref/tgt/src sample lists and method config.
  - Run `sai score` per chromosome/window.
  - Merge scores and run `sai outlier`.
  - Convert SAI outlier windows into the unified candidate-region table.

- SSTAR2 method:
  - Render an `sstar2` config per analysis unit from method config, demes path, sample lists, window settings, phase setting, and training parameters.
  - Train ONNX model from configured demes.
  - Run `sstar2 infer` on the analysis VCF.
  - Convert inferred tract BED into the unified candidate-region table.
  - Do not create artificial Manhattan scores for SSTAR2 in v1.

- Unified downstream:
  - Candidate-region table columns: `method`, `analysis`, `dataset`, `chrom`, `start`, `end`, plus method-specific metadata columns.
  - Annotate candidate regions against ANNOVAR-style multianno outputs or an equivalent derived annotation table.
  - Produce candidate SNP/region lists, candidate gene tables, enrichment input files, GOWINDA enrichment TSVs, enrichment plots, and HTML tables.
  - `rule all` targets candidate region tables, candidate gene HTML/TSV, enrichment TSV/PNG/HTML for every enabled method and analysis unit.

## Interfaces

- Add dedicated environments:
  - `workflow/envs/sai.yaml` based on `sai-analysis` dependencies, with `sai-pg`.
  - `workflow/envs/sstar2.yaml` based on `sstar-qr-analysis`, with `sstar2` dependencies.
  - `workflow/envs/introscape.yaml` for shared preprocessing, annotation, plotting, and Snakemake-side scripts.

- Add schemas:
  - `workflow/schemas/config.schema.yaml`
  - `workflow/schemas/dataset.schema.yaml`
  - `workflow/schemas/sai.schema.yaml`
  - `workflow/schemas/sstar2.schema.yaml`

- Add example configs:
  - `config/main.yaml`
  - `config/datasets/example.yaml`
  - `config/methods/sai.yaml`
  - `config/methods/sstar2.yaml`

## Test Plan

- Run `snakemake -n -p --configfile config/main.yaml` from the repository root.
- Run `snakemake --lint` if supported by the installed Snakemake version.
- Run `snakefmt --check .` if `snakefmt` is available.
- Run focused unit-style checks for Python adapter scripts with small artificial TSV/BED inputs.
- Verify no changes occur under `selscape/`, `sai-analysis/`, or `sstar-qr-analysis/`.
- Verify `git diff` before finalizing implementation.

## Assumptions

- First version analyzes real multisample VCF datasets, not simulated benchmark scenarios.
- Each configured `ref`, `tgt`, and `src` population is present in the dataset metadata and VCF.
- SAI v1 supports only generic 1-source mode.
- SSTAR2 v1 always trains from configured demes before inference.
- Manhattan plots are out of scope for v1; enrichment plots and HTML/TSV tables are in scope.
- Existing top-level GITA/LR workflow and config content should be removed during implementation.

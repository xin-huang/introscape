# Copyright 2026 Xin Huang
#
# GNU General Public License v3.0
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, please see
#
#    https://www.gnu.org/licenses/gpl-3.0.en.html

from typing import Any

import yaml
from snakemake.utils import validate


main_config = config

validate(main_config, schema="../schemas/config.schema.yaml")


def _load_yaml(path: str) -> dict[str, Any]:
    """
    Load a YAML file.

    Parameters
    ----------
    path : str
        YAML file path.

    Returns
    -------
    dict[str, Any]
        Parsed YAML data.
    """
    with open(path, "r") as handle:
        data = yaml.safe_load(handle)
    return data or {}


dataset_configs = {}
for dataset_config_path in main_config["datasets"]:
    dataset_config = _load_yaml(dataset_config_path)
    validate(dataset_config, schema="../schemas/dataset.schema.yaml")
    dataset_configs[dataset_config["dataset"]] = dataset_config


method_configs = {}
for method_name, method_config_path in main_config["methods"].items():
    method_config = _load_yaml(method_config_path)
    if method_name == "sstar2":
        validate(method_config, schema="../schemas/sstar2.schema.yaml")
    method_configs[method_name] = method_config


analysis_configs = {analysis["id"]: analysis for analysis in main_config["analyses"]}
sstar2_analyses = [
    analysis for analysis in main_config["analyses"] if "sstar2" in analysis["methods"]
]


def get_analysis_config(wildcards: Any) -> dict[str, Any]:
    """
    Return one analysis configuration.

    Parameters
    ----------
    wildcards : Any
        Snakemake wildcards with an ``analysis`` attribute.

    Returns
    -------
    dict[str, Any]
        Analysis configuration.
    """
    return analysis_configs[wildcards.analysis]


def get_dataset_config_by_id(dataset_id: str) -> dict[str, Any]:
    """
    Return one dataset configuration.

    Parameters
    ----------
    dataset_id : str
        Dataset identifier.

    Returns
    -------
    dict[str, Any]
        Dataset configuration.
    """
    return dataset_configs[dataset_id]


def get_dataset_config(wildcards: Any) -> dict[str, Any]:
    """
    Return the dataset configuration for wildcards.

    Parameters
    ----------
    wildcards : Any
        Snakemake wildcards with a ``dataset`` attribute.

    Returns
    -------
    dict[str, Any]
        Dataset configuration.
    """
    return get_dataset_config_by_id(wildcards.dataset)


def get_analysis_dataset_config(wildcards: Any) -> dict[str, Any]:
    """
    Return the dataset configuration for an analysis.

    Parameters
    ----------
    wildcards : Any
        Snakemake wildcards with an ``analysis`` attribute.

    Returns
    -------
    dict[str, Any]
        Dataset configuration.
    """
    analysis_config = get_analysis_config(wildcards)
    return get_dataset_config_by_id(analysis_config["dataset"])


def get_dataset_chromosomes(wildcards: Any) -> list[Any]:
    """
    Return chromosomes for a dataset wildcard.

    Parameters
    ----------
    wildcards : Any
        Snakemake wildcards with a ``dataset`` attribute.

    Returns
    -------
    list[Any]
        Chromosome identifiers.
    """
    return get_dataset_config(wildcards)["chromosomes"]


def get_analysis_chromosomes(wildcards: Any) -> list[Any]:
    """
    Return chromosomes for an analysis wildcard.

    Parameters
    ----------
    wildcards : Any
        Snakemake wildcards with an ``analysis`` attribute.

    Returns
    -------
    list[Any]
        Chromosome identifiers.
    """
    return get_analysis_dataset_config(wildcards)["chromosomes"]


def get_vcf_input_path(wildcards: Any) -> str:
    """
    Return a raw VCF path for a dataset chromosome.

    Parameters
    ----------
    wildcards : Any
        Snakemake wildcards with ``dataset`` and ``chrom`` attributes.

    Returns
    -------
    str
        Raw VCF path.
    """
    dataset_config = get_dataset_config(wildcards)
    return (
        f"{dataset_config['data_folder']}/"
        f"{dataset_config['vcf_prefix']}{wildcards.chrom}{dataset_config['vcf_suffix']}"
    )


def get_metadata(wildcards: Any) -> str:
    """
    Return a metadata path for an analysis.

    Parameters
    ----------
    wildcards : Any
        Snakemake wildcards with an ``analysis`` attribute.

    Returns
    -------
    str
        Metadata path.
    """
    return get_analysis_dataset_config(wildcards)["metadata"]


def get_genome_annotation(wildcards: Any) -> str:
    """
    Return a genome annotation path for an analysis.

    Parameters
    ----------
    wildcards : Any
        Snakemake wildcards with an ``analysis`` attribute.

    Returns
    -------
    str
        Genome annotation path.
    """
    return get_analysis_dataset_config(wildcards)["genome_annotation"]


def get_ref_population(wildcards: Any) -> str:
    """
    Return the real-data reference population for an analysis.

    Parameters
    ----------
    wildcards : Any
        Snakemake wildcards with an ``analysis`` attribute.

    Returns
    -------
    str
        Reference population identifier.
    """
    return get_analysis_config(wildcards)["ref_population"]


def get_tgt_population(wildcards: Any) -> str:
    """
    Return the real-data target population for an analysis.

    Parameters
    ----------
    wildcards : Any
        Snakemake wildcards with an ``analysis`` attribute.

    Returns
    -------
    str
        Target population identifier.
    """
    return get_analysis_config(wildcards)["tgt_population"]


def get_src_population(wildcards: Any) -> str:
    """
    Return the real-data source population for an analysis.

    Parameters
    ----------
    wildcards : Any
        Snakemake wildcards with an ``analysis`` attribute.

    Returns
    -------
    str
        Source population identifier, or an empty string when not configured.
    """
    source_population = get_analysis_config(wildcards).get("src_population")
    return source_population or ""


def get_analysis_dataset_id(wildcards: Any) -> str:
    """
    Return the dataset ID for an analysis.

    Parameters
    ----------
    wildcards : Any
        Snakemake wildcards with an ``analysis`` attribute.

    Returns
    -------
    str
        Dataset identifier.
    """
    return get_analysis_config(wildcards)["dataset"]


def get_analysis_demography(wildcards: Any) -> str:
    """
    Return the demography model path for an analysis.

    Parameters
    ----------
    wildcards : Any
        Snakemake wildcards with an ``analysis`` attribute.

    Returns
    -------
    str
        Demography model path.
    """
    return get_analysis_config(wildcards)["demography"]


def get_biallelic_vcf_for_analysis(wildcards: Any) -> str:
    """
    Return a biallelic VCF path for an analysis chromosome.

    Parameters
    ----------
    wildcards : Any
        Snakemake wildcards with ``analysis`` and ``chrom`` attributes.

    Returns
    -------
    str
        Biallelic VCF path.
    """
    dataset_id = get_analysis_dataset_id(wildcards)
    return f"results/processed_data/{dataset_id}/all/{wildcards.chrom}.biallelic.snps.vcf.gz"


def get_biallelic_index_for_analysis(wildcards: Any) -> str:
    """
    Return a biallelic VCF index path for an analysis chromosome.

    Parameters
    ----------
    wildcards : Any
        Snakemake wildcards with ``analysis`` and ``chrom`` attributes.

    Returns
    -------
    str
        Biallelic VCF index path.
    """
    return f"{get_biallelic_vcf_for_analysis(wildcards)}.tbi"


def get_analysis_chrom_vcfs(wildcards: Any) -> list[str]:
    """
    Return per-chromosome analysis VCF paths.

    Parameters
    ----------
    wildcards : Any
        Snakemake wildcards with an ``analysis`` attribute.

    Returns
    -------
    list[str]
        Per-chromosome VCF paths.
    """
    return expand(
        "results/processed_data/{analysis}/chromosomes/{chrom}.biallelic.snps.vcf.gz",
        analysis=wildcards.analysis,
        chrom=get_analysis_chromosomes(wildcards),
    )


def get_analysis_chrom_indexes(wildcards: Any) -> list[str]:
    """
    Return per-chromosome analysis VCF index paths.

    Parameters
    ----------
    wildcards : Any
        Snakemake wildcards with an ``analysis`` attribute.

    Returns
    -------
    list[str]
        Per-chromosome VCF index paths.
    """
    return [f"{path}.tbi" for path in get_analysis_chrom_vcfs(wildcards)]


def get_sstar2_method_config_path(wildcards: Any) -> str:
    """
    Return the configured sstar2 method config path.

    Parameters
    ----------
    wildcards : Any
        Snakemake wildcards. The value is unused but accepted by Snakemake.

    Returns
    -------
    str
        Method configuration path.
    """
    return main_config["methods"]["sstar2"]


def get_method_config_title(wildcards: Any) -> str:
    """
    Return a display title for reports.

    Parameters
    ----------
    wildcards : Any
        Snakemake wildcards with an ``analysis`` attribute.

    Returns
    -------
    str
        Display title.
    """
    return f"{wildcards.analysis} sstar2 candidate genes"


ALL_TARGETS = []
for analysis_config in sstar2_analyses:
    analysis_id = analysis_config["id"]
    ALL_TARGETS.extend(
        [
            f"results/sstar2/{analysis_id}/candidates/{analysis_id}.sstar2.candidate.regions.tsv",
            f"results/sstar2/{analysis_id}/candidates/{analysis_id}.sstar2.annotated.candidate.regions.tsv",
            f"results/sstar2/{analysis_id}/candidates/{analysis_id}.sstar2.candidate.genes.tsv",
            f"results/sstar2/{analysis_id}/candidates/{analysis_id}.sstar2.candidate.genes.html",
        ]
    )

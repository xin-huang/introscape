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

from copy import deepcopy
from pathlib import Path
from typing import Any

import yaml


def load_yaml(path: str) -> dict[str, Any]:
    """
    Load YAML data.

    Parameters
    ----------
    path : str
        YAML file path.

    Returns
    -------
    dict[str, Any]
        Parsed YAML mapping.
    """
    with open(path, "r") as handle:
        data = yaml.safe_load(handle)
    return data or {}


def render_config(
    method_config: dict[str, Any],
    vcf_file: str,
    ref_ind_file: str,
    tgt_ind_file: str,
) -> dict[str, Any]:
    """
    Render an sstar2 configuration for one analysis.

    Parameters
    ----------
    method_config : dict[str, Any]
        Base sstar2 method configuration.
    vcf_file : str
        Processed analysis VCF path.
    ref_ind_file : str
        Reference sample-list path.
    tgt_ind_file : str
        Target sample-list path.

    Returns
    -------
    dict[str, Any]
        Rendered sstar2 configuration.
    """
    rendered = deepcopy(method_config)
    preprocessing = rendered.setdefault("preprocessing", {})
    preprocessing["vcf_file"] = vcf_file
    preprocessing["ref_ind_file"] = ref_ind_file
    preprocessing["tgt_ind_file"] = tgt_ind_file
    return rendered


def write_yaml(data: dict[str, Any], path: str) -> None:
    """
    Write YAML data.

    Parameters
    ----------
    data : dict[str, Any]
        YAML mapping to write.
    path : str
        Output YAML path.
    """
    output_path = Path(path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w") as handle:
        yaml.safe_dump(data, handle, sort_keys=False)


method_config = load_yaml(str(snakemake.input.method_config))
rendered_config = render_config(
    method_config,
    str(snakemake.input.vcf),
    str(snakemake.input.ref),
    str(snakemake.input.tgt),
)
write_yaml(rendered_config, str(snakemake.output.config))

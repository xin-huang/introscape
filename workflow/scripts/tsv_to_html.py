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

import csv
import html
from pathlib import Path


def read_tsv(path: str) -> tuple[list[str], list[dict[str, str]]]:
    """
    Read a TSV file.

    Parameters
    ----------
    path : str
        TSV file path.

    Returns
    -------
    tuple[list[str], list[dict[str, str]]]
        Header and rows.
    """
    with open(path, "r", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return list(reader.fieldnames or []), list(reader)


def write_html(
    title: str,
    fieldnames: list[str],
    rows: list[dict[str, str]],
    output_path: str,
) -> None:
    """
    Write an HTML table.

    Parameters
    ----------
    title : str
        Page title.
    fieldnames : list[str]
        Table field names.
    rows : list[dict[str, str]]
        Table rows.
    output_path : str
        Output HTML path.
    """
    output = Path(output_path)
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w") as handle:
        handle.write("<!doctype html>\n<html lang=\"en\">\n<head>\n")
        handle.write("<meta charset=\"utf-8\">\n")
        handle.write(f"<title>{html.escape(title)}</title>\n")
        handle.write("</head>\n<body>\n")
        handle.write(f"<h1>{html.escape(title)}</h1>\n")
        handle.write("<table>\n<thead><tr>")
        for field in fieldnames:
            handle.write(f"<th>{html.escape(field)}</th>")
        handle.write("</tr></thead>\n<tbody>\n")
        for row in rows:
            handle.write("<tr>")
            for field in fieldnames:
                handle.write(f"<td>{html.escape(row.get(field, ''))}</td>")
            handle.write("</tr>\n")
        handle.write("</tbody>\n</table>\n</body>\n</html>\n")


fields, table_rows = read_tsv(str(snakemake.input.tsv))
write_html(str(snakemake.params.title), fields, table_rows, str(snakemake.output.html))

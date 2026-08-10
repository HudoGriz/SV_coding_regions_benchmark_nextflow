#!/usr/bin/env python3
"""Refresh supplementary result tables and per-simulation sheets from a clean run."""
from __future__ import annotations

import argparse
import csv
import re
import shutil
from pathlib import Path


TARGET_LABELS = {
    "high_confidence": "HCI",
    "gene_panel": "GP",
    "wes_utr": "EX+UTR",
}
TECH_LABELS = {
    "Illumina_WES": "Illumina WES",
    "Illumina_WGS": "Illumina WGS",
    "ONT": "ONT",
    "PacBio": "PacBio",
}
CALLER_LABELS = {"Manta": "Manta", "Sniffles": "Sniffles", "CuteSV": "cuteSV", "Pbsv": "pbsv"}
TARGET_ORDER = {name: index for index, name in enumerate(TARGET_LABELS)}
TECH_ORDER = {name: index for index, name in enumerate(TECH_LABELS)}
CALLER_ORDER = {"Manta": 0, "Sniffles": 1, "CuteSV": 2, "Pbsv": 3}


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def read_simulation_summary(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader)
        rows = list(reader)
    if not rows or any(len(row) != len(header) + 1 for row in rows):
        raise ValueError(f"expected an unnamed row-key column in {path}")
    return [dict(zip(["data_name", *header], row)) for row in rows]


def sort_rows(rows: list[dict[str, str]]) -> list[dict[str, str]]:
    return sorted(rows, key=lambda row: (
        TARGET_ORDER.get(row["range"], 99), TECH_ORDER.get(row["tech"], 99),
        CALLER_ORDER.get(row["caller"], 99)
    ))


def pipeline_label(key: str) -> str:
    tech, caller, target = key.rsplit("-", 2)
    return f"{TECH_LABELS[tech]} {CALLER_LABELS[caller]} / {TARGET_LABELS[target]}"


def real_rows(path: Path) -> list[str]:
    output = []
    for row in sort_rows(read_tsv(path)):
        output.append(
            f'{TARGET_LABELS[row["range"]]} & {TECH_LABELS[row["tech"]]} & {CALLER_LABELS[row["caller"]]} & '
            f'{int(float(row["TP.base"]))} & {int(float(row["FP"]))} & {int(float(row["FN"]))} & '
            f'{float(row["precision"]):.3f} & {float(row["recall"]):.3f} & {float(row["f1"]):.3f} & '
            f'{int(float(row["base.cnt"]))} & {int(float(row["comp.cnt"]))} & {float(row["gt_concordance"]):.3f}'
            + r' \\'
        )
    return output


def simulation_rows(path: Path, summary: bool) -> list[str]:
    rows = read_simulation_summary(path)
    rows.sort(key=lambda row: (
        TARGET_ORDER.get(row["data_name"].rsplit("-", 1)[-1], 99),
        pipeline_label(row["data_name"])
    ))
    output = []
    for row in rows:
        key = row["data_name"]
        label = pipeline_label(key)
        if summary:
            fields = [
                "precision_median", "recall_median", "f1_median", "precision_diff", "recall_diff", "f1_diff",
                "precision_sd", "recall_sd", "f1_sd",
            ]
            values = " & ".join(f'{float(row[field]):.3f}' for field in fields)
        else:
            percentile_fields = ["precision_percentile", "recall_percentile", "f1_percentile"]
            kde_fields = ["precision_kde", "recall_kde", "f1_kde"]
            values = " & ".join(
                [f'{float(row[field]):.2f}' for field in percentile_fields]
                + [f'{float(row[field]):.3f}' for field in kde_fields]
            )
        output.append(f"{label} & {values}" + r" \\")
    return output


def replace_table_rows(text: str, table_number: int, rows: list[str]) -> str:
    heading = rf"(\\subsection\*\{{Supplementary Table {table_number}\..*?\\endhead\n)"
    pattern = heading + r".*?(\n\\bottomrule)"
    updated, count = re.subn(
        pattern,
        lambda match: match.group(1) + "\n".join(rows) + match.group(2),
        text,
        count=1,
        flags=re.DOTALL,
    )
    if count != 1:
        raise RuntimeError(f"could not uniquely replace Supplementary Table {table_number}")
    return updated


def write_csv_from_tsv(source: Path, destination: Path) -> None:
    rows = read_tsv(source)
    if not rows:
        raise ValueError(f"empty source table: {source}")
    with destination.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-root", required=True, type=Path)
    parser.add_argument("--submission", required=True, type=Path)
    args = parser.parse_args()

    tex_path = args.submission / "supplementary" / "supplementary_material.tex"
    text = tex_path.read_text()
    for assembly, numbers in (("GRCh37", (1, 2, 3)), ("GRCh38", (4, 5, 6))):
        table_dir = args.run_root / f"results-{assembly}" / "statistics" / "tables"
        real = table_dir / "truvari_metrics_real_intervals.tsv"
        simulations = table_dir / "truvari_metrics_simulated_intervals.tsv"
        raw = table_dir / "truvari_metrics_simulated_intervals_raw.tsv"
        text = replace_table_rows(text, numbers[0], real_rows(real))
        text = replace_table_rows(text, numbers[1], simulation_rows(simulations, summary=False))
        text = replace_table_rows(text, numbers[2], simulation_rows(simulations, summary=True))
        write_csv_from_tsv(
            raw,
            args.submission / "supplementary" / f"Supplementary_Data_Sheet_{1 if assembly == 'GRCh37' else 2}_{assembly}_per_simulation.csv",
        )
        evidence_dir = args.submission / "evidence" / f"clean_rerun_2026-08-03-v3_{assembly}"
        evidence_dir.mkdir(parents=True, exist_ok=True)
        for source in (real, simulations, raw):
            shutil.copy2(source, evidence_dir / source.name)
    tex_path.write_text(text)


if __name__ == "__main__":
    main()

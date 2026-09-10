#!/usr/bin/env python3

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

from measure_lif_embryo_batch import (
    STAGE_ORDER,
    SPECIES_ORDER,
    exclude_known_mislabeled_rows,
    plot_combined,
    plot_metric,
    write_stats,
    write_summary,
)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--folder", type=Path, required=True)
    args = parser.parse_args()

    table_dir = args.folder / "tables"
    fig_dir = args.folder / "figures"
    code_dir = args.folder / "code"
    fig_dir.mkdir(parents=True, exist_ok=True)
    code_dir.mkdir(parents=True, exist_ok=True)

    current = table_dir / "embryo_dic_size_measurements.csv"
    pre_exclusion = table_dir / "embryo_dic_size_measurements_before_manual_exclusion.csv"
    if pre_exclusion.exists():
        source = pre_exclusion
    else:
        source = current

    raw_df = pd.read_csv(source)
    raw_df["stage"] = pd.Categorical(raw_df["stage"], STAGE_ORDER, ordered=True)
    raw_df["species"] = pd.Categorical(raw_df["species"], SPECIES_ORDER, ordered=True)
    raw_df = raw_df.sort_values(["stage", "species", "series_index"])
    raw_df.to_csv(pre_exclusion, index=False)

    df, excluded = exclude_known_mislabeled_rows(raw_df)
    df["stage"] = pd.Categorical(df["stage"], STAGE_ORDER, ordered=True)
    df["species"] = pd.Categorical(df["species"], SPECIES_ORDER, ordered=True)
    df = df.sort_values(["stage", "species", "series_index"])

    df.to_csv(current, index=False)
    excluded.to_csv(table_dir / "embryo_dic_size_manual_exclusions.csv", index=False)
    write_summary(df, table_dir / "embryo_dic_size_summary_by_stage_species.csv")
    stat_file = table_dir / "embryo_dic_size_wilcoxon_rank_sum_stats_by_stage.csv"
    write_stats(df, stat_file)
    stat_df = pd.read_csv(stat_file)

    plot_metric(df, "length_um", "Length (um)", fig_dir / "Figure1_length_by_stage_species", stat_df)
    plot_metric(df, "width_um", "Width (um)", fig_dir / "Figure2_width_by_stage_species", stat_df)
    plot_metric(df, "depth_um", "Depth (um)", fig_dir / "Figure3_depth_by_stage_species", stat_df)
    plot_metric(df, "ellipsoid_volume_um3", "Ellipsoid volume (um^3)", fig_dir / "Figure4_volume_by_stage_species", stat_df)
    plot_combined(df, fig_dir / "Figure5_size_metrics_combined", stat_df)

    readme_append = """

## Manual revision

This folder was regenerated from the existing smooth-ellipse measurement table without rereading the LIF file.

Revision details:
- Excluded the longest 4-cell C. elegans embryo because it was user-confirmed as mislabeled.
- Added n labels for each species group in every panel.
- Added pairwise Wilcoxon rank-sum significance brackets with BH-adjusted labels.
- Set length, width, and depth y-axis ranges to 0-75 um.
- Set volume y-axis to start from 0.
"""
    readme = args.folder / "README.md"
    old = readme.read_text() if readme.exists() else ""
    if "## Manual revision" not in old:
        readme.write_text(old + readme_append)


if __name__ == "__main__":
    main()

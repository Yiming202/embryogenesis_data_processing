#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import math
import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from PIL import Image, ImageDraw
from readlif.reader import LifFile
from scipy import stats
from skimage import exposure, filters, measure, morphology, segmentation
from skimage.feature import canny
from skimage.draw import ellipse as draw_ellipse


SPECIES_LABELS = {
    "cel": "C. elegans",
    "cbr": "C. briggsae",
    "cni": "C. nigoni",
}
SPECIES_ORDER = ["C. elegans", "C. briggsae", "C. nigoni"]
STAGE_ORDER = ["4_Cell", "350_Cell", "550_Cell"]
COLORS = {
    "C. elegans": "#CD534C",
    "C. briggsae": "#0073C2",
    "C. nigoni": "#EFC000",
}


def parse_stage_species(name: str) -> tuple[str, str, str]:
    parts = name.split("/")
    stage = parts[0] if len(parts) > 0 else "unknown"
    collection = parts[1] if len(parts) > 1 else "unknown"
    prefix = collection.split("_")[0].lower()
    species = SPECIES_LABELS.get(prefix, "unknown")
    return stage, collection, species


def microns_per_pixel(scale_value: float | None) -> float:
    if scale_value is None or not np.isfinite(scale_value) or scale_value <= 0:
        return float("nan")
    return 1.0 / float(scale_value)


def robust_uint8(img: np.ndarray) -> np.ndarray:
    arr = img.astype(float)
    lo, hi = np.nanpercentile(arr, [1, 99])
    if not np.isfinite(lo) or not np.isfinite(hi) or hi <= lo:
        lo, hi = float(np.nanmin(arr)), float(np.nanmax(arr))
    if hi <= lo:
        return np.zeros(arr.shape, dtype=np.uint8)
    out = np.clip((arr - lo) / (hi - lo), 0, 1)
    return (out * 255).astype(np.uint8)


def read_stack(img, channel: int) -> np.ndarray:
    return np.stack([np.asarray(img.get_frame(z=z, c=channel)) for z in range(img.dims.z)], axis=0)


def best_embryo_component(mask: np.ndarray) -> np.ndarray:
    lab = measure.label(mask)
    props = measure.regionprops(lab)
    if not props:
        return mask
    h, w = mask.shape
    center = np.array([h / 2.0, w / 2.0])
    diag = np.hypot(h, w)
    candidates = []
    for p in props:
        area_frac = p.area / float(h * w)
        if area_frac < 0.015 or area_frac > 0.75:
            continue
        centroid = np.array(p.centroid)
        center_penalty = np.linalg.norm(centroid - center) / diag
        # Prefer a large, central, smooth embryo-like object over edge noise or debris.
        score = area_frac - 0.55 * center_penalty - 0.05 * abs(p.eccentricity - 0.75)
        candidates.append((score, p))
    if not candidates:
        candidates = [(p.area, p) for p in props]
    best = max(candidates, key=lambda item: item[0])[1]
    return lab == best.label


def smooth_ellipse_from_mask(mask: np.ndarray, pad_fraction: float = 0.02) -> tuple[np.ndarray, object]:
    lab = measure.label(mask)
    props = measure.regionprops(lab)
    if not props:
        raise RuntimeError("No connected component found for ellipse fitting.")
    prop = max(props, key=lambda p: p.area)
    h, w = mask.shape
    major = prop.major_axis_length * (1.0 + pad_fraction)
    minor = prop.minor_axis_length * (1.0 + pad_fraction)
    rr, cc = draw_ellipse(
        r=prop.centroid[0],
        c=prop.centroid[1],
        r_radius=max(minor / 2.0, 1.0),
        c_radius=max(major / 2.0, 1.0),
        rotation=prop.orientation + np.pi / 2.0,
        shape=mask.shape,
    )
    ellipse_mask = np.zeros(mask.shape, dtype=bool)
    ellipse_mask[rr, cc] = True
    ellipse_props = measure.regionprops(measure.label(ellipse_mask))
    ellipse_prop = max(ellipse_props, key=lambda p: p.area)
    return ellipse_mask, ellipse_prop


def segment_xy(stack: np.ndarray) -> np.ndarray:
    mip = stack.max(axis=0).astype(float)
    lo, hi = np.percentile(mip, [1, 99])
    img = exposure.rescale_intensity(mip, in_range=(float(lo), float(hi)), out_range=(0, 1))
    smooth = filters.gaussian(img, sigma=2)
    edges = canny(smooth, sigma=2)
    dilated = morphology.binary_dilation(edges, morphology.disk(5))
    filled = morphology.remove_small_holes(dilated, area_threshold=20000)
    cleaned = morphology.binary_closing(filled, morphology.disk(9))
    cleaned = morphology.remove_small_objects(cleaned, min_size=2500)
    if cleaned.mean() < 0.01 or cleaned.mean() > 0.75:
        thresh = filters.threshold_otsu(smooth)
        cleaned = smooth > thresh
        cleaned = morphology.binary_closing(cleaned, morphology.disk(9))
        cleaned = morphology.remove_small_objects(cleaned, min_size=2500)
    return best_embryo_component(cleaned)


def estimate_z_extent(stack: np.ndarray, xy_mask: np.ndarray, z_um: float) -> tuple[int, int, float, float, str]:
    """Estimate embryo-containing z-range from DIC contrast inside the xy embryo mask."""
    profile = []
    for z in range(stack.shape[0]):
        plane = stack[z].astype(float)
        inside = plane[xy_mask]
        outside = plane[~xy_mask]
        if inside.size == 0 or outside.size == 0:
            profile.append(0.0)
        else:
            med_delta = abs(np.median(inside) - np.median(outside))
            texture = np.percentile(inside, 90) - np.percentile(inside, 10)
            profile.append(float(med_delta + 0.25 * texture))
    profile = np.asarray(profile)
    if np.nanmax(profile) <= 0:
        return 0, stack.shape[0] - 1, stack.shape[0] * z_um, 0.0, "failed_full_stack"
    threshold = max(np.nanpercentile(profile, 35), 0.25 * np.nanmax(profile))
    valid = np.where(profile >= threshold)[0]
    if valid.size == 0:
        z0, z1, method = 0, stack.shape[0] - 1, "no_valid_full_stack"
    else:
        z0, z1, method = int(valid.min()), int(valid.max()), "contrast_profile"
    return z0, z1, (z1 - z0 + 1) * z_um, float(threshold), method


def save_overlay(stack: np.ndarray, mask: np.ndarray, row: dict, out_png: Path) -> None:
    overlay = Image.fromarray(robust_uint8(stack.max(axis=0))).convert("RGB")
    boundaries = segmentation.find_boundaries(mask, mode="outer")
    arr = np.asarray(overlay).copy()
    arr[boundaries] = [255, 0, 0]
    overlay = Image.fromarray(arr)
    draw = ImageDraw.Draw(overlay)
    text = (
        f"{row['series_name']}\n"
        f"{row['species']}, {row['stage']}; channel={row['dic_channel_used']}\n"
        f"L={row['length_um']:.2f} um; W={row['width_um']:.2f} um; D={row['depth_um']:.2f} um"
    )
    draw.text((10, 10), text, fill=(255, 255, 0))
    overlay.save(out_png)


def measure_one(img, series_index: int, preferred_dic_channel: int, overlay_dir: Path | None) -> dict:
    stage, collection, species = parse_stage_species(img.name)
    dic_channel = min(preferred_dic_channel, img.channels - 1)

    stack = read_stack(img, dic_channel)
    x_um = microns_per_pixel(img.scale[0])
    y_um = microns_per_pixel(img.scale[1])
    z_um = microns_per_pixel(img.scale[2])
    xy_um = float(np.nanmean([x_um, y_um]))

    raw_mask = segment_xy(stack)
    mask, prop = smooth_ellipse_from_mask(raw_mask)

    length_um = float(prop.major_axis_length * xy_um)
    width_um = float(prop.minor_axis_length * xy_um)
    z0, z1, depth_um, z_threshold, z_method = estimate_z_extent(stack, mask, z_um)
    volume = math.pi / 6.0 * length_um * width_um * depth_um
    mask_area_um2 = float(prop.area * xy_um * xy_um)

    row = {
        "series_index": series_index,
        "series_name": img.name,
        "stage": stage,
        "collection": collection,
        "species": species,
        "dic_channel_used": dic_channel,
        "dic_channel_rule": f"min(preferred {preferred_dic_channel}, channels - 1)",
        "channels": img.channels,
        "x_pixels": img.dims.x,
        "y_pixels": img.dims.y,
        "z_slices": img.dims.z,
        "xy_um_per_pixel": xy_um,
        "z_um_per_slice": z_um,
        "z_start_slice": z0,
        "z_end_slice": z1,
        "z_threshold": z_threshold,
        "z_extent_method": z_method,
        "mask_area_um2": mask_area_um2,
        "boundary_type": "smooth fitted ellipse from DIC embryo candidate",
        "length_um": length_um,
        "width_um": width_um,
        "depth_um": depth_um,
        "ellipsoid_volume_um3": volume,
        "qc_flag_full_z_stack": (z0 == 0 and z1 == img.dims.z - 1),
    }
    if overlay_dir is not None:
        safe_name = re.sub(r"[^A-Za-z0-9_.-]+", "_", img.name)
        save_overlay(stack, mask, row, overlay_dir / f"{series_index:03d}_{safe_name}.png")
    return row


def write_summary(df: pd.DataFrame, out_csv: Path) -> None:
    rows = []
    metrics = ["length_um", "width_um", "depth_um", "ellipsoid_volume_um3"]
    for (stage, species), sub in df.groupby(["stage", "species"], sort=False):
        for metric in metrics:
            vals = sub[metric].dropna().to_numpy()
            rows.append(
                {
                    "stage": stage,
                    "species": species,
                    "metric": metric,
                    "n": len(vals),
                    "mean": np.mean(vals) if len(vals) else np.nan,
                    "sd": np.std(vals, ddof=1) if len(vals) > 1 else np.nan,
                    "median": np.median(vals) if len(vals) else np.nan,
                    "q25": np.percentile(vals, 25) if len(vals) else np.nan,
                    "q75": np.percentile(vals, 75) if len(vals) else np.nan,
                }
            )
    pd.DataFrame(rows).to_csv(out_csv, index=False)


def write_stats(df: pd.DataFrame, out_csv: Path) -> None:
    rows = []
    metrics = ["length_um", "width_um", "depth_um", "ellipsoid_volume_um3"]
    pairs = [("C. elegans", "C. briggsae"), ("C. elegans", "C. nigoni"), ("C. briggsae", "C. nigoni")]
    for stage in STAGE_ORDER:
        stage_df = df[df["stage"] == stage]
        for metric in metrics:
            p_values = []
            temp = []
            for a, b in pairs:
                va = stage_df.loc[stage_df["species"] == a, metric].dropna()
                vb = stage_df.loc[stage_df["species"] == b, metric].dropna()
                if len(va) >= 2 and len(vb) >= 2:
                    _, p = stats.ranksums(va, vb)
                else:
                    p = np.nan
                temp.append((a, b, p, len(va), len(vb)))
                if np.isfinite(p):
                    p_values.append(p)
            # BH within stage + metric
            finite = np.array(p_values)
            adj = {}
            if len(finite):
                order = np.argsort(finite)
                ranked = finite[order]
                q_ranked = np.minimum.accumulate((ranked * len(finite) / (np.arange(len(finite)) + 1))[::-1])[::-1]
                q_vals = np.empty_like(q_ranked)
                q_vals[order] = np.minimum(q_ranked, 1)
                idx = 0
                for a, b, p, _, _ in temp:
                    if np.isfinite(p):
                        adj[(a, b)] = q_vals[idx]
                        idx += 1
            for a, b, p, na, nb in temp:
                q = adj.get((a, b), np.nan)
                rows.append(
                    {
                        "stage": stage,
                        "metric": metric,
                        "group1": a,
                        "group2": b,
                        "n1": na,
                        "n2": nb,
                        "wilcoxon_rank_sum_p": p,
                        "BH_q": q,
                        "significance": significance_label(q),
                    }
                )
    pd.DataFrame(rows).to_csv(out_csv, index=False)


def significance_label(q: float) -> str:
    if not np.isfinite(q):
        return "n/a"
    if q < 0.0001:
        return "****"
    if q < 0.001:
        return "***"
    if q < 0.01:
        return "**"
    if q < 0.05:
        return "*"
    return "ns"


def exclude_known_mislabeled_rows(df: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Exclude the longest 4-cell C. elegans row, which was user-confirmed as mislabeled."""
    target = (df["stage"].astype(str) == "4_Cell") & (df["species"].astype(str) == "C. elegans")
    if not target.any():
        return df.copy(), df.iloc[0:0].copy()
    idx = df.loc[target, "length_um"].astype(float).idxmax()
    excluded = df.loc[[idx]].copy()
    excluded["exclusion_reason"] = "User-confirmed mislabeled 4-cell C. elegans embryo; longest length in 4-cell C. elegans group"
    filtered = df.drop(index=idx).copy()
    return filtered, excluded


def add_n_labels(ax, data: list[np.ndarray], positions: np.ndarray, y: float) -> None:
    for pos, vals in zip(positions, data):
        ax.text(pos, y, f"n={len(vals)}", ha="center", va="bottom", fontsize=7)


def add_stat_brackets(ax, stat_df: pd.DataFrame | None, stage: str, metric: str, ylim: tuple[float, float]) -> None:
    if stat_df is None or stat_df.empty:
        return
    pairs = [
        ("C. elegans", "C. briggsae", 0.78),
        ("C. elegans", "C. nigoni", 0.86),
        ("C. briggsae", "C. nigoni", 0.70),
    ]
    species_to_x = {sp: i + 1 for i, sp in enumerate(SPECIES_ORDER)}
    y0, y1 = ylim
    yr = y1 - y0
    for group1, group2, frac in pairs:
        row = stat_df[
            (stat_df["stage"].astype(str) == stage)
            & (stat_df["metric"].astype(str) == metric)
            & (stat_df["group1"].astype(str) == group1)
            & (stat_df["group2"].astype(str) == group2)
        ]
        if row.empty:
            continue
        label = str(row.iloc[0].get("significance", ""))
        x1, x2 = species_to_x[group1], species_to_x[group2]
        y = y0 + yr * frac
        tick = yr * 0.015
        ax.plot([x1, x1, x2, x2], [y - tick, y, y, y - tick], color="black", linewidth=0.7, clip_on=False)
        ax.text((x1 + x2) / 2.0, y + yr * 0.006, label, ha="center", va="bottom", fontsize=7)


def metric_ylim(metric: str, values: pd.Series | np.ndarray) -> tuple[float, float]:
    if metric in {"length_um", "width_um", "depth_um"}:
        return 0.0, 75.0
    vals = np.asarray(values, dtype=float)
    vals = vals[np.isfinite(vals)]
    upper = float(np.nanmax(vals) * 1.18) if vals.size else 1.0
    return 0.0, upper


def plot_metric(df: pd.DataFrame, metric: str, ylabel: str, out_base: Path, stat_df: pd.DataFrame | None = None) -> None:
    fig, axes = plt.subplots(1, 3, figsize=(9.2, 3.2), sharey=True)
    rng = np.random.default_rng(42)
    ylim = metric_ylim(metric, df[metric])
    for ax, stage in zip(axes, STAGE_ORDER):
        sub = df[df["stage"] == stage]
        data = [sub.loc[sub["species"] == sp, metric].dropna().to_numpy() for sp in SPECIES_ORDER]
        positions = np.arange(1, len(SPECIES_ORDER) + 1)
        bp = ax.boxplot(data, positions=positions, widths=0.55, patch_artist=True, showfliers=False)
        for patch, sp in zip(bp["boxes"], SPECIES_ORDER):
            patch.set_facecolor(COLORS[sp])
            patch.set_alpha(0.78)
            patch.set_edgecolor("black")
        for key in ["whiskers", "caps", "medians"]:
            for artist in bp[key]:
                artist.set_color("black")
                artist.set_linewidth(1.0)
        for pos, vals, sp in zip(positions, data, SPECIES_ORDER):
            jitter = rng.normal(0, 0.055, size=len(vals))
            ax.scatter(np.full(len(vals), pos) + jitter, vals, s=12, color=COLORS[sp], edgecolor="black", linewidth=0.2, alpha=0.65)
        add_n_labels(ax, data, positions, ylim[1] * 0.94)
        add_stat_brackets(ax, stat_df, stage, metric, ylim)
        ax.set_ylim(*ylim)
        ax.set_title(stage.replace("_", "-"), fontsize=10, fontweight="bold")
        ax.set_xticks(positions)
        ax.set_xticklabels(["Ce", "Cb", "Cn"])
        ax.set_xlabel("Species")
        ax.spines["top"].set_visible(True)
        ax.spines["right"].set_visible(True)
    axes[0].set_ylabel(ylabel)
    fig.suptitle(ylabel + " by stage and species", fontsize=12, fontweight="bold")
    fig.text(0.5, 0.01, "Pairwise Wilcoxon rank-sum tests; labels show BH-adjusted significance", ha="center", fontsize=7)
    fig.tight_layout()
    fig.savefig(out_base.with_suffix(".png"), dpi=600)
    fig.savefig(out_base.with_suffix(".pdf"))
    plt.close(fig)


def plot_combined(df: pd.DataFrame, out_base: Path, stat_df: pd.DataFrame | None = None) -> None:
    metrics = [
        ("length_um", "Length (um)"),
        ("width_um", "Width (um)"),
        ("depth_um", "Depth (um)"),
        ("ellipsoid_volume_um3", "Ellipsoid volume (um^3)"),
    ]
    fig, axes = plt.subplots(4, 3, figsize=(9.2, 10.6))
    rng = np.random.default_rng(123)
    for r, (metric, ylabel) in enumerate(metrics):
        ylim = metric_ylim(metric, df[metric])
        for c, stage in enumerate(STAGE_ORDER):
            ax = axes[r, c]
            sub = df[df["stage"] == stage]
            data = [sub.loc[sub["species"] == sp, metric].dropna().to_numpy() for sp in SPECIES_ORDER]
            positions = np.arange(1, len(SPECIES_ORDER) + 1)
            bp = ax.boxplot(data, positions=positions, widths=0.55, patch_artist=True, showfliers=False)
            for patch, sp in zip(bp["boxes"], SPECIES_ORDER):
                patch.set_facecolor(COLORS[sp])
                patch.set_alpha(0.78)
                patch.set_edgecolor("black")
            for key in ["whiskers", "caps", "medians"]:
                for artist in bp[key]:
                    artist.set_color("black")
                    artist.set_linewidth(0.9)
            for pos, vals, sp in zip(positions, data, SPECIES_ORDER):
                jitter = rng.normal(0, 0.055, size=len(vals))
                ax.scatter(np.full(len(vals), pos) + jitter, vals, s=9, color=COLORS[sp], edgecolor="black", linewidth=0.15, alpha=0.65)
            add_n_labels(ax, data, positions, ylim[1] * 0.94)
            add_stat_brackets(ax, stat_df, stage, metric, ylim)
            ax.set_ylim(*ylim)
            if r == 0:
                ax.set_title(stage.replace("_", "-"), fontsize=10, fontweight="bold")
            if c == 0:
                ax.set_ylabel(ylabel)
            ax.set_xticks(positions)
            ax.set_xticklabels(["Ce", "Cb", "Cn"])
    fig.suptitle("Embryo DIC size measurements", fontsize=13, fontweight="bold")
    fig.text(0.5, 0.006, "Pairwise Wilcoxon rank-sum tests; labels show BH-adjusted significance", ha="center", fontsize=7)
    fig.tight_layout()
    fig.savefig(out_base.with_suffix(".png"), dpi=600)
    fig.savefig(out_base.with_suffix(".pdf"))
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--lif", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--dic-channel", type=int, default=2, help="Preferred DIC channel. If a series has fewer channels, the last channel is used.")
    parser.add_argument("--save-overlays", action="store_true")
    args = parser.parse_args()

    args.out.mkdir(parents=True, exist_ok=True)
    table_dir = args.out / "tables"
    fig_dir = args.out / "figures"
    code_dir = args.out / "code"
    overlay_dir = args.out / "qc_overlays" if args.save_overlays else None
    for d in [table_dir, fig_dir, code_dir, overlay_dir]:
        if d is not None:
            d.mkdir(parents=True, exist_ok=True)

    rows = []
    errors = []
    lf = LifFile(str(args.lif))
    for i, img in enumerate(lf.get_iter_image()):
        try:
            rows.append(measure_one(img, i, args.dic_channel, overlay_dir))
            print(f"measured {i}: {img.name}", flush=True)
        except Exception as exc:
            errors.append({"series_index": i, "series_name": getattr(img, "name", ""), "error": str(exc)})
            print(f"ERROR {i}: {exc}", flush=True)

    df = pd.DataFrame(rows)
    df["stage"] = pd.Categorical(df["stage"], STAGE_ORDER, ordered=True)
    df["species"] = pd.Categorical(df["species"], SPECIES_ORDER, ordered=True)
    df = df.sort_values(["stage", "species", "series_index"])
    raw_df = df.copy()
    raw_df.to_csv(table_dir / "embryo_dic_size_measurements_before_manual_exclusion.csv", index=False)
    df, excluded_df = exclude_known_mislabeled_rows(raw_df)
    df["stage"] = pd.Categorical(df["stage"], STAGE_ORDER, ordered=True)
    df["species"] = pd.Categorical(df["species"], SPECIES_ORDER, ordered=True)
    df = df.sort_values(["stage", "species", "series_index"])
    df.to_csv(table_dir / "embryo_dic_size_measurements.csv", index=False)
    excluded_df.to_csv(table_dir / "embryo_dic_size_manual_exclusions.csv", index=False)
    pd.DataFrame(errors, columns=["series_index", "series_name", "error"]).to_csv(
        table_dir / "embryo_dic_size_measurement_errors.csv", index=False
    )
    write_summary(df, table_dir / "embryo_dic_size_summary_by_stage_species.csv")
    stat_file = table_dir / "embryo_dic_size_wilcoxon_rank_sum_stats_by_stage.csv"
    write_stats(df, stat_file)
    stat_df = pd.read_csv(stat_file)

    plot_metric(df, "length_um", "Length (um)", fig_dir / "Figure1_length_by_stage_species", stat_df)
    plot_metric(df, "width_um", "Width (um)", fig_dir / "Figure2_width_by_stage_species", stat_df)
    plot_metric(df, "depth_um", "Depth (um)", fig_dir / "Figure3_depth_by_stage_species", stat_df)
    plot_metric(df, "ellipsoid_volume_um3", "Ellipsoid volume (um^3)", fig_dir / "Figure4_volume_by_stage_species", stat_df)
    plot_combined(df, fig_dir / "Figure5_size_metrics_combined", stat_df)

    readme = """# DIC embryo size measurement

Input LIF: `/Volumes/Extreme SSD/Emb_Measure.lif`

Collection parsing:
- stage = first folder in the LIF series name, for example `4_Cell`, `350_Cell`, or `550_Cell`
- species = second folder prefix: `Cel` = C. elegans, `Cbr` = C. briggsae, `Cni` = C. nigoni

Image channel:
- DIC/bright channel was manually set to preferred channel 2 after checking the channel contact-sheet example.
- If a series has only two channels, the script uses the last available channel, channel 1.

Boundary definition:
- For each series, the channel-2 z-stack was max-projected along z.
- The projection was contrast-normalized using the 1st and 99th percentiles.
- Canny edge detection was applied to the smoothed projection.
- Edges were dilated, holes filled, closed, and small objects removed.
- A large, central embryo-like component was selected to avoid small debris/background contaminants.
- A smooth ellipse was then fitted to this embryo candidate.
- The red QC line is this smooth fitted ellipse, not the raw noisy edge contour.
- Length and width are the major and minor axes of this fitted ellipse, converted to microns using LIF calibration.

Depth definition:
- A z-profile was calculated from embryo-mask DIC contrast in each z-slice.
- The embryo-containing z-range was defined as slices above an adaptive contrast threshold.
- Depth = selected z-slice count * calibrated z-step.
- `qc_flag_full_z_stack = TRUE` means the estimated embryo z-range spans the full acquired stack and should be checked manually.

Volume:
- Volume is ellipsoid-based: `pi / 6 * length * width * depth`.

Statistics:
- The stats table uses pairwise Wilcoxon rank-sum tests for species comparisons within each stage and metric, with BH correction within each stage/metric.
- Figure brackets show BH-adjusted significance: ns, *, **, ***.

Manual exclusion:
- The longest 4-cell C. elegans row was excluded because the user confirmed this embryo was mislabeled and should be another species. The excluded row is saved in `tables/embryo_dic_size_manual_exclusions.csv`.

Plot display:
- n is labeled for each stage/species group.
- Length, width, and depth panels use y-axis range 0-75 um.
- Volume panels start from 0.
"""
    (args.out / "README.md").write_text(readme)


if __name__ == "__main__":
    main()

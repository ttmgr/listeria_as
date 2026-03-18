#!/usr/bin/env python3
"""
Build standalone local HTML reports for cohort-level AS-vs-N analysis.

Inputs are read from a downloaded project folder, e.g.:
    python3 build_local_black_report.py /path/to/downloaded_results/codex

Outputs:
    <base_dir>/local_black_report/report.html
    <base_dir>/local_blue_report/report.html
"""

from __future__ import annotations

import html
import math
import os
import re
import sys
from pathlib import Path
from typing import Dict, Iterable, List, Tuple

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


GROUP_SUFFIX_ORDER = ["A", "C", "D"]
GROUP_SUFFIX_COLORS = {
    "A": "#3b82f6",
    "C": "#f59e0b",
    "D": "#8b5cf6",
}
FALLBACK_GROUP_COLORS = ["#64748b", "#0f766e", "#b45309", "#be123c"]
COND_COLORS = {"AS": "#2d8a4e", "N": "#b03060"}
COND_FILLS = {"AS": "#2d8a4e55", "N": "#b0306055"}


def ensure_dir(path: Path) -> None:
    path.mkdir(parents=True, exist_ok=True)


def savefig(fig: plt.Figure, path: Path) -> str:
    fig.savefig(path, dpi=220, bbox_inches="tight")
    plt.close(fig)
    return path.name


def sample_sort_key(sample: str) -> Tuple[int, str]:
    m = re.search(r"(?:r\d+_)?barcode(\d+)_(AS|N)", sample)
    if not m:
        return (999, sample)
    return (int(m.group(1)), m.group(2))


def barcode_num(sample_or_barcode: str) -> int:
    m = re.search(r"(?:r\d+_)?barcode(\d+)", str(sample_or_barcode))
    return int(m.group(1)) if m else -1


def cohort_slug(cohort: str) -> str:
    return re.sub(r"[^a-z0-9]+", "_", str(cohort).strip().lower()).strip("_")


def group_suffix(group_name: str) -> str:
    value = str(group_name)
    return value.split("_")[-1] if "_" in value else value


def group_sort_key(group_name: str) -> Tuple[int, str]:
    suffix = group_suffix(group_name)
    try:
        idx = GROUP_SUFFIX_ORDER.index(suffix)
    except ValueError:
        idx = len(GROUP_SUFFIX_ORDER)
    return (idx, str(group_name))


def build_group_labels(meta: pd.DataFrame) -> Dict[str, str]:
    labels: Dict[str, str] = {}
    for group_name, sub in meta.groupby("group"):
        swabs = [str(v) for v in sub["swab_type"].dropna().unique() if str(v).strip()]
        kits = [str(v) for v in sub["kit"].dropna().unique() if str(v).strip()]
        parts = []
        if swabs:
            parts.append("/".join(swabs))
        if kits:
            parts.append("/".join(kits))
        labels[str(group_name)] = " / ".join(parts) if parts else str(group_name)
    return labels


def build_group_colors(groups: Iterable[str]) -> Dict[str, str]:
    colors: Dict[str, str] = {}
    extra_idx = 0
    for group_name in sorted({str(g) for g in groups}, key=group_sort_key):
        suffix = group_suffix(group_name)
        if suffix in GROUP_SUFFIX_COLORS:
            colors[group_name] = GROUP_SUFFIX_COLORS[suffix]
        else:
            colors[group_name] = FALLBACK_GROUP_COLORS[extra_idx % len(FALLBACK_GROUP_COLORS)]
            extra_idx += 1
    return colors


def species_level_taxon(taxon: str) -> str:
    value = str(taxon).strip()
    if not value:
        return value
    parts = value.split()
    if len(parts) >= 3 and parts[0] == "Candidatus":
        return " ".join(parts[:3])
    if len(parts) >= 2:
        return " ".join(parts[:2])
    return value


def weighted_quantile(lengths: np.ndarray, counts: np.ndarray, q: float) -> float:
    if lengths.size == 0 or counts.sum() == 0:
        return 0.0
    order = np.argsort(lengths)
    lengths = lengths[order]
    counts = counts[order]
    cum = np.cumsum(counts)
    cutoff = q * counts.sum()
    idx = np.searchsorted(cum, cutoff, side="left")
    return float(lengths[min(idx, len(lengths) - 1)])


def weighted_n50(lengths: np.ndarray, counts: np.ndarray) -> float:
    if lengths.size == 0 or counts.sum() == 0:
        return 0.0
    order = np.argsort(lengths)[::-1]
    lengths = lengths[order]
    counts = counts[order]
    bp = lengths * counts
    cum_bp = np.cumsum(bp)
    cutoff = bp.sum() / 2.0
    idx = np.searchsorted(cum_bp, cutoff, side="left")
    return float(lengths[min(idx, len(lengths) - 1)])


def format_int(value: float | int) -> str:
    return f"{int(round(float(value))):,}"


def format_float(value: float, digits: int = 1) -> str:
    if value is None or (isinstance(value, float) and (math.isnan(value) or math.isinf(value))):
        return "NA"
    return f"{float(value):,.{digits}f}"


def df_to_html_table(df: pd.DataFrame, max_rows: int | None = None) -> str:
    if df is None or df.empty:
        return "<p class='note'>No data available.</p>"
    if max_rows is not None:
        df = df.head(max_rows)
    headers = "".join(f"<th onclick=\"sortTable(this)\">{html.escape(str(c))}</th>" for c in df.columns)
    body_rows = []
    for _, row in df.iterrows():
        cells = []
        for val in row:
            if isinstance(val, float):
                cell = format_float(val, 2)
            elif isinstance(val, (int, np.integer)):
                cell = format_int(val)
            else:
                cell = html.escape(str(val))
            cells.append(f"<td>{cell}</td>")
        body_rows.append("<tr>" + "".join(cells) + "</tr>")
    return (
        "<div class='table-wrap'><table class='sortable'><thead><tr>"
        + headers
        + "</tr></thead><tbody>"
        + "".join(body_rows)
        + "</tbody></table></div>"
    )


def load_inputs(base_dir: Path) -> Dict[str, pd.DataFrame]:
    data: Dict[str, pd.DataFrame] = {}

    meta = pd.read_csv(base_dir / "sample_metadata.csv")
    meta = meta.dropna(subset=["barcode"])
    meta["barcode"] = meta["barcode"].astype(int)
    if "basename" in meta.columns:
        meta["sample"] = meta["basename"]
    else:
        meta["sample"] = meta.apply(
            lambda r: f"barcode{int(r['barcode']):02d}_{r['condition']}", axis=1
        )
    data["meta"] = meta

    listeria_path = base_dir / "processing" / "listeria" / "overview" / "listeria_overview.csv"
    if not listeria_path.exists():
        print(f"ERROR: Required file not found: {listeria_path}")
        raise FileNotFoundError(f"Missing {listeria_path}")
    listeria = pd.read_csv(listeria_path)
    listeria = listeria.rename(
        columns={
            "Sample": "sample",
            "Listeria Reads": "listeria_reads",
            "Listeria Bases": "listeria_bases",
            "Mean Read Length": "listeria_mean_len",
            "Median Read Length": "listeria_median_len",
            "Total Reads": "total_reads",
            "Total Bases": "total_bases",
            "Flye Contigs": "flye_contigs",
            "MetaMDBG Contigs": "mdbg_contigs",
            "Myloasm Contigs": "myloasm_contigs",
            "Listeria (%)": "listeria_pct",
            "Barcode": "barcode_label",
            "Type": "condition",
        }
    )
    listeria["barcode"] = listeria["sample"].map(barcode_num)
    data["listeria"] = listeria

    read_lengths_path = base_dir / "processing" / "read_lengths_filtered_agg_rebuilt.tsv"
    if not read_lengths_path.exists():
        read_lengths_path = base_dir / "processing" / "read_lengths_filtered_agg.tsv"
    if read_lengths_path.exists():
        read_lengths = pd.read_csv(
            read_lengths_path,
            sep="\t",
            header=None,
            names=["sample", "length", "state", "count"],
        )
        read_lengths["barcode"] = read_lengths["sample"].map(barcode_num)
        read_lengths["condition"] = read_lengths["sample"].str.extract(r"_(AS|N)$")[0]
        data["read_lengths"] = read_lengths
    else:
        print(f"WARNING: No read lengths file found, skipping read length plots")
        data["read_lengths"] = pd.DataFrame(columns=["sample", "length", "state", "count", "barcode", "condition"])

    amr_reads_path = base_dir / "processing" / "amrfinder" / "overview" / "amr_reads_overview.csv"
    amr_contigs_path = base_dir / "processing" / "amrfinder" / "overview" / "amr_contigs_overview.csv"
    data["amr_reads"] = pd.read_csv(amr_reads_path) if amr_reads_path.exists() else pd.DataFrame()
    data["amr_contigs"] = pd.read_csv(amr_contigs_path) if amr_contigs_path.exists() else pd.DataFrame()

    asm_tables = {}
    for assembler in ["flye", "mdbg", "myloasm"]:
        asm_path = base_dir / "processing" / "stats" / f"assembly_stats_{assembler}.tsv"
        if asm_path.exists():
            df = pd.read_csv(asm_path, sep="\t")
            df["sample"] = df["file"].str.extract(r"((?:r\d+_)?barcode\d+_(?:AS|N))")
            df["assembler"] = assembler
            asm_tables[assembler] = df
        else:
            print(f"WARNING: Assembly stats not found: {asm_path}")
    data["assembly"] = pd.concat(asm_tables.values(), ignore_index=True) if asm_tables else pd.DataFrame()

    contig_frames = []
    for name in ["flye", "mdbg", "myloasm"]:
        path = base_dir / "processing" / "kraken2_csv" / f"{name}_contigs_classification.csv"
        if path.exists():
            contig_frames.append(pd.read_csv(path))
    data["contig_taxa"] = pd.concat(contig_frames, ignore_index=True) if contig_frames else pd.DataFrame()
    return data


def prepare_cohort_tables(data: Dict[str, pd.DataFrame], cohort: str) -> Dict[str, pd.DataFrame]:
    meta = data["meta"][data["meta"]["cohort"].eq(cohort)].copy()
    group_labels = build_group_labels(meta)
    meta["group_label"] = meta["group"].map(group_labels).fillna(meta["group"])
    meta["Barcode"] = meta["barcode"].map(lambda x: f"barcode{int(x):02d}")
    meta["Type"] = meta["condition"]
    keep_samples = set(meta["sample"])

    listeria = data["listeria"][data["listeria"]["sample"].isin(keep_samples)].copy()
    listeria = listeria.merge(
        meta[
            [
                "sample",
                "sample_id",
                "barcode",
                "condition",
                "group",
                "group_label",
                "swab_type",
                "kit",
                "dna_concentration_ng_ul",
            ]
        ],
        on=["sample", "barcode", "condition"],
        how="left",
    )

    rl = data["read_lengths"][data["read_lengths"]["sample"].isin(keep_samples)].copy()
    rl = rl.merge(meta[["sample", "group", "group_label"]], on="sample", how="left")

    amr_reads = data["amr_reads"].copy()
    cohort_barcode_labels = set(meta["Barcode"].unique())
    amr_reads = amr_reads[amr_reads["Barcode"].isin(cohort_barcode_labels)].copy()
    amr_reads = amr_reads.merge(
        meta[["Barcode", "Type", "sample_id", "group_label"]],
        on=["Barcode", "Type"],
        how="left",
    )

    amr_contigs = data["amr_contigs"].copy()
    amr_contigs = amr_contigs[amr_contigs["Barcode"].isin(cohort_barcode_labels)].copy()
    amr_contigs = amr_contigs.merge(
        meta[["Barcode", "Type", "sample_id", "group_label"]],
        on=["Barcode", "Type"],
        how="left",
    )

    asm = data["assembly"].copy()
    asm = asm[asm["sample"].isin(keep_samples)].copy()
    asm = asm.merge(meta[["sample", "barcode", "condition", "group", "group_label"]], on="sample", how="left")

    contig_taxa = data["contig_taxa"].copy()
    if not contig_taxa.empty:
        contig_taxa = contig_taxa[contig_taxa["sample"].isin(keep_samples)].copy()
        contig_taxa = contig_taxa.merge(meta[["sample", "barcode", "condition", "group", "group_label"]], on="sample", how="left")

    return {
        "meta": meta,
        "listeria": listeria,
        "read_lengths": rl,
        "amr_reads": amr_reads,
        "amr_contigs": amr_contigs,
        "assembly": asm,
        "contig_taxa": contig_taxa,
        "cohort": cohort,
        "cohort_slug": cohort_slug(cohort),
        "group_labels": group_labels,
        "group_colors": build_group_colors(meta["group"].dropna().unique()),
        "group_order": sorted(meta["group"].dropna().unique(), key=group_sort_key),
    }


def build_read_metrics(rl: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for sample, g in rl.groupby("sample"):
        lengths = g["length"].to_numpy(dtype=float)
        counts = g["count"].to_numpy(dtype=float)
        total_reads = counts.sum()
        total_bases = (lengths * counts).sum()
        rows.append(
            {
                "sample": sample,
                "barcode": barcode_num(sample),
                "condition": sample.split("_")[-1],
                "total_reads": int(total_reads),
                "total_bases": int(total_bases),
                "mean_read_length": total_bases / total_reads if total_reads else 0.0,
                "median_read_length": weighted_quantile(lengths, counts, 0.5),
                "read_length_n50": weighted_n50(lengths, counts),
            }
        )
    return pd.DataFrame(rows).sort_values(["barcode", "condition"])


def compute_enrichment(listeria: pd.DataFrame) -> pd.DataFrame:
    rows = []
    group_cols = ["round", "barcode"] if "round" in listeria.columns else ["barcode"]
    for group_key, g in listeria.groupby(group_cols):
        barcode = group_key[-1] if isinstance(group_key, tuple) else group_key
        as_row = g[g["condition"] == "AS"]
        n_row = g[g["condition"] == "N"]
        if as_row.empty or n_row.empty:
            continue
        as_row = as_row.iloc[0]
        n_row = n_row.iloc[0]
        abs_enr = (as_row["listeria_reads"] / n_row["listeria_reads"]) if n_row["listeria_reads"] > 0 else np.nan
        rel_enr = (as_row["listeria_pct"] / n_row["listeria_pct"]) if n_row["listeria_pct"] > 0 else np.nan
        rows.append(
            {
                "barcode": int(barcode),
                "sample_id": as_row["sample_id"],
                "group": as_row["group"],
                "method": as_row["group_label"],
                "swab_type": as_row["swab_type"],
                "kit": as_row["kit"],
                "dna_ng_ul": as_row["dna_concentration_ng_ul"],
                "listeria_AS": int(as_row["listeria_reads"]),
                "listeria_N": int(n_row["listeria_reads"]),
                "pct_AS": float(as_row["listeria_pct"]),
                "pct_N": float(n_row["listeria_pct"]),
                "total_reads_AS": int(as_row["total_reads"]),
                "total_reads_N": int(n_row["total_reads"]),
                "absolute_enrichment": abs_enr,
                "relative_enrichment": rel_enr,
            }
        )
    if not rows:
        return pd.DataFrame(columns=[
            "barcode", "sample_id", "group", "method", "swab_type", "kit",
            "dna_ng_ul", "listeria_AS", "listeria_N", "pct_AS", "pct_N",
            "total_reads_AS", "total_reads_N", "absolute_enrichment", "relative_enrichment",
        ])
    return pd.DataFrame(rows).sort_values("barcode")


def build_general_metrics(listeria: pd.DataFrame, read_metrics: pd.DataFrame, cohort: str) -> pd.DataFrame:
    lengths = read_metrics.copy()
    total_reads = int(lengths["total_reads"].sum())
    total_bases = int(lengths["total_bases"].sum())
    mean_len = total_bases / total_reads if total_reads else 0
    median_len = lengths["median_read_length"].median() if not lengths.empty else 0
    mean_n50 = lengths["read_length_n50"].mean() if not lengths.empty else 0
    listeria_reads = int(listeria["listeria_reads"].sum())
    listeria_pct = float((listeria["listeria_reads"].sum() / listeria["total_reads"].sum()) * 100) if listeria["total_reads"].sum() else 0
    return pd.DataFrame(
        [
            (f"Samples in {cohort} cohort", f"{listeria['barcode'].nunique()} barcodes / {len(listeria)} libraries"),
            ("Total filtered reads", format_int(total_reads)),
            ("Total filtered bases", format_int(total_bases)),
            ("Weighted mean read length", f"{format_float(mean_len, 1)} bp"),
            ("Median per-sample read length", f"{format_float(median_len, 1)} bp"),
            ("Mean per-sample read N50", f"{format_float(mean_n50, 1)} bp"),
            ("Total Listeria reads", format_int(listeria_reads)),
            ("Overall Listeria fraction", f"{format_float(listeria_pct, 3)}%"),
        ],
        columns=["Metric", "Value"],
    )


def section_plot_read_lengths_all(rl: pd.DataFrame, out_assets: Path, cohort: str) -> Tuple[List[str], pd.DataFrame]:
    files = []
    max_len = int(rl["length"].max())
    bins = np.logspace(np.log10(50), np.log10(max(max_len, 100000)), 101)

    fig, ax = plt.subplots(figsize=(8.5, 4.2))
    for cond in ["N", "AS"]:
        sub = rl[rl["condition"] == cond]
        ax.hist(
            sub["length"],
            bins=bins,
            weights=sub["count"],
            histtype="stepfilled",
            color=COND_FILLS[cond],
            edgecolor=COND_COLORS[cond],
            linewidth=0.9,
            alpha=0.55,
            label=cond,
        )
    ax.set_xscale("log")
    ax.set_xlabel("Read length (bp)")
    ax.set_ylabel("Read count")
    ax.set_title(f"Filtered read length distribution across {cohort} cohort")
    ax.set_xticks([100, 1000, 10000, 100000])
    ax.set_xticklabels(["100", "1k", "10k", "100k"])
    ax.axvline(400, linestyle="--", color="#6b7280", linewidth=1.0)
    ax.legend(frameon=False)
    files.append(savefig(fig, out_assets / "section2_read_length_distribution.png"))

    fig, axes = plt.subplots(1, 2, figsize=(12.8, 4.2), sharex=True)
    for cond in ["N", "AS"]:
        sub = rl[rl["condition"] == cond]
        axes[0].hist(
            sub["length"],
            bins=bins,
            weights=sub["count"],
            histtype="stepfilled",
            color=COND_FILLS[cond],
            edgecolor=COND_COLORS[cond],
            linewidth=1.2,
            alpha=0.55,
            label=cond,
        )
        axes[1].hist(
            sub["length"],
            bins=bins,
            weights=sub["count"] * sub["length"],
            histtype="stepfilled",
            color=COND_FILLS[cond],
            edgecolor=COND_COLORS[cond],
            linewidth=1.2,
            alpha=0.55,
            label=cond,
        )
    axes[0].set_title(f"All {cohort} barcodes pooled: read counts")
    axes[0].set_ylabel("Read count")
    axes[1].set_title(f"All {cohort} barcodes pooled: base yield")
    axes[1].set_ylabel("Bases per bin")
    for ax in axes:
        ax.set_xscale("log")
        ax.set_xlabel("Read length (bp)")
        ax.set_xticks([100, 1000, 10000, 100000])
        ax.set_xticklabels(["100", "1k", "10k", "100k"])
        ax.axvline(400, linestyle="--", color="#6b7280", linewidth=1.0)
        ax.legend(frameon=False)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
    files.append(savefig(fig, out_assets / "section2_all_barcodes_as_vs_n.png"))

    summary = []
    for cond, g in rl.groupby("condition"):
        lengths = g["length"].to_numpy(dtype=float)
        counts = g["count"].to_numpy(dtype=float)
        total_reads = counts.sum()
        total_bases = (lengths * counts).sum()
        summary.append(
            {
                "Condition": cond,
                "Reads": int(total_reads),
                "Bases": int(total_bases),
                "Mean Length (bp)": round(total_bases / total_reads, 1) if total_reads else 0,
                "Median Length (bp)": round(weighted_quantile(lengths, counts, 0.5), 1),
                "Read N50 (bp)": round(weighted_n50(lengths, counts), 1),
            }
        )
    summary_df = pd.DataFrame(summary).sort_values("Condition")

    return files, summary_df


def section_plot_read_lengths_by_method(
    rl: pd.DataFrame,
    out_assets: Path,
    cohort: str,
    group_order: List[str],
    group_labels: Dict[str, str],
) -> str:
    methods = [g for g in group_order if g in set(rl["group"].dropna())]
    fig, axes = plt.subplots(1, len(methods), figsize=(5.0 * max(len(methods), 1), 4.2), sharey=True)
    axes = np.atleast_1d(axes).ravel()
    max_len = int(rl["length"].max())
    bins = np.logspace(np.log10(50), np.log10(max(max_len, 100000)), 101)

    for ax, group_name in zip(axes, methods):
        sub = rl[rl["group"] == group_name]
        for cond in ["N", "AS"]:
            cond_sub = sub[sub["condition"] == cond]
            ax.hist(
                cond_sub["length"],
                bins=bins,
                weights=cond_sub["count"],
                histtype="stepfilled",
                color=COND_FILLS[cond],
                edgecolor=COND_COLORS[cond],
                linewidth=1.0,
                alpha=0.5,
                label=cond,
            )
        ax.set_xscale("log")
        ax.set_xticks([100, 1000, 10000, 100000])
        ax.set_xticklabels(["100", "1k", "10k", "100k"])
        ax.axvline(400, linestyle="--", color="#6b7280", linewidth=1.0)
        ax.set_xlabel("Read length (bp)")
        ax.set_title(group_labels.get(group_name, group_name))
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
    axes[0].set_ylabel("Read count")
    axes[0].legend(frameon=False)
    fig.suptitle(f"Section 3b. Read Lengths per Method ({cohort})", fontsize=12, y=1.02)
    return savefig(fig, out_assets / "section3b_read_length_by_method.png")


def section_plot_method_comparison(
    rl: pd.DataFrame,
    out_assets: Path,
    cohort: str,
    group_order: List[str],
    group_labels: Dict[str, str],
    group_colors: Dict[str, str],
) -> str:
    methods = [g for g in group_order if g in set(rl["group"].dropna())]
    fig, axes = plt.subplots(1, 2, figsize=(12.8, 4.4), sharex=True, sharey=True)
    max_len = int(rl["length"].max())
    bins = np.logspace(np.log10(50), np.log10(max(max_len, 100000)), 101)

    for ax, cond in zip(axes, ["N", "AS"]):
        for group_name in methods:
            sub = rl[(rl["group"] == group_name) & (rl["condition"] == cond)]
            ax.hist(
                sub["length"],
                bins=bins,
                weights=sub["count"],
                histtype="stepfilled",
                color=group_colors.get(group_name, "#64748b"),
                edgecolor=group_colors.get(group_name, "#64748b"),
                linewidth=1.0,
                alpha=0.28,
                label=group_labels.get(group_name, group_name),
            )
        ax.set_xscale("log")
        ax.set_xticks([100, 1000, 10000, 100000])
        ax.set_xticklabels(["100", "1k", "10k", "100k"])
        ax.axvline(400, linestyle="--", color="#6b7280", linewidth=1.0)
        ax.set_xlabel("Read length (bp)")
        ax.set_title("Normal" if cond == "N" else "Adaptive Sampling")
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
    axes[0].set_ylabel("Read count")
    axes[0].legend(frameon=False)
    fig.suptitle(f"Section 3c. Method Comparison ({cohort})", fontsize=12, y=1.02)
    return savefig(fig, out_assets / "section3c_method_comparison.png")


def section_plot_read_lengths_per_barcode(rl: pd.DataFrame, listeria: pd.DataFrame, out_assets: Path, cohort: str) -> str:
    barcodes = sorted(listeria["barcode"].unique())
    ncols = 4
    nrows = int(math.ceil(len(barcodes) / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(15, 3.0 * nrows), sharex=True, sharey=True)
    axes = np.atleast_1d(axes).ravel()
    max_len = int(rl["length"].max())
    bins = np.logspace(np.log10(50), np.log10(max(max_len, 100000)), 80)

    label_map = (
        listeria.drop_duplicates("barcode")
        .set_index("barcode")[["sample_id", "group_label"]]
        .to_dict("index")
    )

    for ax, barcode in zip(axes, barcodes):
        sub = rl[rl["barcode"] == barcode]
        for cond in ["N", "AS"]:
            g = sub[sub["condition"] == cond]
            ax.hist(
                g["length"],
                bins=bins,
                weights=g["count"],
                histtype="stepfilled",
                color=COND_FILLS[cond],
                edgecolor=COND_COLORS[cond],
                linewidth=0.8,
                alpha=0.5,
                label=cond,
            )
        meta = label_map.get(barcode, {})
        ax.set_title(f"bc{barcode:02d} • {meta.get('group_label', '')}", fontsize=8)
        ax.set_xscale("log")
        ax.set_xticks([100, 1000, 10000, 100000])
        ax.set_xticklabels(["100", "1k", "10k", "100k"], fontsize=7)
        ax.axvline(400, linestyle="--", color="#6b7280", linewidth=0.9)
        ax.tick_params(axis="y", labelsize=7)
    for ax in axes[len(barcodes):]:
        ax.axis("off")
    axes[0].legend(frameon=False, fontsize=8)
    fig.suptitle(f"Read length distributions per barcode ({cohort}: AS vs N)", fontsize=12, y=1.01)
    return savefig(fig, out_assets / "section3_read_lengths_per_barcode.png")


def section_plot_listeria_overview(listeria: pd.DataFrame, out_assets: Path, group_colors: Dict[str, str], cohort: str) -> str:
    df = listeria.sort_values(["group", "barcode", "condition"]).copy()
    x = np.arange(len(df))
    non_listeria = df["total_reads"] - df["listeria_reads"]
    fig, ax1 = plt.subplots(figsize=(14, 4.8))
    colors = [group_colors.get(g, "#64748b") for g in df["group"]]
    ax1.bar(x, non_listeria / df["total_reads"] * 100, color="#dbe2ea", edgecolor="white", linewidth=0.4)
    ax1.bar(x, df["listeria_reads"] / df["total_reads"] * 100, bottom=non_listeria / df["total_reads"] * 100,
            color=colors, edgecolor="white", linewidth=0.4)
    ax1.set_ylabel("Share of total reads (%)")
    ax1.set_ylim(0, 100)
    ax1.set_title(f"Listeria-detection overview: {cohort} cohort")
    ax1.set_xticks(x)
    ax1.set_xticklabels(df["sample"].tolist(), rotation=90, fontsize=7)
    return savefig(fig, out_assets / "section4_listeria_overview.png")


def section_plot_extraction_comparison(
    listeria: pd.DataFrame,
    enrichment: pd.DataFrame,
    out_assets: Path,
    group_order: List[str],
    group_labels: Dict[str, str],
) -> Tuple[List[str], pd.DataFrame]:
    files = []

    method_summary = (
        listeria.groupby(["group", "group_label", "condition"], as_index=False)
        .agg(
            libraries=("sample", "count"),
            total_reads=("total_reads", "sum"),
            listeria_reads=("listeria_reads", "sum"),
        )
    )
    method_summary["pooled_listeria_pct"] = np.where(
        method_summary["total_reads"] > 0,
        method_summary["listeria_reads"] / method_summary["total_reads"] * 100.0,
        0.0,
    )
    method_summary["group"] = pd.Categorical(method_summary["group"], categories=group_order, ordered=True)
    method_summary = method_summary.sort_values(["group", "condition"]).copy()

    x = np.arange(len(group_order))
    width = 0.36
    fig, axes = plt.subplots(1, 3, figsize=(15.0, 4.4), sharex=True)
    specs = [
        ("total_reads", "Total filtered reads", False),
        ("listeria_reads", "Listeria-classified reads", True),
        ("pooled_listeria_pct", "Listeria fraction (%)", False),
    ]
    for ax, (metric, title, log_scale) in zip(axes, specs):
        n_vals = []
        as_vals = []
        for group in group_order:
            sub = method_summary[method_summary["group"] == group]
            n_vals.append(float(sub.loc[sub["condition"] == "N", metric].sum()))
            as_vals.append(float(sub.loc[sub["condition"] == "AS", metric].sum()))
        ax.bar(x - width / 2, n_vals, width, color=COND_FILLS["N"], edgecolor=COND_COLORS["N"], label="N")
        ax.bar(x + width / 2, as_vals, width, color=COND_FILLS["AS"], edgecolor=COND_COLORS["AS"], label="AS")
        ax.set_title(title)
        ax.set_xticks(x)
        ax.set_xticklabels([group_labels[g] for g in group_order], rotation=18, ha="right")
        if log_scale:
            ax.set_yscale("log")
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
    axes[0].legend(frameon=False)
    files.append(savefig(fig, out_assets / "section6_method_comparison.png"))
    method_summary = method_summary.rename(
        columns={
            "group_label": "method",
            "libraries": "libraries",
            "total_reads": "total_reads",
            "listeria_reads": "listeria_reads",
            "pooled_listeria_pct": "pooled_listeria_pct",
        }
    )[
        ["method", "condition", "libraries", "total_reads", "listeria_reads", "pooled_listeria_pct"]
    ]
    return files, method_summary


def section_plot_contigs(listeria: pd.DataFrame, asm: pd.DataFrame, contig_taxa: pd.DataFrame, out_assets: Path) -> Tuple[List[str], pd.DataFrame, pd.DataFrame]:
    files = []
    fig, axes = plt.subplots(2, 2, figsize=(14.5, 9.0), sharex=False)
    axes = axes.ravel()
    assemblers = ["flye", "mdbg", "myloasm"]
    assembler_labels = {"flye": "Flye", "mdbg": "MetaMDBG", "myloasm": "Myloasm"}
    metric_specs = [
        ("num_seqs", "Total contigs"),
        ("sum_len", "Assembly bases"),
        ("N50", "N50"),
        ("Q2", "Median contig length"),
    ]
    for ax, (metric, title) in zip(axes, metric_specs):
        legend_handles = []
        for idx, assembler in enumerate(assemblers):
            sub = asm[asm["assembler"] == assembler].copy()
            for cond, offset, face, edge in [
                ("N", -0.18, COND_FILLS["N"], COND_COLORS["N"]),
                ("AS", 0.18, COND_FILLS["AS"], COND_COLORS["AS"]),
            ]:
                vals = sub.loc[sub["condition"] == cond, metric].dropna().values
                plot_vals = vals[vals > 0]
                if len(plot_vals) == 0:
                    continue
                bp = ax.boxplot(
                    [plot_vals],
                    positions=[idx + offset],
                    widths=0.28,
                    patch_artist=True,
                    showfliers=False,
                    medianprops=dict(color=edge, linewidth=1.6),
                    boxprops=dict(facecolor=face, edgecolor=edge, linewidth=1.2),
                    whiskerprops=dict(color=edge, linewidth=1.0),
                    capprops=dict(color=edge, linewidth=1.0),
                )
                if len(plot_vals) == 1:
                    xvals = np.array([idx + offset])
                else:
                    xvals = np.linspace(idx + offset - 0.045, idx + offset + 0.045, len(plot_vals))
                ax.scatter(
                    xvals,
                    plot_vals,
                    s=18,
                    color=edge,
                    alpha=0.7,
                    linewidths=0.3,
                    edgecolors="white",
                    zorder=3,
                )
                if idx == 0:
                    legend_handles.append(bp["boxes"][0])
        ax.set_xticks(range(len(assemblers)))
        ax.set_xticklabels([assembler_labels[a] for a in assemblers])
        ax.set_title(title)
        ax.set_yscale("log", base=10)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        if legend_handles and metric == "num_seqs":
            ax.legend(legend_handles, ["N", "AS"], frameon=False)
    fig.suptitle("Section 7. Assembly metric comparison across assemblers", fontsize=12, y=1.01)
    files.append(savefig(fig, out_assets / "section7_assembly_metric_comparison.png"))

    # Listeria contig counts per barcode
    df_ctg = listeria[["barcode", "condition", "flye_contigs", "mdbg_contigs", "myloasm_contigs"]].copy()
    barcodes = sorted(df_ctg["barcode"].unique())
    fig, axes = plt.subplots(3, 1, figsize=(13, 8.5), sharex=True)
    specs = [("flye_contigs", "Flye"), ("mdbg_contigs", "MetaMDBG"), ("myloasm_contigs", "Myloasm")]
    for ax, (col, title) in zip(axes, specs):
        x = np.arange(len(barcodes))
        n_vals = []
        as_vals = []
        for barcode in barcodes:
            sub = df_ctg[df_ctg["barcode"] == barcode]
            n_vals.append(float(sub.loc[sub["condition"] == "N", col].iloc[0]))
            as_vals.append(float(sub.loc[sub["condition"] == "AS", col].iloc[0]))
        ax.bar(x - 0.18, n_vals, 0.36, color=COND_FILLS["N"], edgecolor=COND_COLORS["N"], label="N")
        ax.bar(x + 0.18, as_vals, 0.36, color=COND_FILLS["AS"], edgecolor=COND_COLORS["AS"], label="AS")
        ax.set_ylabel("Listeria contigs")
        ax.set_title(title)
    axes[-1].set_xticks(np.arange(len(barcodes)))
    axes[-1].set_xticklabels([f"bc{b:02d}" for b in barcodes], rotation=90, fontsize=7)
    files.append(savefig(fig, out_assets / "section7_listeria_contig_counts.png"))

    asm_summary = (
        asm.groupby(["assembler", "condition"], as_index=False)
        .agg(
            samples=("sample", "count"),
            median_total_contigs=("num_seqs", "median"),
            median_assembly_bases=("sum_len", "median"),
            median_n50=("N50", "median"),
            median_contig_length=("Q2", "median"),
        )
    )

    top_taxa = pd.DataFrame()
    if not contig_taxa.empty:
        taxa = contig_taxa.copy()
        taxa["taxon"] = taxa["taxon"].map(species_level_taxon)
        top_taxa = (
            taxa.groupby(["assembler", "sample", "taxon"], as_index=False)
            .size()
            .rename(columns={"size": "contig_count"})
            .sort_values(["assembler", "sample", "contig_count"], ascending=[True, True, False])
            .groupby(["assembler", "sample"], as_index=False)
            .head(3)
            .drop_duplicates()
        )
    return files, asm_summary, top_taxa


def build_amr_summaries(amr_reads: pd.DataFrame, amr_contigs: pd.DataFrame) -> Tuple[pd.DataFrame, pd.DataFrame]:
    def barcode_sort_series(series: pd.Series) -> pd.Series:
        return series.map(barcode_num)

    reads_summary = (
        amr_reads.groupby(
            ["Barcode", "sample_id", "group_label", "Type", "Gene Symbol", "Class", "Subclass"],
            as_index=False,
        )["Read Count"]
        .sum()
        .sort_values(
            by=["Barcode", "Type", "Read Count", "Gene Symbol"],
            ascending=[True, True, False, True],
            key=lambda col: barcode_sort_series(col) if col.name == "Barcode" else col,
        )
    )
    contigs_summary = (
        amr_contigs.groupby(
            ["Assembler", "Barcode", "sample_id", "group_label", "Type", "Gene Symbol", "Class", "Subclass"],
            as_index=False,
        )["Contig Count"]
        .sum()
        .sort_values(
            by=["Assembler", "Barcode", "Type", "Contig Count", "Gene Symbol"],
            ascending=[True, True, True, False, True],
            key=lambda col: barcode_sort_series(col) if col.name == "Barcode" else col,
        )
    )
    return reads_summary, contigs_summary


def build_barcode_metrics_export(read_metrics: pd.DataFrame, listeria: pd.DataFrame) -> pd.DataFrame:
    merged = read_metrics.merge(
        listeria[
            [
                "sample",
                "barcode",
                "condition",
                "listeria_reads",
                "listeria_pct",
            ]
        ],
        on=["sample", "barcode", "condition"],
        how="left",
    ).copy()
    merged["barcode"] = merged["barcode"].map(lambda x: f"barcode{int(x):02d}")
    export = merged.rename(
        columns={
            "barcode": "Barcode",
            "condition": "Condition",
            "total_bases": "Total Bases",
            "total_reads": "Total Reads",
            "median_read_length": "Median Read Length",
            "read_length_n50": "Read N50",
            "listeria_reads": "Listeria Reads",
            "listeria_pct": "Listeria % of Total Reads",
        }
    )[
        [
            "Barcode",
            "Condition",
            "Total Bases",
            "Total Reads",
            "Median Read Length",
            "Read N50",
            "Listeria Reads",
            "Listeria % of Total Reads",
        ]
    ].sort_values(["Barcode", "Condition"])
    return export


def write_html(
    out_html: Path,
    cohort: str,
    general_metrics: pd.DataFrame,
    section2_summary: pd.DataFrame,
    section2_plots: List[str],
    section3_plot: str,
    section3b_plot: str,
    section3c_plot: str,
    section3_table: pd.DataFrame,
    section4_plot: str,
    section4_table: pd.DataFrame,
    section6_plots: List[str],
    section6_method_summary: pd.DataFrame,
    section7_plots: List[str],
    section7_asm_summary: pd.DataFrame,
    section7_taxa: pd.DataFrame,
    amr_reads_summary: pd.DataFrame,
    amr_contigs_summary: pd.DataFrame,
) -> None:
    def imgs(files: Iterable[str]) -> str:
        return "".join(f"<img src='assets/{html.escape(f)}' alt='{html.escape(f)}'>" for f in files)

    cohort_metric_label = f"Samples in {cohort} cohort"

    body = f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <title>{cohort} Cohort AS vs N Report</title>
  <style>
    body {{ font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif; margin: 0; padding: 0; background: #f8fafc; color: #334155; }}
    .container {{ max-width: 1280px; margin: 0 auto; padding: 20px; background: #fff; box-shadow: 0 4px 6px rgba(0,0,0,0.1); }}
    h1, h2, h3 {{ color: #1e293b; border-bottom: 2px solid #e2e8f0; padding-bottom: 10px; margin-top: 2rem; }}
    .header {{ text-align: center; padding: 2rem; background: linear-gradient(135deg, #1e293b 0%, #334155 100%); color: white; border-radius: 8px; margin-bottom: 2rem; }}
    .header h1 {{ color: white; border: none; margin: 0; }}
    .header p {{ color: #cbd5e1; margin: 0.5rem 0 0; }}
    .summary-stats {{ display: flex; gap: 15px; flex-wrap: wrap; justify-content: center; margin: 1.5rem 0 0; }}
    .stat-card {{ flex: 1; min-width: 130px; max-width: 190px; background: #fff; border: 1px solid #e2e8f0; border-radius: 8px; padding: 1rem; text-align: center; box-shadow: 0 1px 3px rgba(0,0,0,0.1); }}
    .stat-val {{ font-size: 1.6rem; font-weight: bold; color: #2563eb; }}
    .stat-label {{ color: #64748b; font-size: 0.75rem; text-transform: uppercase; letter-spacing: 0.05em; margin-top: 4px; }}
    .section {{ margin: 2rem 0; }}
    .note {{ color: #64748b; }}
    .intro {{ background: #f8fafc; border: 1px solid #e2e8f0; border-radius: 8px; padding: 18px 20px; margin: 0 0 2rem; }}
    .intro p {{ margin: 0 0 0.8rem; line-height: 1.6; }}
    .toc {{ display: flex; flex-wrap: wrap; gap: 10px; margin-top: 0.8rem; }}
    .toc a {{ display: inline-block; text-decoration: none; color: #1d4ed8; background: #eff6ff; border: 1px solid #bfdbfe; border-radius: 999px; padding: 6px 10px; font-size: 0.9rem; }}
    .toc a:hover {{ background: #dbeafe; }}
    .plot-grid {{ display: grid; grid-template-columns: 1fr; gap: 1.4rem; margin: 1.4rem 0; }}
    .plot-grid img {{ max-width: 100%; height: auto; border-radius: 8px; border: 1px solid #e2e8f0; }}
    table {{ width: 100%; border-collapse: collapse; margin: 1rem 0; font-size: 0.85rem; }}
    th, td {{ padding: 10px 12px; text-align: left; border-bottom: 1px solid #e2e8f0; vertical-align: top; }}
    th {{ background: #f1f5f9; font-weight: 600; color: #475569; cursor: pointer; user-select: none; }}
    tr:hover {{ background: #f8fafc; }}
    .table-wrap {{ overflow-x: auto; }}
    @media (max-width: 1100px) {{
      .summary-stats {{ justify-content: stretch; }}
      .stat-card {{ min-width: 45%; }}
    }}
  </style>
</head>
<body>
  <div class="container">
    <div class="header">
      <h1>{cohort} Cohort Report</h1>
      <p>Adaptive Sampling vs Normal sequencing benchmark for the {cohort} cohort.</p>
      <div class="summary-stats">
        <div class="stat-card"><div class="stat-val">{general_metrics.loc[general_metrics['Metric'].eq(cohort_metric_label), 'Value'].iloc[0].split(' barcodes')[0]}</div><div class="stat-label">Barcodes</div></div>
        <div class="stat-card"><div class="stat-val">{general_metrics.loc[general_metrics['Metric'].eq('Total filtered reads'), 'Value'].iloc[0]}</div><div class="stat-label">Filtered Reads</div></div>
        <div class="stat-card"><div class="stat-val">{general_metrics.loc[general_metrics['Metric'].eq('Total Listeria reads'), 'Value'].iloc[0]}</div><div class="stat-label">Listeria Reads</div></div>
        <div class="stat-card"><div class="stat-val">{general_metrics.loc[general_metrics['Metric'].eq('Overall Listeria fraction'), 'Value'].iloc[0]}</div><div class="stat-label">Listeria Fraction</div></div>
      </div>
    </div>

    <div class="intro">
      <p>This report compares Oxford Nanopore <strong>Adaptive Sampling</strong> and <strong>Normal</strong> sequencing within the <strong>{cohort}</strong> cohort. The underlying pipeline converted BAM files to FASTQ, trimmed adapters, filtered reads by minimum length, summarized read metrics, taxonomically classified reads with Kraken2, quantified <em>Listeria</em>-assigned reads, assembled the filtered data with Flye, MetaMDBG, and Myloasm, and summarized AMR hits with AMRFinder.</p>
      <p>The sections below move from global sequencing output to barcode-level read distributions, extraction-method comparisons, <em>Listeria</em> detection summaries, and final assembly/contig outputs. All tables are sortable by clicking the column headers.</p>
      <div class="toc">
        <a href="#section-1">1. Sequencing Overview</a>
        <a href="#section-2">2. Cohort Read Summary</a>
        <a href="#section-3">3. Barcode Read Metrics</a>
        <a href="#section-3b">3b. Read Lengths by Method</a>
        <a href="#section-3c">3c. Method Comparison</a>
        <a href="#section-4">4. Listeria Detection</a>
        <a href="#section-6">6. Extraction Comparison</a>
        <a href="#section-7">7. Assembly and Contigs</a>
      </div>
    </div>

    <div class="section" id="section-1">
      <h2>1. Sequencing Overview</h2>
      <p class="note">This section is the fastest high-level summary of the cohort. It tells you how much filtered data is present overall and how much of it was assigned to <em>Listeria</em>.</p>
      <p class="note">Use this section first to judge scale before looking at per-barcode or per-method differences.</p>
      {df_to_html_table(general_metrics)}
    </div>

    <div class="section" id="section-2">
      <h2>2. Cohort-Level Read Summary</h2>
      <p class="note">All {cohort}-cohort barcodes are pooled here so the AS and N read-length profiles can be compared directly across the full cohort. The dashed line marks 400 bp.</p>
      <p class="note">If one condition has a visibly larger area under the curve, it contributed more reads in that length range. The summary table below gives the same comparison as cohort-level totals.</p>
      <div class="plot-grid">{imgs(section2_plots)}</div>
      <h3>Aggregated Summary</h3>
      {df_to_html_table(section2_summary)}
    </div>

    <div class="section" id="section-3">
      <h2>3. Barcode-Level Read Metrics</h2>
      <p class="note">This section breaks the cohort back down to individual barcode pairs. Each panel compares the N and AS read-length profiles for one barcode only.</p>
      <p class="note">The table is useful if you want exact values rather than the visual distribution, especially for total bases, total reads, median read length, and read N50.</p>
      <div class="plot-grid"><img src="assets/{section3_plot}" alt="Per-barcode read length distributions"></div>
      <h3>Per-Sample Read Metrics</h3>
      {df_to_html_table(section3_table)}
    </div>

    <div class="section" id="section-3b">
      <h2>3b. Read Lengths by Extraction Method</h2>
      <p class="note">Three method panels show Sponge / Mini, Cotton / Mini, and Zymo / Mini separately, each with N and AS overlaid on the same axes.</p>
      <p class="note">This view answers a different question than Section 3: instead of barcode-to-barcode variation, it shows whether one extraction method tends to produce a different read-length profile overall.</p>
      <div class="plot-grid"><img src="assets/{section3b_plot}" alt="Read lengths per method"></div>
    </div>

    <div class="section" id="section-3c">
      <h2>3c. Method-Level Read Length Comparison</h2>
      <p class="note">Two pooled comparison panels split by sequencing mode. Within each panel, methods A, C, and D are overlaid directly.</p>
      <p class="note">Read this section left-to-right: first compare methods within Normal, then compare methods within Adaptive Sampling. This makes method effects easier to see without the AS/N overlay from Section 3b.</p>
      <div class="plot-grid"><img src="assets/{section3c_plot}" alt="Method comparison"></div>
    </div>

    <div class="section" id="section-4">
      <h2>4. Listeria Detection Overview</h2>
      <p class="note">Bars show each library as 100% of its total reads, with the Listeria-classified fraction highlighted.</p>
      <p class="note">This section is about enrichment visibility, not absolute sequencing yield. A small highlighted fraction means <em>Listeria</em> makes up only a small part of that sample even if the total library is large.</p>
      <p class="note">If you need the exact counts and percentages, use the sortable table directly below the figure.</p>
      <div class="plot-grid"><img src="assets/{section4_plot}" alt="Listeria overview"></div>
      <h3>Overview Table</h3>
      {df_to_html_table(section4_table)}
    </div>

    <div class="section" id="section-6">
      <h2>6. Extraction Method Comparison</h2>
      <p class="note">Direct comparison of the three {cohort}-cohort extraction methods: Sponge / Mini, Cotton / Mini, and Zymo / Mini.</p>
      <p class="note">This section collapses barcode-level noise and shows how the methods behave as groups. Use it to ask whether one extraction workflow systematically gives more total reads, more <em>Listeria</em> reads, or different AMR patterns.</p>
      <p class="note">The AMR tables keep barcode-level rows so you can see which sample each call came from, rather than only showing cohort-wide totals.</p>
      <div class="plot-grid">{imgs(section6_plots)}</div>
      <h3>Method-Level Totals</h3>
      {df_to_html_table(section6_method_summary)}
      <h3>AMR Summary in Reads</h3>
      {df_to_html_table(amr_reads_summary)}
      <h3>AMR Summary in Contigs</h3>
      {df_to_html_table(amr_contigs_summary)}
    </div>

    <div class="section" id="section-7">
      <h2>7. Assembly and Contig Analysis</h2>
      <p class="note">This section shifts from reads to assemblies. The first figure compares Flye, MetaMDBG, and Myloasm directly for contig count, total assembled bases, N50, and median contig length.</p>
      <p class="note">The four assembly-metric panels use a log10 y-axis so differences between lower and higher values stay visible in the same figure. Zero-valued assemblies are omitted from those log-scale boxplots.</p>
      <p class="note">The taxa table below shows which organisms dominate the assembled contigs for each sample and assembler, which is useful for checking whether the assemblies are really centered on <em>Listeria</em> or still dominated by background organisms.</p>
      <p class="note">Strain and serotype labels are aggregated to species level in this table, so entries such as <em>Listeria monocytogenes</em> include matching subtype labels assigned by Kraken.</p>
      <div class="plot-grid">{imgs(section7_plots)}</div>
      <h3>Assembly Summary</h3>
      {df_to_html_table(section7_asm_summary)}
      <h3>Top Contig Taxa by Sample and Assembler</h3>
      {df_to_html_table(section7_taxa)}
    </div>
  </div>
  <script>
    function sortTable(th) {{
      const table = th.closest('table');
      const tbody = table.querySelector('tbody');
      const rows = Array.from(tbody.querySelectorAll('tr'));
      const idx = Array.from(th.parentNode.children).indexOf(th);
      const asc = th.dataset.sort !== 'asc';
      rows.sort((a, b) => {{
        const va = a.children[idx].textContent.trim().replace(/,/g, '');
        const vb = b.children[idx].textContent.trim().replace(/,/g, '');
        const na = parseFloat(va);
        const nb = parseFloat(vb);
        if (!Number.isNaN(na) && !Number.isNaN(nb)) {{
          return asc ? na - nb : nb - na;
        }}
        return asc ? va.localeCompare(vb, undefined, {{numeric: true, sensitivity: 'base'}})
                   : vb.localeCompare(va, undefined, {{numeric: true, sensitivity: 'base'}});
      }});
      rows.forEach((row) => tbody.appendChild(row));
      table.querySelectorAll('th').forEach((header) => {{
        if (header !== th) {{
          delete header.dataset.sort;
        }}
      }});
      th.dataset.sort = asc ? 'asc' : 'desc';
    }}
  </script>
</body>
</html>"""
    out_html.write_text(body, encoding="utf-8")


def _filter_data_by_round(raw: Dict[str, pd.DataFrame], rnd) -> Dict[str, pd.DataFrame]:
    """Return a copy of *raw* filtered to a single sequencing round."""
    filtered: Dict[str, pd.DataFrame] = {}
    for key, df in raw.items():
        if "round" in df.columns:
            filtered[key] = df[df["round"] == rnd].copy()
        elif "sample" in df.columns:
            # Match samples that start with the round prefix
            prefix = f"r{int(rnd)}_"
            filtered[key] = df[df["sample"].str.startswith(prefix, na=False)].copy()
        else:
            filtered[key] = df.copy()
    return filtered


def main() -> int:
    if len(sys.argv) < 2:
        print("Usage: python3 build_local_black_report.py <downloaded_base_dir> [cohort]")
        return 1

    base_dir = Path(sys.argv[1]).resolve()
    raw = load_inputs(base_dir)
    available_cohorts = [c for c in raw["meta"]["cohort"].dropna().unique() if str(c) != "Control"]
    if len(sys.argv) > 2:
        requested = sys.argv[2]
        cohorts = [c for c in available_cohorts if str(c).lower() == requested.lower()]
        if not cohorts:
            print(f"Cohort not found: {requested}")
            return 1
    else:
        preferred = {"Black": 0, "Blue": 1}
        cohorts = sorted(available_cohorts, key=lambda c: (preferred.get(str(c), 99), str(c)))

    # Detect rounds
    rounds = sorted(raw["meta"]["round"].dropna().unique()) if "round" in raw["meta"].columns else [None]

    for rnd in rounds:
        round_suffix = f"_r{int(rnd)}" if rnd is not None and len(rounds) > 1 else ""
        if rnd is not None and len(rounds) > 1:
            round_data = _filter_data_by_round(raw, rnd)
        else:
            round_data = raw

        for cohort in cohorts:
            data = prepare_cohort_tables(round_data, cohort)
            out_dir = base_dir / f"local_{data['cohort_slug']}_report{round_suffix}"
            assets_dir = out_dir / "assets"
            data_dir = out_dir / "data"
            ensure_dir(assets_dir)
            ensure_dir(data_dir)

            read_metrics = build_read_metrics(data["read_lengths"])
            general_metrics = build_general_metrics(data["listeria"], read_metrics, cohort)
            enrichment = compute_enrichment(data["listeria"])
            amr_reads_summary, amr_contigs_summary = build_amr_summaries(data["amr_reads"], data["amr_contigs"])

            s2_plots, s2_table = section_plot_read_lengths_all(data["read_lengths"], assets_dir, cohort)
            s3_plot = section_plot_read_lengths_per_barcode(data["read_lengths"], data["listeria"], assets_dir, cohort)
            s3b_plot = section_plot_read_lengths_by_method(
                data["read_lengths"], assets_dir, cohort, data["group_order"], data["group_labels"]
            )
            s3c_plot = section_plot_method_comparison(
                data["read_lengths"], assets_dir, cohort, data["group_order"], data["group_labels"], data["group_colors"]
            )
            s4_plot = section_plot_listeria_overview(data["listeria"], assets_dir, data["group_colors"], cohort)
            s6_plots, s6_method_summary = section_plot_extraction_comparison(
                data["listeria"], enrichment, assets_dir, data["group_order"], data["group_labels"]
            )
            s7_plots, s7_asm_summary, s7_taxa = section_plot_contigs(
                data["listeria"], data["assembly"], data["contig_taxa"], assets_dir
            )

            prefix = data["cohort_slug"]
            read_metrics.to_csv(data_dir / f"{prefix}_read_metrics_by_sample.csv", index=False)
            build_barcode_metrics_export(read_metrics, data["listeria"]).to_csv(
                data_dir / f"{prefix}_barcode_read_metrics.csv", index=False
            )
            s2_table.to_csv(data_dir / f"{prefix}_read_length_summary_all.csv", index=False)
            enrichment.to_csv(data_dir / f"{prefix}_enrichment_summary.csv", index=False)
            s6_method_summary.to_csv(data_dir / f"{prefix}_method_totals.csv", index=False)
            s7_asm_summary.to_csv(data_dir / f"{prefix}_assembly_summary.csv", index=False)
            s7_taxa.to_csv(data_dir / f"{prefix}_top_contig_taxa.csv", index=False)
            data["listeria"].to_csv(data_dir / f"{prefix}_listeria_overview_flat.csv", index=False)

            section4_table = data["listeria"].sort_values(["barcode", "condition"])[
                [
                    "sample",
                    "sample_id",
                    "group_label",
                    "total_reads",
                    "listeria_reads",
                    "listeria_pct",
                    "listeria_mean_len",
                    "listeria_median_len",
                ]
            ]

            section3_table = read_metrics.merge(
                data["meta"][["sample", "sample_id", "group_label"]],
                on="sample",
                how="left",
            )[
                [
                    "sample",
                    "sample_id",
                    "group_label",
                    "condition",
                    "total_reads",
                    "total_bases",
                    "mean_read_length",
                    "median_read_length",
                    "read_length_n50",
                ]
            ].sort_values(["sample"])

            write_html(
                out_dir / "report.html",
                cohort=cohort,
                general_metrics=general_metrics,
                section2_summary=s2_table,
                section2_plots=s2_plots,
                section3_plot=s3_plot,
                section3b_plot=s3b_plot,
                section3c_plot=s3c_plot,
                section3_table=section3_table,
                section4_plot=s4_plot,
                section4_table=section4_table,
                section6_plots=s6_plots,
                section6_method_summary=s6_method_summary,
                section7_plots=s7_plots,
                section7_asm_summary=s7_asm_summary,
                section7_taxa=s7_taxa,
                amr_reads_summary=amr_reads_summary,
                amr_contigs_summary=amr_contigs_summary,
            )

            print(f"Report written to: {out_dir / 'report.html'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

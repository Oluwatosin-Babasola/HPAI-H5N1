#!/usr/bin/env python3
"""
Create host source-sink transmission heatmaps from either:
  1) host_transition_counts.csv   with columns: from_host, into_host, n
  2) host_transition_counts_matrix.csv with first column from_host and remaining columns as sinks
  3) edge_annotations.csv / edge_annotations_simple.csv with columns: from_host, into_host

Outputs:
  - host_transition_heatmap_counts.png
  - host_transition_heatmap_row_normalized.png
  - host_source_sink_summary.csv
  - host_source_sink_barplot.png

Usage examples:
  python host_source_sink_heatmap.py --input host_transition_counts.csv --outdir host_figures
  python host_source_sink_heatmap.py --input edge_annotations.csv --outdir host_figures --exclude-self-summary
"""

from __future__ import annotations

import argparse
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

HOST_ORDER_DEFAULT = ["WB", "DB", "WM", "DM", "H"]


def load_transition_matrix(path: Path, weight_col: str = "n") -> pd.DataFrame:
    df = pd.read_csv(path)
    cols = {c.lower(): c for c in df.columns}
    lower_cols = [c.lower() for c in df.columns]

    # Case 1: long table with counts
    if {"from_host", "into_host"}.issubset(lower_cols):
        from_col = cols["from_host"]
        into_col = cols["into_host"]
        if weight_col in df.columns:
            wcol = weight_col
        elif "n" in df.columns:
            wcol = "n"
        else:
            wcol = None

        if wcol is None:
            out = df.groupby([from_col, into_col], dropna=False).size().unstack(fill_value=0)
        else:
            out = (
                df.groupby([from_col, into_col], dropna=False)[wcol]
                .sum()
                .unstack(fill_value=0)
            )
        out.index.name = "from_host"
        out.columns.name = "into_host"
        return out

    # Case 2: wide matrix with first column from_host
    first_col = df.columns[0]
    if first_col.lower() == "from_host":
        out = df.set_index(first_col)
        out.index.name = "from_host"
        out.columns.name = "into_host"
        return out.apply(pd.to_numeric, errors="coerce").fillna(0)

    raise ValueError(
        "Could not detect input format. Expected from_host / into_host columns "
        "or a matrix file whose first column is from_host."
    )


def reorder_matrix(mat: pd.DataFrame, host_order: list[str] | None = None) -> pd.DataFrame:
    observed = list(pd.Index(mat.index).astype(str).union(pd.Index(mat.columns).astype(str)))
    if host_order:
        ordered = [h for h in host_order if h in observed]
        extras = [h for h in observed if h not in ordered]
        final_order = ordered + sorted(extras)
    else:
        final_order = sorted(observed)

    return mat.reindex(index=final_order, columns=final_order, fill_value=0)


def draw_heatmap(
    mat: pd.DataFrame,
    out_png: Path,
    title: str,
    fmt: str = ".0f",
    cmap: str = "YlOrRd",
    value_label: str = "Count",
) -> None:
    fig, ax = plt.subplots(figsize=(7.5, 6.2))
    im = ax.imshow(mat.values, aspect="auto", cmap=cmap)

    ax.set_xticks(np.arange(mat.shape[1]))
    ax.set_yticks(np.arange(mat.shape[0]))
    ax.set_xticklabels(mat.columns)
    ax.set_yticklabels(mat.index)
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right", rotation_mode="anchor")

    # cell annotations
    for i in range(mat.shape[0]):
        for j in range(mat.shape[1]):
            val = mat.iat[i, j]
            if fmt == ".0f":
                txt = f"{val:.0f}"
            elif fmt == ".2f":
                txt = f"{val:.2f}"
            else:
                txt = format(val, fmt)
            ax.text(j, i, txt, ha="center", va="center", fontsize=10)

    cbar = fig.colorbar(im, ax=ax, shrink=0.88)
    cbar.set_label(value_label)

    ax.set_xlabel("Sink host (into_host)")
    ax.set_ylabel("Source host (from_host)")
    ax.set_title(title)

    fig.tight_layout()
    fig.savefig(out_png, dpi=300, bbox_inches="tight")
    plt.close(fig)


def compute_source_sink_summary(
    mat: pd.DataFrame,
    exclude_self: bool = True,
) -> pd.DataFrame:
    m = mat.copy().astype(float)
    diag = pd.Series(np.diag(m.values), index=m.index, name="self_transitions")

    outgoing_total = m.sum(axis=1)
    incoming_total = m.sum(axis=0)

    if exclude_self:
        outgoing_other = outgoing_total - diag
        incoming_other = incoming_total - diag
    else:
        outgoing_other = outgoing_total.copy()
        incoming_other = incoming_total.copy()

    summary = pd.DataFrame({
        "host": m.index,
        "self_transitions": diag.reindex(m.index).values,
        "outgoing_total": outgoing_total.reindex(m.index).values,
        "incoming_total": incoming_total.reindex(m.index).values,
        "outgoing_other_hosts": outgoing_other.reindex(m.index).values,
        "incoming_other_hosts": incoming_other.reindex(m.index).values,
    })
    summary["net_source_sink"] = summary["outgoing_other_hosts"] - summary["incoming_other_hosts"]
    summary["classification"] = np.where(
        summary["net_source_sink"] > 0,
        "Source",
        np.where(summary["net_source_sink"] < 0, "Sink", "Balanced"),
    )
    summary = summary.sort_values("net_source_sink", ascending=False)
    return summary


def draw_source_sink_barplot(summary: pd.DataFrame, out_png: Path) -> None:
    fig, ax = plt.subplots(figsize=(7.2, 4.8))
    x = np.arange(len(summary))
    vals = summary["net_source_sink"].values
    ax.bar(x, vals)
    ax.axhline(0, linewidth=1)
    ax.set_xticks(x)
    ax.set_xticklabels(summary["host"])
    ax.set_ylabel("Net source-sink score")
    ax.set_xlabel("Host")
    ax.set_title("Net source-sink transmission by host")

    for xi, yi in zip(x, vals):
        ax.text(xi, yi, f"{yi:.0f}", ha="center", va="bottom" if yi >= 0 else "top", fontsize=9)

    fig.tight_layout()
    fig.savefig(out_png, dpi=300, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True, help="Path to input CSV")
    parser.add_argument("--outdir", default="host_source_sink_outputs", help="Output directory")
    parser.add_argument(
        "--host-order",
        nargs="+",
        default=HOST_ORDER_DEFAULT,
        help="Host order for rows/columns, e.g. WB DB WM DM H",
    )
    parser.add_argument(
        "--include-self-summary",
        action="store_true",
        help="Include self-transitions in the source-sink net score. Default excludes them.",
    )
    parser.add_argument(
        "--interhost-only",
        action="store_true",
        help="Zero the diagonal so plots and tables show only host-switch transitions.",
    )
    args = parser.parse_args()

    input_path = Path(args.input)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    mat = load_transition_matrix(input_path)
    mat = reorder_matrix(mat, host_order=args.host_order)

    if args.interhost_only:
        common = [h for h in mat.index if h in mat.columns]
        for h in common:
            mat.loc[h, h] = 0

    mat.to_csv(outdir / "host_transition_matrix_used.csv")

    draw_heatmap(
        mat=mat,
        out_png=outdir / "host_transition_heatmap_counts.png",
        title="Host-to-host transmission heatmap (counts)",
        fmt=".0f",
        cmap="YlOrRd",
        value_label="Transition count",
    )

    row_sums = mat.sum(axis=1).replace(0, np.nan)
    mat_norm = mat.div(row_sums, axis=0).fillna(0)
    mat_norm.to_csv(outdir / "host_transition_matrix_row_normalized.csv")
    draw_heatmap(
        mat=mat_norm,
        out_png=outdir / "host_transition_heatmap_row_normalized.png",
        title="Host-to-host transmission heatmap (row-normalized)",
        fmt=".2f",
        cmap="Blues",
        value_label="Row proportion",
    )

    summary = compute_source_sink_summary(
        mat=mat,
        exclude_self=not args.include_self_summary,
    )
    summary.to_csv(outdir / "host_source_sink_summary.csv", index=False)
    draw_source_sink_barplot(summary, outdir / "host_source_sink_barplot.png")

    print("Saved outputs to:", outdir)
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()

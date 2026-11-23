# Global H5 Host and Region Transition Analysis

Scripts for reconstructing host and geographic transitions on a global influenza
time-calibrated phylogeny, and for plotting the results.

## Scripts

| Script | What it does |
| --- | --- |
| `global_tree_transitions_interhost.R` | Parses tip labels, reconstructs ancestral host and region states, and writes per-edge annotations plus host and region transition tables. |
| `Tree_region_host_heatmaps.R` | Draws the time tree with host-colored tips beside Region and Host heat-map strips. |
| `host_source_sink_heatmap.py` | Builds host-to-host heat maps and a net source-sink summary from the transition tables. |

## Inputs

- `timetree.nwk` — time-calibrated tree. Tip labels are expected as
  `strain|accession|host|region|date`.
- `meta_data.csv` or `meta.csv` — optional metadata, joined on strain name or accession.
  Values here override those parsed from tip labels.

Host classes are `WB` wild bird, `DB` domestic bird, `WM` wild mammal,
`DM` domestic mammal, and `H` human.

## Outputs

The R transition script writes to `global_outputs/`:

- `edge_annotations.csv` — every branch with its parent and child host and region.
- `host_transition_counts.csv` — cross-host transitions only, for source-sink reading.
- `host_transition_counts_all_edges.csv` — all edges, including same-host continuity.
- `region_transition_counts.csv` and matching in/out summaries.

The Python script writes count and row-normalized heat maps plus
`host_source_sink_summary.csv`.

## Usage

```bash
Rscript updated_global_tree_transitions_interhost.R
python host_source_sink_heatmap.py --input global_outputs/host_transition_counts.csv \
  --outdir host_figures --interhost-only
Rscript Tree_region_host_heatmaps_colored_tips.R
```

Edit the file paths at the top of each R script before running.

## Requirements
R with ape, phytools, phangorn, ggtree, ggplot2, patchwork, and the
tidyverse. Python 3 with pandas, numpy, and matplotlib.

#!/usr/bin/env Rscript

# ============================================================
# Cleaner tree + Region and Host heat-map strips
# - tree branches only (no host-colored tips, to avoid duplicate host encoding)
# - one Region legend and one Host legend
# - dynamic region values from the actual data
# ============================================================

suppressPackageStartupMessages({
  library(ape)
  library(ggtree)
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(stringr)
  library(lubridate)
  library(patchwork)
  library(tibble)
})

tree_file <- "~/trees/sample1/timetree/timetree.nwk"
meta_file <- "~/metadata/meta_data.csv"  # optional

out_png  <- "~/Tree_region_host_heatmaps.png"
out_pdf  <- "~/Tree_region_host_heatmaps.pdf"

legend_text_size <- 9
axis_text_size   <- 9
title_text_size  <- 10

host_cols <- c(
  WB = "#1f78b4",
  DB = "#d62728",
  WM = "#6a3d9a",
  DM = "#ff7f00",
  H  = "#7f7f7f"
)

make_region_palette <- function(region_levels) {
  cols <- grDevices::hcl.colors(length(region_levels), palette = "Set 3")
  setNames(cols, region_levels)
}

clean_region <- function(x) {
  x <- as.character(x)
  x <- stringr::str_trim(x)
  x <- stringr::str_replace_all(x, "_", " ")
  x[x %in% c("", "NA", "nan")] <- NA_character_
  x
}

clean_host <- function(x) {
  x <- as.character(x)
  x <- stringr::str_trim(x)
  x <- stringr::str_to_upper(x)
  x[x %in% c("", "NA", "NAN")] <- NA_character_
  dplyr::case_when(
    x %in% c("WB", "WILD_BIRD", "WILD BIRD") ~ "WB",
    x %in% c("DB", "DOMESTIC_BIRD", "DOMESTIC BIRD", "POULTRY") ~ "DB",
    x %in% c("WM", "WILD_MAMMAL", "WILD MAMMAL") ~ "WM",
    x %in% c("DM", "DOMESTIC_MAMMAL", "DOMESTIC MAMMAL") ~ "DM",
    x %in% c("H", "HUMAN") ~ "H",
    TRUE ~ x
  )
}

read_tree_safely <- function(path) {
  lines <- readLines(path, warn = FALSE, n = 50)
  lines_trim <- trimws(lines)
  lines_trim <- lines_trim[lines_trim != ""]
  first_line <- lines_trim[1]

  if (grepl("^#NEXUS", first_line, ignore.case = TRUE) ||
      any(grepl("^begin\\s+trees", lines_trim, ignore.case = TRUE)) ||
      any(grepl("^begin\\s+taxa", lines_trim, ignore.case = TRUE))) {
    tree <- read.nexus(path)
  } else {
    tree <- read.tree(path)
  }
  if (inherits(tree, "multiPhylo")) tree <- tree[[1]]
  tree$tip.label <- gsub("^'|'$", "", tree$tip.label)
  tree
}

make_strip <- function(dat, fill_var, fill_cols, title_text, show_legend = TRUE) {
  ggplot(dat, aes(x = 1, y = y, fill = .data[[fill_var]])) +
    geom_tile(width = 0.98, height = 0.98) +
    scale_fill_manual(values = fill_cols, drop = TRUE, na.translate = FALSE, name = title_text) +
    scale_y_continuous(expand = c(0, 0)) +
    scale_x_continuous(expand = c(0, 0)) +
    labs(title = title_text) +
    theme_void() +
    theme(
      plot.title   = element_text(size = title_text_size, face = "bold", hjust = 0.5),
      legend.title = element_text(size = title_text_size, face = "bold"),
      legend.text  = element_text(size = legend_text_size),
      plot.margin  = margin(5.5, 1, 5.5, 1),
      panel.background = element_rect(fill = "white", color = NA),
      plot.background  = element_rect(fill = "white", color = NA)
    )
}

tree <- read_tree_safely(tree_file)

meta_plot <- tibble(label = tree$tip.label) %>%
  mutate(
    sample_id  = str_split_i(label, "\\|", 1),
    accession  = str_split_i(label, "\\|", 2),
    host_tip   = str_split_i(label, "\\|", 3),
    region_tip = str_split_i(label, "\\|", 4),
    date_tip   = str_split_i(label, "\\|", 5),
    host       = clean_host(host_tip),
    region     = clean_region(region_tip),
    sample_date = suppressWarnings(mdy(date_tip))
  )

if (file.exists(meta_file)) {
  meta_in <- tryCatch(read_csv(meta_file, show_col_types = FALSE), error = function(e) NULL)
  if (!is.null(meta_in) && "name" %in% names(meta_in)) {
    meta_in <- meta_in %>%
      mutate(
        name = as.character(name),
        name = na_if(name, "nan"),
        name = na_if(name, "NA"),
        name = na_if(name, ""),
        host_meta = if ("host" %in% names(.)) clean_host(host) else NA_character_,
        region_meta = if ("region" %in% names(.)) clean_region(region) else NA_character_
      ) %>%
      filter(!is.na(name)) %>%
      distinct(name, .keep_all = TRUE) %>%
      select(name, host_meta, region_meta)

    meta_plot <- meta_plot %>%
      left_join(meta_in, by = c("label" = "name")) %>%
      mutate(
        host = coalesce(host_meta, host),
        region = coalesce(region_meta, region)
      ) %>%
      select(label, accession, sample_id, region, host, sample_date)
  }
}

host_levels_used   <- intersect(names(host_cols), unique(na.omit(meta_plot$host)))
region_levels_used <- sort(unique(na.omit(meta_plot$region)))

host_cols_used   <- host_cols[host_levels_used]
region_cols_used <- make_region_palette(region_levels_used)

meta_plot <- meta_plot %>%
  mutate(
    host   = factor(host, levels = host_levels_used),
    region = factor(region, levels = region_levels_used)
  )

mrsd_date <- max(meta_plot$sample_date, na.rm = TRUE)

p_tree <- ggtree(
  tree,
  layout = "rectangular",
  size   = 0.22,
  color  = "grey75",
  mrsd   = mrsd_date,
  as.Date = TRUE
) %<+% meta_plot +
  geom_tippoint(aes(color = host), size = 0.85, alpha = 0.95, na.rm = TRUE) +
  scale_color_manual(values = host_cols_used, drop = TRUE, name = "Host") +
  theme_tree2() +
  scale_x_date(date_breaks = "6 months", date_labels = "%Y-%m") +
  theme(
    axis.title.x    = element_blank(),
    axis.text.x     = element_text(size = axis_text_size, color = "black"),
    axis.ticks.x    = element_line(color = "black"),
    plot.margin     = margin(5.5, 8, 5.5, 5.5),
    legend.title    = element_text(size = title_text_size, face = "bold"),
    legend.text     = element_text(size = legend_text_size),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background  = element_rect(fill = "white", color = NA)
  )

tip_order <- p_tree$data %>%
  filter(isTip) %>%
  select(label, y)

meta_ord <- tip_order %>%
  left_join(meta_plot, by = "label") %>%
  arrange(y)

p_region <- make_strip(meta_ord, "region", region_cols_used, "Region")
p_host   <- make_strip(meta_ord, "host", host_cols_used, "Host") + guides(fill = "none")

panel_B <- p_tree | p_region | p_host
panel_B <- panel_B +
  plot_layout(widths = c(10.5, 0.5, 0.50), guides = "collect") &
  theme(
    legend.position = "right",
    legend.box = "vertical",
    panel.background = element_rect(fill = "white", color = NA),
    plot.background  = element_rect(fill = "white", color = NA)
  )

ggsave(out_png, panel_B, width = 13.5, height = 8.5, dpi = 600, bg = "white")
ggsave(out_pdf, panel_B, width = 13.5, height = 8.5, bg = "white")

cat("\nSaved:\n")
cat(out_png, "\n")
cat(out_pdf, "\n")

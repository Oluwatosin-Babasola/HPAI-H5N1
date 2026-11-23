suppressPackageStartupMessages({
  library(ape)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(tibble)
  library(stringr)
  library(lubridate)
  library(phytools)
  library(phangorn)
})

# ============================================================
# INPUTS
# ============================================================
TREE_FILE <- "timetree.nwk"
META_FILE <- "meta.csv"   # optional
OUT_DIR   <- "global_outputs"

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# ============================================================
# HELPERS
# ============================================================
normalize_ids <- function(x) {
  x %>%
    as.character() %>%
    trimws() %>%
    gsub("^'+|'+$", "", ., perl = TRUE) %>%
    gsub('^"+|"+$', "", ., perl = TRUE) %>%
    gsub("\\s+", " ", ., perl = TRUE)
}

extract_part <- function(parts, i) {
  sapply(parts, function(x) if (length(x) >= i) x[i] else NA_character_)
}

extract_last <- function(parts) {
  sapply(parts, function(x) if (length(x) >= 1) x[length(x)] else NA_character_)
}

normalize_geo <- function(x) {
  x %>%
    as.character() %>%
    str_replace_all("_", " ") %>%
    str_squish() %>%
    na_if("") %>%
    str_to_title()
}

normalize_host_class <- function(x) {
  y <- toupper(str_squish(as.character(x)))
  y <- case_when(
    y %in% c("WB", "WILD BIRD", "WILD_BIRD", "WATERFOWL", "AVIAN_WILD") ~ "WB",
    y %in% c("DB", "DOMESTIC BIRD", "DOMESTIC_BIRD", "POULTRY", "CHICKEN", "TURKEY", "DUCK DOMESTIC") ~ "DB",
    y %in% c("WM", "WILD MAMMAL", "WILD_MAMMAL") ~ "WM",
    y %in% c("DM", "DOMESTIC MAMMAL", "DOMESTIC_MAMMAL", "CATTLE", "BOVINE", "DAIRY CATTLE") ~ "DM",
    y %in% c("H", "HUMAN", "HOMO SAPIENS") ~ "H",
    TRUE ~ y
  )
  y[y == ""] <- NA_character_
  y
}

# Influenza name parser:
# avian example: A/pheasant/Wales/385129/2021  -> geography = Wales
# human example: A/Texas/37/2024               -> geography = Texas
extract_geo_from_name <- function(name) {
  parts <- str_split(name, "/")
  sapply(parts, function(x) {
    x <- x[nzchar(x)]
    if (length(x) >= 5) return(x[3])  # A/species/location/isolate/year
    if (length(x) == 4) return(x[2])  # A/location/isolate/year
    NA_character_
  })
}

safe_date_parse <- function(x) {
  suppressWarnings(coalesce(
    mdy(x, quiet = TRUE),
    ymd(x, quiet = TRUE),
    dmy(x, quiet = TRUE)
  ))
}

standardize_meta_names <- function(df) {
  names(df) <- names(df) %>%
    str_replace_all("[^A-Za-z0-9]+", "_") %>%
    str_replace_all("(^_+|_+$)", "") %>%
    str_to_lower()
  df
}

choose_first_present <- function(df, candidates) {
  hits <- intersect(candidates, names(df))
  if (length(hits)) hits[1] else NA_character_
}

majority_descendant_states <- function(tr, x_factor) {
  n_tips <- length(tr$tip.label)
  out <- character(tr$Nnode)
  for (node in seq_len(tr$Nnode) + n_tips) {
    desc <- phytools::getDescendants(tr, node)
    tip_desc <- desc[desc <= n_tips]
    vals <- as.character(x_factor[tip_desc])
    vals <- vals[!is.na(vals) & vals != ""]
    out[node - n_tips] <- if (length(vals)) names(sort(table(vals), decreasing = TRUE))[1] else NA_character_
  }
  out
}

infer_node_labels <- function(tr, tip_named_vector, trait_name = "trait") {
  x <- tip_named_vector[tr$tip.label]
  x <- as.character(x)
  x[is.na(x) | x == ""] <- "Unknown"
  x <- factor(x)

  if (nlevels(x) == 1) {
    return(rep(levels(x)[1], tr$Nnode))
  }

  rr <- tryCatch(
    phytools::rerootingMethod(tr, x, model = "ER"),
    error = function(e) NULL
  )

  if (!is.null(rr) && !is.null(rr$marginal.anc) && is.matrix(rr$marginal.anc)) {
    lev <- levels(x)
    return(lev[max.col(rr$marginal.anc, ties.method = "first")])
  }

  ac <- tryCatch(
    ape::ace(x, tr, type = "discrete", model = "ER"),
    error = function(e) NULL
  )

  if (!is.null(ac) && !is.null(ac$lik.anc) && is.matrix(ac$lik.anc)) {
    lev <- colnames(ac$lik.anc)
    return(lev[max.col(ac$lik.anc, ties.method = "first")])
  }

  message("ML ancestral reconstruction failed for ", trait_name,
          "; using majority-descendant fallback.")
  majority_descendant_states(tr, x)
}

# ============================================================
# LOAD TREE
# ============================================================
tr <- read.tree(TREE_FILE)
tr$tip.label <- normalize_ids(tr$tip.label)

# ============================================================
# BUILD METADATA DIRECTLY FROM TREE LABELS
# Expected label pattern:
# STRAIN_NAME|EPI_ISL_xxx|HOST|DATE
# ============================================================
meta_tree <- tibble(label = tr$tip.label) %>%
  mutate(
    parts       = str_split(label, "\\|"),
    name        = extract_part(parts, 1),
    accession   = extract_part(parts, 2),
    host_tree   = normalize_host_class(extract_part(parts, 3)),
    date_raw    = extract_last(parts),
    date_tree   = safe_date_parse(date_raw),
    region_tree = normalize_geo(extract_geo_from_name(name))
  ) %>%
  select(label, name, accession, region_tree, host_tree, date_raw, date_tree)

cat("Preview of parsed tree metadata:\n")
print(head(meta_tree, 10))
cat("Unique parsed hosts in tree:\n")
print(table(meta_tree$host_tree, useNA = "ifany"))
cat("Top parsed regions from tree labels:\n")
print(sort(table(meta_tree$region_tree), decreasing = TRUE)[1:min(15, length(unique(meta_tree$region_tree)))])

# ============================================================
# OPTIONAL EXTERNAL METADATA JOIN
# ============================================================
meta <- meta_tree
if (file.exists(META_FILE)) {
  meta_ext <- read_csv(META_FILE, show_col_types = FALSE) %>%
    standardize_meta_names()

  for (nm in intersect(c("name", "strain", "strain_name", "isolate_id", "label",
                         "accession", "epi_isl", "gisaid_accession"), names(meta_ext))) {
    meta_ext[[nm]] <- normalize_ids(meta_ext[[nm]])
  }

  join_key_tree <- NULL
  join_key_meta <- NULL

  if ("name" %in% names(meta_ext)) {
    join_key_tree <- "name"
    join_key_meta <- "name"
  } else if ("strain" %in% names(meta_ext)) {
    join_key_tree <- "name"
    join_key_meta <- "strain"
  } else if ("strain_name" %in% names(meta_ext)) {
    join_key_tree <- "name"
    join_key_meta <- "strain_name"
  } else if ("accession" %in% names(meta_ext)) {
    join_key_tree <- "accession"
    join_key_meta <- "accession"
  } else if ("epi_isl" %in% names(meta_ext)) {
    join_key_tree <- "accession"
    join_key_meta <- "epi_isl"
  } else if ("gisaid_accession" %in% names(meta_ext)) {
    join_key_tree <- "accession"
    join_key_meta <- "gisaid_accession"
  } else if ("label" %in% names(meta_ext)) {
    join_key_tree <- "label"
    join_key_meta <- "label"
  }

  if (!is.null(join_key_tree)) {
    meta <- meta_tree %>%
      left_join(meta_ext, by = setNames(join_key_meta, join_key_tree), suffix = c("_tree", "_ext"))
    cat("Joined external metadata using:", join_key_tree, "->", join_key_meta, "\n")
  } else {
    cat("meta.csv found, but no usable join key detected; proceeding with tree-derived metadata only.\n")
  }
}

# ============================================================
# FINAL TRAIT SELECTION
# ============================================================
region_col <- choose_first_present(meta, c("region", "hhs", "state", "stateprovince", "country", "location"))
host_col   <- choose_first_present(meta, c("host", "broad_host", "host_group", "host_class"))
date_col   <- choose_first_present(meta, c("date", "collection_date", "sample_date"))

meta_final <- meta %>%
  mutate(
    region = if (!is.na(region_col)) normalize_geo(.data[[region_col]]) else region_tree,
    host   = if (!is.na(host_col)) normalize_host_class(.data[[host_col]]) else host_tree,
    date   = if (!is.na(date_col)) safe_date_parse(.data[[date_col]]) else date_tree
  ) %>%
  mutate(
    region = coalesce(region, region_tree, "Unknown"),
    host   = coalesce(host, host_tree, "Unknown")
  )

# ============================================================
# SANITY CHECKS
# ============================================================
stopifnot(all(tr$tip.label %in% meta_final$label))

missing_region <- sum(is.na(meta_final$region) | meta_final$region == "Unknown")
missing_host   <- sum(is.na(meta_final$host)   | meta_final$host == "Unknown")
cat("Tips with missing/unknown region:", missing_region, "\n")
cat("Tips with missing/unknown host:",   missing_host, "\n")

tip_regions <- setNames(meta_final$region, meta_final$label)
tip_hosts   <- setNames(meta_final$host,   meta_final$label)

# ============================================================
# SANITIZE EDGE LENGTHS
# ============================================================
if (is.null(tr$edge.length)) tr$edge.length <- rep(1, nrow(tr$edge))

bad <- !is.finite(tr$edge.length) | tr$edge.length <= 0
if (any(bad)) {
  pos_min <- suppressWarnings(min(tr$edge.length[tr$edge.length > 0], na.rm = TRUE))
  eps <- if (is.finite(pos_min)) max(pos_min / 1000, 1e-8) else 1e-8
  tr$edge.length[bad] <- eps
}

tr <- di2multi(tr, tol = 1e-8)

if (!is.rooted(tr)) {
  tr <- phangorn::midpoint(tr)
}

# ============================================================
# ANCESTRAL RECONSTRUCTION WITH FALLBACKS
# ============================================================
node_regions_lab <- infer_node_labels(tr, tip_regions, trait_name = "region")
node_hosts_lab   <- infer_node_labels(tr, tip_hosts,   trait_name = "host")

# ============================================================
# BUILD EDGE ANNOTATIONS
# ============================================================
edg <- as.data.frame(tr$edge)
names(edg) <- c("parent", "child")
n_tips <- length(tr$tip.label)

is_child_tip  <- edg$child  <= n_tips
is_parent_tip <- edg$parent <= n_tips

lab_child_region <- character(nrow(edg))
lab_child_region[is_child_tip]  <- as.character(tip_regions[tr$tip.label[edg$child[is_child_tip]]])
lab_child_region[!is_child_tip] <- node_regions_lab[edg$child[!is_child_tip] - n_tips]

lab_parent_region <- character(nrow(edg))
lab_parent_region[is_parent_tip]  <- as.character(tip_regions[tr$tip.label[edg$parent[is_parent_tip]]])
lab_parent_region[!is_parent_tip] <- node_regions_lab[edg$parent[!is_parent_tip] - n_tips]

lab_child_host <- character(nrow(edg))
lab_child_host[is_child_tip]  <- as.character(tip_hosts[tr$tip.label[edg$child[is_child_tip]]])
lab_child_host[!is_child_tip] <- node_hosts_lab[edg$child[!is_child_tip] - n_tips]

lab_parent_host <- character(nrow(edg))
lab_parent_host[is_parent_tip]  <- as.character(tip_hosts[tr$tip.label[edg$parent[is_parent_tip]]])
lab_parent_host[!is_parent_tip] <- node_hosts_lab[edg$parent[!is_parent_tip] - n_tips]

edge_annotations <- tibble(
  parent = edg$parent,
  child = edg$child,
  from_region = lab_parent_region,
  into_region = lab_child_region,
  from_hhs = lab_parent_region,
  into_hhs = lab_child_region,
  from_host = lab_parent_host,
  into_host = lab_child_host,
  branch_length = tr$edge.length,
  interregion = lab_parent_region != lab_child_region,
  interstate = lab_parent_region != lab_child_region,
  interhost = lab_parent_host != lab_child_host,
  same_region = lab_parent_region == lab_child_region,
  same_host = lab_parent_host == lab_child_host
)

write_csv(edge_annotations, file.path(OUT_DIR, "edge_annotations.csv"))

edge_annotations %>%
  select(from_region, into_region, from_hhs, into_hhs,
         from_host, into_host, branch_length, interregion, interstate,
         interhost, same_region, same_host) %>%
  write_csv(file.path(OUT_DIR, "edge_annotations_simple.csv"))

# ============================================================
# DERIVED TABLES
# IMPORTANT:
# - host_transition_counts.csv is now INTERHOST ONLY
# - diagonal same-host edges are written separately for transparency
# ============================================================

# All host-labelled edges, including same-host continuity
host_edge_counts_all <- edge_annotations %>%
  count(from_host, into_host, name = "n", sort = TRUE)
write_csv(host_edge_counts_all, file.path(OUT_DIR, "host_transition_counts_all_edges.csv"))
host_edge_counts_all %>%
  pivot_wider(names_from = into_host, values_from = n, values_fill = 0) %>%
  write_csv(file.path(OUT_DIR, "host_transition_counts_matrix_all_edges.csv"))

# Cross-host transitions only: this is the table to use for source-sink interpretation
host_transition_counts <- edge_annotations %>%
  filter(interhost) %>%
  count(from_host, into_host, name = "n", sort = TRUE)
write_csv(host_transition_counts, file.path(OUT_DIR, "host_transition_counts.csv"))
host_transition_counts %>%
  pivot_wider(names_from = into_host, values_from = n, values_fill = 0) %>%
  write_csv(file.path(OUT_DIR, "host_transition_counts_matrix.csv"))

# Same-host continuity counts for audit trail
host_same_counts <- edge_annotations %>%
  filter(same_host) %>%
  count(host = from_host, name = "n_same_host_edges", sort = TRUE)
write_csv(host_same_counts, file.path(OUT_DIR, "host_same_host_edge_counts.csv"))

# Host in/out summary using cross-host transitions only
host_inout_summary <- bind_rows(
  edge_annotations %>% filter(interhost) %>% count(host = from_host, name = "outbound"),
  edge_annotations %>% filter(interhost) %>% count(host = into_host, name = "inbound")
) %>%
  group_by(host) %>%
  summarise(
    outbound = sum(ifelse(is.na(outbound), 0, outbound), na.rm = TRUE),
    inbound  = sum(ifelse(is.na(inbound), 0, inbound), na.rm = TRUE),
    net_flow = outbound - inbound,
    .groups = "drop"
  ) %>%
  arrange(desc(net_flow))
write_csv(host_inout_summary, file.path(OUT_DIR, "host_inout_summary.csv"))

# Region tables
region_transition_counts <- edge_annotations %>%
  filter(interregion) %>%
  count(from_region, into_region, name = "n", sort = TRUE)
write_csv(region_transition_counts, file.path(OUT_DIR, "region_transition_counts.csv"))
region_transition_counts %>%
  pivot_wider(names_from = into_region, values_from = n, values_fill = 0) %>%
  write_csv(file.path(OUT_DIR, "region_transition_counts_matrix.csv"))

# legacy filenames for old downstream code
write_csv(region_transition_counts, file.path(OUT_DIR, "state_transition_counts.csv"))
region_transition_counts %>%
  pivot_wider(names_from = into_region, values_from = n, values_fill = 0) %>%
  write_csv(file.path(OUT_DIR, "state_transition_counts_matrix.csv"))

region_inout_summary <- bind_rows(
  edge_annotations %>% filter(interregion) %>% count(region = from_region, name = "outbound"),
  edge_annotations %>% filter(interregion) %>% count(region = into_region, name = "inbound")
) %>%
  group_by(region) %>%
  summarise(
    outbound = sum(ifelse(is.na(outbound), 0, outbound), na.rm = TRUE),
    inbound  = sum(ifelse(is.na(inbound), 0, inbound), na.rm = TRUE),
    net_flow = outbound - inbound,
    .groups = "drop"
  ) %>%
  arrange(desc(net_flow))

write_csv(region_inout_summary, file.path(OUT_DIR, "region_inout_summary.csv"))
write_csv(region_inout_summary %>% rename(hhs = region),
          file.path(OUT_DIR, "state_inout_summary.csv"))

cat("Fresh files written to:\n", OUT_DIR, "\n")
cat("\nHost transition tables now mean:\n")
cat("- host_transition_counts.csv               = cross-host transitions only\n")
cat("- host_transition_counts_all_edges.csv     = all edges, including same-host continuity\n")
cat("- host_same_host_edge_counts.csv           = diagonal / same-host edge counts only\n")

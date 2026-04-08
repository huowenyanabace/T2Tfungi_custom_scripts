#!/usr/bin/env Rscript
# ==============================================================
# GO Enrichment Visualization for Coprinellus radians
# Expanded gene families - Bubble Plot + GO DAG
# Figure titles output to separate .txt file (not shown on plots)
#
# Usage:
#   Rscript plot_go_enrichment.R Crad_GO_enrichment_results.csv Crad_results.csv
#
# Required packages: ggplot2, dplyr, tidyr, igraph, ggraph, scales, stringr
# ==============================================================

# ---- Load libraries ----
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(scales)
  library(stringr)
  library(igraph)
  library(ggraph)
})

# ---- Parse command line arguments ----
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  cat("Usage: Rscript plot_go_enrichment.R <enrichment.csv> <go_annotation.csv> [output_prefix]\n")
  cat("Example: Rscript plot_go_enrichment.R Crad_GO_enrichment_results.csv Crad_results.csv Crad_GO\n")
  quit(status = 1)
}

enrichment_file <- args[1]
annotation_file <- args[2]
output_prefix   <- ifelse(length(args) >= 3, args[3], "Crad_GO")

cat("============================================================\n")
cat("GO Enrichment Visualization for C. radians (R version)\n")
cat("============================================================\n\n")

# ---- Global theme for SCI publication ----
theme_pub <- theme_bw(base_size = 10, base_family = "sans") +
  theme(
    plot.title       = element_blank(),
    axis.text        = element_text(color = "black", size = 8.5),
    axis.title       = element_text(size = 10, face = "bold"),
    legend.title     = element_text(size = 8.5, face = "bold"),
    legend.text      = element_text(size = 8),
    legend.key.size  = unit(0.4, "cm"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border     = element_rect(color = "black", linewidth = 0.6),
    plot.margin      = margin(10, 10, 10, 10)
  )

# ============================================================
# 1. Load & merge data
# ============================================================
cat("[1/4] Loading data...\n")

enrich <- read.csv(enrichment_file, stringsAsFactors = FALSE)
annot  <- read.csv(annotation_file, stringsAsFactors = FALSE)

# Build ID -> Name / Namespace mapping
id2name <- setNames(annot$Name, annot$GO_ID)
id2ns   <- setNames(annot$Namespace, annot$GO_ID)

# Replace NOT FOUND with NA
id2name[id2name == "NOT FOUND"] <- NA

enrich$Term <- id2name[enrich$ID]
enrich$Term <- ifelse(is.na(enrich$Term), enrich$ID, enrich$Term)

enrich$Namespace <- id2ns[enrich$ID]
enrich$Namespace <- ifelse(is.na(enrich$Namespace), "unknown", enrich$Namespace)

# Short namespace labels
ns_map <- c(
  "biological_process"  = "BP",
  "cellular_component"  = "CC",
  "molecular_function"  = "MF"
)
enrich$NS_short <- ifelse(
  enrich$Namespace %in% names(ns_map),
  ns_map[enrich$Namespace],
  "Other"
)

cat(sprintf("  Enriched GO terms: %d\n", nrow(enrich)))
cat(sprintf("  Annotated GO terms: %d\n", nrow(annot)))

# ============================================================
# 2. Figure 1: Bubble Plot (Top 20 GO terms)
# ============================================================
cat("\n[2/4] Generating Bubble Plot...\n")

top20 <- enrich %>%
  arrange(p.adjust) %>%
  head(20) %>%
  mutate(
    Label = reorder(str_wrap(Term, width = 42), RichFactor)
  )

p_bubble <- ggplot(top20, aes(x = RichFactor, y = Label)) +
  geom_point(
    aes(size = Count, color = p.adjust, shape = NS_short),
    alpha = 0.9, stroke = 0.4
  ) +
  scale_color_gradientn(
    colors = c("#313695", "#4575B4", "#74ADD1", "#ABD9E9",
               "#FEE090", "#FDAE61", "#F46D43", "#D73027", "#A50026"),
    trans  = "log10",
    name   = "Adjusted\np-value",
    labels = scales::scientific
  ) +
  scale_size_continuous(
    range = c(1.5, 5),
    name  = "Gene Count"
  ) +
  scale_shape_manual(
    values = c("BP" = 16, "CC" = 15, "MF" = 18, "Other" = 17),
    name   = "Category"
  ) +
  labs(x = "Rich Factor", y = NULL) +
  theme_pub +
  theme(
    axis.text.y = element_text(size = 7.5, lineheight = 0.9),
    legend.position = "right",
    legend.box = "vertical"
  ) +
  guides(
    color = guide_colorbar(order = 1, barwidth = 1, barheight = 6),
    size  = guide_legend(order = 2),
    shape = guide_legend(order = 3)
  )

ggsave(
  paste0(output_prefix, "_bubble_plot.pdf"),
  p_bubble, width = 7.5, height = 6.5, dpi = 300, device = cairo_pdf
)
ggsave(
  paste0(output_prefix, "_bubble_plot.png"),
  p_bubble, width = 7.5, height = 6.5, dpi = 300
)
cat("  Bubble plot saved.\n")

# ============================================================
# 3. Figure 2: GO Hierarchy DAG
# ============================================================
cat("\n[3/4] Generating GO Hierarchy DAG...\n")

# ---- Abbreviation dictionary removed: using numbered labels instead ----

abbreviate_name <- function(name) {
  if (is.na(name) || name == "NOT FOUND") return(name)
  return(name)
}

# ---- Build edge list from annotation file ----
padj_map  <- setNames(enrich$p.adjust, enrich$ID)
count_map <- setNames(enrich$Count, enrich$ID)

edges <- data.frame(from = character(), to = character(), stringsAsFactors = FALSE)

for (i in seq_len(nrow(annot))) {
  gid     <- annot$GO_ID[i]
  status  <- annot$Status[i]
  parents <- annot$Parents[i]

  if (status == "Not Found" || is.na(parents) || parents == "") next

  parent_ids <- str_extract_all(parents, "GO:\\d+")[[1]]
  for (pid in parent_ids) {
    edges <- rbind(edges, data.frame(from = pid, to = gid, stringsAsFactors = FALSE))
  }
}

# ---- Select top 25 enriched + direct parents + grandparents ----
top25_ids <- enrich %>% arrange(p.adjust) %>% head(25) %>% pull(ID)

core_nodes <- character()
for (gid in top25_ids) {
  core_nodes <- c(core_nodes, gid)
  # Parents
  parents <- edges$from[edges$to == gid]
  core_nodes <- c(core_nodes, parents)
  # Grandparents
  for (pid in parents) {
    gps <- edges$from[edges$to == pid]
    core_nodes <- c(core_nodes, gps)
  }
}
core_nodes <- unique(core_nodes)

# Filter edges to core nodes
edges_sub <- edges %>%
  filter(from %in% core_nodes & to %in% core_nodes)

# ---- Build igraph and node attributes ----
g <- graph_from_data_frame(edges_sub, directed = TRUE)

all_nodes <- V(g)$name
V(g)$is_enriched <- all_nodes %in% names(padj_map)
V(g)$padj <- ifelse(V(g)$is_enriched, padj_map[all_nodes], NA)
V(g)$node_size <- ifelse(V(g)$is_enriched, 10, 6)

# ---- Assign numeric labels: enriched first (sorted by padj), then parental ----
enriched_in_g  <- all_nodes[all_nodes %in% names(padj_map)]
parental_in_g  <- all_nodes[!(all_nodes %in% names(padj_map))]
enriched_in_g  <- enriched_in_g[order(padj_map[enriched_in_g])]
parental_in_g  <- sort(parental_in_g)
ordered_nodes  <- c(enriched_in_g, parental_in_g)

node_num_map <- setNames(seq_along(ordered_nodes), ordered_nodes)
V(g)$label <- as.character(node_num_map[all_nodes])

# ---- Save DAG node legend to txt ----
raw_names <- id2name[ordered_nodes]
raw_names[is.na(raw_names)] <- ordered_nodes[is.na(raw_names)]
node_types <- ifelse(ordered_nodes %in% names(padj_map), "Enriched", "Parental")

legend_df <- data.frame(
  No.  = seq_along(ordered_nodes),
  GO_ID = ordered_nodes,
  Term  = raw_names,
  Type  = node_types,
  stringsAsFactors = FALSE
)

legend_file <- paste0(output_prefix, "_DAG_node_legend.txt")
sink(legend_file)
cat("DAG Node Legend for C. radians GO Enrichment\n")
cat(strrep("=", 60), "\n")
cat(sprintf("%-6s%-14s%-55s%s\n", "No.", "GO ID", "Term Name", "Type"))
cat(strrep("-", 60), "\n")
for (i in seq_len(nrow(legend_df))) {
  cat(sprintf("%-6d%-14s%-55s%s\n",
              legend_df$No.[i], legend_df$GO_ID[i],
              legend_df$Term[i], legend_df$Type[i]))
}
sink()
cat(sprintf("  DAG node legend saved: %s\n", legend_file))

# ---- Font color: white for enriched (dark bg), dark for parental (light bg) ----
V(g)$label_color <- ifelse(V(g)$is_enriched, "white", "#333333")

# ---- Plot DAG with ggraph ----
set.seed(42)
p_dag <- ggraph(g, layout = "sugiyama") +
  geom_edge_link(
    arrow = arrow(length = unit(2.5, "mm"), type = "closed"),
    edge_color = "#AAAAAA",
    edge_width = 0.4,
    alpha = 0.5,
    end_cap = circle(6, "mm"),
    start_cap = circle(6, "mm")
  ) +
  geom_node_point(
    aes(
      size = node_size,
      fill = padj
    ),
    shape = 21,
    stroke = 0.6,
    color = ifelse(V(g)$is_enriched, "black", "#999999")
  ) +
  geom_node_text(
    aes(label = label, color = label_color),
    size = 4,
    fontface = "bold",
    family = "sans"
  ) +
  scale_fill_gradientn(
    colors = c("#313695", "#4575B4", "#74ADD1", "#ABD9E9",
               "#FEE090", "#FDAE61", "#F46D43", "#D73027", "#A50026"),
    trans  = "log10",
    na.value = "#E0E0E0",
    name   = "Adjusted\np-value",
    labels = scales::scientific
  ) +
  scale_color_identity() +
  scale_size_identity() +
  theme_void(base_family = "sans") +
  theme(
    legend.position  = "right",
    legend.title     = element_text(size = 9, face = "bold"),
    legend.text      = element_text(size = 7.5),
    plot.margin      = margin(10, 10, 10, 10)
  )

ggsave(
  paste0(output_prefix, "_GO_DAG.pdf"),
  p_dag, width = 14, height = 11, dpi = 300, device = cairo_pdf
)
ggsave(
  paste0(output_prefix, "_GO_DAG.png"),
  p_dag, width = 14, height = 11, dpi = 300
)
cat("  GO DAG saved.\n")

# ============================================================
# 4. Save figure titles to txt
# ============================================================
cat("\n[4/4] Saving figure titles...\n")

titles_file <- paste0(output_prefix, "_figure_titles.txt")
writeLines(c(
  "Figure Titles for C. radians GO Enrichment Analysis",
  strrep("=", 60),
  "",
  "Figure 1 (Bubble Plot):",
  paste0(
    "Top 20 enriched GO terms of expanded gene families ",
    "in Coprinellus radians. Bubble size represents gene count, ",
    "color indicates adjusted p-value, and shape denotes GO ",
    "category (BP: biological process; CC: cellular component; ",
    "MF: molecular function)."
  ),
  "",
  "Figure 2 (GO DAG):",
  paste0(
    "Directed acyclic graph (DAG) showing the hierarchical ",
    "relationships among enriched GO terms in expanded gene ",
    "families of Coprinellus radians. Colored nodes represent ",
    "significantly enriched GO terms (color mapped to adjusted ",
    "p-value), and grey nodes represent parental GO terms ",
    "providing hierarchical context."
  )
), con = titles_file)

cat(sprintf("  Titles saved: %s\n", titles_file))

cat("\n============================================================\n")
cat("All figures generated successfully!\n")
cat(sprintf("  - %s_bubble_plot.pdf / .png\n", output_prefix))
cat(sprintf("  - %s_GO_DAG.pdf / .png\n", output_prefix))
cat(sprintf("  - %s_DAG_node_legend.txt\n", output_prefix))
cat(sprintf("  - %s_figure_titles.txt\n", output_prefix))
cat("============================================================\n")

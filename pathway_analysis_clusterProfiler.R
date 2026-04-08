#!/usr/bin/env Rscript
################################################################################
# Genome-wide Pathway Analysis using clusterProfiler
# Input: eggNOG-mapper annotation file
# Output: GO, KEGG, and COG/KOG pathway statistics and visualizations
################################################################################

# Load required libraries
library(tidyverse)
library(clusterProfiler)
library(enrichplot)
library(ggplot2)
library(scales)

################################################################################
# 1. DATA LOADING AND PREPROCESSING
################################################################################

cat("Loading eggNOG-mapper annotation file...\n")

# Read eggNOG-mapper annotation file
# Skip comment lines starting with #
emapper_file <- "20220035.emapper.annotations"
annotations <- read.delim(emapper_file, 
                          header = TRUE, 
                          comment.char = "#", 
                          stringsAsFactors = FALSE,
                          quote = "")

cat(sprintf("Loaded %d gene annotations\n", nrow(annotations)))

# Display column names
cat("\nColumn names in annotation file:\n")
print(colnames(annotations))

################################################################################
# 2. PREPARE GO ANNOTATIONS
################################################################################

cat("\n=====================================\n")
cat("Processing GO annotations...\n")
cat("=====================================\n")

# Extract genes with GO annotations
go_data <- annotations %>%
  filter(GOs != "-" & GOs != "") %>%
  select(query, GOs)

cat(sprintf("Genes with GO annotations: %d\n", nrow(go_data)))

# Parse GO terms (format: GO:0000001,GO:0000002,...)
go_list <- list()
for(i in 1:nrow(go_data)) {
  gene_id <- go_data$query[i]
  go_terms <- unlist(strsplit(go_data$GOs[i], ","))
  go_list[[gene_id]] <- go_terms
}

# Create TERM2GENE data frame for GO analysis
go_term2gene <- data.frame(
  GO = unlist(go_list),
  gene = rep(names(go_list), sapply(go_list, length)),
  stringsAsFactors = FALSE
)

# Remove duplicates
go_term2gene <- unique(go_term2gene)

cat(sprintf("Total GO term-gene associations: %d\n", nrow(go_term2gene)))
cat(sprintf("Unique GO terms: %d\n", length(unique(go_term2gene$GO))))

# Get GO term descriptions
go_terms <- unique(go_term2gene$GO)
go_term2name <- data.frame(
  GO = go_terms,
  Description = go_terms,  # Will be filled with actual names if GO.db is available
  stringsAsFactors = FALSE
)

################################################################################
# 3. PREPARE KEGG ANNOTATIONS
################################################################################

cat("\n=====================================\n")
cat("Processing KEGG annotations...\n")
cat("=====================================\n")

# Extract genes with KEGG annotations
kegg_ko_data <- annotations %>%
  filter(KEGG_ko != "-" & KEGG_ko != "") %>%
  select(query, KEGG_ko)

cat(sprintf("Genes with KEGG KO annotations: %d\n", nrow(kegg_ko_data)))

# Parse KEGG KO terms
kegg_ko_list <- list()
for(i in 1:nrow(kegg_ko_data)) {
  gene_id <- kegg_ko_data$query[i]
  ko_terms <- unlist(strsplit(kegg_ko_data$KEGG_ko[i], ","))
  kegg_ko_list[[gene_id]] <- ko_terms
}

# Create TERM2GENE data frame for KEGG KO
kegg_ko_term2gene <- data.frame(
  KO = unlist(kegg_ko_list),
  gene = rep(names(kegg_ko_list), sapply(kegg_ko_list, length)),
  stringsAsFactors = FALSE
)
kegg_ko_term2gene <- unique(kegg_ko_term2gene)

cat(sprintf("Total KEGG KO term-gene associations: %d\n", nrow(kegg_ko_term2gene)))
cat(sprintf("Unique KEGG KO terms: %d\n", length(unique(kegg_ko_term2gene$KO))))

# Extract KEGG pathway annotations
kegg_pathway_data <- annotations %>%
  filter(KEGG_Pathway != "-" & KEGG_Pathway != "") %>%
  select(query, KEGG_Pathway)

cat(sprintf("Genes with KEGG pathway annotations: %d\n", nrow(kegg_pathway_data)))

# Parse KEGG pathways
kegg_pathway_list <- list()
for(i in 1:nrow(kegg_pathway_data)) {
  gene_id <- kegg_pathway_data$query[i]
  pathways <- unlist(strsplit(kegg_pathway_data$KEGG_Pathway[i], ","))
  kegg_pathway_list[[gene_id]] <- pathways
}

# Create TERM2GENE data frame for KEGG pathways
kegg_pathway_term2gene <- data.frame(
  Pathway = unlist(kegg_pathway_list),
  gene = rep(names(kegg_pathway_list), sapply(kegg_pathway_list, length)),
  stringsAsFactors = FALSE
)
kegg_pathway_term2gene <- unique(kegg_pathway_term2gene)

cat(sprintf("Total KEGG pathway-gene associations: %d\n", nrow(kegg_pathway_term2gene)))
cat(sprintf("Unique KEGG pathways: %d\n", length(unique(kegg_pathway_term2gene$Pathway))))

################################################################################
# 4. PREPARE COG/KOG ANNOTATIONS
################################################################################

cat("\n=====================================\n")
cat("Processing COG/KOG annotations...\n")
cat("=====================================\n")

# Extract genes with COG category annotations
cog_data <- annotations %>%
  filter(COG_category != "-" & COG_category != "") %>%
  select(query, COG_category, Description)

cat(sprintf("Genes with COG category annotations: %d\n", nrow(cog_data)))

# Parse COG categories (single letters or combinations like "KLT")
cog_list <- list()
for(i in 1:nrow(cog_data)) {
  gene_id <- cog_data$query[i]
  cog_cats <- unlist(strsplit(cog_data$COG_category[i], ""))
  cog_list[[gene_id]] <- cog_cats
}

# Create TERM2GENE data frame for COG
cog_term2gene <- data.frame(
  COG = unlist(cog_list),
  gene = rep(names(cog_list), sapply(cog_list, length)),
  stringsAsFactors = FALSE
)
cog_term2gene <- unique(cog_term2gene)

cat(sprintf("Total COG category-gene associations: %d\n", nrow(cog_term2gene)))
cat(sprintf("Unique COG categories: %d\n", length(unique(cog_term2gene$COG))))

# COG category descriptions
cog_descriptions <- data.frame(
  COG = c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", 
          "N", "O", "P", "Q", "R", "S", "T", "U", "V", "W", "X", "Y", "Z"),
  Description = c(
    "RNA processing and modification",
    "Chromatin structure and dynamics",
    "Energy production and conversion",
    "Cell cycle control, cell division, chromosome partitioning",
    "Amino acid transport and metabolism",
    "Nucleotide transport and metabolism",
    "Carbohydrate transport and metabolism",
    "Coenzyme transport and metabolism",
    "Lipid transport and metabolism",
    "Translation, ribosomal structure and biogenesis",
    "Transcription",
    "Replication, recombination and repair",
    "Cell wall/membrane/envelope biogenesis",
    "Cell motility",
    "Posttranslational modification, protein turnover, chaperones",
    "Inorganic ion transport and metabolism",
    "Secondary metabolites biosynthesis, transport and catabolism",
    "General function prediction only",
    "Function unknown",
    "Signal transduction mechanisms",
    "Intracellular trafficking, secretion, and vesicular transport",
    "Defense mechanisms",
    "Extracellular structures",
    "Mobilome: prophages, transposons",
    "Nuclear structure",
    "Cytoskeleton"
  ),
  stringsAsFactors = FALSE
)

################################################################################
# 5. GENOME-WIDE STATISTICS
################################################################################

cat("\n=====================================\n")
cat("GENOME-WIDE ANNOTATION STATISTICS\n")
cat("=====================================\n")

# Overall statistics
total_genes <- nrow(annotations)
genes_with_go <- length(unique(go_term2gene$gene))
genes_with_kegg_ko <- length(unique(kegg_ko_term2gene$gene))
genes_with_kegg_pathway <- length(unique(kegg_pathway_term2gene$gene))
genes_with_cog <- length(unique(cog_term2gene$gene))

cat(sprintf("\nTotal genes: %d\n", total_genes))
cat(sprintf("Genes with GO annotations: %d (%.2f%%)\n", 
            genes_with_go, 100 * genes_with_go / total_genes))
cat(sprintf("Genes with KEGG KO annotations: %d (%.2f%%)\n", 
            genes_with_kegg_ko, 100 * genes_with_kegg_ko / total_genes))
cat(sprintf("Genes with KEGG pathway annotations: %d (%.2f%%)\n", 
            genes_with_kegg_pathway, 100 * genes_with_kegg_pathway / total_genes))
cat(sprintf("Genes with COG annotations: %d (%.2f%%)\n", 
            genes_with_cog, 100 * genes_with_cog / total_genes))

################################################################################
# 6. GO TERM DISTRIBUTION
################################################################################

cat("\n=====================================\n")
cat("GO Term Distribution Analysis\n")
cat("=====================================\n")

# Count genes per GO term
go_counts <- go_term2gene %>%
  group_by(GO) %>%
  summarise(count = n()) %>%
  arrange(desc(count))

cat(sprintf("\nTop 20 GO terms by gene count:\n"))
print(head(go_counts, 20))

# Save GO statistics
write.table(go_counts, 
            "GO_term_statistics.txt", 
            sep = "\t", 
            row.names = FALSE, 
            quote = FALSE)

# Plot GO term distribution
pdf("GO_term_distribution.pdf", width = 10, height = 6)
go_plot_data <- head(go_counts, 30)
ggplot(go_plot_data, aes(x = reorder(GO, count), y = count)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() +
  labs(title = "Top 30 GO Terms by Gene Count",
       x = "GO Term",
       y = "Number of Genes") +
  theme_bw() +
  theme(axis.text.y = element_text(size = 8))
dev.off()

################################################################################
# 7. KEGG PATHWAY DISTRIBUTION
################################################################################

cat("\n=====================================\n")
cat("KEGG Pathway Distribution Analysis\n")
cat("=====================================\n")

# Count genes per KEGG KO
kegg_ko_counts <- kegg_ko_term2gene %>%
  group_by(KO) %>%
  summarise(count = n()) %>%
  arrange(desc(count))

cat(sprintf("\nTop 20 KEGG KO terms by gene count:\n"))
print(head(kegg_ko_counts, 20))

write.table(kegg_ko_counts, 
            "KEGG_KO_statistics.txt", 
            sep = "\t", 
            row.names = FALSE, 
            quote = FALSE)

# Count genes per KEGG pathway
kegg_pathway_counts <- kegg_pathway_term2gene %>%
  group_by(Pathway) %>%
  summarise(count = n()) %>%
  arrange(desc(count))

cat(sprintf("\nTop 20 KEGG pathways by gene count:\n"))
print(head(kegg_pathway_counts, 20))

write.table(kegg_pathway_counts, 
            "KEGG_pathway_statistics.txt", 
            sep = "\t", 
            row.names = FALSE, 
            quote = FALSE)

# Plot KEGG pathway distribution
pdf("KEGG_pathway_distribution.pdf", width = 12, height = 8)
kegg_plot_data <- head(kegg_pathway_counts, 30)
ggplot(kegg_plot_data, aes(x = reorder(Pathway, count), y = count)) +
  geom_bar(stat = "identity", fill = "darkgreen") +
  coord_flip() +
  labs(title = "Top 30 KEGG Pathways by Gene Count",
       x = "KEGG Pathway",
       y = "Number of Genes") +
  theme_bw() +
  theme(axis.text.y = element_text(size = 7))
dev.off()

################################################################################
# 8. COG CATEGORY DISTRIBUTION
################################################################################

cat("\n=====================================\n")
cat("COG Category Distribution Analysis\n")
cat("=====================================\n")

# Count genes per COG category
cog_counts <- cog_term2gene %>%
  group_by(COG) %>%
  summarise(count = n()) %>%
  left_join(cog_descriptions, by = "COG") %>%
  arrange(desc(count))

cat(sprintf("\nCOG category distribution:\n"))
print(cog_counts)

write.table(cog_counts, 
            "COG_category_statistics.txt", 
            sep = "\t", 
            row.names = FALSE, 
            quote = FALSE)

# Plot COG category distribution
pdf("COG_category_distribution.pdf", width = 12, height = 8)
ggplot(cog_counts, aes(x = reorder(paste(COG, "-", Description), count), y = count)) +
  geom_bar(stat = "identity", fill = "coral") +
  coord_flip() +
  labs(title = "COG Category Distribution",
       x = "COG Category",
       y = "Number of Genes") +
  theme_bw() +
  theme(axis.text.y = element_text(size = 9))
dev.off()

# Pie chart for COG categories
pdf("COG_category_pie_chart.pdf", width = 10, height = 10)
cog_counts_filtered <- cog_counts %>% filter(count >= 10)  # Show categories with >= 10 genes
ggplot(cog_counts_filtered, aes(x = "", y = count, fill = paste(COG, "-", Description))) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar("y", start = 0) +
  labs(title = "COG Category Distribution (≥10 genes)",
       fill = "COG Category") +
  theme_void() +
  theme(legend.text = element_text(size = 8))
dev.off()

################################################################################
# 9. ANNOTATION OVERVIEW SUMMARY
################################################################################

cat("\n=====================================\n")
cat("Creating annotation overview summary\n")
cat("=====================================\n")

# Create summary table
summary_df <- data.frame(
  Category = c("Total Genes", 
               "GO Annotated", 
               "KEGG KO Annotated", 
               "KEGG Pathway Annotated", 
               "COG Annotated",
               "Unique GO Terms",
               "Unique KEGG KOs",
               "Unique KEGG Pathways",
               "Unique COG Categories"),
  Count = c(total_genes,
            genes_with_go,
            genes_with_kegg_ko,
            genes_with_kegg_pathway,
            genes_with_cog,
            length(unique(go_term2gene$GO)),
            length(unique(kegg_ko_term2gene$KO)),
            length(unique(kegg_pathway_term2gene$Pathway)),
            length(unique(cog_term2gene$COG))),
  Percentage = c(100,
                 100 * genes_with_go / total_genes,
                 100 * genes_with_kegg_ko / total_genes,
                 100 * genes_with_kegg_pathway / total_genes,
                 100 * genes_with_cog / total_genes,
                 NA, NA, NA, NA)
)

print(summary_df)

write.table(summary_df, 
            "annotation_summary.txt", 
            sep = "\t", 
            row.names = FALSE, 
            quote = FALSE)

# Plot annotation coverage
pdf("annotation_coverage.pdf", width = 10, height = 6)
coverage_data <- summary_df[2:5, ]
ggplot(coverage_data, aes(x = Category, y = Percentage)) +
  geom_bar(stat = "identity", fill = "skyblue", alpha = 0.7) +
  geom_text(aes(label = sprintf("%.1f%%\n(%d genes)", Percentage, Count)), 
            vjust = -0.5, size = 4) +
  ylim(0, 100) +
  labs(title = "Genome Annotation Coverage",
       x = "Annotation Type",
       y = "Percentage of Genes (%)") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

################################################################################
# 10. EXPORT TERM2GENE MAPPINGS FOR ENRICHMENT ANALYSIS
################################################################################

cat("\n=====================================\n")
cat("Exporting TERM2GENE mappings\n")
cat("=====================================\n")

# Save TERM2GENE files for future enrichment analysis
write.table(go_term2gene, 
            "GO_TERM2GENE.txt", 
            sep = "\t", 
            row.names = FALSE, 
            quote = FALSE)

write.table(kegg_ko_term2gene, 
            "KEGG_KO_TERM2GENE.txt", 
            sep = "\t", 
            row.names = FALSE, 
            quote = FALSE)

write.table(kegg_pathway_term2gene, 
            "KEGG_PATHWAY_TERM2GENE.txt", 
            sep = "\t", 
            row.names = FALSE, 
            quote = FALSE)

write.table(cog_term2gene, 
            "COG_TERM2GENE.txt", 
            sep = "\t", 
            row.names = FALSE, 
            quote = FALSE)

cat("\nTERM2GENE mappings saved for:\n")
cat("  - GO_TERM2GENE.txt\n")
cat("  - KEGG_KO_TERM2GENE.txt\n")
cat("  - KEGG_PATHWAY_TERM2GENE.txt\n")
cat("  - COG_TERM2GENE.txt\n")

################################################################################
# 11. EXAMPLE ENRICHMENT ANALYSIS (optional)
################################################################################

cat("\n=====================================\n")
cat("Example: Enrichment Analysis Setup\n")
cat("=====================================\n")

cat("\nTo perform enrichment analysis on a gene list, use:\n\n")
cat("# Example for GO enrichment:\n")
cat("gene_list <- c('gene1', 'gene2', 'gene3')  # Your genes of interest\n")
cat("universe <- unique(go_term2gene$gene)  # Background genes\n")
cat("enrichGO_result <- enricher(gene = gene_list,\n")
cat("                            universe = universe,\n")
cat("                            TERM2GENE = go_term2gene,\n")
cat("                            pvalueCutoff = 0.05,\n")
cat("                            qvalueCutoff = 0.2)\n\n")

cat("# Example for KEGG enrichment:\n")
cat("enrichKEGG_result <- enricher(gene = gene_list,\n")
cat("                              universe = unique(kegg_pathway_term2gene$gene),\n")
cat("                              TERM2GENE = kegg_pathway_term2gene,\n")
cat("                              pvalueCutoff = 0.05,\n")
cat("                              qvalueCutoff = 0.2)\n\n")

cat("# Example for COG enrichment:\n")
cat("enrichCOG_result <- enricher(gene = gene_list,\n")
cat("                             universe = unique(cog_term2gene$gene),\n")
cat("                             TERM2GENE = cog_term2gene,\n")
cat("                             TERM2NAME = cog_descriptions,\n")
cat("                             pvalueCutoff = 0.05,\n")
cat("                             qvalueCutoff = 0.2)\n\n")

################################################################################
# 12. COMPLETION
################################################################################

cat("\n=====================================\n")
cat("ANALYSIS COMPLETE!\n")
cat("=====================================\n")

cat("\nGenerated files:\n")
cat("  1. GO_term_statistics.txt\n")
cat("  2. GO_term_distribution.pdf\n")
cat("  3. KEGG_KO_statistics.txt\n")
cat("  4. KEGG_pathway_statistics.txt\n")
cat("  5. KEGG_pathway_distribution.pdf\n")
cat("  6. COG_category_statistics.txt\n")
cat("  7. COG_category_distribution.pdf\n")
cat("  8. COG_category_pie_chart.pdf\n")
cat("  9. annotation_summary.txt\n")
cat(" 10. annotation_coverage.pdf\n")
cat(" 11. GO_TERM2GENE.txt\n")
cat(" 12. KEGG_KO_TERM2GENE.txt\n")
cat(" 13. KEGG_PATHWAY_TERM2GENE.txt\n")
cat(" 14. COG_TERM2GENE.txt\n")

cat("\nAll outputs saved to current directory.\n")

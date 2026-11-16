#!/usr/bin/env Rscript
# =============================================================================
# ITS phylogenetic validation of Coprinellus domesticus 20220035 and
# Coprinellus radians 20210058
#
# Pipeline:
#   1. Combine the two study strains' Sanger ITS sequences with curated
#      reference sequences (type material / culture-collection / voucher
#      specimens from taxonomic studies; Parasola + Coprinopsis outgroups).
#   2. Multiple sequence alignment with MUSCLE v5.3.
#   3. Trim alignment columns containing > 50% gaps.
#   4. Maximum-likelihood tree with IQ-TREE 3 (ModelFinder + 1000 SH-aLRT
#      + 1000 ultrafast bootstrap replicates).
#   5. Root on the outgroup and produce publication-quality figures
#      (PDF + 300 dpi PNG) and tree files.
#
# Usage: set the working directory to this script's folder, then
#   Rscript ITS_phylogenetic_tree.R
# External tools (bundled): tools/muscle.exe, tools/iqtree3.exe
# R packages required: ape (CRAN)
# =============================================================================

suppressPackageStartupMessages(library(ape))

## ---- 0. Paths --------------------------------------------------------------
# optional args: Rscript ITS_phylogenetic_tree.R [keep_list.txt] [out_suffix]
args <- commandArgs(trailingOnly = TRUE)
keep_file <- if (length(args) >= 1) args[1] else NA_character_
out_suffix <- if (length(args) >= 2) args[2] else ""

study_fasta  <- "Coprinellus.fasta"                 # two study strains (Sanger ITS)
ref_fasta    <- "reference_ITS_sequences.fasta"     # curated NCBI references
ref_info     <- "reference_ITS_sequences_info.tsv"  # accession metadata
muscle_bin   <- file.path("tools", "muscle.exe")
iqtree_bin   <- file.path("tools", "iqtree3.exe")

out_dir      <- paste0("results", out_suffix)
dir.create(out_dir, showWarnings = FALSE)

log_con <- file(file.path(out_dir, "pipeline_log.txt"), open = "wt")
sink(log_con, split = TRUE)
cat("ITS phylogenetic tree pipeline\n")
cat("Date:", as.character(Sys.Date()), "\n\n")

stopifnot(file.exists(study_fasta), file.exists(ref_fasta),
          file.exists(muscle_bin), file.exists(iqtree_bin))

## ---- 1. Combine and rename sequences ---------------------------------------
read_fasta <- function(f) {
  lines <- readLines(f)
  hdr   <- grep("^>", lines)
  ends  <- c(hdr[-1] - 1L, length(lines))
  setNames(Map(function(s, e) paste(lines[(s + 1):e], collapse = ""), hdr, ends),
           sub("^>", "", lines[hdr]))
}

study <- read_fasta(study_fasta)
refs  <- read_fasta(ref_fasta)

# rename study strains to unambiguous tree IDs
study_ids <- c("ThisStudy_Cdomesticus_20220035", "ThisStudy_Cradians_20210058")
names(study) <- study_ids
names(refs)  <- sub("\\s.*", "", names(refs))      # keep accession only

# optional accession keep-list (one accession per line)
if (!is.na(keep_file)) {
  keep <- trimws(readLines(keep_file))
  keep <- keep[nzchar(keep)]
  refs <- refs[names(refs) %in% keep]
  cat("Keep-list applied:", keep_file, "->", length(refs), "references kept\n")
}

all_seqs <- c(study, refs)
all_fasta <- file.path(out_dir, "all_ITS_sequences.fasta")
writeLines(paste0(">", names(all_seqs), "\n", unlist(all_seqs)), all_fasta)
cat("Sequences combined:", length(all_seqs), "\n")

## ---- 2. Alignment (MUSCLE) --------------------------------------------------
# set env ITS_PLOT_ONLY=1 to reuse existing results and only redraw figures
plot_only <- Sys.getenv("ITS_PLOT_ONLY") == "1"
aln_raw <- file.path(out_dir, "ITS_alignment_raw.fasta")
trimmed_fasta <- file.path(out_dir, "ITS_alignment_trimmed.fasta")
prefix <- file.path(normalizePath(out_dir), "ITS_ML")

if (!plot_only) {
  cmd <- sprintf('"%s" -align "%s" -output "%s" -threads 4',
                 normalizePath(muscle_bin), normalizePath(all_fasta),
                 normalizePath(aln_raw))
  cat("MUSCLE:", cmd, "\n")
  status <- system(cmd, intern = TRUE)
  cat(paste(status, collapse = "\n"), "\n")
}

## ---- 3. Trim gappy columns ---------------------------------------------------
aln <- read.dna(aln_raw, format = "fasta")
gap_prop <- apply(as.character(aln) == "-", 2, mean)
keep <- gap_prop <= 0.5
aln_trim <- aln[, keep]
write.dna(aln_trim, trimmed_fasta, format = "fasta", nbcol = -1, colsep = "")
cat(sprintf("Alignment: %d columns -> %d after trimming (gap>50%% removed)\n",
            ncol(aln), ncol(aln_trim)))

## ---- 4. ML tree (IQ-TREE: ModelFinder + 1000 SH-aLRT + 1000 UFBoot) --------
if (!plot_only) {
  cmd <- sprintf('"%s" -s "%s" -m MFP -alrt 1000 -B 1000 -T 2 --prefix "%s" --redo',
                 normalizePath(iqtree_bin), normalizePath(trimmed_fasta), prefix)
  cat("IQ-TREE:", cmd, "\n")
  status <- system(cmd, intern = TRUE)
  cat(tail(status, 25), sep = "\n")
}

## ---- 5. Root tree, parse support -------------------------------------------
tree <- read.tree(paste0(prefix, ".treefile"))
outgroup <- c("MH855977.1", "MK397584.1")  # Parasola plicatilis + P. aporos (sister genus outgroup)
tree <- root(tree, outgroup = outgroup, resolve.root = TRUE)
tree <- ladderize(tree)

# IQ-TREE3 node labels: "SH-aLRT/UFBoot" -> keep UFBoot
split_support <- function(x) {
  parts <- strsplit(x, "/")
  vapply(parts, function(p) if (length(p) >= 2) p[2] else p[1], character(1))
}
ufboot <- suppressWarnings(as.numeric(split_support(tree$node.label)))

## ---- 6. Pretty tip labels (italic species names) -----------------------------
info <- read.delim(ref_info, stringsAsFactors = FALSE)
make_label_txt <- function(acc) {
  note <- info$label[match(acc, info$accession)]
  note <- sub(",.*", "", note)                     # drop source comment
  # keep only "G. epithet [subsp. epithet2]" + accession (genus abbreviated)
  tok <- strsplit(note, " ", fixed = TRUE)[[1]]
  tok[1] <- paste0(substr(tok[1], 1, 1), ".")      # Parasola -> P., etc.
  keep_t <- tok[1:2]
  if (length(tok) >= 4 && tok[3] == "subsp.") keep_t <- c(keep_t, tok[3:4])
  sprintf("%s (%s)", paste(keep_t, collapse = " "), acc)
}
# split "Genus epithet [subsp. epithet2] rest..." -> plotmath expression:
# species name in italics, voucher/accession info in roman
label_expr <- function(lab, study = FALSE) {
  tok <- strsplit(lab, " ", fixed = TRUE)[[1]]
  sp <- tok[1:2]
  rest <- tok[-(1:2)]
  if (length(rest) >= 2 && rest[1] == "subsp.") {
    sp_txt <- sprintf('italic("%s %s")~"subsp."~italic("%s")',
                      sp[1], sp[2], rest[2])
    rest <- rest[-(1:2)]
  } else if (study) {
    sp_txt <- sprintf('bolditalic("%s %s")', sp[1], sp[2])
  } else {
    sp_txt <- sprintf('italic("%s %s")', sp[1], sp[2])
  }
  if (length(rest) > 0) {
    rest_txt <- paste(rest, collapse = " ")
    rest_txt <- if (study) sprintf('bold("%s")', rest_txt)
                else sprintf('"%s"', rest_txt)
    parse(text = paste(sp_txt, rest_txt, sep = "~"))[[1]]
  } else {
    parse(text = sp_txt)[[1]]
  }
}
tip_txt <- vapply(tree$tip.label, function(tip) {
  if (tip == "ThisStudy_Cdomesticus_20220035")
    return("C. domesticus 20220035 (this study)")
  if (tip == "ThisStudy_Cradians_20210058")
    return("C. radians 20210058 (this study)")
  make_label_txt(tip)
}, character(1))
is_study <- grepl("this study", tip_txt)
tip_lab  <- mapply(label_expr, tip_txt, is_study)   # list of expressions

# one colour per species; an edge is coloured only if ALL tips below it
# belong to the same species, otherwise the backbone stays dark grey
sp_key <- sub(" \\(.*", "", tip_txt)                    # e.g. "C. radians"
sp_key[is_study] <- sub("^(\\S+ \\S+).*", "\\1", tip_txt[is_study])
uniq_sp <- unique(sp_key)
fix_cols <- c("C. radians" = "#C0392B", "C. domesticus" = "#1E8449",
              "P. aporos" = "#8A97A5", "P. plicatilis" = "#8A97A5")
other_sp <- setdiff(uniq_sp, names(fix_cols))
palette_sp <- c("#2471A3", "#8E44AD", "#D68910", "#17A589", "#BA4A00",
                "#7D3C98", "#2E86C1", "#9C640C", "#CB4335", "#148F77",
                "#AF601A", "#5D6D7E", "#6C3483", "#0E6655", "#7B241C",
                "#1A5276")
sp_cols <- c(fix_cols, setNames(palette_sp[seq_along(other_sp)], other_sp))
tip_col <- unname(sp_cols[sp_key])

## ---- 7. Identity check (for the response letter) ----------------------------
d <- as.matrix(dist.dna(aln_trim, model = "raw", pairwise.deletion = TRUE))
ident <- function(a, b) 100 * (1 - d[a, b])
cat(sprintf("\nPairwise identity (trimmed alignment):\n"))
report_ident <- function(ours, ref, ref_label) {
  if (ref %in% rownames(d))
    cat(sprintf("  %s vs %s : %.2f%%\n", ours, ref_label, ident(ours, ref)))
}
report_ident("ThisStudy_Cdomesticus_20220035", "MH856480.1", "C. domesticus CBS 175.49")
report_ident("ThisStudy_Cdomesticus_20220035", "OL355086.1", "C. domesticus HMJAU46377")
report_ident("ThisStudy_Cradians_20210058", "MH855978.1", "C. radians  CBS 187.39")
report_ident("ThisStudy_Cradians_20210058", "PP383778.1", "C. radians  TRTC175601")

## ---- 8. Publication-quality plot (circular fan layout) ------------------------
tip_cex <- 0.95          # tip label size
lab_expr <- as.expression(tip_lab)
lab_off <- 0.012         # label offset from tips (data units)

# edge colours: an edge gets the species colour only if ALL tips below it
# belong to that species; the mixed backbone stays dark grey
edge_col <- vapply(seq_len(nrow(tree$edge)), function(i) {
  node <- tree$edge[i, 2]
  tips <- if (node <= Ntip(tree)) tree$tip.label[node]
          else extract.clade(tree, node)$tip.label
  g <- unique(sp_key[match(tips, tree$tip.label)])
  if (length(g) == 1) unname(sp_cols[[g]]) else "#4D5656"
}, character(1))

# plot copy with expression tip labels (rotation handled by plot.phylo)
tree_p <- tree
tree_p$tip.label <- lab_expr

plot_tree <- function(dev_w_in) {
  par(mar = c(0, 0, 0, 0), oma = c(0, 0, 2.4, 0), xpd = NA)
  # pass 1: no labels, just to get fan geometry (tip angles, radius)
  plot(tree_p, type = "fan", show.tip.label = FALSE,
       edge.color = edge_col, edge.width = 1.8, no.margin = TRUE)
  lpp <- get("last_plot.phylo", envir = .PlotPhyloEnv)
  n <- Ntip(tree)
  tipx <- lpp$xx[1:n]; tipy <- lpp$yy[1:n]
  R <- max(sqrt(tipx^2 + tipy^2)) + lab_off
  th <- atan2(tipy, tipx)
  # label sizes in inches (plain-text approximation incl. italic width)
  w_in <- strwidth(tip_txt, cex = tip_cex, font = 3, units = "inches")
  h_in <- strheight(tip_txt, cex = tip_cex, font = 3, units = "inches")
  # rotated bounding box of each label, in inches beyond the tip point
  ext_h <- w_in * abs(cos(th)) + h_in * abs(sin(th))
  ext_v <- w_in * abs(sin(th)) + h_in * abs(cos(th))
  # required half-width K (data units): R_in + ext <= dev_w_in/2
  Kh <- R * dev_w_in / (dev_w_in - 2 * max(ext_h) - 0.15)
  Kv <- R * dev_w_in / (dev_w_in - 2 * max(ext_v) - 0.15)
  K <- max(Kh, Kv) * 1.02
  # pass 2: final centred plot with computed limits
  plot.new()
  plot(tree_p, type = "fan", show.tip.label = TRUE, cex = tip_cex,
       tip.color = tip_col, edge.color = edge_col, edge.width = 1.8,
       label.offset = lab_off, no.margin = TRUE, align.tip.label = TRUE,
       x.lim = c(-K, K), y.lim = c(-K, K))
  # bootstrap support as dots: >=95 solid dark, 70-94 grey
  show95 <- !is.na(ufboot) & ufboot >= 95
  show70 <- !is.na(ufboot) & ufboot >= 70 & ufboot < 95
  nodelabels(pch = 19, col = ifelse(show95, "#222222", NA), cex = 0.8)
  nodelabels(pch = 19, col = ifelse(show70, "#AAAAAA", NA), cex = 0.65)
  # scale bar in the empty top-left corner
  usr <- par("usr")
  sx <- usr[1] + 0.03 * diff(usr[1:2])
  sy <- usr[4] - 0.05 * diff(usr[3:4])
  segments(sx, sy, sx + 0.05, sy, lwd = 2.5)
  text(sx + 0.06, sy, "0.05 substitutions/site", pos = 4, cex = 0.9)
  # legend: bootstrap dots only (species colours are evident at the labels)
  legend("bottomright", inset = c(0.02, 0.02), bty = "n", cex = 0.95, xpd = NA,
         pch = 19, col = c("#222222", "#AAAAAA"),
         legend = c("UFBoot >= 95", "UFBoot 70-94"), pt.cex = c(0.9, 0.75))
  mtext("Maximum-likelihood tree of Coprinellus ITS sequences\n(IQ-TREE 3, ModelFinder + 1000 ultrafast bootstrap)",
        outer = TRUE, cex = 1.0, font = 1)
}

pdf(file.path(out_dir, "ITS_ML_tree.pdf"), width = 10, height = 10)
plot_tree(10)
dev.off()

png(file.path(out_dir, "ITS_ML_tree.png"), width = 3000, height = 3000, res = 300)
plot_tree(10)
dev.off()

# tree files for deposition
write.tree(tree, file.path(out_dir, "ITS_ML_rooted_with_UFBoot.nwk"))

cat("\nDone. Outputs in:", normalizePath(out_dir), "\n")
print(sessionInfo())
sink()
close(log_con)

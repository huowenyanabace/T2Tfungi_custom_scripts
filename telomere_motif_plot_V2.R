#!/usr/bin/env Rscript

# Draw chromosome-scale telomere motif extension plots from a FASTA genome.
#
# V2 uses the legacy approximate terminal-chain statistic by default:
#   * scan both the requested motif and its reverse complement at either end;
#   * allow up to one mismatch per motif-sized word;
#   * chain non-overlapping hits whose starts are no more than 12 bp apart;
#   * choose the longest chain (earliest on the 5' end, latest on the 3' end
#     when chains tie).
# Plot geometry, typography, colours, dimensions, and output devices are kept
# identical to telomere_motif_plot.R.
# Example:
#   Rscript telomere_motif_plot_V2.R --input 020.T2T.fasta --motif AACCCT --scan_len 2000 --output_prefix telomere_motif

usage <- function() {
  cat("
Usage:
  Rscript telomere_motif_plot_V2.R --input genome.fasta --motif AACCCT --scan_len 2000 [options]

Required/exposed parameters:
  --input FILE             Genome FASTA file. If omitted, the first .fasta/.fa/.fna in the working directory is used.
  --motif SEQ              Motif used for the 5' end scan. Default: AACCCT
  --scan_len INT           Number of bases scanned at each chromosome end. Default: 2000

Options:
  --right_motif SEQ        Motif used for the 3' end scan. Default: reverse-complement of --motif
  --output_prefix PREFIX   Output prefix. Default: telomere_motif
  --species TEXT           Species name shown in the title panel. Default: Species name
  --formats LIST           Comma-separated formats: pdf,png,svg. Default: pdf,png,svg
  --dpi INT                PNG resolution. Default: 300
  --scale_mb NUM           Mb per scale tick. Default: 1
  --cap_width NUM          Visual telomere cap width in Mb units. Default: 0.42
  --min_repeats INT        Repeat counts below this value are still tabulated but hidden in the plot. Default: 1
  --scan_mode MODE         approximate, anchored, or longest. Default: approximate
                           approximate = legacy terminal audit: accept the motif or its reverse
                                         complement with bounded mismatches and chain nearby hits.
                           anchored = count perfect tandem repeats touching the chromosome end.
                           longest  = count the longest perfect tandem run inside the terminal scan window.
  --max_mismatch INT       Maximum mismatches allowed per motif word in approximate mode. Default: 1
  --max_start_gap INT      Maximum distance between successive hit starts in approximate mode. Default: 12
  --help                   Show this help message.

Outputs:
  PREFIX_stats.tsv
  PREFIX_main.FORMAT
  PREFIX_legend.FORMAT
  PREFIX_title.FORMAT

In the plot label \"xN, M bp\", N is the number of motif hits in the selected
chain and M is N multiplied by the motif length. The genomic span between the
first and last approximate hits can be longer when the chain contains gaps.
\n")
}

parse_args <- function(args) {
  opts <- list(
    input = NA_character_,
    motif = "AACCCT",
    scan_len = 2000L,
    right_motif = NA_character_,
    output_prefix = "telomere_motif",
    species = "Species name",
    formats = "pdf,png,svg",
    dpi = 300L,
    scale_mb = 1,
    cap_width = 0.42,
    min_repeats = 1L,
    scan_mode = "approximate",
    max_mismatch = 1L,
    max_start_gap = 12L,
    help = FALSE
  )
  if (length(args) == 0) {
    return(opts)
  }
  i <- 1L
  while (i <= length(args)) {
    key <- args[[i]]
    if (key == "--help" || key == "-h") {
      opts$help <- TRUE
      i <- i + 1L
      next
    }
    if (!startsWith(key, "--")) {
      stop("Unexpected argument: ", key)
    }
    name <- sub("^--", "", key)
    if (!name %in% names(opts)) {
      stop("Unknown option: ", key)
    }
    if (i == length(args)) {
      stop("Missing value for option: ", key)
    }
    opts[[name]] <- args[[i + 1L]]
    i <- i + 2L
  }
  opts$motif <- toupper(gsub("\\s+", "", opts$motif))
  opts$right_motif <- if (is.na(opts$right_motif)) NA_character_ else toupper(gsub("\\s+", "", opts$right_motif))
  opts$scan_len <- as.integer(opts$scan_len)
  opts$dpi <- as.integer(opts$dpi)
  opts$scale_mb <- as.numeric(opts$scale_mb)
  opts$cap_width <- as.numeric(opts$cap_width)
  opts$min_repeats <- as.integer(opts$min_repeats)
  opts$scan_mode <- tolower(opts$scan_mode)
  opts$max_mismatch <- as.integer(opts$max_mismatch)
  opts$max_start_gap <- as.integer(opts$max_start_gap)
  opts$formats <- tolower(gsub("\\s+", "", opts$formats))
  opts
}

find_default_fasta <- function() {
  files <- list.files(pattern = "\\.(fasta|fa|fna)$", ignore.case = TRUE)
  if (length(files) == 0) {
    stop("No FASTA file was supplied and no .fasta/.fa/.fna file was found in the working directory.")
  }
  files[[1L]]
}

check_motif <- function(motif, arg_name) {
  if (is.na(motif) || nchar(motif) == 0) {
    stop(arg_name, " must not be empty.")
  }
  if (grepl("[^ACGTRYKMSWBDHVN]", motif)) {
    stop(arg_name, " contains unsupported characters. Use IUPAC DNA letters only.")
  }
  motif
}

rev_comp <- function(seq) {
  comp <- chartr("ACGTRYKMSWBDHVN", "TGCAYRMKSWVHDBN", toupper(seq))
  paste(rev(strsplit(comp, "", fixed = TRUE)[[1L]]), collapse = "")
}

read_fasta <- function(path) {
  lines <- readLines(path, warn = FALSE)
  lines <- trimws(lines)
  lines <- lines[nzchar(lines)]
  header_idx <- grep("^>", lines)
  if (length(header_idx) == 0) {
    stop("No FASTA records were found in: ", path)
  }
  ends <- c(header_idx[-1L] - 1L, length(lines))
  names <- sub("^>", "", lines[header_idx])
  names <- vapply(strsplit(names, "\\s+"), `[`, character(1), 1L)
  seqs <- character(length(header_idx))
  for (i in seq_along(header_idx)) {
    seq_start <- header_idx[[i]] + 1L
    seq_end <- ends[[i]]
    if (seq_start > seq_end) {
      seqs[[i]] <- ""
    } else {
      seqs[[i]] <- toupper(gsub("\\s+", "", paste(lines[seq_start:seq_end], collapse = "")))
    }
  }
  data.frame(name = names, seq = seqs, stringsAsFactors = FALSE)
}

prefix_repeat <- function(window, motif) {
  motif_len <- nchar(motif)
  window_len <- nchar(window)
  count <- 0L
  pos <- 1L
  while ((pos + motif_len - 1L) <= window_len &&
         substr(window, pos, pos + motif_len - 1L) == motif) {
    count <- count + 1L
    pos <- pos + motif_len
  }
  list(count = count, start = if (count > 0L) 1L else NA_integer_,
       end = if (count > 0L) count * motif_len else NA_integer_)
}

suffix_repeat <- function(window, motif) {
  motif_len <- nchar(motif)
  end_pos <- nchar(window)
  count <- 0L
  while ((end_pos - motif_len + 1L) >= 1L &&
         substr(window, end_pos - motif_len + 1L, end_pos) == motif) {
    count <- count + 1L
    end_pos <- end_pos - motif_len
  }
  list(count = count,
       start = if (count > 0L) end_pos + 1L else NA_integer_,
       end = if (count > 0L) nchar(window) else NA_integer_)
}

longest_tandem_run <- function(window, motif) {
  motif_len <- nchar(motif)
  window_len <- nchar(window)
  best_count <- 0L
  best_start <- NA_integer_
  best_end <- NA_integer_
  if (window_len < motif_len) {
    return(list(count = best_count, start = best_start, end = best_end))
  }
  for (start in seq_len(window_len - motif_len + 1L)) {
    count <- 0L
    pos <- start
    while ((pos + motif_len - 1L) <= window_len &&
           substr(window, pos, pos + motif_len - 1L) == motif) {
      count <- count + 1L
      pos <- pos + motif_len
    }
    if (count > best_count) {
      best_count <- count
      best_start <- start
      best_end <- pos - 1L
    }
  }
  list(count = best_count, start = best_start, end = best_end)
}

hamming_distance <- function(a, b) {
  if (nchar(a) != nchar(b)) {
    stop("Cannot calculate Hamming distance for strings of unequal length.")
  }
  a_chars <- strsplit(a, "", fixed = TRUE)[[1L]]
  b_chars <- strsplit(b, "", fixed = TRUE)[[1L]]
  sum(a_chars != b_chars)
}

approximate_terminal_chain <- function(window, motif, end_side,
                                       max_mismatch = 1L,
                                       max_start_gap = 12L) {
  motif <- toupper(motif)
  motif_len <- nchar(motif)
  window_len <- nchar(window)
  queries <- unique(c(motif, rev_comp(motif)))

  empty_result <- list(
    count = 0L,
    start = NA_integer_,
    end = NA_integer_,
    distance = NA_integer_
  )
  if (window_len < motif_len) {
    return(empty_result)
  }

  # Match the legacy audit exactly: after accepting a motif-sized word, move
  # forward by one full motif; otherwise move forward by one base.
  hit_starts <- integer(0)
  hit_ends <- integer(0)
  pos <- 1L
  last_start <- window_len - motif_len + 1L
  while (pos <= last_start) {
    word <- substr(window, pos, pos + motif_len - 1L)
    distance_to_query <- min(vapply(
      queries,
      function(query) hamming_distance(word, query),
      integer(1)
    ))
    if (distance_to_query <= max_mismatch) {
      hit_starts <- c(hit_starts, pos)
      hit_ends <- c(hit_ends, pos + motif_len - 1L)
      pos <- pos + motif_len
    } else {
      pos <- pos + 1L
    }
  }

  if (length(hit_starts) == 0L) {
    return(empty_result)
  }

  chain_id <- integer(length(hit_starts))
  chain_id[[1L]] <- 1L
  if (length(hit_starts) > 1L) {
    for (i in 2:length(hit_starts)) {
      starts_new_chain <- (hit_starts[[i]] - hit_starts[[i - 1L]]) > max_start_gap
      chain_id[[i]] <- chain_id[[i - 1L]] + as.integer(starts_new_chain)
    }
  }

  chain_indices <- split(seq_along(hit_starts), chain_id)
  chain_counts <- vapply(chain_indices, length, integer(1))
  max_count <- max(chain_counts)
  candidates <- which(chain_counts == max_count)

  if (end_side == "left") {
    first_starts <- vapply(chain_indices[candidates], function(idx) hit_starts[idx[[1L]]], integer(1))
    chosen <- candidates[[which.min(first_starts)]]
  } else if (end_side == "right") {
    last_ends <- vapply(chain_indices[candidates], function(idx) hit_ends[idx[[length(idx)]]], integer(1))
    chosen <- candidates[[which.max(last_ends)]]
  } else {
    stop("end_side must be 'left' or 'right'.")
  }

  idx <- chain_indices[[chosen]]
  best_start <- hit_starts[idx[[1L]]]
  best_end <- hit_ends[idx[[length(idx)]]]
  terminal_distance <- if (end_side == "left") best_start - 1L else window_len - best_end

  list(
    count = as.integer(length(idx)),
    start = as.integer(best_start),
    end = as.integer(best_end),
    distance = as.integer(terminal_distance)
  )
}

scan_telomeres <- function(fasta_df, motif_left, motif_right, scan_len, scan_mode,
                           max_mismatch = 1L, max_start_gap = 12L) {
  out <- vector("list", nrow(fasta_df))
  motif_left_len <- nchar(motif_left)
  motif_right_len <- nchar(motif_right)

  for (i in seq_len(nrow(fasta_df))) {
    chr <- fasta_df$name[[i]]
    seq <- fasta_df$seq[[i]]
    chr_len <- nchar(seq)
    this_scan <- min(scan_len, chr_len)
    left_window <- substr(seq, 1L, this_scan)
    right_window <- substr(seq, chr_len - this_scan + 1L, chr_len)

    if (scan_mode == "approximate") {
      left <- approximate_terminal_chain(
        left_window, motif_left, "left", max_mismatch, max_start_gap
      )
      right <- approximate_terminal_chain(
        right_window, motif_right, "right", max_mismatch, max_start_gap
      )
    } else if (scan_mode == "anchored") {
      left <- prefix_repeat(left_window, motif_left)
      right <- suffix_repeat(right_window, motif_right)
    } else if (scan_mode == "longest") {
      left <- longest_tandem_run(left_window, motif_left)
      right <- longest_tandem_run(right_window, motif_right)
    } else {
      stop("--scan_mode must be 'approximate', 'anchored', or 'longest'.")
    }

    left_start <- left$start
    left_end <- left$end
    right_start <- if (is.na(right$start)) NA_integer_ else chr_len - this_scan + right$start
    right_end <- if (is.na(right$end)) NA_integer_ else chr_len - this_scan + right$end

    out[[i]] <- data.frame(
      chromosome = chr,
      length_bp = chr_len,
      length_Mb = chr_len / 1e6,
      left_motif = motif_left,
      left_repeat_count = left$count,
      left_repeat_bp = left$count * motif_left_len,
      left_scan_bp = this_scan,
      left_run_start = left_start,
      left_run_end = left_end,
      right_motif = motif_right,
      right_repeat_count = right$count,
      right_repeat_bp = right$count * motif_right_len,
      right_scan_bp = this_scan,
      right_run_start = right_start,
      right_run_end = right_end,
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, out)
}

open_device <- function(filename, format, width, height, dpi) {
  if (format == "pdf") {
    pdf(filename, width = width, height = height, useDingbats = FALSE)
  } else if (format == "png") {
    png(filename, width = width * dpi, height = height * dpi, res = dpi)
  } else if (format == "svg") {
    svg(filename, width = width, height = height)
  } else {
    stop("Unsupported output format: ", format)
  }
}

semi_capsule <- function(x_center, y, width, height, side, col, border = NA) {
  theta <- if (side == "left") seq(pi / 2, 3 * pi / 2, length.out = 80L) else seq(-pi / 2, pi / 2, length.out = 80L)
  x <- x_center + width * cos(theta)
  yv <- y + (height / 2) * sin(theta)
  polygon(x, yv, col = col, border = border)
}

draw_main_plot <- function(stats, opts) {
  n <- nrow(stats)
  left_col <- "#2B719B"
  right_col <- "#C1534F"
  body_col <- "#B8CFE2"
  text_col <- "#333333"
  grey_col <- "#666666"
  scale_col <- "#333333"
  row_step <- 1
  bar_h <- 0.58
  cap_w <- opts$cap_width
  max_mb <- max(stats$length_Mb)
  tick_max <- ceiling(max_mb / opts$scale_mb) * opts$scale_mb
  y_pos <- rev(seq_len(n)) * row_step
  y_scale <- max(y_pos) + 1.2
  x_left <- -cap_w - max(2.7, tick_max * 0.38)
  x_right <- tick_max + cap_w + max(1.8, tick_max * 0.32)

  par(mar = c(0.6, 0.7, 0.7, 0.7), xaxs = "i", yaxs = "i")
  plot.new()
  plot.window(xlim = c(x_left, x_right), ylim = c(0.25, y_scale + 0.75))

  arrows(0, y_scale, tick_max, y_scale, length = 0.09, angle = 25,
         code = 2, lwd = 2, col = scale_col)
  ticks <- seq(0, tick_max, by = opts$scale_mb)
  segments(ticks, y_scale - 0.13, ticks, y_scale + 0.13, col = scale_col, lwd = 1.8)
  text(ticks, y_scale + 0.3, paste0(format(ticks, trim = TRUE), " Mb"),
       cex = 1.0, font = 2, col = scale_col)

  for (i in seq_len(n)) {
    y <- y_pos[[i]]
    len_mb <- stats$length_Mb[[i]]
    rect(0, y - bar_h / 2, len_mb, y + bar_h / 2, col = body_col, border = NA)
    semi_capsule(0, y, cap_w, bar_h, "left", left_col)
    semi_capsule(len_mb, y, cap_w, bar_h, "right", right_col)

    text(x_left + 0.02, y, stats$chromosome[[i]], adj = c(0, 0.5),
         cex = 1.05, font = 4, col = text_col)
    text(len_mb / 2, y, sprintf("%.3f Mb", stats$length_Mb[[i]]),
         cex = 0.95, font = 2, col = text_col)

    if (stats$left_repeat_count[[i]] >= opts$min_repeats) {
      text(-cap_w - 0.12, y,
           sprintf("x%d, %d bp", stats$left_repeat_count[[i]], stats$left_repeat_bp[[i]]),
           adj = c(1, 0.5), cex = 0.9, font = 2, col = left_col)
    }
    if (stats$right_repeat_count[[i]] >= opts$min_repeats) {
      text(len_mb + cap_w + 0.12, y,
           sprintf("x%d, %d bp", stats$right_repeat_count[[i]], stats$right_repeat_bp[[i]]),
           adj = c(0, 0.5), cex = 0.9, font = 2, col = right_col)
    }
  }
}

draw_legend_plot <- function(opts) {
  left_col <- "#2B719B"
  right_col <- "#C1534F"
  body_col <- "#B8CFE2"
  text_col <- "#333333"
  grey_col <- "#666666"

  par(mar = c(0.2, 0.2, 0.2, 0.2), xaxs = "i", yaxs = "i")
  plot.new()
  plot.window(xlim = c(0, 1), ylim = c(0, 1))
  rect(0.06, 0.06, 0.94, 0.94, border = "#999999", lwd = 1.5)
  text(0.5, 0.86, "Legend", cex = 1.7, font = 2, col = text_col)

  y <- c(0.67, 0.45, 0.23)
  cols <- c(left_col, right_col, body_col)
  labels <- c("5' telomere extension", "3' telomere extension", "Chromosome body")
  sublabels <- c(
    paste0("Motif: (", opts$motif, ")n"),
    paste0("Motif: (", opts$right_motif, ")n"),
    "Drawn to scale"
  )
  for (i in seq_along(y)) {
    rect(0.11, y[[i]] - 0.055, 0.25, y[[i]] + 0.055, col = cols[[i]], border = NA)
    text(0.3, y[[i]] + 0.035, labels[[i]], adj = c(0, 0.5), cex = 1.3, font = 2, col = text_col)
    text(0.3, y[[i]] - 0.05, sublabels[[i]], adj = c(0, 0.5), cex = 1.1,
         font = if (i == 3L) 3 else 1, col = grey_col)
  }
}

draw_title_plot <- function(opts) {
  par(mar = c(0, 0, 0, 0), xaxs = "i", yaxs = "i")
  plot.new()
  plot.window(xlim = c(0, 1), ylim = c(0, 1))
  text(0.5, 0.5,
       paste("Chromosome assembly and telomere extensions of", opts$species),
       cex = 1.15, font = 2, col = "#333333")
}

write_outputs <- function(stats, opts) {
  write.table(stats, paste0(opts$output_prefix, "_stats.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)

  formats <- unlist(strsplit(opts$formats, ",", fixed = TRUE))
  formats <- formats[nzchar(formats)]
  if (length(formats) == 0) {
    stop("No output formats supplied.")
  }

  main_height <- max(4.5, 2.1 + nrow(stats) * 0.72)
  for (fmt in formats) {
    main_file <- paste0(opts$output_prefix, "_main.", fmt)
    legend_file <- paste0(opts$output_prefix, "_legend.", fmt)
    title_file <- paste0(opts$output_prefix, "_title.", fmt)

    open_device(main_file, fmt, width = 10, height = main_height, dpi = opts$dpi)
    draw_main_plot(stats, opts)
    dev.off()

    open_device(legend_file, fmt, width = 4.8, height = 3.4, dpi = opts$dpi)
    draw_legend_plot(opts)
    dev.off()

    open_device(title_file, fmt, width = 8.0, height = 0.8, dpi = opts$dpi)
    draw_title_plot(opts)
    dev.off()
  }
}

main <- function() {
  opts <- parse_args(commandArgs(trailingOnly = TRUE))
  if (isTRUE(opts$help)) {
    usage()
    return(invisible(NULL))
  }
  if (is.na(opts$input)) {
    opts$input <- find_default_fasta()
  }
  if (!file.exists(opts$input)) {
    stop("Input FASTA does not exist: ", opts$input)
  }
  opts$motif <- check_motif(opts$motif, "--motif")
  if (is.na(opts$right_motif)) {
    opts$right_motif <- rev_comp(opts$motif)
  }
  opts$right_motif <- check_motif(opts$right_motif, "--right_motif")
  if (is.na(opts$scan_len) || opts$scan_len <= 0) {
    stop("--scan_len must be a positive integer.")
  }
  if (is.na(opts$dpi) || opts$dpi <= 0) {
    stop("--dpi must be a positive integer.")
  }
  if (is.na(opts$scale_mb) || opts$scale_mb <= 0) {
    stop("--scale_mb must be a positive number.")
  }
  if (is.na(opts$cap_width) || opts$cap_width <= 0) {
    stop("--cap_width must be a positive number.")
  }
  if (is.na(opts$max_mismatch) || opts$max_mismatch < 0) {
    stop("--max_mismatch must be a non-negative integer.")
  }
  if (is.na(opts$max_start_gap) || opts$max_start_gap <= 0) {
    stop("--max_start_gap must be a positive integer.")
  }

  fasta_df <- read_fasta(opts$input)
  stats <- scan_telomeres(
    fasta_df,
    opts$motif,
    opts$right_motif,
    opts$scan_len,
    opts$scan_mode,
    opts$max_mismatch,
    opts$max_start_gap
  )
  write_outputs(stats, opts)
  message("Done. Wrote outputs with prefix: ", opts$output_prefix)
}

main()

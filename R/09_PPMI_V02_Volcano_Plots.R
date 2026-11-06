# Input/output paths are configured via environment variables: PD_TRPV1_ROOT.
if (!nzchar(Sys.getenv("PD_TRPV1_ROOT", unset=""))) stop("Set PD_TRPV1_ROOT before running this script (see README.md).", call. = FALSE)

# REVIEWED static audit 2026-11-05; original: 09_PPMI_V02_Volcano_Plots.R
# Requires original project datasets; full R execution not validated here.
# Supplementary volcano plots with square-shaped panels
# for No-pain / Low-pain / High-pain PD vs HC using age/sex-adjusted DESeq2 results.
#
# Why previous plots looked vertically stretched:
# 1) The saved figure size was rectangular (width != height).
# 2) The legend on the right reduced panel width.
# 3) Even with x-range ~10 and y-range ~10, the plotting panel itself was not square.
#
# This script fixes that by:
# - using x = [-5, 5] and y = [0, 10]
# - forcing a square panel with theme(aspect.ratio = 1)
# - placing the legend at the bottom
# - saving each volcano as a square image

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
  library(grid)
})

if (!requireNamespace("org.Hs.eg.db", quietly = TRUE) ||
    !requireNamespace("AnnotationDbi", quietly = TRUE)) {
  stop(
    "Missing annotation packages. Run once:\n",
    "if (!requireNamespace('BiocManager', quietly=TRUE)) install.packages('BiocManager')\n",
    "BiocManager::install(c('org.Hs.eg.db','AnnotationDbi'))"
  )
}

# ----------------------
# paths and settings
# ----------------------
root <- Sys.getenv("PD_TRPV1_ROOT", unset="")
input_file <- file.path(root, "V02_DESeq2_AgeSex_Adjusted", "01_DESeq2_AgeSex_adjusted_ALL_genes_vs_HC.csv")
out_dir <- file.path(root, "Supplementary_Volcano_AgeSex_Adjusted_WIDE_GRID")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

X_MIN <- -5
X_MAX <- 5
Y_MIN <- 0
Y_MAX <- 10

COL_UP   <- "#FA8072"
COL_DOWN <- "#648FFF"
COL_NS   <- "gray75"

GENES18 <- c(
  "PLCB4", "MAPK10", "MAPK12", "TRPV1", "MAPK13", "IL1R1",
  "MAPK11", "IL1RAP", "PLCB3", "PLCG2", "PLCB1", "IL1B",
  "MAPK14", "MAPK9", "PLCB2", "MAPK8", "NGF", "PLCG1"
)

COMPARISONS <- c("No pain PD vs HC", "Low pain PD vs HC", "High pain PD vs HC")
COMP_LABELS <- c("No-pain PD vs HC", "Low-pain PD vs HC", "High-pain PD vs HC")
names(COMP_LABELS) <- COMPARISONS

# ----------------------
# helpers
# ----------------------
read_and_annotate <- function(path) {
  if (!file.exists(path)) stop("Input file not found: ", path)
  d <- read_csv(path, show_col_types = FALSE)

  req <- c("gene_id", "comparison", "log2FoldChange", "pvalue", "padj")
  if (!all(req %in% names(d))) {
    stop("Unexpected columns in result file. Found:\n", paste(names(d), collapse = ", "))
  }

  d$Ensembl <- sub("\\..*$", "", as.character(d$gene_id))
  d$comparison <- gsub("Control", "HC", d$comparison, fixed = TRUE)
  d$comparison <- gsub("No Pain PD", "No pain PD", d$comparison, fixed = TRUE)
  d$comparison <- gsub("Low Pain PD", "Low pain PD", d$comparison, fixed = TRUE)
  d$comparison <- gsub("High Pain PD", "High pain PD", d$comparison, fixed = TRUE)

  sym_map <- AnnotationDbi::mapIds(
    org.Hs.eg.db::org.Hs.eg.db,
    keys = unique(d$Ensembl),
    keytype = "ENSEMBL",
    column = "SYMBOL",
    multiVals = "first"
  )
  d$Symbol <- as.character(unname(sym_map[d$Ensembl]))
  d$neglog10_padj <- -log10(pmax(d$padj, .Machine$double.xmin))

  d$Status <- ifelse(is.na(d$padj) | is.na(d$log2FoldChange), "No change",
                     ifelse(d$padj < 0.05 & d$log2FoldChange > 0, "Upregulated",
                            ifelse(d$padj < 0.05 & d$log2FoldChange < 0, "Downregulated", "No change")))
  d$Status <- factor(d$Status, levels = c("Downregulated", "No change", "Upregulated"))
  d
}

pick_label_rows <- function(df_sub) {
  target <- df_sub[df_sub$Symbol %in% GENES18, , drop = FALSE]
  if (!nrow(target)) return(target)

  sig <- target[!is.na(target$padj) & target$padj < 0.05, , drop = FALSE]
  keep_syms <- unique(c("TRPV1", sig$Symbol))
  out <- target[target$Symbol %in% keep_syms, , drop = FALSE]

  if (nrow(out)) {
    sp <- split(out, out$Symbol)
    out <- do.call(rbind, lapply(sp, function(z) {
      ix <- if (all(is.na(z$padj))) 1L else which.min(replace(z$padj, is.na(z$padj), Inf))
      z[ix, , drop = FALSE]
    }))
    rownames(out) <- NULL
  }
  out
}

make_single_volcano <- function(d, comp, show_legend = TRUE) {
  sub <- d[d$comparison == comp, , drop = FALSE]
  if (!nrow(sub)) stop("Comparison not found: ", comp)

  label_df <- pick_label_rows(sub)
  # Avoid selecting labels outside the visible axes; these are still in the data.
  label_df <- label_df[!is.na(label_df$log2FoldChange) & !is.na(label_df$neglog10_padj) &
                         label_df$log2FoldChange >= X_MIN & label_df$log2FoldChange <= X_MAX &
                         label_df$neglog10_padj >= Y_MIN & label_df$neglog10_padj <= Y_MAX, , drop = FALSE]

  # Use the same status color for target markers' interior, with black outline.
  # The reference figure has a light grid and leader lines to labeled targets.
  p <- ggplot(sub, aes(x = log2FoldChange, y = neglog10_padj, color = Status)) +
    geom_point(size = 1.15, alpha = 0.78, na.rm = TRUE) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray45", linewidth = 0.5) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray45", linewidth = 0.5) +
    geom_point(
      data = label_df,
      aes(x = log2FoldChange, y = neglog10_padj, fill = Status),
      inherit.aes = FALSE,
      shape = 21, color = "black", stroke = 0.85, size = 3.6
    ) +
    ggrepel::geom_label_repel(
      data = label_df,
      aes(x = log2FoldChange, y = neglog10_padj, label = Symbol),
      inherit.aes = FALSE,
      fill = "white", color = "black",
      size = 3.2, fontface = "plain",
      box.padding = 0.42, point.padding = 0.33,
      label.size = 0.25,
      min.segment.length = 0,
      segment.color = "gray35", segment.size = 0.45,
      force = 2.5, force_pull = 0.1,
      seed = 19,
      max.overlaps = Inf,
      max.time = 3
    ) +
    scale_color_manual(
      values = c("Downregulated" = COL_DOWN, "No change" = COL_NS, "Upregulated" = COL_UP),
      breaks = c("Downregulated", "No change", "Upregulated"), name = NULL
    ) +
    scale_fill_manual(
      values = c("Downregulated" = COL_DOWN, "No change" = COL_NS, "Upregulated" = COL_UP),
      guide = "none"
    ) +
    coord_cartesian(xlim = c(X_MIN, X_MAX), ylim = c(Y_MIN, Y_MAX), expand = FALSE) +
    scale_x_continuous(breaks = seq(-5, 5, 2.5)) +
    scale_y_continuous(breaks = seq(0, 10, 2)) +
    labs(
      title = COMP_LABELS[[comp]],
      x = "log2 Fold Change",
      y = expression(-log[10]~"FDR-adjusted p-value")
    ) +
    theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      axis.title = element_text(face = "bold"),
      axis.text = element_text(color = "black"),
      panel.border = element_blank(),
      panel.grid.major = element_line(color = "grey85", linewidth = 0.32),
      panel.grid.minor = element_blank(),
      axis.line = element_line(color = "black", linewidth = 0.55),
      # Wide plotting area, as in the user's reference image (no square aspect ratio).
      aspect.ratio = 0.72,
      legend.position = if (show_legend) "right" else "none",
      legend.text = element_text(size = 9),
      legend.key = element_blank(),
      plot.margin = margin(8, 12, 8, 8)
    )
  p
}

save_axis_audit <- function(d) {
  out <- do.call(rbind, lapply(COMPARISONS, function(comp) {
    sub <- d[d$comparison == comp, , drop = FALSE]
    data.frame(
      comparison = comp,
      n_total = nrow(sub),
      n_x_below = sum(sub$log2FoldChange < X_MIN, na.rm = TRUE),
      n_x_above = sum(sub$log2FoldChange > X_MAX, na.rm = TRUE),
      n_y_above = sum(sub$neglog10_padj > Y_MAX, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }))
  write_csv(out, file.path(out_dir, "Volcano_axis_audit.csv"))

  out2 <- d[d$log2FoldChange < X_MIN | d$log2FoldChange > X_MAX | d$neglog10_padj > Y_MAX,
            c("gene_id", "Ensembl", "Symbol", "comparison", "log2FoldChange", "pvalue", "padj", "neglog10_padj")]
  write_csv(out2, file.path(out_dir, "Volcano_points_outside_display_limits.csv"))
}

write_18gene_tables <- function(d) {
  tab <- d[d$comparison %in% COMPARISONS & d$Symbol %in% GENES18,
           c("comparison", "Symbol", "log2FoldChange", "pvalue", "padj")]
  tab <- tab[order(match(tab$comparison, COMPARISONS), match(tab$Symbol, GENES18)), , drop = FALSE]
  tab$Significance <- ifelse(is.na(tab$padj), "",
                             ifelse(tab$padj < 0.001, "***",
                                    ifelse(tab$padj < 0.01, "**",
                                           ifelse(tab$padj < 0.05, "*", ""))))
  write_csv(tab, file.path(out_dir, "Supplementary_Table_18genes_ageSex_adjusted_long.csv"))

  # wide format
  wide_lfc <- reshape(tab[, c("Symbol", "comparison", "log2FoldChange")],
                      idvar = "Symbol", timevar = "comparison", direction = "wide")
  wide_padj <- reshape(tab[, c("Symbol", "comparison", "padj")],
                       idvar = "Symbol", timevar = "comparison", direction = "wide")
  wide_sig <- reshape(tab[, c("Symbol", "comparison", "Significance")],
                      idvar = "Symbol", timevar = "comparison", direction = "wide")

  names(wide_lfc) <- sub("log2FoldChange\\.", "log2FC_", names(wide_lfc))
  names(wide_padj) <- sub("padj\\.", "FDR_", names(wide_padj))
  names(wide_sig) <- sub("Significance\\.", "sig_", names(wide_sig))

  wide <- Reduce(function(x, y) merge(x, y, by = "Symbol", all = TRUE), list(wide_lfc, wide_padj, wide_sig))
  hp_col <- "log2FC_High pain PD vs HC"
  if (hp_col %in% names(wide)) {
    ord <- order(-abs(wide[[hp_col]]), wide$Symbol)
    wide <- wide[ord, , drop = FALSE]
  }
  write_csv(wide, file.path(out_dir, "Supplementary_Table_18genes_ageSex_adjusted_wide.csv"))
}

# ----------------------
# run
# ----------------------
res <- read_and_annotate(input_file)
save_axis_audit(res)
write_18gene_tables(res)

p_no   <- make_single_volcano(res, "No pain PD vs HC", show_legend = TRUE)
p_low  <- make_single_volcano(res, "Low pain PD vs HC", show_legend = TRUE)
p_high <- make_single_volcano(res, "High pain PD vs HC", show_legend = TRUE)

# single square panels
# use square device dimensions so the panel looks square overall

ggsave(file.path(out_dir, "Supplementary_Volcano_NoPain_vs_HC_wide_grid.png"), p_no,
       width = 7.6, height = 5.2, dpi = 400, bg = "white")
ggsave(file.path(out_dir, "Supplementary_Volcano_NoPain_vs_HC_wide_grid.pdf"), p_no,
       width = 7.6, height = 5.2, bg = "white")

ggsave(file.path(out_dir, "Supplementary_Volcano_LowPain_vs_HC_wide_grid.png"), p_low,
       width = 7.6, height = 5.2, dpi = 400, bg = "white")
ggsave(file.path(out_dir, "Supplementary_Volcano_LowPain_vs_HC_wide_grid.pdf"), p_low,
       width = 7.6, height = 5.2, bg = "white")

ggsave(file.path(out_dir, "Supplementary_Volcano_HighPain_vs_HC_wide_grid.png"), p_high,
       width = 7.6, height = 5.2, dpi = 400, bg = "white")
ggsave(file.path(out_dir, "Supplementary_Volcano_HighPain_vs_HC_wide_grid.pdf"), p_high,
       width = 7.6, height = 5.2, bg = "white")

# combined panel: remove legends from each panel and keep panels wide
p_no_noleg   <- make_single_volcano(res, "No pain PD vs HC", show_legend = FALSE)
p_low_noleg  <- make_single_volcano(res, "Low pain PD vs HC", show_legend = FALSE)
p_high_noleg <- make_single_volcano(res, "High pain PD vs HC", show_legend = FALSE)

png(file.path(out_dir, "Supplementary_Volcano_Combined_3groups_wide_grid.png"), width = 2400, height = 650, res = 180)
grid::grid.newpage()
print(p_no_noleg,   vp = viewport(x = 0.17, y = 0.53, width = 0.30, height = 0.88))
print(p_low_noleg,  vp = viewport(x = 0.50, y = 0.53, width = 0.30, height = 0.88))
print(p_high_noleg, vp = viewport(x = 0.83, y = 0.53, width = 0.30, height = 0.88))
dev.off()

pdf(file.path(out_dir, "Supplementary_Volcano_Combined_3groups_wide_grid.pdf"), width = 18, height = 5.0)
grid::grid.newpage()
print(p_no_noleg,   vp = viewport(x = 0.17, y = 0.53, width = 0.30, height = 0.88))
print(p_low_noleg,  vp = viewport(x = 0.50, y = 0.53, width = 0.30, height = 0.88))
print(p_high_noleg, vp = viewport(x = 0.83, y = 0.53, width = 0.30, height = 0.88))
dev.off()

# top genes for QC
all_top <- do.call(rbind, lapply(COMPARISONS, function(comp) {
  sub <- res[res$comparison == comp, , drop = FALSE]
  sub <- sub[order(sub$padj, -abs(sub$log2FoldChange)), c("comparison", "Symbol", "log2FoldChange", "pvalue", "padj", "neglog10_padj")]
  head(sub, 20)
}))
write_csv(all_top, file.path(out_dir, "Supplementary_Volcano_top20_by_comparison.csv"))

writeLines(c(
  "Display settings used:",
  paste0("x-axis = [", X_MIN, ", ", X_MAX, "]"),
  paste0("y-axis = [", Y_MIN, ", ", Y_MAX, "]"),
  "Panels are wide (aspect.ratio = 0.72), grid lines are light gray.",
  "Target marker fill matches its Up/Down/NS category; label leader lines are gray."
), file.path(out_dir, "README_square_layout.txt"))

message("Done. Output folder: ", out_dir)

# =========================================================
# corr_matrices_ncomms.R
#
# Between-measure Spearman correlations of ACAM-J-induced change
# (ACAM-J minus combined control, one value per participant),
# Nature Communications figures.
#
#   Main figure:          mean of ACAM-J1-8 minus control (one matrix)
#   Supplementary figure: each ACAM-J state minus control (8 matrices, a-h)
#
# FDR (Benjamini-Hochberg) is applied across the 78 unique pairs within each
# matrix. Asterisks: * q < 0.05, ** q < 0.01, *** q < 0.001.
#
# Figures: sans-serif (Arial on Windows), 5-7 pt text, 8 pt bold lowercase
# panel letters, vector PDF + 600-dpi PNG and TIFF.
# =========================================================

pkgs <- c("dplyr", "tidyr", "ggplot2", "Hmisc", "purrr", "readr", "stringr", "patchwork", "tibble")
to_install <- pkgs[!pkgs %in% rownames(installed.packages())]
if (length(to_install) > 0) install.packages(to_install)
invisible(lapply(pkgs, library, character.only = TRUE))

# =========================================================
# PATHS
# =========================================================
df_path <- "C:/Users/peter/OneDrive/PhD/Farmaceutisk Videnskab/EEG thesis/data/entropy_master_df/entropy_master_df_with_phenomenology_auto.csv"
out_dir <- "C:/Users/peter/OneDrive/PhD/Farmaceutisk Videnskab/EEG thesis/data/corr_matrices_5s_AllACAM_ncomms"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
tab_dir <- file.path(out_dir, "tables")
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)

# =========================================================
# SETTINGS
# =========================================================
ref_state   <- "cm"                    # combined control (mean of memory and counting)
acam_states <- paste0("ACAM", 1:8)

metric_order <- c(
  "total_delta_5s", "total_theta_5s", "total_alpha_5s", "total_low_beta_5s", "total_high_beta_5s",
  "aperiodic_beta_3_35_5s", "aperiodic_offset_3_35_5s",
  "LZsum_c0", "LZsumPN", "PE_1_35_tau4_s2", "SpecEn_1to35Hz", "HFD_kmax200", "DFA_alpha"
)

# Labels match the terms used in the text and the other figures
metric_labels <- c(
  "Total delta", "Total theta", "Total alpha", "Total low beta", "Total high beta",
  "Aperiodic exponent", "Aperiodic offset",
  "LZsum", "LZsumN", "PeEn", "SpecEn", "HFD", "DFA \u03b1"
)
names(metric_labels) <- metric_order

base_pt  <- 7
pt       <- ggplot2::.pt
fill_neg <- "#3B4CC0"
fill_pos <- "#B40426"

# =========================================================
# LOAD DATA
# =========================================================
df <- read_csv(df_path, show_col_types = FALSE)

missing_cols <- setdiff(c("subject", "state", metric_order), names(df))
if (length(missing_cols) > 0) stop("Missing columns: ", paste(missing_cols, collapse = ", "))
if (nrow(filter(count(df, subject, state), n > 1)) > 0) stop("Duplicate subject-state rows in the data file.")

df_sub <- select(df, subject, state, all_of(metric_order))

# =========================================================
# DELTA TABLES
# =========================================================
ref_long <- df_sub %>%
  filter(state == ref_state) %>%
  pivot_longer(all_of(metric_order), names_to = "metric", values_to = "ref_value")

# One ACAM-J state minus control
make_delta_wide <- function(acam_state) {
  df_sub %>%
    filter(state == acam_state) %>%
    pivot_longer(all_of(metric_order), names_to = "metric", values_to = "value") %>%
    inner_join(ref_long, by = c("subject", "metric")) %>%
    mutate(delta = value - ref_value) %>%
    select(subject, metric, delta) %>%
    pivot_wider(names_from = metric, values_from = delta) %>%
    select(subject, all_of(metric_order))
}

# Mean across ACAM-J1-8 (states available for that participant) minus control
make_delta_wide_mean <- function() {
  df_sub %>%
    filter(state %in% acam_states) %>%
    pivot_longer(all_of(metric_order), names_to = "metric", values_to = "value") %>%
    group_by(subject, metric) %>%
    summarise(acam_mean = mean(value, na.rm = TRUE), .groups = "drop") %>%
    inner_join(ref_long, by = c("subject", "metric")) %>%
    mutate(delta = acam_mean - ref_value) %>%
    select(subject, metric, delta) %>%
    pivot_wider(names_from = metric, values_from = delta) %>%
    select(subject, all_of(metric_order))
}

# =========================================================
# SPEARMAN + FDR
# =========================================================
corr_with_fdr <- function(dat_wide) {
  rc <- Hmisc::rcorr(as.matrix(dat_wide[, metric_order]), type = "spearman")
  lt <- lower.tri(rc$P)
  q  <- matrix(NA_real_, nrow(rc$P), ncol(rc$P), dimnames = dimnames(rc$P))
  q[lt] <- p.adjust(rc$P[lt], method = "BH")
  q[upper.tri(q)] <- t(q)[upper.tri(q)]
  list(rho = rc$r, p = rc$P, q = q, n = rc$n)
}

save_corr_tables <- function(res, delta_wide, tag) {
  write_csv(delta_wide, file.path(tab_dir, paste0("delta_table_", tag, "_minus_control.csv")))
  for (nm in c("rho", "p", "q", "n")) {
    write_csv(rownames_to_column(as.data.frame(res[[nm]]), "metric"),
              file.path(tab_dir, paste0("corr_", nm, "_", tag, ".csv")))
  }
}

make_plot_df <- function(res) {
  nvar <- length(metric_order)
  expand.grid(i = seq_len(nvar), j = seq_len(nvar)) %>%
    as_tibble() %>%
    filter(i > j) %>%
    mutate(
      # lower triangle only: the first measure has no row and the last no column
      label_y = factor(metric_labels[metric_order[i]], levels = rev(metric_labels[-1])),
      label_x = factor(metric_labels[metric_order[j]], levels = metric_labels[-length(metric_labels)]),
      rho     = map2_dbl(i, j, ~ res$rho[.x, .y]),
      q       = map2_dbl(i, j, ~ res$q[.x, .y]),
      stars   = case_when(q < 0.001 ~ "***", q < 0.01 ~ "**", q < 0.05 ~ "*", TRUE ~ ""),
      txt     = sub("^-0\\.00$", "0.00", sprintf("%.2f", rho)),
      txt_col = ifelse(abs(rho) > 0.6, "white", "black")
    )
}

# =========================================================
# PLOT ONE MATRIX
# =========================================================
plot_corr_matrix <- function(plot_df, title_text = NULL, value_pt = 5, show_values = TRUE,
                             show_x = TRUE, show_y = TRUE) {
  p <- ggplot(plot_df, aes(x = label_x, y = label_y, fill = rho)) +
    geom_tile(colour = "white", linewidth = 0.4)
  if (show_values) {
    p <- p +
      geom_text(aes(label = txt, colour = txt_col), size = value_pt / pt,
                nudge_y = ifelse(plot_df$stars == "", 0, 0.14)) +
      geom_text(aes(label = stars, colour = txt_col), size = (value_pt + 0.5) / pt,
                nudge_y = -0.26)
  } else {
    # values are in the CSV tables; only the FDR asterisks are drawn
    p <- p + geom_text(aes(label = stars, colour = txt_col), size = (value_pt + 1) / pt,
                       vjust = 0.75)
  }
  p <- p +
    scale_colour_identity() +
    scale_fill_gradient2(low = fill_neg, mid = "white", high = fill_pos, midpoint = 0,
                         limits = c(-1, 1), breaks = c(-1, -0.5, 0, 0.5, 1),
                         name = "Spearman's \u03c1",
                         guide = guide_colourbar(barheight = unit(30, "mm"), barwidth = unit(2.5, "mm"))) +
    scale_x_discrete(drop = FALSE) +
    scale_y_discrete(drop = FALSE) +
    labs(title = title_text, x = NULL, y = NULL) +
    theme_minimal(base_size = base_pt, base_family = "sans") +
    theme(
      panel.grid   = element_blank(),
      axis.text.x  = element_text(angle = 90, vjust = 0.5, hjust = 1, size = base_pt - 1, colour = "black"),
      axis.text.y  = element_text(size = base_pt - 1, colour = "black"),
      plot.title   = element_text(size = base_pt, hjust = 0.5, face = "plain"),
      plot.tag     = element_text(size = base_pt + 1, face = "bold"),
      legend.title = element_text(size = base_pt - 1),
      legend.text  = element_text(size = base_pt - 1),
      plot.background  = element_rect(fill = "white", colour = NA),
      panel.background = element_rect(fill = "white", colour = NA),
      plot.margin = margin(1, 1, 1, 1, "mm")
    )
  if (!show_x) p <- p + theme(axis.text.x = element_blank())
  if (!show_y) p <- p + theme(axis.text.y = element_blank())
  p
}

save_fig <- function(plot, stem, width_mm, height_mm) {
  pdf_device <- if (capabilities("cairo")) cairo_pdf else "pdf"
  ggsave(file.path(out_dir, paste0(stem, ".pdf")), plot, device = pdf_device,
         width = width_mm, height = height_mm, units = "mm", bg = "white")
  ggsave(file.path(out_dir, paste0(stem, ".png")), plot,
         width = width_mm, height = height_mm, units = "mm", dpi = 600, bg = "white")
  ggsave(file.path(out_dir, paste0(stem, ".tiff")), plot,
         width = width_mm, height = height_mm, units = "mm", dpi = 600, bg = "white",
         compression = "lzw")
  message("Saved ", stem, " (", width_mm, " x ", height_mm, " mm)")
}

# =========================================================
# MAIN FIGURE: mean ACAM-J minus control
# =========================================================
delta_mean <- make_delta_wide_mean()
res_mean   <- corr_with_fdr(delta_mean)
save_corr_tables(res_mean, delta_mean, "mean_ACAM")

p_main <- plot_corr_matrix(make_plot_df(res_mean), value_pt = 5.5) +
  coord_fixed() +
  theme(legend.position = c(0.97, 0.97), legend.justification = c(1, 1))
save_fig(p_main, "FIG_corr_matrix_mean_ACAM_minus_control", width_mm = 120, height_mm = 108)

# =========================================================
# SUPPLEMENTARY FIGURE: each ACAM-J state minus control (a-h)
# =========================================================
n_text <- function(n) {
  lo <- min(n[n > 0]); hi <- max(n)
  if (lo == hi) as.character(lo) else paste0(lo, "\u2013", hi)
}

plot_list <- list()
for (k in seq_along(acam_states)) {
  acam <- acam_states[k]
  message("Processing ", acam)
  delta_k <- make_delta_wide(acam)
  res_k   <- corr_with_fdr(delta_k)
  save_corr_tables(res_k, delta_k, acam)
  col <- (k - 1) %% 2 + 1; row <- (k - 1) %/% 2 + 1
  plot_list[[acam]] <- plot_corr_matrix(
    make_plot_df(res_k),
    title_text = paste0("ACAM-J", k, " \u2212 control (n = ", n_text(res_k$n), ")"),
    value_pt = 5,
    show_values = FALSE,    # cells are too small for numbers; values are in tables/
    show_x = row == 4,      # x labels on the bottom row only
    show_y = col == 1       # y labels on the left column only
  )
}

p_supp <- wrap_plots(plot_list, ncol = 2, byrow = TRUE, guides = "collect") +
  plot_annotation(tag_levels = "a") &
  theme(legend.position = "right")
save_fig(p_supp, "SUPP_FIG_corr_matrices_per_ACAM_state", width_mm = 183, height_mm = 247)

message("\nDone. Outputs saved to:\n", out_dir)

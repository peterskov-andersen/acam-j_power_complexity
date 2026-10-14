# corr_matrix_script.R
#
# Spearman correlations between EEG measures, computed on ACAM-J minus
# control (cm) differences, for each ACAM-J state and for the mean of all
# ACAM-J states. BH-FDR across the unique pairs of each matrix.

library(grid)
library(dplyr)
library(tidyr)
library(ggplot2)
library(Hmisc)
library(purrr)
library(readr)
library(stringr)
library(patchwork)
library(tibble)

# =========================================================
# PATHS
# =========================================================
df_path  <- "C:/Users/peter/OneDrive/PhD/Farmaceutisk Videnskab/EEG thesis/data/entropy_master_df/entropy_master_df_with_phenomenology_auto.csv"
out_dir  <- "C:/Users/peter/OneDrive/PhD/Farmaceutisk Videnskab/EEG thesis/data/corr_matrices_5s_AllACAM"

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# =========================================================
# SETTINGS
# =========================================================
ref_state   <- "cm"
acam_states <- paste0("ACAM", 1:8)

metric_order <- c(
  "total_delta_5s",
  "total_theta_5s",
  "total_alpha_5s",
  "total_low_beta_5s",
  "total_high_beta_5s",
  "aperiodic_beta_3_35_5s",
  "aperiodic_offset_3_35_5s",
  "LZsum_c0",
  "LZsumPN",
  "PE_1_35_tau4_s2",
  "SpecEn_1to35Hz",
  "HFD_kmax200",
  "DFA_alpha"
)

# Note: aperiodic_beta_3_35_5s is the aperiodic exponent (beta = -slope);
# it is labelled "Aperiodic slope" in the figures.
metric_labels <- c(
  "PSD mixed delta",
  "PSD mixed theta",
  "PSD mixed alpha",
  "PSD mixed low beta",
  "PSD mixed high beta",
  "Aperiodic slope",
  "Aperiodic offset",
  "LZsum",
  "LZsumN",
  "PeEn",
  "SpecEn",
  "HFD",
  "DFA"
)
names(metric_labels) <- metric_order

# =========================================================
# LOAD DATA
# =========================================================
df <- read_csv(df_path, show_col_types = FALSE)

missing_cols <- setdiff(c("subject", "state", metric_order), names(df))
if (length(missing_cols) > 0) {
  stop("These required columns are missing from df:\n", paste(missing_cols, collapse = ", "))
}

dup_check <- df %>%
  count(subject, state) %>%
  filter(n > 1)

if (nrow(dup_check) > 0) {
  stop("Your dataframe has duplicate subject-state rows. Fix that first.")
}

df_sub <- df %>%
  select(subject, state, all_of(metric_order))

# =========================================================
# HELPERS
# =========================================================
pretty_acam <- function(acam_state) {
  paste0("ACAM-J", str_extract(acam_state, "\\d+"))
}

# =========================================================
# FUNCTION: build delta dataframe for one ACAM
# =========================================================
make_delta_wide <- function(dat, acam_state, ref_state = "cm", metric_order) {
  dat %>%
    filter(state %in% c(ref_state, acam_state)) %>%
    pivot_longer(cols = all_of(metric_order), names_to = "metric", values_to = "value") %>%
    pivot_wider(names_from = state, values_from = value) %>%
    mutate(delta = .data[[acam_state]] - .data[[ref_state]]) %>%
    select(subject, metric, delta) %>%
    pivot_wider(names_from = metric, values_from = delta) %>%
    select(subject, all_of(metric_order))
}

# =========================================================
# FUNCTION: build delta dataframe for mean across all ACAMs
# =========================================================
make_delta_wide_all_acam <- function(dat, acam_states, ref_state = "cm", metric_order) {
  
  acam_mean_long <- dat %>%
    filter(state %in% acam_states) %>%
    pivot_longer(cols = all_of(metric_order), names_to = "metric", values_to = "value") %>%
    group_by(subject, metric) %>%
    summarise(acam_mean = mean(value, na.rm = TRUE), .groups = "drop")
  
  ref_long <- dat %>%
    filter(state == ref_state) %>%
    pivot_longer(cols = all_of(metric_order), names_to = "metric", values_to = "ref_value")
  
  acam_mean_long %>%
    inner_join(ref_long, by = c("subject", "metric")) %>%
    mutate(delta = acam_mean - ref_value) %>%
    select(subject, metric, delta) %>%
    pivot_wider(names_from = metric, values_from = delta) %>%
    select(subject, all_of(metric_order))
}

# =========================================================
# FUNCTION: Spearman correlation + BH-FDR (lower triangle)
# =========================================================
corr_with_fdr <- function(dat_wide, metric_order) {
  x <- as.matrix(dat_wide[, metric_order])
  
  rc <- Hmisc::rcorr(x, type = "spearman")
  
  rho_mat <- rc$r
  p_mat   <- rc$P
  n_mat   <- rc$n
  
  lt_idx <- lower.tri(p_mat)
  q_vals <- p.adjust(p_mat[lt_idx], method = "fdr")
  
  q_mat <- matrix(
    NA_real_,
    nrow = nrow(p_mat),
    ncol = ncol(p_mat),
    dimnames = dimnames(p_mat)
  )
  q_mat[lt_idx] <- q_vals
  q_mat[upper.tri(q_mat)] <- t(q_mat)[upper.tri(q_mat)]
  diag(q_mat) <- NA_real_
  
  list(rho = rho_mat, p = p_mat, q = q_mat, n = n_mat)
}

# =========================================================
# FUNCTION: make plotting dataframe
# =========================================================
make_plot_df <- function(rho_mat, q_mat, metric_order, metric_labels) {
  nvar <- length(metric_order)
  
  expand.grid(i = seq_len(nvar), j = seq_len(nvar)) %>%
    as_tibble() %>%
    filter(i > j) %>%
    mutate(
      var_y   = metric_order[i],
      var_x   = metric_order[j],
      label_y = metric_labels[var_y],
      label_x = metric_labels[var_x],
      rho     = map2_dbl(i, j, ~ rho_mat[.x, .y]),
      q       = map2_dbl(i, j, ~ q_mat[.x, .y]),
      stars   = case_when(
        !is.na(q) & q < 0.01 ~ "**",
        !is.na(q) & q < 0.05 ~ "*",
        TRUE ~ ""
      ),
      text_lab = sprintf("%.2f", rho)
    )
}

# =========================================================
# FUNCTION: plot one matrix
# =========================================================
plot_corr_matrix <- function(plot_df, title_text) {
  
  ggplot(
    plot_df,
    aes(
      x = factor(label_x, levels = metric_labels[metric_order]),
      y = factor(label_y, levels = rev(metric_labels[metric_order])),
      fill = rho
    )
  ) +
    geom_tile(color = "white", linewidth = 0.7) +
    geom_text(
      aes(label = text_lab),
      size = 4.5,
      color = "black"
    ) +
    geom_text(
      aes(label = stars),
      size = 6.0,
      nudge_y = 0.23,
      color = "black",
      fontface = "bold"
    ) +
    scale_fill_gradient2(
      low = "#3B4CC0",
      mid = "white",
      high = "#B40426",
      midpoint = 0,
      limits = c(-1, 1),
      name = "Spearman rho",
      guide = guide_colorbar(
        barheight = unit(13, "cm"),
        barwidth  = unit(0.8, "cm")
      )
    ) +
    coord_fixed() +
    labs(
      title = title_text,
      x = NULL,
      y = NULL
    ) +
    theme_minimal(base_size = 12) +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 11),
      axis.text.y = element_text(size = 11),
      plot.title = element_text(face = "bold", hjust = 0.5),
      legend.justification = c(1, 1),
      legend.direction = "vertical",
      legend.background = element_rect(fill = scales::alpha("white", 0.85), color = "white"),
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA)
    )
}

# =========================================================
# RUN FOR ALL INDIVIDUAL ACAM STATES
# =========================================================
plot_list <- list()

for (acam in acam_states) {
  
  message("Processing ", acam, "...")
  
  delta_wide <- make_delta_wide(
    dat = df_sub,
    acam_state = acam,
    ref_state = ref_state,
    metric_order = metric_order
  )
  
  write_csv(
    delta_wide,
    file.path(out_dir, paste0("delta_table_", acam, "_minus_", ref_state, ".csv"))
  )
  
  corr_res <- corr_with_fdr(delta_wide, metric_order)
  
  write_csv(
    as.data.frame(corr_res$rho) %>% rownames_to_column("metric"),
    file.path(out_dir, paste0("corr_rho_", acam, ".csv"))
  )
  
  write_csv(
    as.data.frame(corr_res$p) %>% rownames_to_column("metric"),
    file.path(out_dir, paste0("corr_p_raw_", acam, ".csv"))
  )
  
  write_csv(
    as.data.frame(corr_res$q) %>% rownames_to_column("metric"),
    file.path(out_dir, paste0("corr_q_fdr_", acam, ".csv"))
  )
  
  write_csv(
    as.data.frame(corr_res$n) %>% rownames_to_column("metric"),
    file.path(out_dir, paste0("corr_n_", acam, ".csv"))
  )
  
  plot_df <- make_plot_df(
    rho_mat = corr_res$rho,
    q_mat = corr_res$q,
    metric_order = metric_order,
    metric_labels = metric_labels
  )
  
  p <- plot_corr_matrix(
    plot_df,
    title_text = paste0("\u0394 ", pretty_acam(acam), " - Control")
  )
  
  plot_list[[acam]] <- p
  
  ggsave(
    filename = file.path(out_dir, paste0("corr_matrix_spearman_", acam, "_delta_vs_", ref_state, ".png")),
    plot = p,
    width = 8.6,
    height = 7.3,
    dpi = 300
  )
  
  ggsave(
    filename = file.path(out_dir, paste0("corr_matrix_spearman_", acam, "_delta_vs_", ref_state, ".pdf")),
    plot = p,
    width = 8.6,
    height = 7.3
  )
}

# =========================================================
# RUN FOR POOLED MEAN ACAM-J VS CONTROL
# =========================================================
message("Processing pooled All ACAM-J mean vs control...")

delta_wide_all_acam <- make_delta_wide_all_acam(
  dat = df_sub,
  acam_states = acam_states,
  ref_state = ref_state,
  metric_order = metric_order
)

write_csv(
  delta_wide_all_acam,
  file.path(out_dir, paste0("delta_table_ALL_ACAM_mean_minus_", ref_state, ".csv"))
)

corr_res_all_acam <- corr_with_fdr(delta_wide_all_acam, metric_order)

write_csv(
  as.data.frame(corr_res_all_acam$rho) %>% rownames_to_column("metric"),
  file.path(out_dir, "corr_rho_ALL_ACAM_mean.csv")
)

write_csv(
  as.data.frame(corr_res_all_acam$p) %>% rownames_to_column("metric"),
  file.path(out_dir, "corr_p_raw_ALL_ACAM_mean.csv")
)

write_csv(
  as.data.frame(corr_res_all_acam$q) %>% rownames_to_column("metric"),
  file.path(out_dir, "corr_q_fdr_ALL_ACAM_mean.csv")
)

write_csv(
  as.data.frame(corr_res_all_acam$n) %>% rownames_to_column("metric"),
  file.path(out_dir, "corr_n_ALL_ACAM_mean.csv")
)

plot_df_all_acam <- make_plot_df(
  rho_mat = corr_res_all_acam$rho,
  q_mat = corr_res_all_acam$q,
  metric_order = metric_order,
  metric_labels = metric_labels
)

p_all_acam <- plot_corr_matrix(
  plot_df_all_acam,
  title_text = "\u0394 Mean ACAM-J - Control"
)

ggsave(
  filename = file.path(out_dir, "corr_matrix_spearman_ALL_ACAM_mean_delta_vs_cm.png"),
  plot = p_all_acam,
  width = 8.6,
  height = 7.3,
  dpi = 300
)

ggsave(
  filename = file.path(out_dir, "corr_matrix_spearman_ALL_ACAM_mean_delta_vs_cm.pdf"),
  plot = p_all_acam,
  width = 8.6,
  height = 7.3
)

# =========================================================
# COMBINED 2 x 4 FIGURE (INDIVIDUAL ACAMs ONLY)
# =========================================================
combined_plot <- wrap_plots(
  plot_list,
  ncol = 4,
  nrow = 2,
  byrow = TRUE,
  guides = "collect"
) &
  theme(
    legend.position = c(0.955, 0.86),
    legend.justification = c(1, 1),
    legend.background = element_rect(fill = scales::alpha("white", 0.85), color = "grey80")
  )

ggsave(
  filename = file.path(out_dir, "corr_matrices_all_ACAM_spearman_delta_vs_cm.png"),
  plot = combined_plot,
  width = 26,
  height = 11,
  dpi = 300
)

ggsave(
  filename = file.path(out_dir, "corr_matrices_all_ACAM_spearman_delta_vs_cm.pdf"),
  plot = combined_plot,
  width = 24,
  height = 12
)

# =========================================================
# Print two of the figures
# =========================================================
print(plot_list[["ACAM1"]])
print(p_all_acam)
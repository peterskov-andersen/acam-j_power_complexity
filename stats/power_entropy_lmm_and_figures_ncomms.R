# ============================================================
# entropy_lmm_and_figures_ncomms.R
#
# Linear mixed models, supplementary statistics tables and
# Nature Communications figures for the global EEG measures
# (complexity measures and IRASA power, 5-s epochs).
#
# Models (per measure), each data point entering once:
#   Main model:    measure ~ state + (1 | subject)
#                  state = combined control (cm = mean of memory and counting,
#                  one value per participant) and ACAM-J1..8
#                  -> ACAM-J1..8 vs combined control, BH-FDR across the 8 contrasts
#                  -> omnibus F test of state
#   Control model: measure ~ state + (1 | subject), state = memory (m), counting (c)
#                  -> counting - memory, uncorrected (single contrast)
#   Memory and counting rows are not in the main model, and cm is not in the
#   control model, so no participant's control data are counted twice.
#   Satterthwaite df (lmerTest / emmeans). Cohen's d = estimate / residual SD
#   of the model the contrast comes from.
#
# Outputs (in out_dir):
#   per_measure/LMM_<metric>_*.csv, DESC_<metric>_by_state.csv, LMM_<metric>_model.rds
#   SUPP_TABLE_*.csv and Supplementary_Data_1_LMM_statistics.xlsx
#   FIG_complexity_panel.pdf/.png/.tiff          (main-text figure, panels a-f)
#   FIG_power_total_panel.pdf/.png/.tiff         (supplementary)
#   FIG_power_oscillatory_panel.pdf/.png/.tiff   (supplementary)
#   FIG_power_aperiodic_panel.pdf/.png/.tiff     (supplementary)
#
# Figures: 183 mm wide (double column), sans-serif (Arial on Windows),
# 6-7 pt text, 8 pt bold lowercase panel letters, vector PDF + 600-dpi PNG/TIFF.
# ============================================================

# ---------------------------
# 0) Packages
# ---------------------------
pkgs <- c("readr", "dplyr", "tidyr", "tibble", "lme4", "lmerTest", "emmeans",
          "ggplot2", "ggdist", "patchwork", "writexl")
to_install <- pkgs[!pkgs %in% rownames(installed.packages())]
if (length(to_install) > 0) install.packages(to_install)
invisible(lapply(pkgs, library, character.only = TRUE))

# ---------------------------
# 1) User settings
# ---------------------------
# Paths (for the public code release, change these to relative paths,
# e.g. "data/entropy_master_df_with_phenomenology_auto.csv" and "outputs")
csv_path <- "C:/Users/peter/OneDrive/PhD/Farmaceutisk Videnskab/EEG thesis/data/entropy_master_df/entropy_master_df_with_phenomenology_auto.csv"
out_dir  <- "C:/Users/peter/OneDrive/PhD/Farmaceutisk Videnskab/EEG thesis/data/entropy_master_df/entropy_outputs_ncomms"

# Measures: column name, colour, label used in tables and figures
run_specs <- tibble::tribble(
  ~metric,                     ~measure_color,   ~y_label,
  # ---- Complexity ----
  "LZsum_c0",                  "coral",          "LZsum",
  "LZsumPN",                   "coral2",         "LZsumN",
  "PE_1_35_tau4_s2",           "orangered",      "PeEn",
  "SpecEn_1to35Hz",            "turquoise4",     "SpecEn",
  "HFD_kmax200",               "seagreen",       "HFD",
  "DFA_alpha",                 "seagreen2",      "DFA alpha",
  # ---- Total (mixed) power, 5-s epochs ----
  "total_delta_5s",            "firebrick3",     "Total power delta",
  "total_theta_5s",            "tomato3",        "Total power theta",
  "total_alpha_5s",            "orangered3",     "Total power alpha",
  "total_low_beta_5s",         "darkorange3",    "Total power low beta",
  "total_high_beta_5s",        "goldenrod3",     "Total power high beta",
  # ---- Aperiodic power, 5-s epochs ----
  "aperiodic_delta_5s",        "darkorchid4",    "Aperiodic power delta",
  "aperiodic_theta_5s",        "darkorchid3",    "Aperiodic power theta",
  "aperiodic_alpha_5s",        "mediumorchid3",  "Aperiodic power alpha",
  "aperiodic_low_beta_5s",     "orchid4",        "Aperiodic power low beta",
  "aperiodic_high_beta_5s",    "plum4",          "Aperiodic power high beta",
  # ---- Oscillatory power, 5-s epochs ----
  "oscillatory_delta_5s",      "forestgreen",    "Oscillatory power delta",
  "oscillatory_theta_5s",      "seagreen3",      "Oscillatory power theta",
  "oscillatory_alpha_5s",      "springgreen4",   "Oscillatory power alpha",
  "oscillatory_low_beta_5s",   "olivedrab4",     "Oscillatory power low beta",
  "oscillatory_high_beta_5s",  "yellowgreen",    "Oscillatory power high beta"
)

# Which measures go into which multi-panel figure (in panel order a, b, c, ...)
figure_sets <- list(
  FIG_complexity_panel        = c("LZsum_c0", "LZsumPN", "PE_1_35_tau4_s2",
                                  "SpecEn_1to35Hz", "HFD_kmax200", "DFA_alpha"),
  FIG_power_total_panel       = paste0("total_",       c("delta", "theta", "alpha", "low_beta", "high_beta"), "_5s"),
  FIG_power_oscillatory_panel = paste0("oscillatory_", c("delta", "theta", "alpha", "low_beta", "high_beta"), "_5s"),
  FIG_power_aperiodic_panel   = paste0("aperiodic_",   c("delta", "theta", "alpha", "low_beta", "high_beta"), "_5s")
)

# Figure settings (Nature Communications: 183 mm double column, 5-8 pt text)
fig_width_mm  <- 183
fig_row_mm    <- 55
fig_ncol      <- 2
base_pt       <- 7
control_color <- "grey70"

# State coding in the data file
state_levels_main  <- c("cm", paste0("ACAM", 1:8))      # rows used in the main model
state_levels_ctrl  <- c("m", "c")                         # rows used in the control model
state_levels_plot  <- c("cm", paste0("ACAM", 1:8))       # rows shown in the figures
acam_levels        <- paste0("ACAM", 1:8)
x_labels_fig       <- c("cm" = "C", setNames(paste0("J", 1:8), acam_levels))  # legend: C, combined control; J1-J8, ACAM-J1-8

# ---------------------------
# 2) Helpers
# ---------------------------
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
per_dir <- file.path(out_dir, "per_measure")
dir.create(per_dir, recursive = TRUE, showWarnings = FALSE)

safe_name <- function(x) gsub("[^A-Za-z0-9_]+", "_", x)

# Exact P values; anything below 0.0001 written as "<0.0001"
fmt_p <- function(p) {
  ifelse(is.na(p), NA_character_,
         ifelse(p < 1e-4, "<0.0001", formatC(p, format = "f", digits = 4)))
}

pick_ci <- function(d) {
  if (all(c("lower.CL", "upper.CL") %in% names(d))) {
    d <- rename(d, d_lower = lower.CL, d_upper = upper.CL)
  } else if (all(c("asymp.LCL", "asymp.UCL") %in% names(d))) {
    d <- rename(d, d_lower = asymp.LCL, d_upper = asymp.UCL)
  } else {
    d$d_lower <- NA_real_; d$d_upper <- NA_real_
  }
  select(d, contrast, d = effect.size, d_SE = SE, d_df = df, d_lower, d_upper)
}

std_ci_names <- function(x) {
  if (!all(c("lower.CL", "upper.CL") %in% names(x)) &&
      all(c("asymp.LCL", "asymp.UCL") %in% names(x))) {
    x <- rename(x, lower.CL = asymp.LCL, upper.CL = asymp.UCL)
  }
  x
}

make_desc_table <- function(df_in, metric) {
  df_in %>%
    filter(!is.na(.data[[metric]])) %>%
    group_by(state) %>%
    summarise(n_subjects = n_distinct(subject), n_rows = n(),
              mean = mean(.data[[metric]]), sd = sd(.data[[metric]]),
              median = median(.data[[metric]]),
              min = min(.data[[metric]]), max = max(.data[[metric]]),
              .groups = "drop")
}

metric_annotation <- function(metric) {
  domain <- NA_character_; measure_type <- NA_character_
  band <- NA_character_;   epoch_length <- NA_character_
  if (grepl("_5s$", metric)) epoch_length <- "5s"
  if (grepl("^(total|aperiodic|oscillatory)_", metric)) {
    domain <- "power"
    measure_type <- paste0(sub("_.*$", "", metric), "_power")
    band <- sub("_5s$", "", sub("^(total|aperiodic|oscillatory)_", "", metric))
  } else if (grepl("^PE_", metric)) {
    domain <- "complexity"; measure_type <- "permutation_entropy"; band <- "1-35Hz"; epoch_length <- "2s"
  } else if (grepl("^SpecEn_", metric)) {
    domain <- "complexity"; measure_type <- "spectral_entropy"; band <- "1-35Hz"; epoch_length <- "2s"
  } else if (grepl("^LZ", metric)) {
    domain <- "complexity"; measure_type <- "lempel_ziv"; epoch_length <- "2s"
  } else if (grepl("^HFD", metric)) {
    domain <- "complexity"; measure_type <- "higuchi_fd"; epoch_length <- "2s"
  } else if (grepl("^DFA", metric)) {
    domain <- "complexity"; measure_type <- "dfa"; epoch_length <- "10s"
  }
  tibble(metric = metric, domain = domain, measure_type = measure_type,
         band = band, epoch_length = epoch_length)
}

# ---------------------------
# 3) Model for one measure
# ---------------------------
run_lmm <- function(metric, df_in) {
  if (!metric %in% names(df_in)) stop("Measure '", metric, "' not found in the data file.")
  message("Running LMM for: ", metric)

  fit_lmm <- function(levels_used, ref) {
    d <- df_in %>%
      filter(!is.na(.data[[metric]]), !is.na(subject), state %in% levels_used) %>%
      mutate(subject = factor(subject),
             state   = droplevels(factor(state, levels = levels_used)),
             state   = relevel(state, ref = ref))
    fit <- lmer(as.formula(paste0(metric, " ~ state + (1 | subject)")), data = d,
                REML = TRUE, control = lmerControl(optimizer = "bobyqa"))
    list(data = d, fit = fit, emm = emmeans(fit, ~ state, lmer.df = "satterthwaite"))
  }

  # ---- Main model: combined control + ACAM-J1..8 ----
  main <- fit_lmm(state_levels_main, ref = "cm")
  m_fit <- main$fit
  null_fit <- update(m_fit, . ~ 1 + (1 | subject))
  anova_fixed <- anova(m_fit)              # Satterthwaite F test of state
  lrt_vs_null <- anova(null_fit, m_fit)    # likelihood-ratio test (refitted with ML)
  emm <- main$emm

  res_acam <- as.data.frame(summary(contrast(emm, method = "trt.vs.ctrl", ref = "cm", adjust = "none"),
                                    infer = TRUE)) %>%
    std_ci_names() %>%
    filter(grepl("^ACAM[1-8]\\s*-\\s*cm$", contrast)) %>%
    mutate(p_fdr_bh = p.adjust(p.value, method = "BH"))
  es_acam <- as.data.frame(eff_size(emm, sigma = sigma(m_fit), edf = df.residual(m_fit),
                                    method = "trt.vs.ctrl", ref = "cm")) %>%
    mutate(contrast = gsub("[()]", "", contrast))
  res_acam <- left_join(res_acam, pick_ci(es_acam), by = "contrast")

  # ---- Control model: counting vs memory ----
  ctrl <- fit_lmm(state_levels_ctrl, ref = "m")
  w_cm <- list("c - m" = c(m = -1, c = 1)[levels(ctrl$data$state)])
  res_cm <- as.data.frame(summary(contrast(ctrl$emm, method = w_cm, adjust = "none"), infer = TRUE)) %>%
    std_ci_names()
  es_cm <- as.data.frame(eff_size(ctrl$emm, sigma = sigma(ctrl$fit), edf = df.residual(ctrl$fit),
                                  method = w_cm))
  res_cm <- left_join(res_cm, pick_ci(es_cm), by = "contrast")

  # Per-measure files
  ms <- safe_name(metric)
  write_csv(res_acam, file.path(per_dir, paste0("LMM_", ms, "_ACAM_vs_control_FDR_BH_with_d.csv")))
  write_csv(res_cm,   file.path(per_dir, paste0("LMM_", ms, "_counting_minus_memory_with_d.csv")))
  write_csv(as.data.frame(emm), file.path(per_dir, paste0("LMM_", ms, "_EMMeans_by_state.csv")))
  write_csv(rownames_to_column(as.data.frame(anova_fixed), "term"), file.path(per_dir, paste0("LMM_", ms, "_anova_fixed_effect.csv")))
  write_csv(rownames_to_column(as.data.frame(lrt_vs_null), "model"), file.path(per_dir, paste0("LMM_", ms, "_LRT_vs_null.csv")))
  write_csv(make_desc_table(df_in, metric), file.path(per_dir, paste0("DESC_", ms, "_by_state.csv")))
  saveRDS(m_fit, file.path(per_dir, paste0("LMM_", ms, "_model.rds")))

  list(metric = metric, model = m_fit, emm = emm, res_acam = res_acam, res_cm = res_cm,
       anova_fixed = anova_fixed, lrt_vs_null = lrt_vs_null)
}

# ---------------------------
# 4) Load data and run all models
# ---------------------------
df <- read_csv(csv_path, show_col_types = FALSE) %>%
  mutate(subject = factor(subject), state = as.character(state))

all_results <- setNames(lapply(run_specs$metric, run_lmm, df_in = df), run_specs$metric)

# ---------------------------
# 5) Supplementary statistics tables
# ---------------------------
combine_supp_table <- function(which_table = c("acam", "cm"), metric_pattern) {
  which_table <- match.arg(which_table)
  res_name <- if (which_table == "acam") "res_acam" else "res_cm"
  metrics_to_use <- run_specs$metric[grepl(metric_pattern, run_specs$metric)]

  bind_rows(lapply(metrics_to_use, function(metric_i) {
    x <- all_results[[metric_i]][[res_name]]
    if (!"p_fdr_bh" %in% names(x)) x$p_fdr_bh <- NA_real_
    x %>%
      mutate(metric = metric_i,
             label  = run_specs$y_label[run_specs$metric == metric_i]) %>%
      left_join(metric_annotation(metric_i), by = "metric") %>%
      transmute(
        measure      = label,
        metric,
        measure_type,
        band,
        epoch_length,
        contrast     = case_when(contrast == "c - m" ~ "counting - memory",
                                 TRUE ~ sub("ACAM([1-8]) - cm", "ACAM-J\\1 - combined control", contrast)),
        estimate     = signif(estimate, 4),
        SE           = signif(SE, 4),
        lower_95CI   = signif(lower.CL, 4),
        upper_95CI   = signif(upper.CL, 4),
        df           = round(df, 1),
        t            = round(t.ratio, 2),
        p            = fmt_p(p.value),
        p_FDR        = fmt_p(p_fdr_bh),
        cohens_d     = round(d, 2),
        d_lower_95CI = round(d_lower, 2),
        d_upper_95CI = round(d_upper, 2)
      )
  }))
}

power_pattern      <- "_5s$"
complexity_pattern <- "^(LZ|PE_|SpecEn|HFD|DFA)"

tables <- list(
  SUPP_TABLE_complexity_ACAM_vs_control    = combine_supp_table("acam", complexity_pattern),
  SUPP_TABLE_power_ACAM_vs_control         = combine_supp_table("acam", power_pattern),
  SUPP_TABLE_complexity_counting_vs_memory = combine_supp_table("cm",   complexity_pattern),
  SUPP_TABLE_power_counting_vs_memory      = combine_supp_table("cm",   power_pattern)
)

tables$SUPP_TABLE_omnibus_F_all_measures <- bind_rows(lapply(run_specs$metric, function(m) {
  a <- as.data.frame(all_results[[m]]$anova_fixed)
  tibble(measure = run_specs$y_label[run_specs$metric == m], metric = m,
         F = round(a[["F value"]][1], 2), num_df = a[["NumDF"]][1],
         den_df = round(a[["DenDF"]][1], 1), p = fmt_p(a[["Pr(>F)"]][1]))
}))

for (nm in names(tables)) {
  write_csv(tables[[nm]], file.path(out_dir, paste0(nm, ".csv")), na = "")
  message("Saved ", nm, ".csv (", nrow(tables[[nm]]), " rows)")
}

readme <- tibble(
  column = c("measure", "contrast", "estimate", "SE", "lower_95CI / upper_95CI", "df", "t",
             "p", "p_FDR", "cohens_d (95% CI)", "F, num_df, den_df"),
  description = c(
    "EEG measure (global value averaged across all channels)",
    "ACAM-J state minus the combined control (mean of the counting and memory conditions), or counting minus memory",
    "Estimated difference from the linear mixed-effects model measure ~ state + (1 | subject), in the units of the measure",
    "Standard error of the estimate",
    "95% confidence interval of the estimate",
    "Satterthwaite degrees of freedom",
    "t statistic",
    "Exact two-sided P value, uncorrected",
    "Benjamini-Hochberg FDR-corrected P across the eight ACAM-J contrasts within each measure (not applied to counting - memory)",
    "Estimated difference divided by the model residual standard deviation, with 95% CI",
    "Omnibus F test of state (Satterthwaite) and its degrees of freedom"))

write_xlsx(list(
  "README"                        = readme,
  "Complexity ACAM-J vs control"  = tables$SUPP_TABLE_complexity_ACAM_vs_control,
  "Power ACAM-J vs control"       = tables$SUPP_TABLE_power_ACAM_vs_control,
  "Complexity counting vs memory" = tables$SUPP_TABLE_complexity_counting_vs_memory,
  "Power counting vs memory"      = tables$SUPP_TABLE_power_counting_vs_memory,
  "Omnibus F test of state"       = tables$SUPP_TABLE_omnibus_F_all_measures
), file.path(out_dir, "Supplementary_Data_1_LMM_statistics.xlsx"))

# ---------------------------
# 6) Figures
# ---------------------------
pt <- ggplot2::.pt   # converts pt to ggplot text size units

# Tick labels with as many decimals as the range needs
axis_labels <- function(rng) {
  digits <- max(0, ceiling(-log10(diff(rng))) + 1)
  function(x) formatC(x, format = "f", digits = digits)
}

y_axis_label <- function(metric, label) {
  if (grepl("^DFA", metric)) return(quote(DFA ~ alpha))
  sub("^(Total|Aperiodic|Oscillatory) power ", "\\1 ", label)   # shorter labels inside panels
}

make_panel_plot <- function(metric) {
  spec <- run_specs[run_specs$metric == metric, ]
  res  <- all_results[[metric]]$res_acam

  dfp <- df %>%
    filter(state %in% state_levels_plot, !is.na(.data[[metric]])) %>%
    mutate(state = factor(state, levels = state_levels_plot),
           fill_group = ifelse(state == "cm", "Control", "ACAM"))

  lab <- res %>%
    mutate(state = sub("\\s*-\\s*cm$", "", contrast),
           stars = case_when(p_fdr_bh < 0.001 ~ "***", p_fdr_bh < 0.01 ~ "**",
                             p_fdr_bh < 0.05 ~ "*", TRUE ~ ""),
           d_txt = ifelse(is.na(d), "", sub("^-0\\.00$", "0.00", sprintf("%.2f", d)))) %>%
    select(state, stars, d_txt)
  lab <- tibble(state = acam_levels) %>%
    left_join(lab, by = "state") %>%
    mutate(stars = coalesce(stars, ""), d_txt = coalesce(d_txt, ""),
           state = factor(state, levels = state_levels_plot))

  y_range <- range(dfp[[metric]], na.rm = TRUE)
  y_diff  <- diff(y_range); if (!is.finite(y_diff) || y_diff == 0) y_diff <- 1
  lab$y_d     <- y_range[2] + 0.05 * y_diff
  lab$y_stars <- y_range[2] + 0.14 * y_diff
  y_lim_top   <- y_range[2] + 0.26 * y_diff

  mean_df <- dfp %>% group_by(state) %>%
    summarise(mean_val = mean(.data[[metric]]), .groups = "drop")

  ggplot(dfp, aes(x = state, y = .data[[metric]], fill = fill_group)) +
    stat_halfeye(adjust = 0.8, width = 0.7, justification = -0.15,
                 point_colour = NA, interval_colour = NA, alpha = 0.7) +
    geom_boxplot(inherit.aes = FALSE, aes(x = state, y = .data[[metric]]),
                 width = 0.22, outlier.shape = NA, fill = "white", colour = "black", linewidth = 0.25) +
    geom_point(position = position_jitter(width = 0.08, height = 0, seed = 1),
               alpha = 0.6, size = 0.5, stroke = 0) +
    geom_line(data = mean_df, aes(x = state, y = mean_val, group = 1),
              inherit.aes = FALSE, colour = "grey30", linewidth = 0.4) +
    geom_text(data = lab, aes(x = state, y = y_stars, label = stars), inherit.aes = FALSE,
              vjust = 0, size = (base_pt + 1) / pt, fontface = "bold") +
    geom_text(data = lab, aes(x = state, y = y_d, label = d_txt), inherit.aes = FALSE,
              vjust = 0, size = (base_pt - 1.5) / pt) +
    scale_fill_manual(values = c("Control" = control_color, "ACAM" = spec$measure_color)) +
    scale_x_discrete(labels = x_labels_fig, drop = FALSE) +
    scale_y_continuous(labels = axis_labels(y_range)) +
    coord_cartesian(ylim = c(y_range[1], y_lim_top), clip = "off") +
    guides(fill = "none") +
    labs(x = NULL, y = y_axis_label(metric, spec$y_label)) +
    theme_classic(base_size = base_pt, base_family = "sans") +
    theme(axis.text    = element_text(colour = "black", size = base_pt - 1),
          axis.title.y = element_text(colour = "black", size = base_pt),
          axis.line    = element_line(linewidth = 0.3),
          axis.ticks   = element_line(linewidth = 0.3),
          plot.margin  = margin(2, 3, 2, 2, "mm"),
          plot.tag     = element_text(face = "bold", size = base_pt + 1))
}

save_figure <- function(fig, file_stem, n_panels) {
  height_mm  <- fig_row_mm * ceiling(n_panels / fig_ncol)
  out_base   <- file.path(out_dir, file_stem)
  pdf_device <- if (capabilities("cairo")) cairo_pdf else "pdf"
  ggsave(paste0(out_base, ".pdf"),  fig, device = pdf_device,
         width = fig_width_mm, height = height_mm, units = "mm", bg = "white")
  ggsave(paste0(out_base, ".png"),  fig,
         width = fig_width_mm, height = height_mm, units = "mm", dpi = 600, bg = "white")
  ggsave(paste0(out_base, ".tiff"), fig,
         width = fig_width_mm, height = height_mm, units = "mm", dpi = 600, bg = "white",
         compression = "lzw")
  message("Saved ", file_stem, " (", fig_width_mm, " x ", height_mm, " mm)")
}

for (file_stem in names(figure_sets)) {
  metrics <- figure_sets[[file_stem]]
  fig <- wrap_plots(lapply(metrics, make_panel_plot), ncol = fig_ncol) +
    plot_annotation(tag_levels = "a")
  save_figure(fig, file_stem, length(metrics))
}

message("\nDone. All outputs saved to:\n", out_dir)

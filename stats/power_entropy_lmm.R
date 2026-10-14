# ============================================================
# power_entropy_lmm.R
#
# Linear mixed models and raincloud figures for complexity and power measures
#
# Model:
#   metric ~ state + (1 | subject)
#   States in the model: c, m, cm and ACAM1..ACAM8 (cm is the mean of c and m)
#
# Contrasts:
#   - ACAM1..ACAM8 vs cm  -> BH-FDR across 8 contrasts
#   - c - m              -> no FDR
#
# Saves for each metric:
#   1) LMM_<metric>_ACAM_vs_cm_FDR_BH_with_d.csv
#   2) LMM_<metric>_c_minus_m_noFDR_with_d.csv
#   3) LMM_<metric>_EMMeans_by_state.csv
#   4) LMM_<metric>_anova_fixed_effect.csv
#   5) LMM_<metric>_LRT_vs_null.csv
#   6) DESC_<metric>_by_state.csv
#   7) FIG_<metric>_raincloud.png
#   8) FIG_<metric>_raincloud.pdf
#   9) LMM_<metric>_model.rds
#
# Additional combined supplementary outputs:
#   10) SUPP_TABLE_5s_ACAM_vs_cm_FDR_BH_with_d_combined.csv
#   11) SUPP_TABLE_5s_ACAM_vs_cm_FDR_BH_with_d_combined.tsv
#   12) SUPP_TABLE_5s_c_minus_m_noFDR_with_d_combined.csv
#   13) SUPP_TABLE_5s_c_minus_m_noFDR_with_d_combined.tsv
#
# ============================================================

# ---------------------------
# 0) Packages
# ---------------------------
pkgs <- c(
  "readr", "dplyr", "lme4", "lmerTest", "emmeans",
  "tidyr", "ggplot2", "ggdist", "tibble"
)

to_install <- pkgs[!pkgs %in% rownames(installed.packages())]
if (length(to_install) > 0) install.packages(to_install)

library(readr)
library(dplyr)
library(lme4)
library(lmerTest)
library(emmeans)
library(tidyr)
library(ggplot2)
library(ggdist)
library(tibble)

# ---------------------------
# 1) USER SETTINGS
# ---------------------------

# Input table (one row per subject x state)
csv_path <- "C:/Users/peter/OneDrive/PhD/Farmaceutisk Videnskab/EEG thesis/data/entropy_master_df/entropy_master_df_with_phenomenology_auto.csv"

# Output directory
out_dir <- "C:/Users/peter/OneDrive/PhD/Farmaceutisk Videnskab/EEG thesis/data/entropy_master_df/entropy_outputs_w_power_5s_comb"

# Measures to analyse, with plot colour and y-axis label
run_specs <- tibble::tribble(
  ~metric,                     ~measure_color,   ~y_label,
  
  # ---- Lempel-Ziv / entropy-like ----
  "LZsum_c0",                  "coral",          "LZsum",
  "LZsumPN",                   "coral2",         "LZsumN",
  
  # ---- Higuchi / DFA ----
  "HFD_kmax200",               "seagreen",       "HFD",
  "DFA_alpha",                 "seagreen2",      "DFA α",
  
  # ---- Permutation entropy ----
  "PE_1_35_tau4_s2",           "orangered",      "PeEn",
  
  # ---- Spectral entropy ----
  "SpecEn_1to35Hz",            "turquoise4",     "SpecEn",
  
  # ---- Total power (5 s) ----
  "total_delta_5s",            "firebrick3",     "Total power delta",
  "total_theta_5s",            "tomato3",        "Total power theta",
  "total_alpha_5s",            "orangered3",     "Total power alpha",
  "total_low_beta_5s",         "darkorange3",    "Total power low beta",
  "total_high_beta_5s",        "goldenrod3",     "Total power high beta",
  
  # ---- Aperiodic power (5 s) ----
  "aperiodic_delta_5s",        "darkorchid4",    "Aperiodic power delta",
  "aperiodic_theta_5s",        "darkorchid3",    "Aperiodic power theta",
  "aperiodic_alpha_5s",        "mediumorchid3",  "Aperiodic power alpha",
  "aperiodic_low_beta_5s",     "orchid4",        "Aperiodic power low beta",
  "aperiodic_high_beta_5s",    "plum4",          "Aperiodic power high beta",
  
  # ---- Oscillatory power (5 s) ----
  "oscillatory_delta_5s",      "forestgreen",    "Oscillatory power delta",
  "oscillatory_theta_5s",      "seagreen3",      "Oscillatory power theta",
  "oscillatory_alpha_5s",      "springgreen4",   "Oscillatory power alpha",
  "oscillatory_low_beta_5s",   "olivedrab4",     "Oscillatory power low beta",
  "oscillatory_high_beta_5s",  "yellowgreen",    "Oscillatory power high beta"
)

# Fixed settings
control_color <- "grey70"
add_spaghetti <- FALSE
save_width    <- 11
save_height   <- 4
save_dpi      <- 600

# State ordering used in the model
state_levels_model <- c("m", "c", "cm", paste0("ACAM", 1:8))

# State ordering used in the figure
state_levels_plot  <- c("cm", paste0("ACAM", 1:8))
acam_levels        <- paste0("ACAM", 1:8)

# Pretty x-axis labels for plot
x_labels <- c(
  "cm" = "Control",
  setNames(paste0("ACAM-J", 1:8), paste0("ACAM", 1:8))
)

# ---------------------------
# 2) Helpers
# ---------------------------

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

safe_name <- function(x) {
  gsub("[^A-Za-z0-9_]+", "_", x)
}

pick_ci <- function(d) {
  if (all(c("lower.CL", "upper.CL") %in% names(d))) {
    d <- dplyr::rename(d, d_lower = lower.CL, d_upper = upper.CL)
  } else if (all(c("asymp.LCL", "asymp.UCL") %in% names(d))) {
    d <- dplyr::rename(d, d_lower = asymp.LCL, d_upper = asymp.UCL)
  } else {
    d$d_lower <- NA_real_
    d$d_upper <- NA_real_
  }
  
  d %>%
    dplyr::select(
      contrast,
      d = effect.size,
      d_SE = SE,
      d_df = df,
      d_lower,
      d_upper
    )
}

make_desc_table <- function(df_in, metric) {
  df_in %>%
    group_by(state) %>%
    summarise(
      n_subjects = n_distinct(subject),
      n_rows     = n(),
      mean       = mean(.data[[metric]], na.rm = TRUE),
      sd         = sd(.data[[metric]], na.rm = TRUE),
      median     = median(.data[[metric]], na.rm = TRUE),
      min        = min(.data[[metric]], na.rm = TRUE),
      max        = max(.data[[metric]], na.rm = TRUE),
      .groups    = "drop"
    )
}

metric_annotation <- function(metric) {
  domain <- NA_character_
  measure_type <- NA_character_
  band <- NA_character_
  epoch_length <- NA_character_
  
  if (grepl("_5s$", metric)) {
    epoch_length <- "5s"
  }
  
  if (grepl("^total_", metric)) {
    domain <- "power"
    measure_type <- "total_power"
    band <- sub("^total_", "", metric)
    band <- sub("_5s$", "", band)
    
  } else if (grepl("^aperiodic_", metric)) {
    domain <- "power"
    measure_type <- "aperiodic_power"
    band <- sub("^aperiodic_", "", metric)
    band <- sub("_5s$", "", band)
    
  } else if (grepl("^oscillatory_", metric)) {
    domain <- "power"
    measure_type <- "oscillatory_power"
    band <- sub("^oscillatory_", "", metric)
    band <- sub("_5s$", "", band)
    
  } else if (grepl("^PE_", metric)) {
    domain <- "entropy"
    measure_type <- "permutation_entropy"
    epoch_length <- ifelse(is.na(epoch_length), "2s", epoch_length)
    
    band <- dplyr::case_when(
      metric == "PE_1_4_tau64_s32" ~ "1-4Hz",
      metric == "PE_1_8_tau32_s16" ~ "1-8Hz",
      metric == "PE_1_13_tau16_s8" ~ "1-13Hz",
      metric == "PE_1_20_tau8_s4"  ~ "1-20Hz",
      metric == "PE_1_35_tau4_s2"  ~ "1-35Hz",
      TRUE ~ NA_character_
    )
    
  } else if (grepl("^SpecEn_", metric)) {
    domain <- "entropy"
    measure_type <- "spectral_entropy"
    band <- "1-35Hz"
    epoch_length <- ifelse(is.na(epoch_length), "2s", epoch_length)
    
  } else if (grepl("^LZ", metric)) {
    domain <- "entropy"
    measure_type <- "lempel_ziv"
    band <- NA_character_
    epoch_length <- ifelse(is.na(epoch_length), "2s", epoch_length)
    
  } else if (grepl("^HFD", metric)) {
    domain <- "entropy"
    measure_type <- "higuchi_fd"
    band <- NA_character_
    epoch_length <- ifelse(is.na(epoch_length), "2s", epoch_length)
    
  } else if (grepl("^DFA", metric)) {
    domain <- "entropy"
    measure_type <- "dfa"
    band <- NA_character_
    epoch_length <- ifelse(is.na(epoch_length), "10s", epoch_length)
  }
  
  tibble(
    metric = metric,
    domain = domain,
    measure_type = measure_type,
    band = band,
    epoch_length = epoch_length
  )
}

combine_supp_table <- function(all_results,
                               run_specs,
                               which_table = c("acam", "cm"),
                               metric_pattern = "5s") {
  
  which_table <- match.arg(which_table)
  res_name <- if (which_table == "acam") "res_acam" else "res_cm"
  
  metrics_to_use <- run_specs$metric[grepl(metric_pattern, run_specs$metric, ignore.case = TRUE)]
  
  if (length(metrics_to_use) == 0) {
    warning("No metrics matched pattern: ", metric_pattern)
    return(NULL)
  }
  
  out <- bind_rows(lapply(metrics_to_use, function(metric_i) {
    
    x <- all_results[[metric_i]][[res_name]]
    
    if (is.null(x) || nrow(x) == 0) return(NULL)
    
    # Standardize CI column names if needed
    if (!all(c("lower.CL", "upper.CL") %in% names(x))) {
      if (all(c("asymp.LCL", "asymp.UCL") %in% names(x))) {
        x <- x %>%
          rename(
            lower.CL = asymp.LCL,
            upper.CL = asymp.UCL
          )
      }
    }
    
    # c - m table does not have BH-FDR, so create an empty column for consistency
    if (!"p_fdr_bh" %in% names(x)) {
      x$p_fdr_bh <- NA_real_
    }
    
    x %>%
      mutate(metric = metric_i) %>%
      left_join(metric_annotation(metric_i), by = "metric") %>%
      select(
        metric,
        domain,
        measure_type,
        band,
        epoch_length,
        contrast,
        estimate,
        lower.CL,
        upper.CL,
        p_fdr_bh,
        d,
        d_lower,
        d_upper
      )
  }))
  
  if (is.null(out) || nrow(out) == 0) {
    warning("No rows available for supplementary table: ", which_table)
    return(NULL)
  }
  
  out %>%
    mutate(metric = factor(metric, levels = run_specs$metric)) %>%
    arrange(metric, contrast) %>%
    mutate(metric = as.character(metric))
}

run_entropy_lmm <- function(metric, df_in, out_dir) {
  
  if (!metric %in% names(df_in)) {
    stop(
      paste0(
        "Metric '", metric, "' not found in dataframe.\n",
        "Available columns:\n",
        paste(names(df_in), collapse = ", ")
      )
    )
  }
  
  message("\n==============================")
  message("Running LMM for: ", metric)
  message("==============================")
  
  # Keep needed rows only
  df0 <- df_in %>%
    filter(!is.na(.data[[metric]]), !is.na(subject), !is.na(state)) %>%
    mutate(
      subject = factor(subject),
      state   = factor(state, levels = state_levels_model)
    ) %>%
    filter(!is.na(state)) %>%
    mutate(state = relevel(state, ref = "cm"))
  
  if (nrow(df0) == 0) {
    stop("No non-missing rows left for metric: ", metric)
  }
  
  # Fit model
  fml  <- as.formula(paste0(metric, " ~ state + (1 | subject)"))
  fml0 <- as.formula(paste0(metric, " ~ 1 + (1 | subject)"))
  
  m_fit <- lmer(
    fml,
    data = df0,
    REML = TRUE,
    control = lmerControl(optimizer = "bobyqa")
  )
  
  null_fit <- lmer(
    fml0,
    data = df0,
    REML = TRUE,
    control = lmerControl(optimizer = "bobyqa")
  )
  
  # Fixed effect ANOVA + LRT vs null
  anova_fixed <- anova(m_fit)
  lrt_vs_null <- anova(null_fit, m_fit)
  
  # EMMs
  emm <- emmeans(m_fit, ~ state, lmer.df = "satterthwaite")
  
  # ACAM vs cm
  con_trt <- contrast(emm, method = "trt.vs.ctrl", ref = "cm", adjust = "none")
  res_all <- as.data.frame(summary(con_trt, infer = TRUE))
  
  res_acam <- res_all %>%
    filter(grepl("^ACAM[1-8]\\s*-\\s*cm$", contrast)) %>%
    mutate(p_fdr_bh = p.adjust(p.value, method = "BH"))
  
  # c - m
  states <- levels(df0$state)
  w_cm <- setNames(rep(0, length(states)), states)
  w_cm["c"] <-  1
  w_cm["m"] <- -1
  
  con_cm <- contrast(emm, method = list("c - m" = w_cm), adjust = "none")
  res_cm <- as.data.frame(summary(con_cm, infer = TRUE))
  
  # Effect sizes (Cohen's d: contrast / residual SD)
  sig <- sigma(m_fit)
  edf <- df.residual(m_fit)
  
  es_trt <- as.data.frame(
    eff_size(
      emm,
      sigma = sig,
      edf   = edf,
      method = "trt.vs.ctrl",
      ref = "cm"
    )
  )
  
  es_cm <- as.data.frame(
    eff_size(
      emm,
      sigma = sig,
      edf   = edf,
      method = list("c - m" = w_cm)
    )
  )
  
  es_acam <- es_trt %>%
    filter(grepl("^ACAM[1-8]\\s*-\\s*cm$", contrast))
  
  es_acam2 <- pick_ci(es_acam)
  es_cm2   <- pick_ci(es_cm)
  
  res_acam <- res_acam %>% left_join(es_acam2, by = "contrast")
  res_cm   <- res_cm   %>% left_join(es_cm2, by = "contrast")
  
  # Descriptives
  desc_tab <- make_desc_table(df0, metric)
  
  # Save files
  metric_safe <- safe_name(metric)
  
  out_acam <- file.path(out_dir, paste0("LMM_", metric_safe, "_ACAM_vs_cm_FDR_BH_with_d.csv"))
  out_cm   <- file.path(out_dir, paste0("LMM_", metric_safe, "_c_minus_m_noFDR_with_d.csv"))
  out_emm  <- file.path(out_dir, paste0("LMM_", metric_safe, "_EMMeans_by_state.csv"))
  out_anv  <- file.path(out_dir, paste0("LMM_", metric_safe, "_anova_fixed_effect.csv"))
  out_lrt  <- file.path(out_dir, paste0("LMM_", metric_safe, "_LRT_vs_null.csv"))
  out_desc <- file.path(out_dir, paste0("DESC_", metric_safe, "_by_state.csv"))
  out_rds  <- file.path(out_dir, paste0("LMM_", metric_safe, "_model.rds"))
  
  write_csv(res_acam, out_acam)
  write_csv(res_cm, out_cm)
  write_csv(as.data.frame(emm), out_emm)
  
  write_csv(
    tibble::rownames_to_column(as.data.frame(anova_fixed), var = "term"),
    out_anv
  )
  
  write_csv(
    tibble::rownames_to_column(as.data.frame(lrt_vs_null), var = "model"),
    out_lrt
  )
  
  write_csv(desc_tab, out_desc)
  saveRDS(m_fit, out_rds)
  
  message("Saved model outputs for ", metric)
  
  list(
    metric       = metric,
    model        = m_fit,
    emm          = emm,
    res_acam     = res_acam,
    res_cm       = res_cm,
    anova_fixed  = anova_fixed,
    lrt_vs_null  = lrt_vs_null,
    desc         = desc_tab
  )
}

make_entropy_plot <- function(metric,
                              df_in,
                              res_acam,
                              measure_color,
                              y_label = metric,
                              out_dir,
                              control_color = "grey70",
                              add_spaghetti = FALSE,
                              save_width = 11,
                              save_height = 6,
                              save_dpi = 300) {
  
  metric_safe <- safe_name(metric)
  
  # Keep only cm + ACAM1:8 for plotting
  dfp <- df_in %>%
    mutate(
      subject = factor(subject),
      state   = factor(state, levels = state_levels_plot)
    ) %>%
    filter(state %in% state_levels_plot) %>%
    filter(!is.na(.data[[metric]])) %>%
    mutate(fill_group = ifelse(state == "cm", "Control", "ACAM"))
  
  if (nrow(dfp) == 0) {
    warning("No rows available for plotting metric: ", metric)
    return(NULL)
  }
  
  # Stars from FDR-adjusted p-values
  sig <- res_acam %>%
    mutate(
      state = sub("\\s*-\\s*cm$", "", contrast),
      stars = case_when(
        p_fdr_bh < 0.001 ~ "***",
        p_fdr_bh < 0.01  ~ "**",
        p_fdr_bh < 0.05  ~ "*",
        TRUE ~ ""
      )
    ) %>%
    filter(grepl("^ACAM[1-8]$", state)) %>%
    select(state, stars)
  
  sig <- tibble(state = acam_levels) %>%
    left_join(sig, by = "state") %>%
    mutate(
      stars = ifelse(is.na(stars), "", stars),
      state = factor(state, levels = acam_levels)
    )
  
  # d labels
  d_lab <- res_acam %>%
    mutate(
      state = sub("\\s*-\\s*cm$", "", contrast),
      d_txt = ifelse(is.na(d), "", paste0("d=", formatC(d, format = "f", digits = 2)))
    ) %>%
    filter(grepl("^ACAM[1-8]$", state)) %>%
    select(state, d_txt)
  
  d_lab <- tibble(state = acam_levels) %>%
    left_join(d_lab, by = "state") %>%
    mutate(
      d_txt = ifelse(is.na(d_txt), "", d_txt),
      state = factor(state, levels = acam_levels)
    )
  
  # y placements
  y_range <- range(dfp[[metric]], na.rm = TRUE)
  y_diff  <- diff(y_range)
  
  if (is.na(y_diff) || y_diff == 0) {
    y_diff <- max(abs(y_range), na.rm = TRUE) * 0.15
    if (!is.finite(y_diff) || y_diff == 0) y_diff <- 1
  }
  
  y_top     <- y_range[2] + 0.08 * y_diff
  y_mid     <- y_range[2] + 0.03 * y_diff
  y_lim_top <- y_range[2] + 0.22 * y_diff
  
  sig$y   <- y_top
  d_lab$y <- y_mid
  
  # Mean per state
  mean_df <- dfp %>%
    group_by(state) %>%
    summarise(mean_val = mean(.data[[metric]], na.rm = TRUE), .groups = "drop") %>%
    mutate(state = factor(state, levels = state_levels_plot))
  
  p <- ggplot(dfp, aes(x = state, y = .data[[metric]], fill = fill_group)) +
    
    {if (add_spaghetti)
      geom_line(
        aes(group = subject),
        color = "grey85",
        linewidth = 0.35,
        alpha = 0.6
      )
    } +
    
    stat_halfeye(
      adjust = 0.8,
      width = 0.75,
      justification = 0,
      point_colour = NA,
      alpha = 0.7
    ) +
    
    geom_boxplot(
      inherit.aes = FALSE,
      aes(x = state, y = .data[[metric]]),
      width = 0.18,
      outlier.shape = NA,
      fill = "white",
      color = "black",
      alpha = 1
    ) +
    
    geom_point(
      position = position_jitter(width = 0.12, height = 0),
      alpha = 0.6,
      size = 1.4
    ) +
    
    geom_line(
      data = mean_df,
      aes(x = state, y = mean_val, group = 1),
      inherit.aes = FALSE,
      color = "grey30",
      linewidth = 0.8
    ) +
    
    geom_text(
      data = sig,
      aes(x = state, y = y, label = stars),
      inherit.aes = FALSE,
      vjust = 0,
      size = 9,
      fontface = "bold"
    ) +
    
    geom_text(
      data = d_lab,
      aes(x = state, y = y, label = d_txt),
      inherit.aes = FALSE,
      vjust = 0,
      size = 4.5
    ) +
    
    scale_fill_manual(values = c("Control" = control_color, "ACAM" = measure_color)) +
    scale_x_discrete(labels = x_labels) +
    scale_y_continuous(labels = function(x) sprintf("%.2f", x)) +
    guides(fill = "none", color = "none") +
    coord_cartesian(ylim = c(y_range[1], y_lim_top)) +
    theme_classic(base_size = 16) +
    theme(
      axis.text.x = element_text(angle = 0, hjust = 0.5, face = "bold", colour = "black"),
      axis.text.y = element_text(face = "bold", colour = "black"),
      axis.title.x = element_text(face = "bold", colour = "black"),
      axis.title.y = element_text(face = "bold", colour = "black")
    ) +
    labs(x = NULL, y = y_label)
  
  # Save
  out_png <- file.path(out_dir, paste0("FIG_", metric_safe, "_raincloud.png"))
  out_pdf <- file.path(out_dir, paste0("FIG_", metric_safe, "_raincloud.pdf"))
  
  ggsave(out_png, plot = p, width = save_width, height = save_height, dpi = save_dpi, bg = "white")
  ggsave(out_pdf, plot = p, width = save_width, height = save_height, bg = "white")
  
  message("Saved figure for ", metric)
  print(p)
  
  invisible(p)
}

# ---------------------------
# 3) Load data
# ---------------------------
df <- read_csv(csv_path, show_col_types = FALSE) %>%
  mutate(
    subject = factor(subject),
    state   = factor(state)
  )

message("\nLoaded dataframe with columns:")
print(names(df))

# Print the numeric columns
candidate_metrics <- names(df)[sapply(df, is.numeric)]
candidate_metrics <- setdiff(candidate_metrics, c("n_epochs", "DFA_n_epochs"))
message("\nDetected numeric metric columns:")
print(candidate_metrics)

# ---------------------------
# 4) Run everything
# ---------------------------
all_results <- vector("list", nrow(run_specs))
names(all_results) <- run_specs$metric

for (i in seq_len(nrow(run_specs))) {
  
  metric_i <- run_specs$metric[i]
  color_i  <- run_specs$measure_color[i]
  ylab_i   <- run_specs$y_label[i]
  
  # run LMM
  res_i <- run_entropy_lmm(
    metric = metric_i,
    df_in = df,
    out_dir = out_dir
  )
  
  # make plot
  make_entropy_plot(
    metric         = metric_i,
    df_in          = df,
    res_acam       = res_i$res_acam,
    measure_color  = color_i,
    y_label        = ylab_i,
    out_dir        = out_dir,
    control_color  = control_color,
    add_spaghetti  = add_spaghetti,
    save_width     = save_width,
    save_height    = save_height,
    save_dpi       = save_dpi
  )
  
  all_results[[metric_i]] <- res_i
}

# ---------------------------
# 5) Combined supplementary tables for 5s metrics
# ---------------------------
supp_5s_acam <- combine_supp_table(
  all_results = all_results,
  run_specs = run_specs,
  which_table = "acam",
  metric_pattern = "5s"
)

supp_5s_cm <- combine_supp_table(
  all_results = all_results,
  run_specs = run_specs,
  which_table = "cm",
  metric_pattern = "5s"
)

if (!is.null(supp_5s_acam) && nrow(supp_5s_acam) > 0) {
  write_csv(
    supp_5s_acam,
    file.path(out_dir, "SUPP_TABLE_5s_ACAM_vs_cm_FDR_BH_with_d_combined.csv")
  )
  
  write_tsv(
    supp_5s_acam,
    file.path(out_dir, "SUPP_TABLE_5s_ACAM_vs_cm_FDR_BH_with_d_combined.tsv")
  )
}

if (!is.null(supp_5s_cm) && nrow(supp_5s_cm) > 0) {
  write_csv(
    supp_5s_cm,
    file.path(out_dir, "SUPP_TABLE_5s_c_minus_m_noFDR_with_d_combined.csv")
  )
  
  write_tsv(
    supp_5s_cm,
    file.path(out_dir, "SUPP_TABLE_5s_c_minus_m_noFDR_with_d_combined.tsv")
  )
}

message("\nSaved combined 5s supplementary tables.")
message("\nDone. All outputs saved to:\n", out_dir)
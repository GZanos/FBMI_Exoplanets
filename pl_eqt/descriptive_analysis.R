################################################################################
# Descriptive Analysis: Exoplanet Equilibrium Temperature Dataset
# Target: pl_eqt | Predictors: all other columns except pl_name
# Note: pl_name identifies planets; some planets have multiple measurements
################################################################################

# ---- Packages ----
required_pkgs <- c(
  "tidyverse", "ggplot2", "dplyr", "tidyr", "readr",
  "corrplot", "GGally", "scales", "gridExtra", "skimr"
)

install_if_missing <- function(pkgs) {
  missing <- pkgs[!pkgs %in% rownames(installed.packages())]
  if (length(missing) > 0) {
    options(repos = c(CRAN = "https://cloud.r-project.org"))
    install.packages(missing, dependencies = TRUE)
  }
}
install_if_missing(required_pkgs)

library(tidyverse)
library(corrplot)
library(GGally)
library(scales)
library(gridExtra)
library(skimr)

# ---- Paths ----
data_path  <- "planet_eqt_train_all_features.csv"
plot_dir   <- "descriptive_plots"
dir.create(plot_dir, showWarnings = FALSE)

save_plot <- function(p, filename, width = 10, height = 7) {
  ggsave(
    filename = file.path(plot_dir, filename),
    plot = p, width = width, height = height, dpi = 300, bg = "white"
  )
  message("Saved: ", file.path(plot_dir, filename))
}

theme_set(theme_bw(base_size = 12) +
            theme(
              plot.title = element_text(face = "bold", size = 14),
              plot.subtitle = element_text(color = "grey40"),
              legend.position = "bottom"
            ))

# ---- Load data ----
raw <- read_csv(data_path, show_col_types = FALSE)

id_col     <- "pl_name"
target_col <- "pl_eqt"
pred_cols  <- setdiff(names(raw), c(id_col, target_col))

# Unique-planet view: keep first row per planet for planet-level summaries
# (duplicates = multiple measurements of the same planet)
unique_planets <- raw %>%
  group_by(.data[[id_col]]) %>%
  slice(1) %>%
  ungroup()

# Measurement multiplicity
meas_counts <- raw %>%
  count(.data[[id_col]], name = "n_measurements")

cat("\n========== DATASET OVERVIEW ==========\n")
cat("Total rows (measurements):     ", nrow(raw), "\n")
cat("Unique planets (pl_name):      ", n_distinct(raw[[id_col]]), "\n")
cat("Planets with >1 measurement:   ", sum(meas_counts$n_measurements > 1), "\n")
cat("Max measurements per planet:   ", max(meas_counts$n_measurements), "\n")
cat("Predictors:                    ", length(pred_cols), "\n")
cat("Target:                        ", target_col, "\n\n")

# ============================================================================
# 1. MISSINGNESS
# ============================================================================

missing_row <- tibble(
  variable = names(raw),
  n_missing_rows = colSums(is.na(raw)),
  pct_missing_rows = 100 * colSums(is.na(raw)) / nrow(raw)
) %>%
  arrange(desc(pct_missing_rows))

# Planet-level missingness (unique planets)
missing_planet <- tibble(
  variable = names(unique_planets),
  n_missing_planets = colSums(is.na(unique_planets)),
  pct_missing_planets = 100 * colSums(is.na(unique_planets)) / nrow(unique_planets)
) %>%
  arrange(desc(pct_missing_planets))

missing_summary <- missing_row %>%
  left_join(missing_planet, by = "variable")

cat("========== MISSINGNESS SUMMARY ==========\n")
print(as.data.frame(missing_summary), row.names = FALSE)
cat("\n")

# Target missingness among unique planets
n_unique <- nrow(unique_planets)
n_target_obs   <- sum(!is.na(unique_planets[[target_col]]))
n_target_miss  <- sum(is.na(unique_planets[[target_col]]))

cat("Target (pl_eqt) among UNIQUE planets:\n")
cat("  Observed: ", n_target_obs,  " (", round(100 * n_target_obs / n_unique, 1), "%)\n", sep = "")
cat("  Missing:  ", n_target_miss, " (", round(100 * n_target_miss / n_unique, 1), "%)\n\n", sep = "")

# Missingness bar plot
p_miss <- missing_summary %>%
  filter(n_missing_rows > 0 | n_missing_planets > 0) %>%
  pivot_longer(
    cols = c(pct_missing_rows, pct_missing_planets),
    names_to = "level", values_to = "pct"
  ) %>%
  mutate(level = recode(level,
                        pct_missing_rows = "Rows (measurements)",
                        pct_missing_planets = "Unique planets")) %>%
  ggplot(aes(x = reorder(variable, pct), y = pct, fill = level)) +
  geom_col(position = "dodge", width = 0.7) +
  coord_flip() +
  labs(
    title = "Missing Values by Variable",
    subtitle = paste0(
      "Only variables with missing data shown | ",
      n_distinct(raw[[id_col]]), " unique planets, ", nrow(raw), " rows"
    ),
    x = NULL, y = "% Missing", fill = NULL
  ) +
  scale_fill_manual(values = c("#2C7BB6", "#D7191C"))

if (nrow(missing_summary %>% filter(n_missing_rows > 0)) == 0) {
  # Still plot target if somehow filtered; otherwise create informative empty note
  p_miss <- missing_summary %>%
    filter(variable == target_col) %>%
    pivot_longer(cols = c(pct_missing_rows, pct_missing_planets),
                 names_to = "level", values_to = "pct") %>%
    mutate(level = recode(level,
                          pct_missing_rows = "Rows (measurements)",
                          pct_missing_planets = "Unique planets")) %>%
    ggplot(aes(x = variable, y = pct, fill = level)) +
    geom_col(position = "dodge", width = 0.6) +
    labs(
      title = "Missing Values — Target Variable (pl_eqt)",
      subtitle = "All predictors are complete; only the target has missing values",
      x = NULL, y = "% Missing", fill = NULL
    ) +
    scale_fill_manual(values = c("#2C7BB6", "#D7191C")) +
    ylim(0, 100)
}
save_plot(p_miss, "01_missingness.png", width = 9, height = 5)

# ============================================================================
# 2. DUPLICATE / MEASUREMENT MULTIPLICITY
# ============================================================================

p_meas <- meas_counts %>%
  count(n_measurements) %>%
  ggplot(aes(x = factor(n_measurements), y = n)) +
  geom_col(fill = "#4575B4", width = 0.6) +
  geom_text(aes(label = n), vjust = -0.4, size = 3.5) +
  labs(
    title = "Number of Measurements per Planet",
    subtitle = paste0(
      sum(meas_counts$n_measurements > 1), " planets have multiple rows | ",
      n_distinct(raw[[id_col]]), " unique planets total"
    ),
    x = "Measurements per planet", y = "Number of unique planets"
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12)))

save_plot(p_meas, "02_measurements_per_planet.png", width = 8, height = 5)

# Among duplicated planets, how often does pl_eqt differ?
dup_mass_var <- raw %>%
  group_by(.data[[id_col]]) %>%
  filter(n() > 1) %>%
  summarise(
    n_meas = n(),
    n_unique_teq = n_distinct(.data[[target_col]], na.rm = FALSE),
    teq_range = if (all(is.na(.data[[target_col]]))) NA_real_
                 else max(.data[[target_col]], na.rm = TRUE) -
                      min(.data[[target_col]], na.rm = TRUE),
    .groups = "drop"
  )

cat("========== DUPLICATE TEQ VARIATION ==========\n")
cat("Duplicated planets:                    ", nrow(dup_mass_var), "\n")
cat("With differing pl_eqt across rows:  ",
    sum(dup_mass_var$n_unique_teq > 1), "\n")
cat("Teq range summary (where Teq varies):\n")
print(summary(dup_mass_var$teq_range[dup_mass_var$n_unique_teq > 1]))
cat("\n")

# ============================================================================
# 3. TARGET VARIABLE — DESCRIPTIVE STATS & DISTRIBUTIONS
# ============================================================================

# Planet-level target (unique planets)
target_unique <- unique_planets[[target_col]]
target_obs    <- target_unique[!is.na(target_unique)]

desc_stats <- function(x, label) {
  tibble(
    level = label,
    n = length(x),
    n_missing = sum(is.na(x)),
    mean = mean(x, na.rm = TRUE),
    sd = sd(x, na.rm = TRUE),
    min = min(x, na.rm = TRUE),
    q25 = quantile(x, 0.25, na.rm = TRUE),
    median = median(x, na.rm = TRUE),
    q75 = quantile(x, 0.75, na.rm = TRUE),
    max = max(x, na.rm = TRUE),
    skewness = {
      m <- mean(x, na.rm = TRUE); s <- sd(x, na.rm = TRUE)
      mean(((x - m) / s)^3, na.rm = TRUE)
    },
    kurtosis = {
      m <- mean(x, na.rm = TRUE); s <- sd(x, na.rm = TRUE)
      mean(((x - m) / s)^4, na.rm = TRUE) - 3
    }
  )
}

target_stats <- bind_rows(
  desc_stats(raw[[target_col]], "All rows (measurements)"),
  desc_stats(target_unique, "Unique planets")
)

cat("========== TARGET: pl_eqt ==========\n")
print(as.data.frame(target_stats), row.names = FALSE)
cat("\n")

# Histogram + density (linear)
p_target_hist <- unique_planets %>%
  filter(!is.na(.data[[target_col]])) %>%
  ggplot(aes(x = .data[[target_col]])) +
  geom_histogram(aes(y = after_stat(density)), bins = 40,
                 fill = "#74ADD1", color = "white", alpha = 0.85) +
  geom_density(color = "#D73027", linewidth = 1) +
  labs(
    title = "Distribution of Equilibrium Temperature (pl_eqt)",
    subtitle = paste0("Unique planets with observed Teq: n = ", length(target_obs)),
    x = expression(T[eq]~(K)), y = "Density"
  )
save_plot(p_target_hist, "03_target_histogram.png", width = 9, height = 5)

# Log10 histogram (Teq is right-skewed)
p_target_log <- unique_planets %>%
  filter(!is.na(.data[[target_col]]), .data[[target_col]] > 0) %>%
  ggplot(aes(x = .data[[target_col]])) +
  geom_histogram(aes(y = after_stat(density)), bins = 40,
                 fill = "#FDAE61", color = "white", alpha = 0.85) +
  geom_density(color = "#D73027", linewidth = 1) +
  scale_x_log10(labels = label_comma()) +
  annotation_logticks(sides = "b") +
  labs(
    title = "Distribution of Equilibrium Temperature (log10 scale)",
    subtitle = "Right-skewed Teq values are clearer on a log scale",
    x = expression(T[eq]~(K)~"[log10]"), y = "Density"
  )
save_plot(p_target_log, "04_target_histogram_log10.png", width = 9, height = 5)

# Boxplot + violin
p_target_box <- unique_planets %>%
  filter(!is.na(.data[[target_col]]), .data[[target_col]] > 0) %>%
  ggplot(aes(x = "", y = .data[[target_col]])) +
  geom_violin(fill = "#ABD9E9", alpha = 0.7, color = NA) +
  geom_boxplot(width = 0.15, fill = "white", outlier.alpha = 0.4) +
  scale_y_log10(labels = label_comma()) +
  labs(
    title = "Equilibrium Temperature — Violin & Boxplot (log10)",
    subtitle = paste0("Unique planets, n = ", length(target_obs)),
    x = NULL, y = expression(T[eq]~(K))
  )
save_plot(p_target_box, "05_target_violin_boxplot.png", width = 6, height = 7)

# QQ plot (log Teq)
p_qq <- unique_planets %>%
  filter(!is.na(.data[[target_col]]), .data[[target_col]] > 0) %>%
  ggplot(aes(sample = log10(.data[[target_col]]))) +
  stat_qq(alpha = 0.5, color = "#4575B4") +
  stat_qq_line(color = "#D73027", linewidth = 0.9) +
  labs(
    title = "Normal Q–Q Plot of log10(pl_eqt)",
    subtitle = "Unique planets with observed Teq",
    x = "Theoretical quantiles", y = "Sample quantiles (log10 Teq)"
  )
save_plot(p_qq, "06_target_qq_log10.png", width = 7, height = 6)

# Observed vs missing target — compare key predictors
unique_planets <- unique_planets %>%
  mutate(teq_status = ifelse(is.na(.data[[target_col]]), "Missing Teq", "Observed Teq"))

cat("Unique planets by Teq availability:\n")
print(table(unique_planets$teq_status))
cat("\n")

# ============================================================================
# 4. PREDICTOR DESCRIPTIVE STATISTICS
# ============================================================================

pred_stats_unique <- unique_planets %>%
  select(all_of(pred_cols)) %>%
  summarise(across(everything(), list(
    n = ~sum(!is.na(.)),
    mean = ~mean(., na.rm = TRUE),
    sd = ~sd(., na.rm = TRUE),
    min = ~min(., na.rm = TRUE),
    q25 = ~quantile(., 0.25, na.rm = TRUE),
    median = ~median(., na.rm = TRUE),
    q75 = ~quantile(., 0.75, na.rm = TRUE),
    max = ~max(., na.rm = TRUE)
  ), .names = "{.col}__{.fn}")) %>%
  pivot_longer(everything(), names_to = "key", values_to = "value") %>%
  separate(key, into = c("variable", "stat"), sep = "__") %>%
  pivot_wider(names_from = stat, values_from = value) %>%
  arrange(variable)

cat("========== PREDICTOR SUMMARY (unique planets) ==========\n")
print(as.data.frame(pred_stats_unique), row.names = FALSE)
write_csv(pred_stats_unique, file.path(plot_dir, "predictor_summary_unique_planets.csv"))
write_csv(target_stats, file.path(plot_dir, "target_summary.csv"))
write_csv(missing_summary, file.path(plot_dir, "missingness_summary.csv"))
cat("Summary tables written to", plot_dir, "\n\n")

# Skim overview
cat("========== SKIMR OVERVIEW (unique planets) ==========\n")
print(skim(unique_planets %>% select(all_of(c(target_col, pred_cols)))))
cat("\n")

# ============================================================================
# 5. PREDICTOR DISTRIBUTIONS
# ============================================================================

# Group predictors for readable faceted plots
planet_phys <- c("pl_rade", "pl_orbper", "pl_orbsmax")
stellar     <- c("st_teff", "st_mass", "st_rad", "st_logg", "st_age")
system_cnt  <- c("sy_snum", "sy_pnum")
coords      <- c("ra", "dec", "glat", "glon", "elat", "elon")
mags        <- c("sy_vmag", "sy_jmag", "sy_hmag", "sy_kmag", "sy_tmag")
mag_errs    <- c("sy_vmagerr1", "sy_vmagerr2", "sy_tmagerr1", "sy_tmagerr2")
counts_meta <- c("pl_nnotes", "st_nphot", "st_nrvc", "st_nspec", "pl_nespec", "pl_ntranspec")

plot_hist_facet <- function(data, vars, title, filename,
                            log_vars = character(0), ncol = 2) {
  long <- data %>%
    select(all_of(vars)) %>%
    pivot_longer(everything(), names_to = "variable", values_to = "value") %>%
    filter(!is.na(value))

  p <- ggplot(long, aes(x = value)) +
    geom_histogram(bins = 30, fill = "#74ADD1", color = "white", alpha = 0.9) +
    facet_wrap(~ variable, scales = "free", ncol = ncol) +
    labs(title = title, x = NULL, y = "Count",
         subtitle = paste0("Unique planets, n = ", nrow(data)))

  # Optional log scales via annotation in title for heavily skewed vars
  if (length(log_vars) > 0) {
    long_log <- data %>%
      select(all_of(intersect(vars, log_vars))) %>%
      pivot_longer(everything(), names_to = "variable", values_to = "value") %>%
      filter(!is.na(value), value > 0)

    p_log <- ggplot(long_log, aes(x = value)) +
      geom_histogram(bins = 30, fill = "#FDAE61", color = "white", alpha = 0.9) +
      scale_x_log10(labels = label_comma()) +
      facet_wrap(~ variable, scales = "free", ncol = ncol) +
      labs(
        title = paste0(title, " (log10)"),
        subtitle = paste0("Unique planets, n = ", nrow(data)),
        x = NULL, y = "Count"
      )
    save_plot(p_log, gsub("\\.png$", "_log10.png", filename),
              width = 10, height = 2.8 * ceiling(length(intersect(vars, log_vars)) / ncol))
  }

  save_plot(p, filename,
            width = 10,
            height = 2.8 * ceiling(length(vars) / ncol))
}

plot_hist_facet(unique_planets, planet_phys,
                "Planet Physical Parameters",
                "07_dist_planet_physical.png",
                log_vars = c("pl_orbper", "pl_orbsmax", "pl_rade"))

plot_hist_facet(unique_planets, stellar,
                "Stellar Parameters",
                "08_dist_stellar.png",
                log_vars = c("st_teff", "st_mass", "st_rad", "st_age"),
                ncol = 3)

plot_hist_facet(unique_planets, mags,
                "System Magnitudes",
                "09_dist_magnitudes.png", ncol = 3)

plot_hist_facet(unique_planets, mag_errs,
                "Magnitude Errors",
                "10_dist_magnitude_errors.png", ncol = 2)

plot_hist_facet(unique_planets, coords,
                "Sky Coordinates",
                "11_dist_coordinates.png", ncol = 3)

plot_hist_facet(unique_planets, counts_meta,
                "Count / Metadata Features",
                "12_dist_count_features.png", ncol = 3)

# Discrete system counts
p_sy <- unique_planets %>%
  select(all_of(system_cnt)) %>%
  pivot_longer(everything(), names_to = "variable", values_to = "value") %>%
  count(variable, value) %>%
  ggplot(aes(x = factor(value), y = n)) +
  geom_col(fill = "#4575B4", width = 0.7) +
  geom_text(aes(label = n), vjust = -0.3, size = 3) +
  facet_wrap(~ variable, scales = "free_x") +
  labs(
    title = "System Architecture Counts",
    subtitle = paste0("Unique planets, n = ", nrow(unique_planets)),
    x = "Value", y = "Number of planets"
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15)))
save_plot(p_sy, "13_system_counts.png", width = 9, height = 5)

# ============================================================================
# 6. TARGET vs KEY PREDICTORS
# ============================================================================

# Temperature–radius
p_mr <- unique_planets %>%
  filter(!is.na(.data[[target_col]]), !is.na(pl_rade),
         .data[[target_col]] > 0, pl_rade > 0) %>%
  ggplot(aes(x = pl_rade, y = .data[[target_col]], color = teq_status)) +
  geom_point(alpha = 0.55, size = 1.8) +
  scale_x_log10() + scale_y_log10(labels = label_comma()) +
  scale_color_manual(values = c("Observed Teq" = "#2C7BB6")) +
  labs(
    title = "Temperature–Radius Relation",
    subtitle = "Unique planets with observed Teq",
    x = expression(Planet~radius~(R[Earth])),
    y = expression(T[eq]~(K)),
    color = NULL
  ) +
  guides(color = "none")
save_plot(p_mr, "14_teq_radius.png", width = 8, height = 6)

# Teq vs stellar effective temperature
p_meqt <- unique_planets %>%
  filter(!is.na(.data[[target_col]]), .data[[target_col]] > 0) %>%
  ggplot(aes(x = st_teff, y = .data[[target_col]])) +
  geom_point(alpha = 0.45, color = "#2C7BB6", size = 1.6) +
  scale_y_log10(labels = label_comma()) +
  labs(
    title = "Teq vs Stellar Effective Temperature",
    subtitle = "Unique planets with observed Teq",
    x = "Stellar Teff (K)",
    y = expression(T[eq]~(K)~"[log10]")
  )
save_plot(p_meqt, "15_teq_vs_st_teff.png", width = 8, height = 6)

# Teq vs orbital period
p_mper <- unique_planets %>%
  filter(!is.na(.data[[target_col]]), .data[[target_col]] > 0, pl_orbper > 0) %>%
  ggplot(aes(x = pl_orbper, y = .data[[target_col]])) +
  geom_point(alpha = 0.45, color = "#D7191C", size = 1.6) +
  scale_x_log10() + scale_y_log10(labels = label_comma()) +
  labs(
    title = "Teq vs Orbital Period",
    subtitle = "Unique planets with observed Teq",
    x = "Orbital period (days) [log10]",
    y = expression(T[eq]~(K)~"[log10]")
  )
save_plot(p_mper, "16_teq_vs_orbper.png", width = 8, height = 6)

# Teq vs stellar parameters
p_mst <- unique_planets %>%
  filter(!is.na(.data[[target_col]]), .data[[target_col]] > 0) %>%
  select(pl_eqt, st_mass, st_teff, st_rad, st_age) %>%
  pivot_longer(-pl_eqt, names_to = "variable", values_to = "value") %>%
  filter(!is.na(value), value > 0) %>%
  ggplot(aes(x = value, y = pl_eqt)) +
  geom_point(alpha = 0.35, color = "#4575B4", size = 1.2) +
  scale_y_log10(labels = label_comma()) +
  facet_wrap(~ variable, scales = "free_x") +
  labs(
    title = "Teq vs Stellar Parameters",
    subtitle = "Unique planets with observed Teq",
    x = NULL, y = expression(T[eq]~(K)~"[log10]")
  )
save_plot(p_mst, "17_teq_vs_stellar.png", width = 10, height = 7)

# Compare predictors: observed vs missing Teq
compare_vars <- c("pl_rade", "pl_orbsmax", "pl_orbper", "st_teff", "st_mass", "sy_vmag")
p_obs_miss <- unique_planets %>%
  select(teq_status, all_of(compare_vars)) %>%
  pivot_longer(-teq_status, names_to = "variable", values_to = "value") %>%
  filter(!is.na(value)) %>%
  ggplot(aes(x = teq_status, y = value, fill = teq_status)) +
  geom_boxplot(outlier.alpha = 0.3, alpha = 0.75) +
  facet_wrap(~ variable, scales = "free_y", ncol = 3) +
  scale_fill_manual(values = c("Observed Teq" = "#74ADD1", "Missing Teq" = "#F46D43")) +
  labs(
    title = "Predictor Distributions by Target Availability",
    subtitle = paste0(
      "Unique planets | Observed: ", n_target_obs,
      " | Missing: ", n_target_miss
    ),
    x = NULL, y = NULL, fill = NULL
  ) +
  theme(axis.text.x = element_text(angle = 15, hjust = 1))
save_plot(p_obs_miss, "18_predictors_by_teq_status.png", width = 11, height = 8)

# ============================================================================
# 7. CORRELATION ANALYSIS
# ============================================================================

# Correlations involving Teq
cor_data <- unique_planets %>%
  select(all_of(c(target_col, pred_cols))) %>%
  select(where(is.numeric))

# Pairwise complete correlations (all unique planets; NA target handled pairwise)
cor_data <- cor_data[, vapply(cor_data, function(x) is.numeric(x) && sd(x, na.rm = TRUE) > 0, logical(1)), drop = FALSE]
cor_mat <- cor(cor_data, use = "pairwise.complete.obs")
cor_mat[!is.finite(cor_mat)] <- 0
diag(cor_mat) <- 1

# Target correlations (sorted)
target_cor <- sort(cor_mat[target_col, setdiff(colnames(cor_mat), target_col)],
                   decreasing = TRUE)

cat("========== CORRELATION WITH pl_eqt (unique planets) ==========\n")
print(round(target_cor, 3))
cat("\n")

target_cor_df <- tibble(
  predictor = names(target_cor),
  correlation = as.numeric(target_cor)
) %>%
  arrange(desc(abs(correlation)))

write_csv(target_cor_df, file.path(plot_dir, "target_correlations.csv"))

p_tcor <- target_cor_df %>%
  mutate(predictor = fct_reorder(predictor, correlation)) %>%
  ggplot(aes(x = predictor, y = correlation, fill = correlation > 0)) +
  geom_col(width = 0.7) +
  coord_flip() +
  scale_fill_manual(values = c("TRUE" = "#2C7BB6", "FALSE" = "#D7191C"),
                    guide = "none") +
  geom_hline(yintercept = 0, linewidth = 0.3) +
  labs(
    title = "Pearson Correlation with Equilibrium Temperature (pl_eqt)",
    subtitle = "Unique planets; pairwise complete observations",
    x = NULL, y = "Pearson r"
  )
save_plot(p_tcor, "19_target_correlations.png", width = 9, height = 10)

# Full correlation heatmap (png via corrplot device)
png(file.path(plot_dir, "20_correlation_heatmap.png"),
    width = 1400, height = 1400, res = 150)
corrplot(
  cor_mat,
  method = "color",
  type = "upper",
  order = "hclust",
  tl.cex = 0.55,
  tl.col = "black",
  col = colorRampPalette(c("#D73027", "#FFFFBF", "#4575B4"))(200),
  title = "Feature Correlation Matrix (unique planets)",
  mar = c(0, 0, 2, 0)
)
dev.off()
message("Saved: ", file.path(plot_dir, "20_correlation_heatmap.png"))

# Focused correlation: physically relevant subset
focus_vars <- c(target_col, planet_phys, stellar, "sy_vmag", "sy_pnum")
cor_focus <- cor(unique_planets %>% select(all_of(focus_vars)),
                 use = "pairwise.complete.obs")
cor_focus[!is.finite(cor_focus)] <- 0
diag(cor_focus) <- 1

png(file.path(plot_dir, "21_correlation_heatmap_focus.png"),
    width = 1000, height = 1000, res = 150)
corrplot(
  cor_focus,
  method = "number",
  type = "upper",
  order = "hclust",
  tl.cex = 0.8,
  tl.col = "black",
  number.cex = 0.7,
  col = colorRampPalette(c("#D73027", "#FFFFBF", "#4575B4"))(200),
  title = "Focused Correlations (planet / stellar / target)",
  mar = c(0, 0, 2, 0)
)
dev.off()
message("Saved: ", file.path(plot_dir, "21_correlation_heatmap_focus.png"))

# ============================================================================
# 8. SKY DISTRIBUTION
# ============================================================================

p_sky <- unique_planets %>%
  mutate(teq_status = factor(teq_status,
                              levels = c("Observed Teq", "Missing Teq"))) %>%
  ggplot(aes(x = ra, y = dec, color = teq_status)) +
  geom_point(alpha = 0.55, size = 1.5) +
  scale_color_manual(values = c("Observed Teq" = "#2C7BB6",
                                "Missing Teq" = "#F46D43")) +
  scale_x_reverse() +
  labs(
    title = "Sky Distribution of Planets (Equatorial)",
    subtitle = paste0("Unique planets, n = ", n_unique),
    x = "Right Ascension (deg)", y = "Declination (deg)", color = NULL
  )
save_plot(p_sky, "22_sky_equatorial.png", width = 10, height = 6)

p_gal <- unique_planets %>%
  ggplot(aes(x = glon, y = glat, color = teq_status)) +
  geom_point(alpha = 0.55, size = 1.5) +
  scale_color_manual(values = c("Observed Teq" = "#2C7BB6",
                                "Missing Teq" = "#F46D43")) +
  labs(
    title = "Sky Distribution of Planets (Galactic)",
    subtitle = paste0("Unique planets, n = ", n_unique),
    x = "Galactic longitude (deg)", y = "Galactic latitude (deg)", color = NULL
  )
save_plot(p_gal, "23_sky_galactic.png", width = 10, height = 6)

# ============================================================================
# 9. PAIRWISE RELATIONSHIPS (key features)
# ============================================================================

pair_vars <- c(target_col, "pl_rade", "pl_orbsmax", "pl_orbper", "st_teff", "st_mass")
pair_data <- unique_planets %>%
  filter(!is.na(.data[[target_col]])) %>%
  select(all_of(pair_vars)) %>%
  mutate(
    log_teq = log10(pl_eqt),
    log_rade = log10(pl_rade),
    log_orbper = log10(pl_orbper),
    log_orbsmax = log10(pmax(pl_orbsmax, 1e-8))
  ) %>%
  select(log_teq, log_rade, log_orbsmax, log_orbper, st_teff, st_mass)

# Use a sample if too large for ggpairs speed (here n is fine)
p_pairs <- ggpairs(
  pair_data,
  upper = list(continuous = wrap("cor", size = 3)),
  lower = list(continuous = wrap("points", alpha = 0.25, size = 0.6)),
  diag  = list(continuous = wrap("densityDiag", alpha = 0.5)),
  title = "Pairwise Relationships (unique planets with observed Teq)"
) + theme_bw(base_size = 9)

ggsave(file.path(plot_dir, "24_pairs_key_features.png"),
       p_pairs, width = 12, height = 12, dpi = 250, bg = "white")
message("Saved: ", file.path(plot_dir, "24_pairs_key_features.png"))

# ============================================================================
# 10. OUTLIER SCREENING (IQR rule on unique planets)
# ============================================================================

iqr_outliers <- function(x) {
  q <- quantile(x, c(0.25, 0.75), na.rm = TRUE)
  iqr <- q[2] - q[1]
  lo <- q[1] - 1.5 * iqr
  hi <- q[2] + 1.5 * iqr
  sum(x < lo | x > hi, na.rm = TRUE)
}

outlier_tbl <- tibble(
  variable = c(target_col, pred_cols),
  n_iqr_outliers = map_dbl(
    c(target_col, pred_cols),
    ~ iqr_outliers(unique_planets[[.x]])
  ),
  pct_iqr_outliers = 100 * n_iqr_outliers / n_unique
) %>%
  arrange(desc(n_iqr_outliers))

cat("========== IQR OUTLIERS (unique planets) ==========\n")
print(as.data.frame(outlier_tbl), row.names = FALSE)
write_csv(outlier_tbl, file.path(plot_dir, "outlier_iqr_summary.csv"))
cat("\n")

p_out <- outlier_tbl %>%
  filter(n_iqr_outliers > 0) %>%
  mutate(variable = fct_reorder(variable, n_iqr_outliers)) %>%
  ggplot(aes(x = variable, y = n_iqr_outliers)) +
  geom_col(fill = "#F46D43", width = 0.7) +
  coord_flip() +
  labs(
    title = "IQR Outlier Counts by Variable",
    subtitle = paste0("Unique planets, n = ", n_unique,
                      " | Beyond Q1−1.5·IQR or Q3+1.5·IQR"),
    x = NULL, y = "Number of outlier planets"
  )
save_plot(p_out, "25_outlier_counts.png", width = 9, height = 9)

# ============================================================================
# 11. CONSOLE SUMMARY REPORT
# ============================================================================

cat("\n============================================================\n")
cat("DESCRIPTIVE ANALYSIS COMPLETE\n")
cat("============================================================\n")
cat("Rows (measurements):     ", nrow(raw), "\n")
cat("Unique planets:          ", n_unique, "\n")
cat("Duplicated planets:      ", sum(meas_counts$n_measurements > 1), "\n")
cat("Target observed/missing: ", n_target_obs, " / ", n_target_miss,
    " (unique planets)\n", sep = "")
cat("Predictors:              ", length(pred_cols), " (all complete)\n")
cat("Strongest |r| with Teq: ",
    target_cor_df$predictor[1], " (r = ",
    round(target_cor_df$correlation[1], 3), ")\n", sep = "")
cat("Plots & tables saved to: ", normalizePath(plot_dir), "\n")
cat("============================================================\n")

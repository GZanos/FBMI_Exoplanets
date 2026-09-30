required_pkgs <- c("readxl", "writexl", "dplyr", "tidyr", "tibble")
install_if_missing <- function(pkgs) {
  missing <- pkgs[!pkgs %in% rownames(installed.packages())]
  if (length(missing) > 0) {
    options(repos = c(CRAN = "https://cloud.r-project.org"))
    install.packages(missing, dependencies = TRUE)
  }
}
install_if_missing(required_pkgs)

library(readxl)
library(writexl)
library(dplyr)
library(tidyr)
library(tibble)

data_dir <- "/Users/georgecharizanos/Desktop/Current Papers/Exoplanets/03.09.2026"
work_dir <- "/Users/georgecharizanos/Desktop/Current Papers/Exoplanets/10.09.2026"
results_xlsx <- file.path(data_dir, "Results_Final_V1.xlsx")
applications <- c("pl_bmasse", "pl_eqt")

# Robust RMSE: squared errors, drop the lowest and highest 10%, RMSE on the rest.
# The ASOC 115833 paper used a 5% trim; 10% is used here for stronger tail protection
# on the wide-range exoplanet targets.
# Charizanos, G., Demirhan, H., & İçen, D. (2026). A new fuzzy-Bayesian
# multiple imputation approach for missing data. Applied Soft Computing, 202,
# 115833. doi:10.1016/j.asoc.2026.115833
robust_rmse_trim <- 0.10

na_metrics <- function() {
  tibble(
    n = 0L,
    MAE = NA_real_,
    MdAE = NA_real_,
    RMSE = NA_real_,
    Robust_RMSE = NA_real_,
    RMSLRE = NA_real_,
    NRMSE = NA_real_
  )
}

# Median Absolute Error (MdAE): median(|e|).
# Hyndman, R. J., & Koehler, A. B. (2006). Another look at measures of forecast
# accuracy. International Journal of Forecasting, 22(4), 679-688.
# doi:10.1016/j.ijforecast.2006.03.001
#
# RMSLRE (Root Mean Squared Logarithmic Ratio Error):
#   RMSLRE = sqrt( mean( [ln(y_i / yhat_i)]^2 ) )
# Tasker, E. J., Laneuville, M., & Guttenberg, N. (2020). Estimating planetary
# mass with deep learning. The Astronomical Journal, 159, 41.
# Lalande, F., Tasker, E., & Doya, K. (2024). Estimating exoplanet mass using
# machine learning on incomplete datasets. arXiv:2410.06922.
# Requires strictly positive observed and imputed values.
#
# NRMSE (Normalised RMSE): RMSE / s_y, where s_y is the standard deviation of
# the target (actual) values in the evaluation set. The same s_y is used for
# every method in that set.
# Stekhoven, D. J., & Bühlmann, P. (2012). MissForest: Non-parametric missing
# value imputation for mixed-type data. Bioinformatics, 28(1), 112-118.
# Keerin, P., & Boongoen, T. (2022). Estimation of missing values in astronomical
# survey data. Information Processing & Management, 59, 102881.
compute_metrics <- function(pred, actual, trim = robust_rmse_trim, sd_actual = NULL) {
  pred   <- as.numeric(pred)
  actual <- as.numeric(actual)
  ok <- is.finite(pred) & is.finite(actual)
  n  <- sum(ok)

  if (n == 0L) {
    return(na_metrics())
  }

  e  <- pred[ok] - actual[ok]
  se <- e^2
  rmse <- sqrt(mean(se))

  if (is.null(sd_actual)) {
    sd_actual <- stats::sd(actual[ok])
  }
  nrmse <- NA_real_
  if (is.finite(sd_actual) && sd_actual > 0) {
    nrmse <- rmse / sd_actual
  }

  robust_rmse <- NA_real_
  if (n >= 3L) {
    q_lo <- stats::quantile(se, probs = trim)
    q_hi <- stats::quantile(se, probs = 1 - trim)
    keep <- (se >= q_lo) & (se <= q_hi)
    if (sum(keep) >= 1L) {
      robust_rmse <- sqrt(mean(se[keep]))
    }
  }

  # RMSLRE uses only strictly positive pairs (log-ratio is otherwise undefined).
  pos <- ok & (pred > 0) & (actual > 0)
  rmslre <- NA_real_
  if (sum(pos) > 0L) {
    log_ratio <- log(actual[pos] / pred[pos])
    rmslre <- sqrt(mean(log_ratio^2))
  }

  tibble(
    n = n,
    MAE = mean(abs(e)),
    MdAE = stats::median(abs(e)),
    RMSE = rmse,
    Robust_RMSE = robust_rmse,
    RMSLRE = rmslre,
    NRMSE = nrmse
  )
}

compute_from_sheet <- function(
    xlsx_path,
    sheet,
    actual_col = "actuals",
    id_col = "pl_name"
) {
  pred_df <- read_excel(xlsx_path, sheet = sheet)

  if (!actual_col %in% names(pred_df)) {
    stop("Column '", actual_col, "' not found in sheet '", sheet, "'.")
  }

  method_cols <- setdiff(names(pred_df), c(id_col, actual_col))
  actual <- as.numeric(pred_df[[actual_col]])
  # Shared NRMSE denominator: SD of all actuals, same for every method.
  sy <- stats::sd(actual[is.finite(actual)])

  bind_rows(lapply(method_cols, function(method) {
    compute_metrics(pred_df[[method]], actual, sd_actual = sy) %>%
      mutate(method = method, .before = 1)
  }))
}

metric_cols <- c("MAE", "MdAE", "RMSE", "Robust_RMSE", "RMSLRE", "NRMSE")
pct_impr <- function(fbmi, other) (other - fbmi) / other * 100

all_computed <- list()
vs_avg <- list()
vs_2nd <- list()

for (label in applications) {
  cat("\n============================================================\n")
  cat(label, ":", basename(results_xlsx), "\n")
  cat("============================================================\n")

  computed <- compute_from_sheet(results_xlsx, sheet = label)
  all_computed[[label]] <- computed %>% mutate(application = label, .before = 1)

  print(as.data.frame(computed), row.names = FALSE, digits = 6)

  fbmi_row <- computed[computed$method == "fbmi", , drop = FALSE]
  others <- computed[computed$method != "fbmi", , drop = FALSE]

  if (nrow(fbmi_row) != 1L) {
    warning(
      "Expected exactly one 'fbmi' row in ", label,
      "; found ", nrow(fbmi_row), "."
    )
    next
  }

  for (m in metric_cols) {
    fbmi_v <- fbmi_row[[m]]
    avg_other <- mean(others[[m]], na.rm = TRUE)
    j <- which.min(others[[m]])

    vs_avg[[length(vs_avg) + 1L]] <- tibble(
      application = label,
      metric = m,
      FBMI = fbmi_v,
      others_average = avg_other,
      pct_improvement = pct_impr(fbmi_v, avg_other)
    )
    vs_2nd[[length(vs_2nd) + 1L]] <- tibble(
      application = label,
      metric = m,
      FBMI = fbmi_v,
      second_best_method = others$method[j],
      second_best = others[[m]][j],
      pct_improvement = pct_impr(fbmi_v, others[[m]][j])
    )
  }
}

vs_avg_df <- as.data.frame(bind_rows(vs_avg)) %>% arrange(application, metric)
vs_2nd_df <- as.data.frame(bind_rows(vs_2nd)) %>% arrange(application, metric)
metrics_df <- as.data.frame(bind_rows(all_computed)) %>% arrange(application, method)

cat("\n============================================================\n")
cat("FBMI % improvement vs average of all other methods\n")
cat("============================================================\n")
print(vs_avg_df, row.names = FALSE)

cat("\n============================================================\n")
cat("FBMI % improvement vs 2nd-best method (best non-FBMI, per metric)\n")
cat("============================================================\n")
print(vs_2nd_df, row.names = FALSE)

out_xlsx <- file.path(work_dir, "benchmark_metrics_v3.xlsx")
write_xlsx(
  list(
    metrics = metrics_df,
    vs_average = vs_avg_df,
    vs_second_best = vs_2nd_df
  ),
  path = out_xlsx
)

cat("\nSaved Excel:\n", out_xlsx, "\n", sep = "")
cat("\nDone.\n")

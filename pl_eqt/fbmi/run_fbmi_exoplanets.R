# ==============================================================================
# FBMI (Fuzzy Bayesian Multiple Imputation) — optimized exoplanet runner
# Source: https://github.com/GZanos/Fuzzy_Multiple_Imputation
# ==============================================================================

args <- commandArgs(trailingOnly = TRUE)
work_dir <- if (length(args) >= 1) args[[1]] else getwd()
setwd(work_dir)
feature_list_file <- if (length(args) >= 2 && nzchar(args[[2]])) args[[2]] else NA_character_
run_tag <- if (length(args) >= 3 && nzchar(args[[3]])) args[[3]] else "default"

util_candidates <- c(
  file.path(work_dir, "DBDA2E-utilities.R"),
  file.path(dirname(work_dir), "10pc", "DBDA2E-utilities.R"),
  file.path(dirname(work_dir), "30pc", "DBDA2E-utilities.R")
)
util_path <- util_candidates[file.exists(util_candidates)][1]
if (!is.na(util_path)) suppressWarnings(try(source(util_path), silent = TRUE))

req <- c("dplyr", "rjags", "runjags", "coda")
to_install <- setdiff(req, rownames(installed.packages()))
if (length(to_install)) {
  options(repos = c(CRAN = "https://cloud.r-project.org"))
  install.packages(to_install, dependencies = TRUE)
}
library(dplyr)
library(rjags)
library(runjags)
library(coda)

# ---- User-specified FBLIR Hyperparameter Configuration ----
JAGS_N_CHAINS <- 4
JAGS_ADAPT <- 1500
JAGS_BURNIN <- 2000
JAGS_THIN <- 15
JAGS_SAMPLE <- 10000

M_VALUES <- seq(-4, 4, by = 0.2)
SYMMETRY_THRESHOLDS <- c(0.1, 0.5, 1, 1.5, 2.5)
K_VALUES <- seq(-4, 4, by = 0.2)
UNCERTAINTY_WEIGHTS <- c(0, 0.25, 0.75, 1)
FUZZIFY_SCALES <- c(0.01, 0.05, 0.1, 0.2, 0.4)

TUNE_HOLDOUT_FRAC <- 0.20

id_col <- "pl_name"
target_col <- "pl_eqt"
set.seed(123)

train <- read.csv("train_fbmi.csv", stringsAsFactors = FALSE)
train[[target_col]] <- suppressWarnings(as.numeric(train[[target_col]]))
all_feature_cols <- setdiff(colnames(train), c(id_col, target_col))
for (fc in all_feature_cols) train[[fc]] <- suppressWarnings(as.numeric(train[[fc]]))

if (!is.na(feature_list_file) && file.exists(feature_list_file)) {
  feature_cols <- intersect(readLines(feature_list_file), all_feature_cols)
  feature_cols <- feature_cols[nzchar(feature_cols)]
  if (length(feature_cols) < 2) stop("Feature list too small")
} else {
  feature_cols <- all_feature_cols
}

obs_idx <- which(!is.na(train[[target_col]]))
miss_idx <- which(is.na(train[[target_col]]))
cat("FBMI tag=", run_tag, "| P=", length(feature_cols),
    "| obs=", length(obs_idx), "| miss=", length(miss_idx), "\n")
cat("Features:", paste(feature_cols, collapse = ", "), "\n")

GFN.multi_vec <- function(beta_gfn, X_mean, X_var, symmetry.threshold = 4) {
  b_m <- beta_gfn[1]; b_v <- beta_gfn[2]
  mean <- b_m * X_mean
  variance <- (X_var * b_m^2) + (b_v * X_mean^2) + (X_var * b_v)
  for (iter in 1:20) {
    ok <- is.finite(variance) & (variance > 0) & (mean != 0) &
      (abs(mean / sqrt(pmax(variance, 1e-300))) < symmetry.threshold)
    if (!any(ok)) break
    variance[ok] <- variance[ok] * 0.1
  }
  cbind(Mean = mean, Variance = variance)
}

defuzzify_vec <- function(means, vars, k, m, symmetry.threshold) {
  out <- means
  valid <- is.finite(means) & is.finite(vars) & (vars > 0)
  if (!any(valid)) return(out)
  delta <- abs(means[valid] / sqrt(vars[valid]))
  adj <- m / (1 + exp(-k * (delta - symmetry.threshold)))
  use_adj <- delta <= symmetry.threshold
  tmp <- means[valid]
  tmp[use_adj] <- means[valid][use_adj] + adj[use_adj] * vars[valid][use_adj]
  out[valid] <- tmp
  out
}

gfn_predict_matrix <- function(Xz, beta0_gfn, beta_gfn_w, sigma_gfn, fuzz_var, sym_thr) {
  n <- nrow(Xz); p <- ncol(Xz)
  Y_mean <- rep(beta0_gfn[1], n)
  Y_var <- rep(beta0_gfn[2], n)
  for (j in seq_len(p)) {
    prod <- GFN.multi_vec(beta_gfn_w[j, ], Xz[, j], rep(fuzz_var, n), symmetry.threshold = sym_thr)
    Y_mean <- Y_mean + prod[, 1]
    Y_var <- Y_var + prod[, 2]
  }
  cbind(Mean = Y_mean + sigma_gfn[1], Variance = Y_var + sigma_gfn[2])
}

X_full_df <- train[, feature_cols, drop = FALSE]
Y_full <- train[, target_col]

X_obs <- as.matrix(X_full_df[obs_idx, , drop = FALSE])
x_means <- colMeans(X_obs, na.rm = TRUE)
x_sds <- apply(X_obs, 2, sd, na.rm = TRUE)
x_sds[!is.finite(x_sds) | x_sds == 0] <- 1e-6
X_obs_z <- scale(X_obs, center = x_means, scale = x_sds)
Y_obs_raw <- Y_full[obs_idx]
y_mean <- mean(Y_obs_raw)
y_sd <- sd(Y_obs_raw)
if (!is.finite(y_sd) || y_sd == 0) y_sd <- 1e-6
Y_obs_z <- (Y_obs_raw - y_mean) / y_sd
P <- ncol(X_obs_z)

modelString <- "
model {
  for (i in 1:N_obs) {
    y[i] ~ dnorm(mu[i], tau)
    mu[i] <- beta0 + inprod(beta[1:P], X[i,1:P])
  }
  beta0 ~ dnorm(0, tau_beta0)
  tau_beta0 ~ dgamma(1, 1)
  for (j in 1:P) {
    beta[j] ~ dnorm(0, tau_beta[j])
    tau_beta[j] ~ dgamma(alpha_beta, beta_beta)
  }
  alpha_beta ~ dgamma(1, 1)
  beta_beta ~ dgamma(1, 1)
  tau ~ dgamma(1, 1)
  sigma <- 1 / sqrt(tau)
  sigma2_resid <- 1/tau
}
"
model_file <- paste0("exoplanet_FBMI_model_", run_tag, ".txt")
writeLines(modelString, con = model_file)

cat("JAGS (", run_tag, "): N=", nrow(X_obs_z), " P=", P,
    " chains=", JAGS_N_CHAINS, " adapt=", JAGS_ADAPT,
    " burnin=", JAGS_BURNIN, " thin=", JAGS_THIN, " sample=", JAGS_SAMPLE, "\n")
t0 <- proc.time()
runJagsOut <- run.jags(
  method = "parallel",
  model = model_file,
  data = list(X = as.matrix(X_obs_z), y = as.numeric(Y_obs_z), N_obs = nrow(X_obs_z), P = P),
  n.chains = JAGS_N_CHAINS,
  adapt = JAGS_ADAPT,
  burnin = JAGS_BURNIN,
  thin = JAGS_THIN,
  sample = JAGS_SAMPLE,
  monitor = c("beta0", "beta", "sigma", "tau_beta0", "tau_beta", "sigma2_resid")
)
mcmc_mat <- as.matrix(as.mcmc.list(runJagsOut))
cat("JAGS done in", round((proc.time() - t0)[3] / 60, 2), "min | samples=", nrow(mcmc_mat), "\n")

beta0_gfn <- c(mean(mcmc_mat[, "beta0"]), var(mcmc_mat[, "beta0"]))
beta_gfn_matrix <- matrix(NA_real_, nrow = P, ncol = 2,
                          dimnames = list(colnames(X_obs_z), c("Mean", "Variance")))
for (j in seq_len(P)) {
  beta_samples <- mcmc_mat[, paste0("beta[", j, "]")]
  tau_beta_samples <- mcmc_mat[, paste0("tau_beta[", j, "]")]
  beta_gfn_matrix[j, ] <- c(mean(beta_samples), var(beta_samples) + 0.5 * mean(1 / tau_beta_samples))
}
sigma_gfn <- c(0, mean(mcmc_mat[, "sigma2_resid"]))

X_full_z <- sweep(as.matrix(X_full_df), 2, x_means, FUN = "-")
X_full_z <- sweep(X_full_z, 2, x_sds, FUN = "/")
X_full_z[is.na(X_full_z)] <- 0

# Holdout indices among observed rows (positions in obs_idx)
n_tune <- max(30, floor(TUNE_HOLDOUT_FRAC * length(obs_idx)))
tune_pos <- sample.int(length(obs_idx), n_tune)
tune_rows <- obs_idx[tune_pos]
X_tune_z <- X_full_z[tune_rows, , drop = FALSE]
Y_tune_raw <- Y_full[tune_rows]

cache_keys <- expand.grid(
  uncertainty_weight = UNCERTAINTY_WEIGHTS,
  fuzzify_variance = FUZZIFY_SCALES,
  symmetry.threshold = SYMMETRY_THRESHOLDS
)
cache <- new.env(parent = emptyenv())
cat("Caching", nrow(cache_keys), "GFN configs...\n")
for (i in seq_len(nrow(cache_keys))) {
  uw <- cache_keys$uncertainty_weight[i]
  fv <- cache_keys$fuzzify_variance[i]
  st <- cache_keys$symmetry.threshold[i]
  beta_w <- beta_gfn_matrix
  beta_w[, 2] <- beta_gfn_matrix[, 2] * uw + (1 - uw) * mean(beta_gfn_matrix[, 2])
  key <- paste(c(uw, fv, st), collapse = "|")
  cache[[key]] <- gfn_predict_matrix(X_tune_z, beta0_gfn, beta_w, sigma_gfn, fv, st)
}

fblr_grid <- expand.grid(
  m = M_VALUES,
  symmetry.threshold = SYMMETRY_THRESHOLDS,
  k = K_VALUES,
  uncertainty_weight = UNCERTAINTY_WEIGHTS,
  fuzzify_variance = FUZZIFY_SCALES
)
cat("Scanning", nrow(fblr_grid), "defuzzify combinations...\n")
best_mae <- Inf
best_params <- NULL
for (i in seq_len(nrow(fblr_grid))) {
  params <- fblr_grid[i, ]
  key <- paste(c(params$uncertainty_weight, params$fuzzify_variance, params$symmetry.threshold), collapse = "|")
  Y_est <- cache[[key]]
  pred_z <- defuzzify_vec(Y_est[, 1], Y_est[, 2], params$k, params$m, params$symmetry.threshold)
  pred <- pred_z * y_sd + y_mean
  mae <- mean(abs(pred - Y_tune_raw), na.rm = TRUE)
  if (is.finite(mae) && mae < best_mae) {
    best_mae <- mae
    best_params <- params
  }
}
cat("Best holdout MAE (log10):", best_mae, "\n")
print(best_params)

beta_w_best <- beta_gfn_matrix
beta_w_best[, 2] <- beta_gfn_matrix[, 2] * best_params$uncertainty_weight +
  (1 - best_params$uncertainty_weight) * mean(beta_gfn_matrix[, 2])
Y_miss <- gfn_predict_matrix(
  X_full_z[miss_idx, , drop = FALSE],
  beta0_gfn, beta_w_best, sigma_gfn,
  best_params$fuzzify_variance, best_params$symmetry.threshold
)
pred_z <- defuzzify_vec(Y_miss[, 1], Y_miss[, 2], best_params$k, best_params$m, best_params$symmetry.threshold)
pred_log <- as.numeric(pred_z) * y_sd + y_mean

out <- data.frame(
  pl_name = train[[id_col]][miss_idx],
  predicted_mass_log10 = pred_log,
  stringsAsFactors = FALSE
)
write.csv(out, paste0("fbmi_predictions_log_", run_tag, ".csv"), row.names = FALSE)
write.csv(
  cbind(best_params, holdout_MAE_log10 = best_mae, n_features = length(feature_cols),
        features = paste(feature_cols, collapse = ";")),
  paste0("fbmi_best_params_", run_tag, ".csv"),
  row.names = FALSE
)
if (run_tag == "default") {
  write.csv(out, "fbmi_predictions_log.csv", row.names = FALSE)
  write.csv(best_params, "fbmi_best_params.csv", row.names = FALSE)
}
cat("Wrote predictions for", run_tag, "(", nrow(out), "rows)\n")
cat("=== FBMI COMPLETE:", run_tag, "===\n")

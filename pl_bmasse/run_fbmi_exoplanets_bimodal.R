# ==============================================================================
# FBMI (Fuzzy Bayesian Multiple Imputation) — optimized exoplanet runner
# Source: https://github.com/GZanos/Fuzzy_Multiple_Imputation
# Args: work_dir, optional feature_list_file, optional run_tag
# ==============================================================================

args <- commandArgs(trailingOnly = TRUE)
work_dir <- if (length(args) >= 1) args[[1]] else getwd()
setwd(work_dir)

setwd("~/Documents/makaleler/George/paper11/fbmi")
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
library(readxl)
source("~/Documents/makaleler/George/paper11/fbmi/DBDA2E-utilities.R")
# ---- User-specified FBLIR Hyperparameter Configuration ----


M_VALUES <- seq(-4, 4, by = 0.2)
SYMMETRY_THRESHOLDS <- c(0.1, 0.5, 1, 1.5, 2.5)
K_VALUES <- seq(-4, 4, by = 0.2)
UNCERTAINTY_WEIGHTS <- c(0, 0.25, 0.75, 1)
FUZZIFY_SCALES <- c(0.01, 0.05, 0.1, 0.2, 0.4) # HD: Modified

TUNE_HOLDOUT_FRAC <- 0.20

id_col <- "pl_name"
target_col <- "pl_bmasse"
set.seed(123)

train <- read.csv("~/Documents/makaleler/George/paper11/fbmi/_exoplanet_run/train_fbmi.csv", stringsAsFactors = FALSE)
train[[target_col]] <- suppressWarnings(as.numeric(train[[target_col]]))

# HD: Validation data is not bimodal?
validate <- read.csv("~/Documents/makaleler/George/paper11/fbmi/_exoplanet_run/validate_fbmi.csv")
hist(log10(validate[,'true_mass']))

# HD: I couldn't replicate Figure 3:
Y_full <- train[, target_col]
hist(Y_full)

full_mass_data <- cbind(Y_full, validate[,'true_mass'])
hist(full_mass_data)

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



# gfn_predict_matrix <- function(Xz, beta0_gfn, beta02_gfn, beta_gfn_w, sigma_gfn, fuzz_var, sym_thr, w1, delta_gfn) {
#   n <- nrow(Xz); p <- ncol(Xz)
#   Y_mean <- rep(beta0_gfn[1], n)
#   Y_var <- rep(beta0_gfn[2], n)
#   for (j in seq_len(p)) {
#     prod <- GFN.multi_vec(beta_gfn_w[j, ], Xz[, j], rep(fuzz_var, n), symmetry.threshold = sym_thr)
#     Y_mean <- Y_mean + prod[, 1]
#     Y_var <- Y_var + prod[, 2]
#   }
#   cbind(Mean = Y_mean + sigma_gfn[1], Variance = Y_var + sigma_gfn[2])
# }

X_full_df <- train[, feature_cols, drop = FALSE]
Y_full <- train[, target_col]

hist(Y_full)

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

hist(Y_obs_raw)
summary(Y_obs_z)
hist(Y_obs_z)

#HD: Assign observed samples to two modes.
trainClass <- array(NA, length(Y_obs_raw))
trainClass[(Y_obs_z < 0)] <- 1
trainClass[!(Y_obs_z < 0)] <- 2 

# classThreshold <- 0.27
# z_ind <- ifelse(Y_obs_z < classThreshold, 1, 2) # HD: Fixed classification based on the the data. Here 2 is from the histogram where the two distributions meet.

modelString <- "
model {
  for (i in 1:N_obs) {
    z[i] ~ dcat(catprob[1:2])
    y[i] ~ dnorm(mu[i,z[i]], tau[z[i]])
    mu[i,1] <- beta01 + inprod(beta1[1:P], X[i,1:P])
    mu[i,2] <- beta02 + inprod(beta1[1:P], X[i,1:P])#mu[i,1] + delta #
    yrep[i] ~ dnorm(mu[i,z[i]], tau[z[i]])
  }
  catprob[1:2] ~ ddirch(alpha[1:2])
  alpha[1] ~ dunif(0,1) #dnorm(0, 1)T(0,) #<- 5
  alpha[2] ~ dunif(0,1) #dnorm(0, 1)T(0,) #<- 1

  eta_delta ~ dnorm(0, tau_delta)
  tau_delta ~ dgamma(1, 0.01)
  delta <- exp(eta_delta)

  beta01 ~ dnorm(0, tau_beta01)
  tau_beta01 ~ dgamma(1, 0.01)
  beta02 ~ dnorm(0, tau_beta02)
  tau_beta02 ~ dgamma(1, 0.01)
  for (j in 1:P) {
    beta1[j] ~ dnorm(0, tau_beta1[j])
    tau_beta1[j] ~ dgamma(1, 0.01) 
  }
  for (k in 1:2) {
    tau[k] ~ dgamma(1, 0.01)
    sigma[k] <- 1 / sqrt(tau[k])
  }
  sigma2_resid <- 1/(tau[1]+tau[2])
}
"

# modelString <- "
# model {
#   for (i in 1:N_obs) {
#     y[i] ~ dnorm(mu[i,trainClass[i]], tau[trainClass[i]])
#     mu[i,1] <- beta01 + inprod(beta1[1:P], X[i,1:P])
#     mu[i,2] <- mu[i,1] + delta #beta02 + inprod(beta2[1:P], X[i,1:P])#
# 
#     yrep[i] ~ dnorm(mu[i,trainClass[i]], tau[trainClass[i]])
#   }
# 
#   eta_delta ~ dnorm(0, tau_delta)
#   tau_delta ~ dgamma(1, 0.01)
#   delta <- exp(eta_delta)
# 
#   beta01 ~ dnorm(0, tau_beta01)
#   tau_beta01 ~ dgamma(1, 0.01)
#   beta02 ~ dnorm(0, tau_beta02)
#   tau_beta02 ~ dgamma(1, 0.01)
#   for (j in 1:P) {
#     beta1[j] ~ dnorm(0, tau_beta1[j])
#     tau_beta1[j] ~ dgamma(1, 0.01)
#     beta2[j] ~ dnorm(0, tau_beta2[j])
#     tau_beta2[j] ~ dgamma(1, 0.01) 
#   }
# 
#   for (k in 1:2) {
#     tau[k] ~ dgamma(1, 0.01)
#     sigma[k] <- 1 / sqrt(tau[k])
#   }
#   sigma2_resid <- 1/(tau[1]+tau[2])
# }
# "
model_file <- paste0("exoplanet_FBMI_model_", run_tag, ".txt")
writeLines(modelString, con = model_file)
  
JAGS_N_CHAINS <- 2
JAGS_ADAPT <- 100
JAGS_BURNIN <- 2000
JAGS_THIN <- 3
JAGS_SAMPLE <- 1000

cat("JAGS (", run_tag, "): N=", nrow(X_obs_z), " P=", P,
    " chains=", JAGS_N_CHAINS, " adapt=", JAGS_ADAPT,
    " burnin=", JAGS_BURNIN, " thin=", JAGS_THIN, " sample=", JAGS_SAMPLE, "\n")
t0 <- proc.time()
runJagsOut <- run.jags(
  method = "parallel",
  model = model_file,
  data = list(X = as.matrix(X_obs_z), y = as.numeric(Y_obs_z), N_obs = nrow(X_obs_z), 
              P = P, trainClass = trainClass),
  n.chains = JAGS_N_CHAINS,
  adapt = JAGS_ADAPT,
  burnin = JAGS_BURNIN,
  thin = JAGS_THIN,
  sample = JAGS_SAMPLE,
  monitor = c("beta01","beta02", "beta1", "sigma", "tau_beta01", "tau_beta02", "tau_beta1", "sigma2_resid", 
              "alpha", "delta", "tau_delta", "yrep")
)
mcmc_mat <- as.matrix(as.mcmc.list(runJagsOut))
cat("JAGS done in", round((proc.time() - t0)[3] / 60, 2), "min | samples=", nrow(mcmc_mat), "\n")

# summary(runJagsOut)

#HD: To check if we really fit a bimodal distribution
yFitted <- array(NA, nrow(X_obs_z))
for (i in 1:nrow(X_obs_z)){
  yFitted[i] <- mean(mcmc_mat[, paste0("yrep[", i ,"]")])
}
hist(Y_obs_z, freq = FALSE,ylim = c(0,1))
lines(density(yFitted),col = "red", lwd = 2)
#HD: To check if we really fit a bimodal distribution

# diagMCMC( codaObject=as.mcmc.list(runJagsOut) , parName="yrep" )
# alpha1 <- mean(mcmc_mat[, "alpha[1]"])
# alpha2 <- mean(mcmc_mat[, "alpha[2]"])
# alpha0 <-  alpha1 + alpha2
# 
# w1 <- alpha1/alpha0
# w2 <- alpha2/alpha0
# 
# delta <- mean(mcmc_mat[, "delta"])
# varDelta <-  var(mcmc_mat[, "delta"]) # 1/mean(mcmc_mat[, "tau_delta"]) #

delta_gfn <- c(0,0)#c(delta, varDelta) #
  
beta01_gfn <- c(mean(mcmc_mat[, "beta01"]), mean(1/mcmc_mat[, "tau_beta01"])) # HD: I changed variances here and use taubeta0
beta02_gfn <- c(mean(mcmc_mat[, "beta02"]), mean(1/mcmc_mat[, "tau_beta02"]))
beta1_gfn_matrix <- matrix(NA_real_, nrow = P, ncol = 2,
                          dimnames = list(colnames(X_obs_z), c("Mean", "Variance")))
# beta2_gfn_matrix <- matrix(NA_real_, nrow = P, ncol = 2,
#                            dimnames = list(colnames(X_obs_z), c("Mean", "Variance")))
for (j in seq_len(P)) {
  beta1_samples <- mcmc_mat[, paste0("beta1[", j, "]")]
  tau_beta1_samples <- mcmc_mat[, paste0("tau_beta1[", j, "]")]
  beta1_gfn_matrix[j, ] <- c(mean(beta1_samples), var(beta1_samples) + 0.5 * mean(1 / tau_beta1_samples))
  
  # beta2_samples <- mcmc_mat[, paste0("beta2[", j, "]")]
  # tau_beta2_samples <- mcmc_mat[, paste0("tau_beta2[", j, "]")]
  # beta2_gfn_matrix[j, ] <- c(mean(beta2_samples), var(beta2_samples) + 0.5 * mean(1 / tau_beta2_samples))
}
#HD: Why sigma_gfn has a mean of 0? Zero variance is not suitable.
sigma_gfn <- c(mean(mcmc_mat[, "sigma2_resid"]), (mean(mcmc_mat[, "sigma[1]"]) + mean(mcmc_mat[, "sigma[2]"]))/2) # (mean(mcmc_mat[, "sigma[1]"]) + mean(mcmc_mat[, "sigma[2]"]))) # sigma2resid was not out from jags 


X_full_z <- sweep(as.matrix(X_full_df), 2, x_means, FUN = "-")
X_full_z <- sweep(X_full_z, 2, x_sds, FUN = "/")
X_full_z[is.na(X_full_z)] <- 0

# Holdout indices among observed rows (positions in obs_idx)
n_tune <- max(30, floor(TUNE_HOLDOUT_FRAC * length(obs_idx)))
set.seed(98942) #243242
tune_pos <- sample.int(length(obs_idx), n_tune)
# write.csv(tune_pos, "tune_pos.csv")
# tune_pos <- as.vector(read.csv("tune_pos.csv")[-1])
# tune_pos <- tune_pos$x
tune_rows <- obs_idx[tune_pos]
X_tune_z <- X_full_z[tune_rows, , drop = FALSE]
Y_tune_raw <- Y_full[tune_rows] #HD: Model is fitted on scaled data but tuning is done on non-scaled data!
Y_tune_z <- scale(Y_tune_raw)

y_tune_mean <- mean(Y_tune_raw)
y_tune_var <- var(Y_tune_raw)

hist(Y_tune_z, freq = FALSE,ylim = c(0,1))
lines(density(yFitted),col = "red", lwd = 2)

cache_keys <- expand.grid(
  uncertainty_weight = UNCERTAINTY_WEIGHTS,
  fuzzify_variance = FUZZIFY_SCALES,
  symmetry.threshold = SYMMETRY_THRESHOLDS
)
cache <- new.env(parent = emptyenv())
cat("Caching", nrow(cache_keys), "GFN configs...\n")

kmeans_result <- kmeans(X_tune_z, centers = 2, nstart = 10)
clustersTune <- kmeans_result$cluster
# clustersTuneSave <- clustersTune
# clustersTune <- ifelse(clustersTune == 1, 2, ifelse(clustersTune == 2, 1, clustersTune))

# write.csv(clustersTune, "clustersTune.csv")
# clustersTune <- as.vector(read.csv("clustersTune.csv")[-1])
# clustersTune <- clustersTune$x
# library(mclust)
# gmm <- Mclust(X_tune_z, G = 2)
# clustersTune <- gmm$classification


# regData <- data.frame(y = as.numeric(trainClass)-1, x = as.matrix(X_obs_z))
# head(regData)
# model <- glm(y ~ ., data = regData, family = "binomial")
# 
# predData <- data.frame(x = as.matrix(X_tune_z))
# clustersTune <- array(NA, nrow(predData) )
# clustersTune[(predict(model, predData, type = "response") < 0.5)] <- 1
# clustersTune[!(predict(model, predData, type = "response") < 0.5)] <- 2

for (i in seq_len(nrow(cache_keys))) {
  uw <- cache_keys$uncertainty_weight[i]
  fv <- cache_keys$fuzzify_variance[i]
  st <- cache_keys$symmetry.threshold[i]
  beta1_w <- beta1_gfn_matrix
  beta1_w[, 2] <- beta1_gfn_matrix[, 2] * uw + (1 - uw) * mean(beta1_gfn_matrix[, 2])
  # beta2_w <- beta2_gfn_matrix
  # beta2_w[, 2] <- beta2_gfn_matrix[, 2] * uw + (1 - uw) * mean(beta2_gfn_matrix[, 2])
  key <- paste(c(uw, fv, st), collapse = "|")
  cache[[key]] <- gfn_predict_matrix(X_tune_z, beta01_gfn, beta02_gfn, beta1_w, beta1_w, sigma_gfn, fv, st, 
                                     w1=w1, delta_gfn=delta_gfn, clusters = clustersTune) # HD modified
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
  # pred <- pred_z * y_sd + y_mean
  pred <- pred_z * sqrt(y_tune_var) + y_tune_mean
  mae <- mean(abs(pred - Y_tune_raw), na.rm = TRUE) 
  if (is.finite(mae) && mae < best_mae) {
    best_mae <- mae
    best_params <- params
  }
}
cat("Best holdout MAE (log10):", best_mae, "\n")
print(best_params)

kmeans_result <- kmeans(X_full_z[miss_idx, , drop = FALSE], centers = 2, nstart = 10)
clusters <- kmeans_result$cluster
# clustersSave <- clusters
# clusters <- ifelse(clusters == 1, 2, ifelse(clusters == 2, 1, clusters))
# write.csv(clusters, "clusters.csv")
# clusters <- as.vector(read.csv("clusters.csv")[-1])
# clusters <- clusters$x

# dist_matrix <- dist(X_full_z[miss_idx, , drop = FALSE], method = "euclidean")
# hc <- hclust(dist_matrix, method = "ward.D2")
# clusters <- cutree(hc, k = 2)
# 
# library(mclust)
# gmm <- Mclust(X_full_z[miss_idx, , drop = FALSE], G = 2)
# clusters <- gmm$classification



#HD: Use the fitted model to the observed data to predict the class of the new data
# predData2 <- data.frame(x = as.matrix(X_full_z[miss_idx, , drop = FALSE]))
# predX_full_z_miss <- array(NA, nrow(predData2) )
# predX_full_z_miss[(predict(model, predData2, type = "response") < 0.5)] <- 1
# predX_full_z_miss[!(predict(model, predData2, type = "response") < 0.5)] <- 2


beta1_w_best <- beta1_gfn_matrix
beta1_w_best[, 2] <- beta1_gfn_matrix[, 2] * best_params$uncertainty_weight +
  (1 - best_params$uncertainty_weight) * mean(beta1_gfn_matrix[, 2])

# beta2_w_best <- beta2_gfn_matrix
# beta2_w_best[, 2] <- beta2_gfn_matrix[, 2] * best_params$uncertainty_weight +
#   (1 - best_params$uncertainty_weight) * mean(beta2_gfn_matrix[, 2])

Y_miss <- gfn_predict_matrix(X_full_z[miss_idx, , drop = FALSE],
                             beta01_gfn, beta02_gfn, beta1_w_best, beta1_w_best, sigma_gfn,
                             best_params$fuzzify_variance, best_params$symmetry.threshold,
                             w1=w1, delta_gfn=delta_gfn, clusters = clusters) # HD: Modified.
                              
pred_z <- defuzzify_vec(Y_miss[, 1], Y_miss[, 2], best_params$k, best_params$m, best_params$symmetry.threshold)
pred_log <- as.numeric(pred_z) * y_sd + y_mean

out <- data.frame(
  pl_name = train[[id_col]][miss_idx],
  predicted_mass_log10 = pred_log,
  stringsAsFactors = FALSE
)


out_raw <- 10^out[,'predicted_mass_log10']
pl_bmasse <- read_excel("~/Documents/makaleler/George/paper11/pl_bmasse_benchmark_summary.xlsx", sheet = "predictions")
pl_bmasse$fbmi <- out_raw

outSave <- out

library(openxlsx)
wb = loadWorkbook("~/Documents/makaleler/George/paper11/pl_bmasse_benchmark_summary.xlsx")
removeWorksheet(wb, sheet = "predictions")
saveWorkbook(wb, "~/Documents/makaleler/George/paper11/pl_bmasse_benchmark_summary.xlsx", overwrite = TRUE)
pl_bmasse <- as.data.frame(pl_bmasse)
addWorksheet(wb, sheetName = "predictions")
writeData(wb, sheet = "predictions", x = pl_bmasse)

saveWorkbook(wb, "~/Documents/makaleler/George/paper11/pl_bmasse_benchmark_summary.xlsx", overwrite = TRUE)

source("~/Documents/makaleler/George/paper11/calculate_benchmark_metrics.R")

cat("Best holdout MAE (log10):", best_mae, "\n")
print(best_params)



out <- outSave










# write.csv(out, paste0("fbmi_predictions_log_", run_tag, "HD4.csv"), row.names = FALSE)
# write.csv(
#   cbind(best_params, holdout_MAE_log10 = best_mae, n_features = length(feature_cols),
#         features = paste(feature_cols, collapse = ";")),
#   paste0("fbmi_best_params_", run_tag, "HD4.csv"),
#   row.names = FALSE
# )
# if (run_tag == "default") {
#   write.csv(out, "fbmi_predictions_logHD4.csv", row.names = FALSE)
#   write.csv(best_params, "fbmi_best_paramsHD4.csv", row.names = FALSE)
# }
# cat("Wrote predictions for", run_tag, "(", nrow(out), "rows)\n")
# cat("=== FBMI COMPLETE:", run_tag, "===\n")


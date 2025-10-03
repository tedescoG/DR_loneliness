# Cross-Fitted Bootstrap Doubly Robust ATT Estimation Analysis
# Comparison: Increase vs Decrease
# This script runs cf_DR_att() with 4 different outcome models
# Using stratified K-fold cross-fitting to reduce overfitting bias

# Setup
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
run = 123
source("utils.R")

# Load data
d = readRDS("data/data_iptw.rds")

# Define propensity score formula
ps_formula = as.formula(
  remote_contact ~
    female +
      age_cat +
      edu +
      emp_status +
      income +
      marital +
      coliving +
      change_res +
      kinless +
      health_pre +
      chronic +
      death_due_covid +
      ppl_infected +
      income_loss +
      job_loss +
      neighborhood +
      baseline_depr +
      baseline_lone
)

# Define 4 different outcome model formulas (RHS only, as character strings)
outcome_formula_1 = "1"

outcome_formula_2 = "baseline_lone"

outcome_formula_3 = "baseline_lone + baseline_depr"

outcome_formula_4 = "baseline_lone + baseline_depr + female + age_cat + edu +
                     emp_status + income + marital + coliving + change_res +
                     kinless + health_pre + chronic + death_due_covid +
                     ppl_infected + income_loss + job_loss + neighborhood"

# Tuned parameters for GBM - Increase vs Decrease
gbm_inc_dec = readRDS("results/weighting/inc_dec_params.rds")

# Store results in a list
results_cf = list()

# =============================================================================
# Model 1: Cross-fitted Bootstrap
# =============================================================================

cat("\n\n")
cat("=======================================================================\n")
cat("MODEL 1: Cross-Fitted Bootstrap (Intercept-only outcome model)\n")
cat(
  "=======================================================================\n\n"
)

results_cf$model1 = cf_DR_att(
  outcome = "severe_loneliness",
  treatment = "remote_contact",
  treated_level = "increase",
  control_level = "decrease",
  f.ps = ps_formula,
  f.out = outcome_formula_1,
  data = d,
  ps_params = gbm_inc_dec,
  k = 2,
  stratified = TRUE,
  n_boot = 50,
  seed = run,
  verbose = TRUE,
  parallel = F,
  n_cores = NULL,
  plot_diagnostics = TRUE
)

# =============================================================================
# Model 2: Cross-fitted Bootstrap
# =============================================================================

cat("\n\n")
cat("=======================================================================\n")
cat("MODEL 2: Cross-Fitted Bootstrap (Baseline loneliness)\n")
cat(
  "=======================================================================\n\n"
)

results_cf$model2 = cf_DR_att(
  outcome = "severe_loneliness",
  treatment = "remote_contact",
  treated_level = "increase",
  control_level = "decrease",
  f.ps = ps_formula,
  f.out = outcome_formula_2,
  data = d,
  ps_params = gbm_inc_dec,
  k = 2,
  stratified = TRUE,
  n_boot = 500,
  seed = run,
  verbose = TRUE,
  parallel = TRUE,
  n_cores = NULL,
  plot_diagnostics = TRUE
)

# =============================================================================
# Model 3: Cross-fitted Bootstrap
# =============================================================================

cat("\n\n")
cat("=======================================================================\n")
cat("MODEL 3: Cross-Fitted Bootstrap (Baseline loneliness + depression)\n")
cat(
  "=======================================================================\n\n"
)

results_cf$model3 = cf_DR_att(
  outcome = "severe_loneliness",
  treatment = "remote_contact",
  treated_level = "increase",
  control_level = "decrease",
  f.ps = ps_formula,
  f.out = outcome_formula_3,
  data = d,
  ps_params = gbm_inc_dec,
  k = 2,
  stratified = TRUE,
  n_boot = 500,
  seed = run,
  verbose = TRUE,
  parallel = TRUE,
  n_cores = NULL,
  plot_diagnostics = TRUE
)

# =============================================================================
# Model 4: Cross-fitted Bootstrap
# =============================================================================

cat("\n\n")
cat("=======================================================================\n")
cat("MODEL 4: Cross-Fitted Bootstrap (Fully saturated)\n")
cat(
  "=======================================================================\n\n"
)

results_cf$model4 = cf_DR_att(
  outcome = "severe_loneliness",
  treatment = "remote_contact",
  treated_level = "increase",
  control_level = "decrease",
  f.ps = ps_formula,
  f.out = outcome_formula_4,
  data = d,
  ps_params = gbm_inc_dec,
  k = 2,
  stratified = TRUE,
  n_boot = 500,
  seed = run,
  verbose = TRUE,
  parallel = TRUE,
  n_cores = NULL,
  plot_diagnostics = TRUE
)

# =============================================================================
# Diagnostic plots - Cross-fitted estimates across models
# =============================================================================

cat("\n\n")
cat("=======================================================================\n")
cat("CREATING COMPARATIVE DIAGNOSTIC PLOTS\n")
cat(
  "=======================================================================\n\n"
)

# Comparison across models for AIPW
par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))

hist(
  results_cf$model1$aipw$bootstrap_samples,
  main = "CF-AIPW Model 1 - Inc vs Dec",
  xlab = "ATT",
  col = "lightblue",
  breaks = 30
)
abline(v = results_cf$model1$aipw$att, col = "red", lwd = 2, lty = 2)

hist(
  results_cf$model2$aipw$bootstrap_samples,
  main = "CF-AIPW Model 2 - Inc vs Dec",
  xlab = "ATT",
  col = "lightblue",
  breaks = 30
)
abline(v = results_cf$model2$aipw$att, col = "red", lwd = 2, lty = 2)

hist(
  results_cf$model3$aipw$bootstrap_samples,
  main = "CF-AIPW Model 3 - Inc vs Dec",
  xlab = "ATT",
  col = "lightblue",
  breaks = 30
)
abline(v = results_cf$model3$aipw$att, col = "red", lwd = 2, lty = 2)

hist(
  results_cf$model4$aipw$bootstrap_samples,
  main = "CF-AIPW Model 4 - Inc vs Dec",
  xlab = "ATT",
  col = "lightblue",
  breaks = 30
)
abline(v = results_cf$model4$aipw$att, col = "red", lwd = 2, lty = 2)

par(mfrow = c(1, 1))

# Comparison across models for DRS
par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))

hist(
  results_cf$model1$drs$bootstrap_samples,
  main = "CF-DRS Model 1 - Inc vs Dec",
  xlab = "ATT",
  col = "lightgreen",
  breaks = 30
)
abline(v = results_cf$model1$drs$att, col = "red", lwd = 2, lty = 2)

hist(
  results_cf$model2$drs$bootstrap_samples,
  main = "CF-DRS Model 2 - Inc vs Dec",
  xlab = "ATT",
  col = "lightgreen",
  breaks = 30
)
abline(v = results_cf$model2$drs$att, col = "red", lwd = 2, lty = 2)

hist(
  results_cf$model3$drs$bootstrap_samples,
  main = "CF-DRS Model 3 - Inc vs Dec",
  xlab = "ATT",
  col = "lightgreen",
  breaks = 30
)
abline(v = results_cf$model3$drs$att, col = "red", lwd = 2, lty = 2)

hist(
  results_cf$model4$drs$bootstrap_samples,
  main = "CF-DRS Model 4 - Inc vs Dec",
  xlab = "ATT",
  col = "lightgreen",
  breaks = 30
)
abline(v = results_cf$model4$drs$att, col = "red", lwd = 2, lty = 2)

par(mfrow = c(1, 1))

# Overlay densities for AIPW across models
plot(
  density(results_cf$model1$aipw$bootstrap_samples),
  main = "CF-AIPW Estimates Across Models - Inc vs Dec",
  xlab = "ATT",
  col = "blue",
  lwd = 2,
  xlim = range(c(
    results_cf$model1$aipw$bootstrap_samples,
    results_cf$model2$aipw$bootstrap_samples,
    results_cf$model3$aipw$bootstrap_samples,
    results_cf$model4$aipw$bootstrap_samples
  ))
)
lines(
  density(results_cf$model2$aipw$bootstrap_samples),
  col = "red",
  lwd = 2
)
lines(
  density(results_cf$model3$aipw$bootstrap_samples),
  col = "green",
  lwd = 2
)
lines(
  density(results_cf$model4$aipw$bootstrap_samples),
  col = "purple",
  lwd = 2
)
legend(
  "topleft",
  legend = c("Model 1", "Model 2", "Model 3", "Model 4"),
  col = c("blue", "red", "green", "purple"),
  lwd = 2
)

# Overlay densities for DRS across models
plot(
  density(results_cf$model1$drs$bootstrap_samples),
  main = "CF-DRS Estimates Across Models - Inc vs Dec",
  xlab = "ATT",
  col = "blue",
  lwd = 2,
  xlim = range(c(
    results_cf$model1$drs$bootstrap_samples,
    results_cf$model2$drs$bootstrap_samples,
    results_cf$model3$drs$bootstrap_samples,
    results_cf$model4$drs$bootstrap_samples
  ))
)
lines(
  density(results_cf$model2$drs$bootstrap_samples),
  col = "red",
  lwd = 2
)
lines(
  density(results_cf$model3$drs$bootstrap_samples),
  col = "green",
  lwd = 2
)
lines(
  density(results_cf$model4$drs$bootstrap_samples),
  col = "purple",
  lwd = 2
)
legend(
  "topleft",
  legend = c("Model 1", "Model 2", "Model 3", "Model 4"),
  col = c("blue", "red", "green", "purple"),
  lwd = 2
)

# =============================================================================
# Save results
# =============================================================================

cat("\n\n")
cat("=======================================================================\n")
cat("SAVING RESULTS\n")
cat(
  "=======================================================================\n\n"
)

# Create output directory if it doesn't exist
if (!dir.exists("results/outcome/cross_fitted")) {
  dir.create("results/outcome/cross_fitted", recursive = TRUE)
}

# Save all individual results
saveRDS(results_cf$model1, "results/outcome/cross_fitted/model1_cf.rds")
saveRDS(results_cf$model2, "results/outcome/cross_fitted/model2_cf.rds")
saveRDS(results_cf$model3, "results/outcome/cross_fitted/model3_cf.rds")
saveRDS(results_cf$model4, "results/outcome/cross_fitted/model4_cf.rds")

# Save entire results list
saveRDS(results_cf, "results/outcome/cross_fitted/all_results_cf.rds")

# Create comprehensive summary table
summary_table_cf = data.frame(
  Model = paste0("Model ", 1:4),
  Method = rep("Cross-Fitted", 4),
  K_folds = rep(5, 4),
  AIPW_ATT = c(
    results_cf$model1$aipw$att,
    results_cf$model2$aipw$att,
    results_cf$model3$aipw$att,
    results_cf$model4$aipw$att
  ),
  AIPW_SE = c(
    results_cf$model1$aipw$se,
    results_cf$model2$aipw$se,
    results_cf$model3$aipw$se,
    results_cf$model4$aipw$se
  ),
  AIPW_CI_lower = c(
    results_cf$model1$aipw$ci_percentile[1],
    results_cf$model2$aipw$ci_percentile[1],
    results_cf$model3$aipw$ci_percentile[1],
    results_cf$model4$aipw$ci_percentile[1]
  ),
  AIPW_CI_upper = c(
    results_cf$model1$aipw$ci_percentile[2],
    results_cf$model2$aipw$ci_percentile[2],
    results_cf$model3$aipw$ci_percentile[2],
    results_cf$model4$aipw$ci_percentile[2]
  ),
  AIPW_pval = c(
    results_cf$model1$aipw$pval,
    results_cf$model2$aipw$pval,
    results_cf$model3$aipw$pval,
    results_cf$model4$aipw$pval
  ),
  DRS_ATT = c(
    results_cf$model1$drs$att,
    results_cf$model2$drs$att,
    results_cf$model3$drs$att,
    results_cf$model4$drs$att
  ),
  DRS_SE = c(
    results_cf$model1$drs$se,
    results_cf$model2$drs$se,
    results_cf$model3$drs$se,
    results_cf$model4$drs$se
  ),
  DRS_CI_lower = c(
    results_cf$model1$drs$ci_percentile[1],
    results_cf$model2$drs$ci_percentile[1],
    results_cf$model3$drs$ci_percentile[1],
    results_cf$model4$drs$ci_percentile[1]
  ),
  DRS_CI_upper = c(
    results_cf$model1$drs$ci_percentile[2],
    results_cf$model2$drs$ci_percentile[2],
    results_cf$model3$drs$ci_percentile[2],
    results_cf$model4$drs$ci_percentile[2]
  ),
  DRS_pval = c(
    results_cf$model1$drs$pval,
    results_cf$model2$drs$pval,
    results_cf$model3$drs$pval,
    results_cf$model4$drs$pval
  ),
  # Add balance diagnostics
  Avg_ASD_mean = c(
    results_cf$model1$balance_diagnostics$avg_asd$mean,
    results_cf$model2$balance_diagnostics$avg_asd$mean,
    results_cf$model3$balance_diagnostics$avg_asd$mean,
    results_cf$model4$balance_diagnostics$avg_asd$mean
  ),
  Max_ASD_mean = c(
    results_cf$model1$balance_diagnostics$max_asd$mean,
    results_cf$model2$balance_diagnostics$max_asd$mean,
    results_cf$model3$balance_diagnostics$max_asd$mean,
    results_cf$model4$balance_diagnostics$max_asd$mean
  ),
  Avg_ESS = c(
    results_cf$model1$balance_diagnostics$ess$mean,
    results_cf$model2$balance_diagnostics$ess$mean,
    results_cf$model3$balance_diagnostics$ess$mean,
    results_cf$model4$balance_diagnostics$ess$mean
  ),
  N_boot_successful = c(
    results_cf$model1$n_boot,
    results_cf$model2$n_boot,
    results_cf$model3$n_boot,
    results_cf$model4$n_boot
  ),
  N_boot_failed = c(
    results_cf$model1$n_failed,
    results_cf$model2$n_failed,
    results_cf$model3$n_failed,
    results_cf$model4$n_failed
  )
)

# Round numeric columns for display
summary_table_cf[, 4:ncol(summary_table_cf)] =
  round(summary_table_cf[, 4:ncol(summary_table_cf)], 4)

cat("\n\nCross-Fitted Summary Table:\n")
print(summary_table_cf)

# Save as CSV
write.csv(
  summary_table_cf,
  "results/outcome/cross_fitted/bootstrap_CF_summary.csv",
  row.names = FALSE
)

cat("\n\nResults saved successfully!\n")
cat("Output directory: results/outcome/cross_fitted/\n")
cat("  - Individual model results: model1_cf.rds - model4_cf.rds\n")
cat("  - Combined results: all_results_cf.rds\n")
cat("  - Summary table: bootstrap_CF_summary.csv\n\n")

cat("=======================================================================\n")
cat("ANALYSIS COMPLETE\n")
cat(
  "=======================================================================\n\n"
)

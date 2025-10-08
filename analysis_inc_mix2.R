# Bootstrap Doubly Robust ATT Estimation Analysis - Complete Comparison
# Comparison: Increase vs Mix
# Using refactored utils2.R with unified DR_att() function
#
# This script compares:
# - 4 outcome model specifications (Model 1-4)
# - 2 estimation approaches (reweight vs cross-fitted k=2)
# - Balanced bootstrap sampling
# - All available confidence intervals (norm, basic, perc)

# Setup
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
run = 123
source("utils2.R")

# Load data
d = readRDS("data/data_iptw2.rds")

# Create output directories
comparison = "inc_mix"
dir.create("results/outcome/inc_mix", showWarnings = FALSE, recursive = TRUE)
dir.create(
  "results/outcome/inc_mix/figures",
  showWarnings = FALSE,
  recursive = TRUE
)

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
      health_pre +
      chronic +
      death_due_covid +
      ppl_infected +
      income_loss +
      neighborhood +
      baseline_depr +
      baseline_lone
)

# Define 4 different outcome model formulas (RHS only, as character strings)
outcome_formula_1 = "1" # No covariates

outcome_formula_2 = "baseline_lone" # Just baseline loneliness

outcome_formula_3 = "baseline_lone + baseline_depr" # Baseline mental health

outcome_formula_4 = "baseline_lone + baseline_depr + female + age_cat + edu +
                     emp_status + income + marital + coliving +
                     health_pre + chronic + death_due_covid +
                     ppl_infected + income_loss + neighborhood" # Full adjustment

# Load tuned GBM parameters for this comparison
gbm_params = readRDS("results/weighting/inc_mix_params.rds")

# Bootstrap parameters
n_boot = 2000
n_cores = NULL # Auto-detect

cat(
  "==========================================================================\n"
)
cat("BOOTSTRAP DOUBLY ROBUST ATT ESTIMATION\n")
cat("Comparison: Increase vs Mix\n")
cat(
  "==========================================================================\n\n"
)

cat("Analysis configuration:\n")
cat("  - 4 outcome model specifications\n")
cat("  - 2 estimation approaches (reweight, cross-fitted k=2)\n")
cat("  - Bootstrap: balanced sampling, n =", n_boot, "\n")
cat("  - Confidence intervals: normal, basic, percentile\n\n")

# Store results
results = list()

# =============================================================================
# PART 1: REWEIGHT METHOD (No Cross-Fitting)
# =============================================================================

cat("\n\n")
cat(
  "==========================================================================\n"
)
cat("PART 1: REWEIGHT METHOD (No Cross-Fitting)\n")
cat(
  "==========================================================================\n\n"
)

# -----------------------------------------------------------------------------
# Model 1 - Reweight
# -----------------------------------------------------------------------------
cat("\n--- Model 1 - Reweight ---\n")
cat("Outcome formula: No covariates (intercept only)\n\n")

results$model1_reweight = DR_att(
  outcome = "severe_loneliness",
  treatment = "remote_contact",
  treated_level = "increase",
  control_level = "mix",
  f.ps = ps_formula,
  f.out = outcome_formula_1,
  data = d,
  gbm_params = gbm_params,
  bootstrap_method = "reweight",
  cross_fitting = FALSE,
  stratified = TRUE,
  n_boot = n_boot,
  seed = run,
  verbose = TRUE,
  parallel = TRUE,
  n_cores = n_cores,
  plot_diagnostics = TRUE,
  ci_type = "all",
  sim = "balanced",
  save_to = "results/outcome/inc_mix/figures"
)

# -----------------------------------------------------------------------------
# Model 2 - Reweight
# -----------------------------------------------------------------------------
cat("\n--- Model 2 - Reweight ---\n")
cat("Outcome formula: baseline_lone\n\n")

results$model2_reweight = DR_att(
  outcome = "severe_loneliness",
  treatment = "remote_contact",
  treated_level = "increase",
  control_level = "mix",
  f.ps = ps_formula,
  f.out = outcome_formula_2,
  data = d,
  gbm_params = gbm_params,
  bootstrap_method = "reweight",
  cross_fitting = FALSE,
  stratified = TRUE,
  n_boot = n_boot,
  seed = run,
  verbose = TRUE,
  parallel = TRUE,
  n_cores = n_cores,
  plot_diagnostics = TRUE,
  ci_type = "all",
  sim = "balanced",
  save_to = "results/outcome/inc_mix/figures"
)

# -----------------------------------------------------------------------------
# Model 3 - Reweight
# -----------------------------------------------------------------------------
cat("\n--- Model 3 - Reweight ---\n")
cat("Outcome formula: baseline_lone + baseline_depr\n\n")

results$model3_reweight = DR_att(
  outcome = "severe_loneliness",
  treatment = "remote_contact",
  treated_level = "increase",
  control_level = "mix",
  f.ps = ps_formula,
  f.out = outcome_formula_3,
  data = d,
  gbm_params = gbm_params,
  bootstrap_method = "reweight",
  cross_fitting = FALSE,
  stratified = TRUE,
  n_boot = n_boot,
  seed = run,
  verbose = TRUE,
  parallel = TRUE,
  n_cores = n_cores,
  plot_diagnostics = TRUE,
  ci_type = "all",
  sim = "balanced",
  save_to = "results/outcome/inc_mix/figures"
)

# -----------------------------------------------------------------------------
# Model 4 - Reweight
# -----------------------------------------------------------------------------
cat("\n--- Model 4 - Reweight ---\n")
cat("Outcome formula: Full adjustment\n\n")

results$model4_reweight = DR_att(
  outcome = "severe_loneliness",
  treatment = "remote_contact",
  treated_level = "increase",
  control_level = "mix",
  f.ps = ps_formula,
  f.out = outcome_formula_4,
  data = d,
  gbm_params = gbm_params,
  bootstrap_method = "reweight",
  cross_fitting = FALSE,
  stratified = TRUE,
  n_boot = n_boot,
  seed = run,
  verbose = TRUE,
  parallel = TRUE,
  n_cores = n_cores,
  plot_diagnostics = TRUE,
  ci_type = "all",
  sim = "balanced",
  save_to = "results/outcome/inc_mix/figures"
)

# =============================================================================
# PART 2: CROSS-FITTED METHOD (k=2)
# =============================================================================

cat("\n\n")
cat(
  "==========================================================================\n"
)
cat("PART 2: CROSS-FITTED METHOD (k=2)\n")
cat(
  "==========================================================================\n\n"
)

# -----------------------------------------------------------------------------
# Model 1 - Cross-Fitted
# -----------------------------------------------------------------------------
cat("\n--- Model 1 - Cross-Fitted (k=2) ---\n")
cat("Outcome formula: No covariates (intercept only)\n\n")

results$model1_cf2 = DR_att(
  outcome = "severe_loneliness",
  treatment = "remote_contact",
  treated_level = "increase",
  control_level = "mix",
  f.ps = ps_formula,
  f.out = outcome_formula_1,
  data = d,
  gbm_params = gbm_params,
  bootstrap_method = "reweight",
  cross_fitting = TRUE,
  k = 5,
  stratify_folds = TRUE,
  stratified = TRUE,
  n_boot = n_boot,
  seed = run,
  verbose = TRUE,
  parallel = TRUE,
  n_cores = n_cores,
  plot_diagnostics = TRUE,
  ci_type = "all",
  sim = "balanced",
  save_to = "results/outcome/inc_mix/figures"
)

# -----------------------------------------------------------------------------
# Model 2 - Cross-Fitted
# -----------------------------------------------------------------------------
cat("\n--- Model 2 - Cross-Fitted (k=2) ---\n")
cat("Outcome formula: baseline_lone\n\n")

results$model2_cf2 = DR_att(
  outcome = "severe_loneliness",
  treatment = "remote_contact",
  treated_level = "increase",
  control_level = "mix",
  f.ps = ps_formula,
  f.out = outcome_formula_2,
  data = d,
  gbm_params = gbm_params,
  bootstrap_method = "reweight",
  cross_fitting = TRUE,
  k = 5,
  stratify_folds = TRUE,
  stratified = TRUE,
  n_boot = n_boot,
  seed = run,
  verbose = TRUE,
  parallel = TRUE,
  n_cores = n_cores,
  plot_diagnostics = TRUE,
  ci_type = "all",
  sim = "balanced",
  save_to = "results/outcome/inc_mix/figures"
)

# -----------------------------------------------------------------------------
# Model 3 - Cross-Fitted
# -----------------------------------------------------------------------------
cat("\n--- Model 3 - Cross-Fitted (k=2) ---\n")
cat("Outcome formula: baseline_lone + baseline_depr\n\n")

results$model3_cf2 = DR_att(
  outcome = "severe_loneliness",
  treatment = "remote_contact",
  treated_level = "increase",
  control_level = "mix",
  f.ps = ps_formula,
  f.out = outcome_formula_3,
  data = d,
  gbm_params = gbm_params,
  bootstrap_method = "reweight",
  cross_fitting = TRUE,
  k = 5,
  stratify_folds = TRUE,
  stratified = TRUE,
  n_boot = n_boot,
  seed = run,
  verbose = TRUE,
  parallel = TRUE,
  n_cores = n_cores,
  plot_diagnostics = TRUE,
  ci_type = "all",
  sim = "balanced",
  save_to = "results/outcome/inc_mix/figures"
)

# -----------------------------------------------------------------------------
# Model 4 - Cross-Fitted
# -----------------------------------------------------------------------------
cat("\n--- Model 4 - Cross-Fitted (k=2) ---\n")
cat("Outcome formula: Full adjustment\n\n")

results$model4_cf2 = DR_att(
  outcome = "severe_loneliness",
  treatment = "remote_contact",
  treated_level = "increase",
  control_level = "mix",
  f.ps = ps_formula,
  f.out = outcome_formula_4,
  data = d,
  gbm_params = gbm_params,
  bootstrap_method = "reweight",
  cross_fitting = TRUE,
  k = 5,
  stratify_folds = TRUE,
  stratified = TRUE,
  n_boot = n_boot,
  seed = run,
  verbose = TRUE,
  parallel = TRUE,
  n_cores = n_cores,
  plot_diagnostics = TRUE,
  ci_type = "all",
  sim = "balanced",
  save_to = "results/outcome/inc_mix/figures"
)

# =============================================================================
# COMPARISON PLOTS
# =============================================================================

cat("\n\n")
cat(
  "==========================================================================\n"
)
cat("GENERATING COMPARISON PLOTS\n")
cat(
  "==========================================================================\n\n"
)

# Create high-resolution comparison plot
png(
  "results/outcome/inc_mix/figures/model_comparison_distributions.png",
  width = 4000,
  height = 3000,
  res = 300
)
par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))

# AIPW - Reweight Method
plot(
  density(results$model1_reweight$aipw$bootstrap_samples),
  main = "AIPW: Reweight Method",
  xlab = "ATT",
  col = "blue",
  lwd = 2,
  xlim = range(c(
    results$model1_reweight$aipw$bootstrap_samples,
    results$model2_reweight$aipw$bootstrap_samples,
    results$model3_reweight$aipw$bootstrap_samples,
    results$model4_reweight$aipw$bootstrap_samples
  ))
)
lines(
  density(results$model2_reweight$aipw$bootstrap_samples),
  col = "red",
  lwd = 2
)
lines(
  density(results$model3_reweight$aipw$bootstrap_samples),
  col = "green",
  lwd = 2
)
lines(
  density(results$model4_reweight$aipw$bootstrap_samples),
  col = "purple",
  lwd = 2
)
legend(
  "topleft",
  legend = c("Model 1", "Model 2", "Model 3", "Model 4"),
  col = c("blue", "red", "green", "purple"),
  lwd = 2,
  cex = 0.8
)

# AIPW - Cross-Fitted Method
plot(
  density(results$model1_cf2$aipw$bootstrap_samples),
  main = "AIPW: Cross-Fitted (k=2)",
  xlab = "ATT",
  col = "blue",
  lwd = 2,
  xlim = range(c(
    results$model1_cf2$aipw$bootstrap_samples,
    results$model2_cf2$aipw$bootstrap_samples,
    results$model3_cf2$aipw$bootstrap_samples,
    results$model4_cf2$aipw$bootstrap_samples
  ))
)
lines(density(results$model2_cf2$aipw$bootstrap_samples), col = "red", lwd = 2)
lines(
  density(results$model3_cf2$aipw$bootstrap_samples),
  col = "green",
  lwd = 2
)
lines(
  density(results$model4_cf2$aipw$bootstrap_samples),
  col = "purple",
  lwd = 2
)
legend(
  "topleft",
  legend = c("Model 1", "Model 2", "Model 3", "Model 4"),
  col = c("blue", "red", "green", "purple"),
  lwd = 2,
  cex = 0.8
)

# DRS - Reweight Method
plot(
  density(results$model1_reweight$drs$bootstrap_samples),
  main = "DRS: Reweight Method",
  xlab = "ATT",
  col = "blue",
  lwd = 2,
  xlim = range(c(
    results$model1_reweight$drs$bootstrap_samples,
    results$model2_reweight$drs$bootstrap_samples,
    results$model3_reweight$drs$bootstrap_samples,
    results$model4_reweight$drs$bootstrap_samples
  ))
)
lines(
  density(results$model2_reweight$drs$bootstrap_samples),
  col = "red",
  lwd = 2
)
lines(
  density(results$model3_reweight$drs$bootstrap_samples),
  col = "green",
  lwd = 2
)
lines(
  density(results$model4_reweight$drs$bootstrap_samples),
  col = "purple",
  lwd = 2
)
legend(
  "topleft",
  legend = c("Model 1", "Model 2", "Model 3", "Model 4"),
  col = c("blue", "red", "green", "purple"),
  lwd = 2,
  cex = 0.8
)

# DRS - Cross-Fitted Method
plot(
  density(results$model1_cf2$drs$bootstrap_samples),
  main = "DRS: Cross-Fitted (k=2)",
  xlab = "ATT",
  col = "blue",
  lwd = 2,
  xlim = range(c(
    results$model1_cf2$drs$bootstrap_samples,
    results$model2_cf2$drs$bootstrap_samples,
    results$model3_cf2$drs$bootstrap_samples,
    results$model4_cf2$drs$bootstrap_samples
  ))
)
lines(density(results$model2_cf2$drs$bootstrap_samples), col = "red", lwd = 2)
lines(density(results$model3_cf2$drs$bootstrap_samples), col = "green", lwd = 2)
lines(
  density(results$model4_cf2$drs$bootstrap_samples),
  col = "purple",
  lwd = 2
)
legend(
  "topleft",
  legend = c("Model 1", "Model 2", "Model 3", "Model 4"),
  col = c("blue", "red", "green", "purple"),
  lwd = 2,
  cex = 0.8
)

par(mfrow = c(1, 1))
dev.off()

cat(
  "Comparison plot saved to: results/outcome/inc_mix/figures/model_comparison_distributions.png\n"
)

# =============================================================================
# CREATE COMPREHENSIVE SUMMARY TABLE
# =============================================================================

cat("\n\n")
cat(
  "==========================================================================\n"
)
cat("CREATING SUMMARY TABLES\n")
cat(
  "==========================================================================\n\n"
)

# Helper function to extract results
extract_results = function(result, model_num, method, estimator) {
  est = result[[tolower(estimator)]]

  data.frame(
    Model = model_num,
    Method = method,
    Estimator = estimator,
    ATT = est$att,
    SE = est$se,
    CI_Normal_Lower = est$ci_normal[1],
    CI_Normal_Upper = est$ci_normal[2],
    CI_Basic_Lower = ifelse(!is.null(est$ci_basic), est$ci_basic[1], NA),
    CI_Basic_Upper = ifelse(!is.null(est$ci_basic), est$ci_basic[2], NA),
    CI_Perc_Lower = est$ci_percentile[1],
    CI_Perc_Upper = est$ci_percentile[2],
    p_value = est$pval,
    n_boot = result$n_boot,
    n_failed = result$n_failed
  )
}

# Create comprehensive summary
summary_rows = list()

# Model 1
summary_rows[[1]] = extract_results(
  results$model1_reweight,
  1,
  "Reweight",
  "AIPW"
)
summary_rows[[2]] = extract_results(
  results$model1_reweight,
  1,
  "Reweight",
  "DRS"
)
summary_rows[[3]] = extract_results(
  results$model1_cf2,
  1,
  "CrossFit_k2",
  "AIPW"
)
summary_rows[[4]] = extract_results(results$model1_cf2, 1, "CrossFit_k2", "DRS")

# Model 2
summary_rows[[5]] = extract_results(
  results$model2_reweight,
  2,
  "Reweight",
  "AIPW"
)
summary_rows[[6]] = extract_results(
  results$model2_reweight,
  2,
  "Reweight",
  "DRS"
)
summary_rows[[7]] = extract_results(
  results$model2_cf2,
  2,
  "CrossFit_k2",
  "AIPW"
)
summary_rows[[8]] = extract_results(results$model2_cf2, 2, "CrossFit_k2", "DRS")

# Model 3
summary_rows[[9]] = extract_results(
  results$model3_reweight,
  3,
  "Reweight",
  "AIPW"
)
summary_rows[[10]] = extract_results(
  results$model3_reweight,
  3,
  "Reweight",
  "DRS"
)
summary_rows[[11]] = extract_results(
  results$model3_cf2,
  3,
  "CrossFit_k2",
  "AIPW"
)
summary_rows[[12]] = extract_results(
  results$model3_cf2,
  3,
  "CrossFit_k2",
  "DRS"
)

# Model 4
summary_rows[[13]] = extract_results(
  results$model4_reweight,
  4,
  "Reweight",
  "AIPW"
)
summary_rows[[14]] = extract_results(
  results$model4_reweight,
  4,
  "Reweight",
  "DRS"
)
summary_rows[[15]] = extract_results(
  results$model4_cf2,
  4,
  "CrossFit_k2",
  "AIPW"
)
summary_rows[[16]] = extract_results(
  results$model4_cf2,
  4,
  "CrossFit_k2",
  "DRS"
)

summary_table = do.call(rbind, summary_rows)
rownames(summary_table) = NULL

# Add balance diagnostics
extract_balance = function(result, model_num, method) {
  if (!is.null(result$balance_diagnostics)) {
    data.frame(
      Model = model_num,
      Method = method,
      Avg_ASD_Mean = result$balance_diagnostics$avg_asd$mean,
      Avg_ASD_SD = result$balance_diagnostics$avg_asd$sd,
      Max_ASD_Mean = result$balance_diagnostics$max_asd$mean,
      Max_ASD_SD = result$balance_diagnostics$max_asd$sd,
      ESS_Mean = result$balance_diagnostics$ess$mean,
      ESS_SD = result$balance_diagnostics$ess$sd
    )
  } else {
    NULL
  }
}

balance_rows = list(
  extract_balance(results$model1_reweight, 1, "Reweight"),
  extract_balance(results$model1_cf2, 1, "CrossFit_k2"),
  extract_balance(results$model2_reweight, 2, "Reweight"),
  extract_balance(results$model2_cf2, 2, "CrossFit_k2"),
  extract_balance(results$model3_reweight, 3, "Reweight"),
  extract_balance(results$model3_cf2, 3, "CrossFit_k2"),
  extract_balance(results$model4_reweight, 4, "Reweight"),
  extract_balance(results$model4_cf2, 4, "CrossFit_k2")
)

balance_table = do.call(rbind, balance_rows)
rownames(balance_table) = NULL

# Print summaries
cat("\n=== MAIN RESULTS SUMMARY ===\n")
print(summary_table, digits = 4)

cat("\n\n=== BALANCE DIAGNOSTICS SUMMARY ===\n")
print(balance_table, digits = 4)

# =============================================================================
# SAVE RESULTS
# =============================================================================

cat("\n\n")
cat(
  "==========================================================================\n"
)
cat("SAVING RESULTS\n")
cat(
  "==========================================================================\n\n"
)

# Save complete results object
saveRDS(results, "results/outcome/inc_mix/complete_results.rds")
cat("Complete results saved to: results/outcome/inc_mix/complete_results.rds\n")

# Save summary tables
write.csv(
  summary_table,
  "results/outcome/inc_mix/summary_table.csv",
  row.names = FALSE
)
cat("Summary table saved to: results/outcome/inc_mix/summary_table.csv\n")

saveRDS(summary_table, "results/outcome/inc_mix/summary_table.rds")
cat("Summary table saved to: results/outcome/inc_mix/summary_table.rds\n")

write.csv(
  balance_table,
  "results/outcome/inc_mix/balance_diagnostics.csv",
  row.names = FALSE
)
cat(
  "Balance diagnostics saved to: results/outcome/inc_mix/balance_diagnostics.csv\n"
)

saveRDS(balance_table, "results/outcome/inc_mix/balance_diagnostics.rds")
cat(
  "Balance diagnostics saved to: results/outcome/inc_mix/balance_diagnostics.rds\n"
)

# Save individual model results
saveRDS(results$model1_reweight, "results/outcome/inc_mix/model1_reweight.rds")
saveRDS(results$model2_reweight, "results/outcome/inc_mix/model2_reweight.rds")
saveRDS(results$model3_reweight, "results/outcome/inc_mix/model3_reweight.rds")
saveRDS(results$model4_reweight, "results/outcome/inc_mix/model4_reweight.rds")
saveRDS(results$model1_cf2, "results/outcome/inc_mix/model1_cf2.rds")
saveRDS(results$model2_cf2, "results/outcome/inc_mix/model2_cf2.rds")
saveRDS(results$model3_cf2, "results/outcome/inc_mix/model3_cf2.rds")
saveRDS(results$model4_cf2, "results/outcome/inc_mix/model4_cf2.rds")

cat("Individual model results saved.\n")

cat(
  "\n==========================================================================\n"
)
cat("ANALYSIS COMPLETE!\n")
cat(
  "==========================================================================\n"
)
cat("\nAll results saved to: results/outcome/inc_mix/\n")
cat("All plots saved to: results/outcome/inc_mix/figures/\n\n")

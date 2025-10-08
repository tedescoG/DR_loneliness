# Bootstrap Doubly Robust ATT Estimation Analysis - Cross-Fitted Method
# Both Comparisons: Increase vs Decrease AND Increase vs Mix
# Using refactored utils2.R with unified DR_att() function
#
# This script runs cross-fitted analysis (k=5) for:
# - 4 outcome model specifications (Model 1-4)
# - 2 pairwise comparisons (inc_dec and inc_mix)
# - Balanced bootstrap sampling
# - All available confidence intervals (norm, basic, perc)

# Setup
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
run = 123
source("utils2.R")

# Load data
d = readRDS("data/data_iptw2.rds")

# Define propensity score formula (same for both comparisons)
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

# Bootstrap parameters
n_boot = 2000
n_cores = NULL # Auto-detect
k_folds = 5 # Number of folds for cross-fitting

# =============================================================================
# COMPARISON 1: INCREASE VS DECREASE
# =============================================================================

cat(
  "\n=============================================================================\n"
)
cat("STARTING ANALYSIS: INCREASE VS DECREASE (CROSS-FITTED)\n")
cat(
  "=============================================================================\n\n"
)

# Create output directories
dir.create(
  "results/outcome/crossfitting/inc_dec",
  showWarnings = FALSE,
  recursive = TRUE
)
dir.create(
  "results/outcome/crossfitting/inc_dec/figures",
  showWarnings = FALSE,
  recursive = TRUE
)

# Load tuned GBM parameters
gbm_params_inc_dec = readRDS("results/weighting/inc_dec_params.rds")

# Store results
results_inc_dec = list()

# -----------------------------------------------------------------------------
# Model 1
# -----------------------------------------------------------------------------

results_inc_dec$model1_cf = DR_att(
  outcome = "loneliness",
  treatment = "remote_contact",
  treated_level = "increase",
  control_level = "decrease",
  f.ps = ps_formula,
  f.out = outcome_formula_1,
  data = d,
  gbm_params = gbm_params_inc_dec,
  bootstrap_method = "reweight",
  cross_fitting = TRUE,
  k = k_folds,
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
  save_to = "results/outcome/crossfitting/inc_dec/figures"
)

# -----------------------------------------------------------------------------
# Model 2
# -----------------------------------------------------------------------------

results_inc_dec$model2_cf = DR_att(
  outcome = "loneliness",
  treatment = "remote_contact",
  treated_level = "increase",
  control_level = "decrease",
  f.ps = ps_formula,
  f.out = outcome_formula_2,
  data = d,
  gbm_params = gbm_params_inc_dec,
  bootstrap_method = "reweight",
  cross_fitting = TRUE,
  k = k_folds,
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
  save_to = "results/outcome/crossfitting/inc_dec/figures"
)

# -----------------------------------------------------------------------------
# Model 3
# -----------------------------------------------------------------------------

results_inc_dec$model3_cf = DR_att(
  outcome = "loneliness",
  treatment = "remote_contact",
  treated_level = "increase",
  control_level = "decrease",
  f.ps = ps_formula,
  f.out = outcome_formula_3,
  data = d,
  gbm_params = gbm_params_inc_dec,
  bootstrap_method = "reweight",
  cross_fitting = TRUE,
  k = k_folds,
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
  save_to = "results/outcome/crossfitting/inc_dec/figures"
)

# -----------------------------------------------------------------------------
# Model 4
# -----------------------------------------------------------------------------

results_inc_dec$model4_cf = DR_att(
  outcome = "loneliness",
  treatment = "remote_contact",
  treated_level = "increase",
  control_level = "decrease",
  f.ps = ps_formula,
  f.out = outcome_formula_4,
  data = d,
  gbm_params = gbm_params_inc_dec,
  bootstrap_method = "reweight",
  cross_fitting = TRUE,
  k = k_folds,
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
  save_to = "results/outcome/crossfitting/inc_dec/figures"
)

# -----------------------------------------------------------------------------
# COMPARISON PLOTS - INC_DEC
# -----------------------------------------------------------------------------

png(
  "results/outcome/crossfitting/inc_dec/figures/model_comparison_distributions.png",
  width = 3000,
  height = 2000,
  res = 300
)
par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))

# AIPW - Cross-Fitted Method
plot(
  density(results_inc_dec$model1_cf$aipw$bootstrap_samples),
  main = paste0("AIPW: Cross-Fitted (k=", k_folds, ")"),
  xlab = "ATT",
  col = "blue",
  lwd = 2,
  xlim = range(c(
    results_inc_dec$model1_cf$aipw$bootstrap_samples,
    results_inc_dec$model2_cf$aipw$bootstrap_samples,
    results_inc_dec$model3_cf$aipw$bootstrap_samples,
    results_inc_dec$model4_cf$aipw$bootstrap_samples
  ))
)
lines(
  density(results_inc_dec$model2_cf$aipw$bootstrap_samples),
  col = "red",
  lwd = 2
)
lines(
  density(results_inc_dec$model3_cf$aipw$bootstrap_samples),
  col = "green",
  lwd = 2
)
lines(
  density(results_inc_dec$model4_cf$aipw$bootstrap_samples),
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
  density(results_inc_dec$model1_cf$drs$bootstrap_samples),
  main = paste0("DRS: Cross-Fitted (k=", k_folds, ")"),
  xlab = "ATT",
  col = "blue",
  lwd = 2,
  xlim = range(c(
    results_inc_dec$model1_cf$drs$bootstrap_samples,
    results_inc_dec$model2_cf$drs$bootstrap_samples,
    results_inc_dec$model3_cf$drs$bootstrap_samples,
    results_inc_dec$model4_cf$drs$bootstrap_samples
  ))
)
lines(
  density(results_inc_dec$model2_cf$drs$bootstrap_samples),
  col = "red",
  lwd = 2
)
lines(
  density(results_inc_dec$model3_cf$drs$bootstrap_samples),
  col = "green",
  lwd = 2
)
lines(
  density(results_inc_dec$model4_cf$drs$bootstrap_samples),
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
  "Comparison plot saved to: results/outcome/crossfitting/inc_dec/figures/model_comparison_distributions.png\n"
)

# -----------------------------------------------------------------------------
# SUMMARY TABLES - INC_DEC
# -----------------------------------------------------------------------------

# Create comprehensive summary
summary_rows_inc_dec = list()

# Model 1
summary_rows_inc_dec[[1]] = extract_results(
  results_inc_dec$model1_cf,
  1,
  paste0("CrossFit_k", k_folds),
  "AIPW"
)
summary_rows_inc_dec[[2]] = extract_results(
  results_inc_dec$model1_cf,
  1,
  paste0("CrossFit_k", k_folds),
  "DRS"
)

# Model 2
summary_rows_inc_dec[[3]] = extract_results(
  results_inc_dec$model2_cf,
  2,
  paste0("CrossFit_k", k_folds),
  "AIPW"
)
summary_rows_inc_dec[[4]] = extract_results(
  results_inc_dec$model2_cf,
  2,
  paste0("CrossFit_k", k_folds),
  "DRS"
)

# Model 3
summary_rows_inc_dec[[5]] = extract_results(
  results_inc_dec$model3_cf,
  3,
  paste0("CrossFit_k", k_folds),
  "AIPW"
)
summary_rows_inc_dec[[6]] = extract_results(
  results_inc_dec$model3_cf,
  3,
  paste0("CrossFit_k", k_folds),
  "DRS"
)

# Model 4
summary_rows_inc_dec[[7]] = extract_results(
  results_inc_dec$model4_cf,
  4,
  paste0("CrossFit_k", k_folds),
  "AIPW"
)
summary_rows_inc_dec[[8]] = extract_results(
  results_inc_dec$model4_cf,
  4,
  paste0("CrossFit_k", k_folds),
  "DRS"
)

summary_table_inc_dec = do.call(rbind, summary_rows_inc_dec)
rownames(summary_table_inc_dec) = NULL

balance_rows_inc_dec = list(
  extract_balance(results_inc_dec$model1_cf, 1, paste0("CrossFit_k", k_folds)),
  extract_balance(results_inc_dec$model2_cf, 2, paste0("CrossFit_k", k_folds)),
  extract_balance(results_inc_dec$model3_cf, 3, paste0("CrossFit_k", k_folds)),
  extract_balance(results_inc_dec$model4_cf, 4, paste0("CrossFit_k", k_folds))
)

balance_table_inc_dec = do.call(rbind, balance_rows_inc_dec)
rownames(balance_table_inc_dec) = NULL

# Print summaries
cat("\n=== INC_DEC SUMMARY TABLE ===\n")
print(summary_table_inc_dec, digits = 4)

cat("\n=== INC_DEC BALANCE TABLE ===\n")
print(balance_table_inc_dec, digits = 4)

# Save results
saveRDS(
  results_inc_dec,
  "results/outcome/crossfitting/inc_dec/complete_results.rds"
)

write.csv(
  summary_table_inc_dec,
  "results/outcome/crossfitting/inc_dec/summary_table.csv",
  row.names = FALSE
)

write.csv(
  balance_table_inc_dec,
  "results/outcome/crossfitting/inc_dec/balance_diagnostics.csv",
  row.names = FALSE
)

# Save individual model results
saveRDS(
  results_inc_dec$model1_cf,
  "results/outcome/crossfitting/inc_dec/model1_cf.rds"
)
saveRDS(
  results_inc_dec$model2_cf,
  "results/outcome/crossfitting/inc_dec/model2_cf.rds"
)
saveRDS(
  results_inc_dec$model3_cf,
  "results/outcome/crossfitting/inc_dec/model3_cf.rds"
)
saveRDS(
  results_inc_dec$model4_cf,
  "results/outcome/crossfitting/inc_dec/model4_cf.rds"
)

cat(
  "\nINC_DEC analysis complete. Results saved to: results/outcome/crossfitting/inc_dec/\n\n"
)

# =============================================================================
# COMPARISON 2: INCREASE VS MIX
# =============================================================================

cat(
  "\n=============================================================================\n"
)
cat("STARTING ANALYSIS: INCREASE VS MIX (CROSS-FITTED)\n")
cat(
  "=============================================================================\n\n"
)

# Create output directories
dir.create(
  "results/outcome/crossfitting/inc_mix",
  showWarnings = FALSE,
  recursive = TRUE
)
dir.create(
  "results/outcome/crossfitting/inc_mix/figures",
  showWarnings = FALSE,
  recursive = TRUE
)

# Load tuned GBM parameters
gbm_params_inc_mix = readRDS("results/weighting/inc_mix_params.rds")

# Store results
results_inc_mix = list()

# -----------------------------------------------------------------------------
# Model 1
# -----------------------------------------------------------------------------

results_inc_mix$model1_cf = DR_att(
  outcome = "loneliness",
  treatment = "remote_contact",
  treated_level = "increase",
  control_level = "mix",
  f.ps = ps_formula,
  f.out = outcome_formula_1,
  data = d,
  gbm_params = gbm_params_inc_mix,
  bootstrap_method = "reweight",
  cross_fitting = TRUE,
  k = k_folds,
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
  save_to = "results/outcome/crossfitting/inc_mix/figures"
)

# -----------------------------------------------------------------------------
# Model 2
# -----------------------------------------------------------------------------

results_inc_mix$model2_cf = DR_att(
  outcome = "loneliness",
  treatment = "remote_contact",
  treated_level = "increase",
  control_level = "mix",
  f.ps = ps_formula,
  f.out = outcome_formula_2,
  data = d,
  gbm_params = gbm_params_inc_mix,
  bootstrap_method = "reweight",
  cross_fitting = TRUE,
  k = k_folds,
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
  save_to = "results/outcome/crossfitting/inc_mix/figures"
)

# -----------------------------------------------------------------------------
# Model 3
# -----------------------------------------------------------------------------

results_inc_mix$model3_cf = DR_att(
  outcome = "loneliness",
  treatment = "remote_contact",
  treated_level = "increase",
  control_level = "mix",
  f.ps = ps_formula,
  f.out = outcome_formula_3,
  data = d,
  gbm_params = gbm_params_inc_mix,
  bootstrap_method = "reweight",
  cross_fitting = TRUE,
  k = k_folds,
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
  save_to = "results/outcome/crossfitting/inc_mix/figures"
)

# -----------------------------------------------------------------------------
# Model 4
# -----------------------------------------------------------------------------

results_inc_mix$model4_cf = DR_att(
  outcome = "loneliness",
  treatment = "remote_contact",
  treated_level = "increase",
  control_level = "mix",
  f.ps = ps_formula,
  f.out = outcome_formula_4,
  data = d,
  gbm_params = gbm_params_inc_mix,
  bootstrap_method = "reweight",
  cross_fitting = TRUE,
  k = k_folds,
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
  save_to = "results/outcome/crossfitting/inc_mix/figures"
)

# -----------------------------------------------------------------------------
# COMPARISON PLOTS - INC_MIX
# -----------------------------------------------------------------------------

png(
  "results/outcome/crossfitting/inc_mix/figures/model_comparison_distributions.png",
  width = 3000,
  height = 2000,
  res = 300
)
par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))

# AIPW - Cross-Fitted Method
plot(
  density(results_inc_mix$model1_cf$aipw$bootstrap_samples),
  main = paste0("AIPW: Cross-Fitted (k=", k_folds, ")"),
  xlab = "ATT",
  col = "blue",
  lwd = 2,
  xlim = range(c(
    results_inc_mix$model1_cf$aipw$bootstrap_samples,
    results_inc_mix$model2_cf$aipw$bootstrap_samples,
    results_inc_mix$model3_cf$aipw$bootstrap_samples,
    results_inc_mix$model4_cf$aipw$bootstrap_samples
  ))
)
lines(
  density(results_inc_mix$model2_cf$aipw$bootstrap_samples),
  col = "red",
  lwd = 2
)
lines(
  density(results_inc_mix$model3_cf$aipw$bootstrap_samples),
  col = "green",
  lwd = 2
)
lines(
  density(results_inc_mix$model4_cf$aipw$bootstrap_samples),
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
  density(results_inc_mix$model1_cf$drs$bootstrap_samples),
  main = paste0("DRS: Cross-Fitted (k=", k_folds, ")"),
  xlab = "ATT",
  col = "blue",
  lwd = 2,
  xlim = range(c(
    results_inc_mix$model1_cf$drs$bootstrap_samples,
    results_inc_mix$model2_cf$drs$bootstrap_samples,
    results_inc_mix$model3_cf$drs$bootstrap_samples,
    results_inc_mix$model4_cf$drs$bootstrap_samples
  ))
)
lines(
  density(results_inc_mix$model2_cf$drs$bootstrap_samples),
  col = "red",
  lwd = 2
)
lines(
  density(results_inc_mix$model3_cf$drs$bootstrap_samples),
  col = "green",
  lwd = 2
)
lines(
  density(results_inc_mix$model4_cf$drs$bootstrap_samples),
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
  "Comparison plot saved to: results/outcome/crossfitting/inc_mix/figures/model_comparison_distributions.png\n"
)

# -----------------------------------------------------------------------------
# SUMMARY TABLES - INC_MIX
# -----------------------------------------------------------------------------

# Create comprehensive summary
summary_rows_inc_mix = list()

# Model 1
summary_rows_inc_mix[[1]] = extract_results(
  results_inc_mix$model1_cf,
  1,
  paste0("CrossFit_k", k_folds),
  "AIPW"
)
summary_rows_inc_mix[[2]] = extract_results(
  results_inc_mix$model1_cf,
  1,
  paste0("CrossFit_k", k_folds),
  "DRS"
)

# Model 2
summary_rows_inc_mix[[3]] = extract_results(
  results_inc_mix$model2_cf,
  2,
  paste0("CrossFit_k", k_folds),
  "AIPW"
)
summary_rows_inc_mix[[4]] = extract_results(
  results_inc_mix$model2_cf,
  2,
  paste0("CrossFit_k", k_folds),
  "DRS"
)

# Model 3
summary_rows_inc_mix[[5]] = extract_results(
  results_inc_mix$model3_cf,
  3,
  paste0("CrossFit_k", k_folds),
  "AIPW"
)
summary_rows_inc_mix[[6]] = extract_results(
  results_inc_mix$model3_cf,
  3,
  paste0("CrossFit_k", k_folds),
  "DRS"
)

# Model 4
summary_rows_inc_mix[[7]] = extract_results(
  results_inc_mix$model4_cf,
  4,
  paste0("CrossFit_k", k_folds),
  "AIPW"
)
summary_rows_inc_mix[[8]] = extract_results(
  results_inc_mix$model4_cf,
  4,
  paste0("CrossFit_k", k_folds),
  "DRS"
)

summary_table_inc_mix = do.call(rbind, summary_rows_inc_mix)
rownames(summary_table_inc_mix) = NULL

balance_rows_inc_mix = list(
  extract_balance(results_inc_mix$model1_cf, 1, paste0("CrossFit_k", k_folds)),
  extract_balance(results_inc_mix$model2_cf, 2, paste0("CrossFit_k", k_folds)),
  extract_balance(results_inc_mix$model3_cf, 3, paste0("CrossFit_k", k_folds)),
  extract_balance(results_inc_mix$model4_cf, 4, paste0("CrossFit_k", k_folds))
)

balance_table_inc_mix = do.call(rbind, balance_rows_inc_mix)
rownames(balance_table_inc_mix) = NULL

# Print summaries
cat("\n=== INC_MIX SUMMARY TABLE ===\n")
print(summary_table_inc_mix, digits = 4)

cat("\n=== INC_MIX BALANCE TABLE ===\n")
print(balance_table_inc_mix, digits = 4)

# Save results
saveRDS(
  results_inc_mix,
  "results/outcome/crossfitting/inc_mix/complete_results.rds"
)

write.csv(
  summary_table_inc_mix,
  "results/outcome/crossfitting/inc_mix/summary_table.csv",
  row.names = FALSE
)

write.csv(
  balance_table_inc_mix,
  "results/outcome/crossfitting/inc_mix/balance_diagnostics.csv",
  row.names = FALSE
)

# Save individual model results
saveRDS(
  results_inc_mix$model1_cf,
  "results/outcome/crossfitting/inc_mix/model1_cf.rds"
)
saveRDS(
  results_inc_mix$model2_cf,
  "results/outcome/crossfitting/inc_mix/model2_cf.rds"
)
saveRDS(
  results_inc_mix$model3_cf,
  "results/outcome/crossfitting/inc_mix/model3_cf.rds"
)
saveRDS(
  results_inc_mix$model4_cf,
  "results/outcome/crossfitting/inc_mix/model4_cf.rds"
)

cat(
  "\nINC_MIX analysis complete. Results saved to: results/outcome/crossfitting/inc_mix/\n\n"
)

# =============================================================================
# ALL ANALYSES COMPLETE
# =============================================================================

cat(
  "\n=============================================================================\n"
)
cat("ALL CROSS-FITTED ANALYSES COMPLETE\n")
cat(
  "=============================================================================\n"
)
cat("\nResults saved to:\n")
cat("  - results/outcome/crossfitting/inc_dec/\n")
cat("  - results/outcome/crossfitting/inc_mix/\n\n")

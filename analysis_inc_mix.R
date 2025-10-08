# Bootstrap Doubly Robust ATT Estimation Analysis - Complete Comparison
# Comparison: Increase vs Mix
# Using refactored utils2.R with unified DR_att() function
#
# This script compares:
# - 4 outcome model specifications (Model 1-4)
# - Reweight estimation approach only
# - Balanced bootstrap sampling
# - All available confidence intervals (norm, basic, perc)

# Setup
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
run = 123
source("utils.R")

# Load data
d = readRDS("data/data_iptw.rds")

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

# Store results
results = list()

# =============================================================================
# REWEIGHT METHOD
# =============================================================================

# -----------------------------------------------------------------------------
# Model 1
# -----------------------------------------------------------------------------

results$model1_reweight = DR_att(
  outcome = "loneliness",
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
  save_to = "results/outcome/inc_mix/figures/model1"
)

# -----------------------------------------------------------------------------
# Model 2
# -----------------------------------------------------------------------------

results$model2_reweight = DR_att(
  outcome = "loneliness",
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
  save_to = "results/outcome/inc_mix/figures/model2"
)

# -----------------------------------------------------------------------------
# Model 3
# -----------------------------------------------------------------------------

results$model3_reweight = DR_att(
  outcome = "loneliness",
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
  save_to = "results/outcome/inc_mix/figures/model3"
)

# -----------------------------------------------------------------------------
# Model 4
# -----------------------------------------------------------------------------

results$model4_reweight = DR_att(
  outcome = "loneliness",
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
  save_to = "results/outcome/inc_mix/figures/model4"
)

# =============================================================================
# COMPARISON PLOTS
# =============================================================================

# Create high-resolution comparison plot
png(
  "results/outcome/inc_mix/figures/model_comparison_distributions.png",
  width = 3000,
  height = 2000,
  res = 300
)
par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))

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

par(mfrow = c(1, 1))
dev.off()

cat(
  "Comparison plot saved to: results/outcome/inc_mix/figures/model_comparison_distributions.png\n"
)

# =============================================================================
# CREATE COMPREHENSIVE SUMMARY TABLE
# =============================================================================

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

# Model 2
summary_rows[[3]] = extract_results(
  results$model2_reweight,
  2,
  "Reweight",
  "AIPW"
)
summary_rows[[4]] = extract_results(
  results$model2_reweight,
  2,
  "Reweight",
  "DRS"
)

# Model 3
summary_rows[[5]] = extract_results(
  results$model3_reweight,
  3,
  "Reweight",
  "AIPW"
)
summary_rows[[6]] = extract_results(
  results$model3_reweight,
  3,
  "Reweight",
  "DRS"
)

# Model 4
summary_rows[[7]] = extract_results(
  results$model4_reweight,
  4,
  "Reweight",
  "AIPW"
)
summary_rows[[8]] = extract_results(
  results$model4_reweight,
  4,
  "Reweight",
  "DRS"
)

summary_table = do.call(rbind, summary_rows)
rownames(summary_table) = NULL


balance_rows = list(
  extract_balance(results$model1_reweight, 1, "Reweight"),
  extract_balance(results$model2_reweight, 2, "Reweight"),
  extract_balance(results$model3_reweight, 3, "Reweight"),
  extract_balance(results$model4_reweight, 4, "Reweight")
)

balance_table = do.call(rbind, balance_rows)
rownames(balance_table) = NULL

# Print summaries

print(summary_table, digits = 4)


print(balance_table, digits = 4)

# =============================================================================
# SAVE RESULTS
# =============================================================================

# Save complete results object
saveRDS(results, "results/outcome/inc_mix/complete_results.rds")
cat("Complete results saved to: results/outcome/inc_mix/complete_results.rds\n")

# Save summary tables
write.csv(
  summary_table,
  "results/outcome/inc_mix/summary_table.csv",
  row.names = FALSE
)


write.csv(
  balance_table,
  "results/outcome/inc_mix/balance_diagnostics.csv",
  row.names = FALSE
)


# Save individual model results
saveRDS(results$model1_reweight, "results/outcome/inc_mix/model1_reweight.rds")
saveRDS(results$model2_reweight, "results/outcome/inc_mix/model2_reweight.rds")
saveRDS(results$model3_reweight, "results/outcome/inc_mix/model3_reweight.rds")
saveRDS(results$model4_reweight, "results/outcome/inc_mix/model4_reweight.rds")

cat("Analysis complete. Results saved to: results/outcome/inc_mix/\n")

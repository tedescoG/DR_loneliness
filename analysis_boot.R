# Bootstrap Doubly Robust ATT Estimation using boot package
# Comparison: Increase vs Decrease
# Model 4 (fully adjusted) - Reweight method

# Setup
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
run = 123
source("utils.R")

# Load data
d = readRDS("data/data_iptw2.rds")

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
outcome_formula_1 = "1"

outcome_formula_2 = "baseline_lone"

outcome_formula_3 = "baseline_lone + baseline_depr"

outcome_formula_4 = "baseline_lone + baseline_depr + female + age_cat + edu +
                     emp_status + income + marital + coliving +
                     health_pre + chronic + death_due_covid +
                     ppl_infected + income_loss  + neighborhood"

# Tuned parameters for GBM - Increase vs Decrease
gbm_inc_dec = readRDS("results/weighting/inc_dec_params.rds")

# Store results in a list
results_boot = list()

# =============================================================================
# Model 1: Ordinary Bootstrap
# =============================================================================

results_boot$model1_ord = DR_att_boot(
  outcome = "severe_loneliness",
  treatment = "remote_contact",
  treated_level = "increase",
  control_level = "decrease",
  f.ps = ps_formula,
  f.out = outcome_formula_1,
  data = d,
  gbm_params = gbm_inc_dec,
  n_boot = 1000,
  seed = run,
  verbose = TRUE,
  parallel = TRUE,
  n_cores = NULL,
  plot_diagnostics = TRUE,
  bootstrap_method = "reweight",
  ci_type = "perc",
  sim = "ordinary"
)

# =============================================================================
# Model 1: Balanced Bootstrap
# =============================================================================

results_boot$model1_bal = DR_att_boot(
  outcome = "severe_loneliness",
  treatment = "remote_contact",
  treated_level = "increase",
  control_level = "decrease",
  f.ps = ps_formula,
  f.out = outcome_formula_1,
  data = d,
  gbm_params = gbm_inc_dec,
  n_boot = 1000,
  seed = run,
  verbose = TRUE,
  parallel = TRUE,
  n_cores = NULL,
  plot_diagnostics = TRUE,
  bootstrap_method = "reweight",
  ci_type = "perc",
  sim = "balanced"
)

# =============================================================================
# Model 2: Ordinary Bootstrap
# =============================================================================

results_boot$model2_ord = DR_att_boot(
  outcome = "severe_loneliness",
  treatment = "remote_contact",
  treated_level = "increase",
  control_level = "decrease",
  f.ps = ps_formula,
  f.out = outcome_formula_2,
  data = d,
  gbm_params = gbm_inc_dec,
  n_boot = 1000,
  seed = run,
  verbose = TRUE,
  parallel = TRUE,
  n_cores = NULL,
  plot_diagnostics = TRUE,
  bootstrap_method = "reweight",
  ci_type = "perc",
  sim = "ordinary"
)

# =============================================================================
# Model 2: Balanced Bootstrap
# =============================================================================

results_boot$model2_bal = DR_att_boot(
  outcome = "severe_loneliness",
  treatment = "remote_contact",
  treated_level = "increase",
  control_level = "decrease",
  f.ps = ps_formula,
  f.out = outcome_formula_2,
  data = d,
  gbm_params = gbm_inc_dec,
  n_boot = 1000,
  seed = run,
  verbose = TRUE,
  parallel = TRUE,
  n_cores = NULL,
  plot_diagnostics = TRUE,
  bootstrap_method = "reweight",
  ci_type = "perc",
  sim = "balanced"
)

# =============================================================================
# Model 3: Ordinary Bootstrap
# =============================================================================

results_boot$model3_ord = DR_att_boot(
  outcome = "severe_loneliness",
  treatment = "remote_contact",
  treated_level = "increase",
  control_level = "decrease",
  f.ps = ps_formula,
  f.out = outcome_formula_3,
  data = d,
  gbm_params = gbm_inc_dec,
  n_boot = 1000,
  seed = run,
  verbose = TRUE,
  parallel = TRUE,
  n_cores = NULL,
  plot_diagnostics = TRUE,
  bootstrap_method = "reweight",
  ci_type = "perc",
  sim = "ordinary"
)

# =============================================================================
# Model 3: Balanced Bootstrap
# =============================================================================

results_boot$model3_bal = DR_att_boot(
  outcome = "severe_loneliness",
  treatment = "remote_contact",
  treated_level = "increase",
  control_level = "decrease",
  f.ps = ps_formula,
  f.out = outcome_formula_3,
  data = d,
  gbm_params = gbm_inc_dec,
  n_boot = 1000,
  seed = run,
  verbose = TRUE,
  parallel = TRUE,
  n_cores = NULL,
  plot_diagnostics = TRUE,
  bootstrap_method = "reweight",
  ci_type = "perc",
  sim = "balanced"
)

# =============================================================================
# Model 4: Ordinary Bootstrap
# =============================================================================

results_boot$model4_ord = DR_att_boot(
  outcome = "severe_loneliness",
  treatment = "remote_contact",
  treated_level = "increase",
  control_level = "decrease",
  f.ps = ps_formula,
  f.out = outcome_formula_4,
  data = d,
  gbm_params = gbm_inc_dec,
  n_boot = 1000,
  seed = run,
  verbose = TRUE,
  parallel = TRUE,
  n_cores = NULL,
  plot_diagnostics = TRUE,
  bootstrap_method = "reweight",
  ci_type = "perc",
  sim = "ordinary"
)

# =============================================================================
# Model 4: Balanced Bootstrap
# =============================================================================

results_boot$model4_bal = DR_att_boot(
  outcome = "severe_loneliness",
  treatment = "remote_contact",
  treated_level = "increase",
  control_level = "decrease",
  f.ps = ps_formula,
  f.out = outcome_formula_4,
  data = d,
  gbm_params = gbm_inc_dec,
  n_boot = 1000,
  seed = run,
  verbose = TRUE,
  parallel = TRUE,
  n_cores = NULL,
  plot_diagnostics = TRUE,
  bootstrap_method = "reweight",
  ci_type = "perc",
  sim = "balanced"
)

# =============================================================================
# Diagnostic plots - Bootstrap estimates across models and schemes
# =============================================================================

# Comparison across models for AIPW - Ordinary vs Balanced
par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))

hist(
  results_boot$model1_ord$aipw$bootstrap_samples,
  main = "AIPW Model 1 - Ordinary",
  xlab = "ATT",
  col = "lightblue",
  breaks = 30
)
abline(v = results_boot$model1_ord$aipw$att, col = "red", lwd = 2, lty = 2)

hist(
  results_boot$model1_bal$aipw$bootstrap_samples,
  main = "AIPW Model 1 - Balanced",
  xlab = "ATT",
  col = "lightcoral",
  breaks = 30
)
abline(v = results_boot$model1_bal$aipw$att, col = "red", lwd = 2, lty = 2)

hist(
  results_boot$model2_ord$aipw$bootstrap_samples,
  main = "AIPW Model 2 - Ordinary",
  xlab = "ATT",
  col = "lightblue",
  breaks = 30
)
abline(v = results_boot$model2_ord$aipw$att, col = "red", lwd = 2, lty = 2)

hist(
  results_boot$model2_bal$aipw$bootstrap_samples,
  main = "AIPW Model 2 - Balanced",
  xlab = "ATT",
  col = "lightcoral",
  breaks = 30
)
abline(v = results_boot$model2_bal$aipw$att, col = "red", lwd = 2, lty = 2)

par(mfrow = c(1, 1))

par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))

hist(
  results_boot$model3_ord$aipw$bootstrap_samples,
  main = "AIPW Model 3 - Ordinary",
  xlab = "ATT",
  col = "lightblue",
  breaks = 30
)
abline(v = results_boot$model3_ord$aipw$att, col = "red", lwd = 2, lty = 2)

hist(
  results_boot$model3_bal$aipw$bootstrap_samples,
  main = "AIPW Model 3 - Balanced",
  xlab = "ATT",
  col = "lightcoral",
  breaks = 30
)
abline(v = results_boot$model3_bal$aipw$att, col = "red", lwd = 2, lty = 2)

hist(
  results_boot$model4_ord$aipw$bootstrap_samples,
  main = "AIPW Model 4 - Ordinary",
  xlab = "ATT",
  col = "lightblue",
  breaks = 30
)
abline(v = results_boot$model4_ord$aipw$att, col = "red", lwd = 2, lty = 2)

hist(
  results_boot$model4_bal$aipw$bootstrap_samples,
  main = "AIPW Model 4 - Balanced",
  xlab = "ATT",
  col = "lightcoral",
  breaks = 30
)
abline(v = results_boot$model4_bal$aipw$att, col = "red", lwd = 2, lty = 2)

par(mfrow = c(1, 1))

# Comparison across models for DRS - Ordinary vs Balanced
par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))

hist(
  results_boot$model1_ord$drs$bootstrap_samples,
  main = "DRS Model 1 - Ordinary",
  xlab = "ATT",
  col = "lightgreen",
  breaks = 30
)
abline(v = results_boot$model1_ord$drs$att, col = "red", lwd = 2, lty = 2)

hist(
  results_boot$model1_bal$drs$bootstrap_samples,
  main = "DRS Model 1 - Balanced",
  xlab = "ATT",
  col = "lightyellow",
  breaks = 30
)
abline(v = results_boot$model1_bal$drs$att, col = "red", lwd = 2, lty = 2)

hist(
  results_boot$model2_ord$drs$bootstrap_samples,
  main = "DRS Model 2 - Ordinary",
  xlab = "ATT",
  col = "lightgreen",
  breaks = 30
)
abline(v = results_boot$model2_ord$drs$att, col = "red", lwd = 2, lty = 2)

hist(
  results_boot$model2_bal$drs$bootstrap_samples,
  main = "DRS Model 2 - Balanced",
  xlab = "ATT",
  col = "lightyellow",
  breaks = 30
)
abline(v = results_boot$model2_bal$drs$att, col = "red", lwd = 2, lty = 2)

par(mfrow = c(1, 1))

par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))

hist(
  results_boot$model3_ord$drs$bootstrap_samples,
  main = "DRS Model 3 - Ordinary",
  xlab = "ATT",
  col = "lightgreen",
  breaks = 30
)
abline(v = results_boot$model3_ord$drs$att, col = "red", lwd = 2, lty = 2)

hist(
  results_boot$model3_bal$drs$bootstrap_samples,
  main = "DRS Model 3 - Balanced",
  xlab = "ATT",
  col = "lightyellow",
  breaks = 30
)
abline(v = results_boot$model3_bal$drs$att, col = "red", lwd = 2, lty = 2)

hist(
  results_boot$model4_ord$drs$bootstrap_samples,
  main = "DRS Model 4 - Ordinary",
  xlab = "ATT",
  col = "lightgreen",
  breaks = 30
)
abline(v = results_boot$model4_ord$drs$att, col = "red", lwd = 2, lty = 2)

hist(
  results_boot$model4_bal$drs$bootstrap_samples,
  main = "DRS Model 4 - Balanced",
  xlab = "ATT",
  col = "lightyellow",
  breaks = 30
)
abline(v = results_boot$model4_bal$drs$att, col = "red", lwd = 2, lty = 2)

par(mfrow = c(1, 1))

# Overlay densities for AIPW across models - Ordinary Bootstrap
plot(
  density(results_boot$model1_ord$aipw$bootstrap_samples),
  main = "AIPW Estimates Across Models - Ordinary Bootstrap",
  xlab = "ATT",
  col = "blue",
  lwd = 2,
  xlim = range(c(
    results_boot$model1_ord$aipw$bootstrap_samples,
    results_boot$model2_ord$aipw$bootstrap_samples,
    results_boot$model3_ord$aipw$bootstrap_samples,
    results_boot$model4_ord$aipw$bootstrap_samples
  ))
)
lines(
  density(results_boot$model2_ord$aipw$bootstrap_samples),
  col = "red",
  lwd = 2
)
lines(
  density(results_boot$model3_ord$aipw$bootstrap_samples),
  col = "green",
  lwd = 2
)
lines(
  density(results_boot$model4_ord$aipw$bootstrap_samples),
  col = "purple",
  lwd = 2
)
legend(
  "topleft",
  legend = c("Model 1", "Model 2", "Model 3", "Model 4"),
  col = c("blue", "red", "green", "purple"),
  lwd = 2
)

# Overlay densities for AIPW across models - Balanced Bootstrap
plot(
  density(results_boot$model1_bal$aipw$bootstrap_samples),
  main = "AIPW Estimates Across Models - Balanced Bootstrap",
  xlab = "ATT",
  col = "blue",
  lwd = 2,
  xlim = range(c(
    results_boot$model1_bal$aipw$bootstrap_samples,
    results_boot$model2_bal$aipw$bootstrap_samples,
    results_boot$model3_bal$aipw$bootstrap_samples,
    results_boot$model4_bal$aipw$bootstrap_samples
  ))
)
lines(
  density(results_boot$model2_bal$aipw$bootstrap_samples),
  col = "red",
  lwd = 2
)
lines(
  density(results_boot$model3_bal$aipw$bootstrap_samples),
  col = "green",
  lwd = 2
)
lines(
  density(results_boot$model4_bal$aipw$bootstrap_samples),
  col = "purple",
  lwd = 2
)
legend(
  "topleft",
  legend = c("Model 1", "Model 2", "Model 3", "Model 4"),
  col = c("blue", "red", "green", "purple"),
  lwd = 2
)

# Overlay densities for DRS across models - Ordinary Bootstrap
plot(
  density(results_boot$model1_ord$drs$bootstrap_samples),
  main = "DRS Estimates Across Models - Ordinary Bootstrap",
  xlab = "ATT",
  col = "blue",
  lwd = 2,
  xlim = range(c(
    results_boot$model1_ord$drs$bootstrap_samples,
    results_boot$model2_ord$drs$bootstrap_samples,
    results_boot$model3_ord$drs$bootstrap_samples,
    results_boot$model4_ord$drs$bootstrap_samples
  ))
)
lines(
  density(results_boot$model2_ord$drs$bootstrap_samples),
  col = "red",
  lwd = 2
)
lines(
  density(results_boot$model3_ord$drs$bootstrap_samples),
  col = "green",
  lwd = 2
)
lines(
  density(results_boot$model4_ord$drs$bootstrap_samples),
  col = "purple",
  lwd = 2
)
legend(
  "topleft",
  legend = c("Model 1", "Model 2", "Model 3", "Model 4"),
  col = c("blue", "red", "green", "purple"),
  lwd = 2
)

# Overlay densities for DRS across models - Balanced Bootstrap
plot(
  density(results_boot$model1_bal$drs$bootstrap_samples),
  main = "DRS Estimates Across Models - Balanced Bootstrap",
  xlab = "ATT",
  col = "blue",
  lwd = 2,
  xlim = range(c(
    results_boot$model1_bal$drs$bootstrap_samples,
    results_boot$model2_bal$drs$bootstrap_samples,
    results_boot$model3_bal$drs$bootstrap_samples,
    results_boot$model4_bal$drs$bootstrap_samples
  ))
)
lines(
  density(results_boot$model2_bal$drs$bootstrap_samples),
  col = "red",
  lwd = 2
)
lines(
  density(results_boot$model3_bal$drs$bootstrap_samples),
  col = "green",
  lwd = 2
)
lines(
  density(results_boot$model4_bal$drs$bootstrap_samples),
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

# Create output directory if it doesn't exist
if (!dir.exists("results/outcome/boot")) {
  dir.create("results/outcome/boot", recursive = TRUE)
}

# Save all individual results
saveRDS(results_boot$model1_ord, "results/outcome/boot/model1_ord.rds")
saveRDS(results_boot$model1_bal, "results/outcome/boot/model1_bal.rds")
saveRDS(results_boot$model2_ord, "results/outcome/boot/model2_ord.rds")
saveRDS(results_boot$model2_bal, "results/outcome/boot/model2_bal.rds")
saveRDS(results_boot$model3_ord, "results/outcome/boot/model3_ord.rds")
saveRDS(results_boot$model3_bal, "results/outcome/boot/model3_bal.rds")
saveRDS(results_boot$model4_ord, "results/outcome/boot/model4_ord.rds")
saveRDS(results_boot$model4_bal, "results/outcome/boot/model4_bal.rds")

# Save entire results list
saveRDS(results_boot, "results/outcome/boot/all_results_boot.rds")

# Create comprehensive summary table
summary_table_boot = data.frame(
  Model = rep(paste0("Model ", 1:4), each = 2),
  Bootstrap_Type = rep(c("Ordinary", "Balanced"), 4),
  AIPW_ATT = c(
    results_boot$model1_ord$aipw$att,
    results_boot$model1_bal$aipw$att,
    results_boot$model2_ord$aipw$att,
    results_boot$model2_bal$aipw$att,
    results_boot$model3_ord$aipw$att,
    results_boot$model3_bal$aipw$att,
    results_boot$model4_ord$aipw$att,
    results_boot$model4_bal$aipw$att
  ),
  AIPW_SE = c(
    results_boot$model1_ord$aipw$se,
    results_boot$model1_bal$aipw$se,
    results_boot$model2_ord$aipw$se,
    results_boot$model2_bal$aipw$se,
    results_boot$model3_ord$aipw$se,
    results_boot$model3_bal$aipw$se,
    results_boot$model4_ord$aipw$se,
    results_boot$model4_bal$aipw$se
  ),
  AIPW_CI_lower = c(
    results_boot$model1_ord$aipw$ci_percentile[1],
    results_boot$model1_bal$aipw$ci_percentile[1],
    results_boot$model2_ord$aipw$ci_percentile[1],
    results_boot$model2_bal$aipw$ci_percentile[1],
    results_boot$model3_ord$aipw$ci_percentile[1],
    results_boot$model3_bal$aipw$ci_percentile[1],
    results_boot$model4_ord$aipw$ci_percentile[1],
    results_boot$model4_bal$aipw$ci_percentile[1]
  ),
  AIPW_CI_upper = c(
    results_boot$model1_ord$aipw$ci_percentile[2],
    results_boot$model1_bal$aipw$ci_percentile[2],
    results_boot$model2_ord$aipw$ci_percentile[2],
    results_boot$model2_bal$aipw$ci_percentile[2],
    results_boot$model3_ord$aipw$ci_percentile[2],
    results_boot$model3_bal$aipw$ci_percentile[2],
    results_boot$model4_ord$aipw$ci_percentile[2],
    results_boot$model4_bal$aipw$ci_percentile[2]
  ),
  AIPW_pval = c(
    results_boot$model1_ord$aipw$pval,
    results_boot$model1_bal$aipw$pval,
    results_boot$model2_ord$aipw$pval,
    results_boot$model2_bal$aipw$pval,
    results_boot$model3_ord$aipw$pval,
    results_boot$model3_bal$aipw$pval,
    results_boot$model4_ord$aipw$pval,
    results_boot$model4_bal$aipw$pval
  ),
  DRS_ATT = c(
    results_boot$model1_ord$drs$att,
    results_boot$model1_bal$drs$att,
    results_boot$model2_ord$drs$att,
    results_boot$model2_bal$drs$att,
    results_boot$model3_ord$drs$att,
    results_boot$model3_bal$drs$att,
    results_boot$model4_ord$drs$att,
    results_boot$model4_bal$drs$att
  ),
  DRS_SE = c(
    results_boot$model1_ord$drs$se,
    results_boot$model1_bal$drs$se,
    results_boot$model2_ord$drs$se,
    results_boot$model2_bal$drs$se,
    results_boot$model3_ord$drs$se,
    results_boot$model3_bal$drs$se,
    results_boot$model4_ord$drs$se,
    results_boot$model4_bal$drs$se
  ),
  DRS_CI_lower = c(
    results_boot$model1_ord$drs$ci_percentile[1],
    results_boot$model1_bal$drs$ci_percentile[1],
    results_boot$model2_ord$drs$ci_percentile[1],
    results_boot$model2_bal$drs$ci_percentile[1],
    results_boot$model3_ord$drs$ci_percentile[1],
    results_boot$model3_bal$drs$ci_percentile[1],
    results_boot$model4_ord$drs$ci_percentile[1],
    results_boot$model4_bal$drs$ci_percentile[1]
  ),
  DRS_CI_upper = c(
    results_boot$model1_ord$drs$ci_percentile[2],
    results_boot$model1_bal$drs$ci_percentile[2],
    results_boot$model2_ord$drs$ci_percentile[2],
    results_boot$model2_bal$drs$ci_percentile[2],
    results_boot$model3_ord$drs$ci_percentile[2],
    results_boot$model3_bal$drs$ci_percentile[2],
    results_boot$model4_ord$drs$ci_percentile[2],
    results_boot$model4_bal$drs$ci_percentile[2]
  ),
  DRS_pval = c(
    results_boot$model1_ord$drs$pval,
    results_boot$model1_bal$drs$pval,
    results_boot$model2_ord$drs$pval,
    results_boot$model2_bal$drs$pval,
    results_boot$model3_ord$drs$pval,
    results_boot$model3_bal$drs$pval,
    results_boot$model4_ord$drs$pval,
    results_boot$model4_bal$drs$pval
  ),
  N_boot = c(
    results_boot$model1_ord$n_boot,
    results_boot$model1_bal$n_boot,
    results_boot$model2_ord$n_boot,
    results_boot$model2_bal$n_boot,
    results_boot$model3_ord$n_boot,
    results_boot$model3_bal$n_boot,
    results_boot$model4_ord$n_boot,
    results_boot$model4_bal$n_boot
  )
)

# Round numeric columns for display
summary_table_boot[, 3:ncol(summary_table_boot)] =
  round(summary_table_boot[, 3:ncol(summary_table_boot)], 4)

cat("\n\nBootstrap Summary Table:\n")
print(summary_table_boot)

# Save as CSV
write.csv(
  summary_table_boot,
  "results/outcome/boot/bootstrap_summary.csv",
  row.names = FALSE
)

cat("\n\nResults saved successfully!\n")

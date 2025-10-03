# Bootstrap Doubly Robust ATT Estimation Analysis
# This script runs DR_att() with 4 different outcome models

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

outcome_formula_2 = "baseline_lone "

outcome_formula_3 = "baseline_lone + baseline_depr "

outcome_formula_4 = "baseline_lone + baseline_depr + female + age_cat + edu +
                     emp_status + income + marital + coliving + change_res +
                     kinless + health_pre + chronic + death_due_covid +
                     ppl_infected + income_loss + job_loss + neighborhood"


# Tuned parameters for GBM
gbm_inc_dec = readRDS("results/weighting/inc_dec_params.rds")
gbm_inc_mix = readRDS("results/weighting/inc_mix_params.rds")

# Store results in a list
results = list()

# =============================================================================
# Model 1: Increase vs Decrease
# =============================================================================

results$model1_inc_dec = DR_att(
  outcome = "severe_loneliness",
  treatment = "remote_contact",
  treated_level = "increase",
  control_level = "decrease",
  f.ps = ps_formula,
  f.out = outcome_formula_1,
  data = d,
  ps_params = gbm_inc_dec,
  n_boot = 500,
  seed = run,
  verbose = TRUE,
  parallel = TRUE,
  n_cores = NULL,
  plot_diagnostics = TRUE
)

# =============================================================================
# Model 1: Increase vs Mix
# =============================================================================

results$model1_inc_mix = DR_att(
  outcome = "severe_loneliness",
  treatment = "remote_contact",
  treated_level = "increase",
  control_level = "mix",
  f.ps = ps_formula,
  f.out = outcome_formula_1,
  data = d,
  ps_params = gbm_inc_mix,
  n_boot = 500,
  seed = run,
  verbose = TRUE,
  parallel = TRUE,
  n_cores = NULL,
  plot_diagnostics = TRUE
)

# =============================================================================
# Model 2: Increase vs Decrease
# =============================================================================

results$model2_inc_dec = DR_att(
  outcome = "severe_loneliness",
  treatment = "remote_contact",
  treated_level = "increase",
  control_level = "decrease",
  f.ps = ps_formula,
  f.out = outcome_formula_2,
  data = d,
  ps_params = gbm_inc_dec,
  n_boot = 500,
  seed = run,
  verbose = TRUE,
  parallel = TRUE,
  n_cores = NULL,
  plot_diagnostics = TRUE
)

# =============================================================================
# Model 2: Increase vs Mix
# =============================================================================

results$model2_inc_mix = DR_att(
  outcome = "severe_loneliness",
  treatment = "remote_contact",
  treated_level = "increase",
  control_level = "mix",
  f.ps = ps_formula,
  f.out = outcome_formula_2,
  data = d,
  ps_params = gbm_inc_mix,
  n_boot = 500,
  seed = run,
  verbose = TRUE,
  parallel = TRUE,
  n_cores = NULL,
  plot_diagnostics = TRUE
)

# =============================================================================
# Model 3: Increase vs Decrease
# =============================================================================

results$model3_inc_dec = DR_att(
  outcome = "severe_loneliness",
  treatment = "remote_contact",
  treated_level = "increase",
  control_level = "decrease",
  f.ps = ps_formula,
  f.out = outcome_formula_3,
  data = d,
  ps_params = gbm_inc_dec,
  n_boot = 500,
  seed = run,
  verbose = TRUE,
  parallel = TRUE,
  n_cores = NULL,
  plot_diagnostics = TRUE
)

# =============================================================================
# Model 3: Increase vs Mix
# =============================================================================

results$model3_inc_mix = DR_att(
  outcome = "severe_loneliness",
  treatment = "remote_contact",
  treated_level = "increase",
  control_level = "mix",
  f.ps = ps_formula,
  f.out = outcome_formula_3,
  data = d,
  ps_params = gbm_inc_mix,
  n_boot = 500,
  seed = run,
  verbose = TRUE,
  parallel = TRUE,
  n_cores = NULL,
  plot_diagnostics = TRUE
)

# =============================================================================
# Model 4: Increase vs Decrease
# =============================================================================

results$model4_inc_dec = DR_att(
  outcome = "severe_loneliness",
  treatment = "remote_contact",
  treated_level = "increase",
  control_level = "decrease",
  f.ps = ps_formula,
  f.out = outcome_formula_4,
  data = d,
  ps_params = gbm_inc_dec,
  n_boot = 500,
  seed = run,
  verbose = TRUE,
  parallel = TRUE,
  n_cores = NULL,
  plot_diagnostics = TRUE
)

# =============================================================================
# Model 4: Increase vs Mix
# =============================================================================

results$model4_inc_mix = DR_att(
  outcome = "severe_loneliness",
  treatment = "remote_contact",
  treated_level = "increase",
  control_level = "mix",
  f.ps = ps_formula,
  f.out = outcome_formula_4,
  data = d,
  ps_params = gbm_inc_mix,
  n_boot = 500,
  seed = run,
  verbose = TRUE,
  parallel = TRUE,
  n_cores = NULL,
  plot_diagnostics = TRUE
)

# =============================================================================
# Diagnostic plots
# =============================================================================

# Comparison plots for Model 1 - Increase vs Decrease
par(mfrow = c(1, 2), mar = c(4, 4, 2, 1))

hist(
  results$model1_inc_dec$aipw$bootstrap_samples,
  main = "AIPW - Model 1",
  xlab = "ATT",
  col = "lightblue",
  breaks = 30,
  xlim = range(c(
    results$model1_inc_dec$aipw$bootstrap_samples,
    results$model1_inc_dec$drs$bootstrap_samples
  ))
)

hist(
  results$model1_inc_dec$drs$bootstrap_samples,
  main = "DRS - Model 1",
  xlab = "ATT",
  col = "lightgreen",
  breaks = 30,
  xlim = range(c(
    results$model1_inc_dec$aipw$bootstrap_samples,
    results$model1_inc_dec$drs$bootstrap_samples
  ))
)

par(mfrow = c(1, 1))

# Overlay distributions for Model 1
plot(
  density(results$model1_inc_dec$aipw$bootstrap_samples),
  main = "Model 1: AIPW vs DRS Bootstrap Distributions",
  xlab = "ATT",
  col = "blue",
  lwd = 2,
  xlim = range(c(
    results$model1_inc_dec$aipw$bootstrap_samples,
    results$model1_inc_dec$drs$bootstrap_samples
  ))
)
lines(
  density(results$model1_inc_dec$drs$bootstrap_samples),
  col = "green",
  lwd = 2
)
legend("topright", legend = c("AIPW", "DRS"), col = c("blue", "green"), lwd = 2)

# Overlay distributions for Model 2
plot(
  density(results$model2_inc_dec$aipw$bootstrap_samples),
  main = "Model 2: AIPW vs DRS Bootstrap Distributions",
  xlab = "ATT",
  col = "blue",
  lwd = 2,
  xlim = range(c(
    results$model2_inc_dec$aipw$bootstrap_samples,
    results$model2_inc_dec$drs$bootstrap_samples
  ))
)
lines(
  density(results$model2_inc_dec$drs$bootstrap_samples),
  col = "green",
  lwd = 2
)
legend("topright", legend = c("AIPW", "DRS"), col = c("blue", "green"), lwd = 2)

# Overlay distributions for Model 3
plot(
  density(results$model3_inc_dec$aipw$bootstrap_samples),
  main = "Model 3: AIPW vs DRS Bootstrap Distributions",
  xlab = "ATT",
  col = "blue",
  lwd = 2,
  xlim = range(c(
    results$model3_inc_dec$aipw$bootstrap_samples,
    results$model3_inc_dec$drs$bootstrap_samples
  ))
)
lines(
  density(results$model3_inc_dec$drs$bootstrap_samples),
  col = "green",
  lwd = 2
)
legend("topright", legend = c("AIPW", "DRS"), col = c("blue", "green"), lwd = 2)

# Overlay distributions for Model 4
plot(
  density(results$model4_inc_dec$aipw$bootstrap_samples),
  main = "Model 4: AIPW vs DRS Bootstrap Distributions",
  xlab = "ATT",
  col = "blue",
  lwd = 2,
  xlim = range(c(
    results$model4_inc_dec$aipw$bootstrap_samples,
    results$model4_inc_dec$drs$bootstrap_samples
  ))
)
lines(
  density(results$model4_inc_dec$drs$bootstrap_samples),
  col = "green",
  lwd = 2
)
legend("topright", legend = c("AIPW", "DRS"), col = c("blue", "green"), lwd = 2)

# Comparison across models for AIPW (Increase vs Decrease)
plot(
  density(results$model1_inc_dec$aipw$bootstrap_samples),
  main = "AIPW Estimates Across Models (Inc vs Dec)",
  xlab = "ATT",
  col = "blue",
  lwd = 2,
  xlim = range(c(
    results$model1_inc_dec$aipw$bootstrap_samples,
    results$model2_inc_dec$aipw$bootstrap_samples,
    results$model3_inc_dec$aipw$bootstrap_samples,
    results$model4_inc_dec$aipw$bootstrap_samples
  ))
)
lines(
  density(results$model2_inc_dec$aipw$bootstrap_samples),
  col = "red",
  lwd = 2
)
lines(
  density(results$model3_inc_dec$aipw$bootstrap_samples),
  col = "green",
  lwd = 2
)
lines(
  density(results$model4_inc_dec$aipw$bootstrap_samples),
  col = "purple",
  lwd = 2
)
legend(
  "topleft",
  legend = c("Model 1", "Model 2", "Model 3", "Model 4"),
  col = c("blue", "red", "green", "purple"),
  lwd = 2
)

# Comparison across models for DRS (Increase vs Decrease)
plot(
  density(results$model1_inc_dec$drs$bootstrap_samples),
  main = "DRS Estimates Across Models (Inc vs Dec)",
  xlab = "ATT",
  col = "blue",
  lwd = 2,
  xlim = range(c(
    results$model1_inc_dec$drs$bootstrap_samples,
    results$model2_inc_dec$drs$bootstrap_samples,
    results$model3_inc_dec$drs$bootstrap_samples,
    results$model4_inc_dec$drs$bootstrap_samples
  ))
)
lines(
  density(results$model2_inc_dec$drs$bootstrap_samples),
  col = "red",
  lwd = 2
)
lines(
  density(results$model3_inc_dec$drs$bootstrap_samples),
  col = "green",
  lwd = 2
)
lines(
  density(results$model4_inc_dec$drs$bootstrap_samples),
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

cat("\n\nSAVING RESULTS:\n")
cat("---------------\n")

# Save all individual results
saveRDS(results$model1_inc_dec, "results/model1_inc_dec.rds")
saveRDS(results$model1_inc_mix, "results/model1_inc_mix.rds")
saveRDS(results$model2_inc_dec, "results/model2_inc_dec.rds")
saveRDS(results$model2_inc_mix, "results/model2_inc_mix.rds")
saveRDS(results$model3_inc_dec, "results/model3_inc_dec.rds")
saveRDS(results$model3_inc_mix, "results/model3_inc_mix.rds")
saveRDS(results$model4_inc_dec, "results/model4_inc_dec.rds")
saveRDS(results$model4_inc_mix, "results/model4_inc_mix.rds")

# Save entire results list
saveRDS(results, "results/all_results.rds")

# Create comprehensive summary table
summary_table = data.frame(
  Model = rep(paste0("Model ", 1:4), each = 2),
  Comparison = rep(c("Increase vs Decrease", "Increase vs Mix"), 4),
  AIPW_ATT = c(
    results$model1_inc_dec$aipw$att,
    results$model1_inc_mix$aipw$att,
    results$model2_inc_dec$aipw$att,
    results$model2_inc_mix$aipw$att,
    results$model3_inc_dec$aipw$att,
    results$model3_inc_mix$aipw$att,
    results$model4_inc_dec$aipw$att,
    results$model4_inc_mix$aipw$att
  ),
  AIPW_SE = c(
    results$model1_inc_dec$aipw$se,
    results$model1_inc_mix$aipw$se,
    results$model2_inc_dec$aipw$se,
    results$model2_inc_mix$aipw$se,
    results$model3_inc_dec$aipw$se,
    results$model3_inc_mix$aipw$se,
    results$model4_inc_dec$aipw$se,
    results$model4_inc_mix$aipw$se
  ),
  AIPW_CI_lower = c(
    results$model1_inc_dec$aipw$ci_percentile[1],
    results$model1_inc_mix$aipw$ci_percentile[1],
    results$model2_inc_dec$aipw$ci_percentile[1],
    results$model2_inc_mix$aipw$ci_percentile[1],
    results$model3_inc_dec$aipw$ci_percentile[1],
    results$model3_inc_mix$aipw$ci_percentile[1],
    results$model4_inc_dec$aipw$ci_percentile[1],
    results$model4_inc_mix$aipw$ci_percentile[1]
  ),
  AIPW_CI_upper = c(
    results$model1_inc_dec$aipw$ci_percentile[2],
    results$model1_inc_mix$aipw$ci_percentile[2],
    results$model2_inc_dec$aipw$ci_percentile[2],
    results$model2_inc_mix$aipw$ci_percentile[2],
    results$model3_inc_dec$aipw$ci_percentile[2],
    results$model3_inc_mix$aipw$ci_percentile[2],
    results$model4_inc_dec$aipw$ci_percentile[2],
    results$model4_inc_mix$aipw$ci_percentile[2]
  ),
  AIPW_pval = c(
    results$model1_inc_dec$aipw$pval,
    results$model1_inc_mix$aipw$pval,
    results$model2_inc_dec$aipw$pval,
    results$model2_inc_mix$aipw$pval,
    results$model3_inc_dec$aipw$pval,
    results$model3_inc_mix$aipw$pval,
    results$model4_inc_dec$aipw$pval,
    results$model4_inc_mix$aipw$pval
  ),
  DRS_ATT = c(
    results$model1_inc_dec$drs$att,
    results$model1_inc_mix$drs$att,
    results$model2_inc_dec$drs$att,
    results$model2_inc_mix$drs$att,
    results$model3_inc_dec$drs$att,
    results$model3_inc_mix$drs$att,
    results$model4_inc_dec$drs$att,
    results$model4_inc_mix$drs$att
  ),
  DRS_SE = c(
    results$model1_inc_dec$drs$se,
    results$model1_inc_mix$drs$se,
    results$model2_inc_dec$drs$se,
    results$model2_inc_mix$drs$se,
    results$model3_inc_dec$drs$se,
    results$model3_inc_mix$drs$se,
    results$model4_inc_dec$drs$se,
    results$model4_inc_mix$drs$se
  ),
  DRS_CI_lower = c(
    results$model1_inc_dec$drs$ci_percentile[1],
    results$model1_inc_mix$drs$ci_percentile[1],
    results$model2_inc_dec$drs$ci_percentile[1],
    results$model2_inc_mix$drs$ci_percentile[1],
    results$model3_inc_dec$drs$ci_percentile[1],
    results$model3_inc_mix$drs$ci_percentile[1],
    results$model4_inc_dec$drs$ci_percentile[1],
    results$model4_inc_mix$drs$ci_percentile[1]
  ),
  DRS_CI_upper = c(
    results$model1_inc_dec$drs$ci_percentile[2],
    results$model1_inc_mix$drs$ci_percentile[2],
    results$model2_inc_dec$drs$ci_percentile[2],
    results$model2_inc_mix$drs$ci_percentile[2],
    results$model3_inc_dec$drs$ci_percentile[2],
    results$model3_inc_mix$drs$ci_percentile[2],
    results$model4_inc_dec$drs$ci_percentile[2],
    results$model4_inc_mix$drs$ci_percentile[2]
  ),
  DRS_pval = c(
    results$model1_inc_dec$drs$pval,
    results$model1_inc_mix$drs$pval,
    results$model2_inc_dec$drs$pval,
    results$model2_inc_mix$drs$pval,
    results$model3_inc_dec$drs$pval,
    results$model3_inc_mix$drs$pval,
    results$model4_inc_dec$drs$pval,
    results$model4_inc_mix$drs$pval
  )
)

print(summary_table)

# Save as CSV
write.csv(summary_table, "results/bootstrap_DR_summary.csv", row.names = FALSE)

cat("\nResults saved successfully!\n")
cat("- Individual results: results/model*_inc_*.rds\n")
cat("- All results: results/all_results.rds\n")
cat("- Summary table: results/bootstrap_DR_summary.csv\n")

# Bootstrap Doubly Robust ATT Estimation Analysis
# Comparison: Increase vs Decrease
# This script runs DR_att() with 4 different outcome models
# Each model is tested with both "resample" and "reweight" bootstrap methods

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


# Tuned parameters for GBM - Increase vs Decrease
gbm_inc_dec = readRDS("results/weighting/inc_dec_params.rds")

# Store results in a list
results = list()

# =============================================================================
# Model 1: Resample Method
# =============================================================================

results$model1_resample = DR_att(
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
  parallel = T,
  n_cores = NULL,
  plot_diagnostics = TRUE,
  bootstrap_method = "resample"
)

# =============================================================================
# Model 1: Reweight Method
# =============================================================================

results$model1_reweight = DR_att(
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
  plot_diagnostics = TRUE,
  bootstrap_method = "reweight"
)

# =============================================================================
# Model 2: Resample Method
# =============================================================================

results$model2_resample = DR_att(
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
  plot_diagnostics = TRUE,
  bootstrap_method = "resample"
)

# =============================================================================
# Model 2: Reweight Method
# =============================================================================

results$model2_reweight = DR_att(
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
  plot_diagnostics = TRUE,
  bootstrap_method = "reweight"
)

# =============================================================================
# Model 3: Resample Method
# =============================================================================

results$model3_resample = DR_att(
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
  plot_diagnostics = TRUE,
  bootstrap_method = "resample"
)

# =============================================================================
# Model 3: Reweight Method
# =============================================================================

results$model3_reweight = DR_att(
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
  plot_diagnostics = TRUE,
  bootstrap_method = "reweight"
)

# =============================================================================
# Model 4: Resample Method
# =============================================================================

results$model4_resample = DR_att(
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
  plot_diagnostics = TRUE,
  bootstrap_method = "resample"
)

# =============================================================================
# Model 4: Reweight Method
# =============================================================================

results$model4_reweight = DR_att(
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
  plot_diagnostics = TRUE,
  bootstrap_method = "reweight"
)

# =============================================================================
# Diagnostic plots - Increase vs Decrease
# =============================================================================

# Comparison plots for Model 1 - Resample vs Reweight
par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))

hist(
  results$model1_resample$aipw$bootstrap_samples,
  main = "AIPW Model 1 (Resample) - Inc vs Dec",
  xlab = "ATT",
  col = "lightblue",
  breaks = 30
)

hist(
  results$model1_reweight$aipw$bootstrap_samples,
  main = "AIPW Model 1 (Reweight) - Inc vs Dec",
  xlab = "ATT",
  col = "lightcoral",
  breaks = 30
)

hist(
  results$model1_resample$drs$bootstrap_samples,
  main = "DRS Model 1 (Resample) - Inc vs Dec",
  xlab = "ATT",
  col = "lightgreen",
  breaks = 30
)

hist(
  results$model1_reweight$drs$bootstrap_samples,
  main = "DRS Model 1 (Reweight) - Inc vs Dec",
  xlab = "ATT",
  col = "lightyellow",
  breaks = 30
)

par(mfrow = c(1, 1))

# Overlay distributions for Model 1: AIPW Resample vs Reweight
plot(
  density(results$model1_resample$aipw$bootstrap_samples),
  main = "Model 1 AIPW: Resample vs Reweight - Inc vs Dec",
  xlab = "ATT",
  col = "blue",
  lwd = 2,
  xlim = range(c(
    results$model1_resample$aipw$bootstrap_samples,
    results$model1_reweight$aipw$bootstrap_samples
  ))
)
lines(
  density(results$model1_reweight$aipw$bootstrap_samples),
  col = "red",
  lwd = 2
)
legend(
  "topright",
  legend = c("Resample", "Reweight"),
  col = c("blue", "red"),
  lwd = 2
)

# Overlay distributions for Model 1: DRS Resample vs Reweight
plot(
  density(results$model1_resample$drs$bootstrap_samples),
  main = "Model 1 DRS: Resample vs Reweight - Inc vs Dec",
  xlab = "ATT",
  col = "blue",
  lwd = 2,
  xlim = range(c(
    results$model1_resample$drs$bootstrap_samples,
    results$model1_reweight$drs$bootstrap_samples
  ))
)
lines(
  density(results$model1_reweight$drs$bootstrap_samples),
  col = "red",
  lwd = 2
)
legend(
  "topright",
  legend = c("Resample", "Reweight"),
  col = c("blue", "red"),
  lwd = 2
)

# Comparison across models for AIPW (Resample) - Increase vs Decrease
plot(
  density(results$model1_resample$aipw$bootstrap_samples),
  main = "AIPW Estimates Across Models (Resample) - Inc vs Dec",
  xlab = "ATT",
  col = "blue",
  lwd = 2,
  xlim = range(c(
    results$model1_resample$aipw$bootstrap_samples,
    results$model2_resample$aipw$bootstrap_samples,
    results$model3_resample$aipw$bootstrap_samples,
    results$model4_resample$aipw$bootstrap_samples
  ))
)
lines(
  density(results$model2_resample$aipw$bootstrap_samples),
  col = "red",
  lwd = 2
)
lines(
  density(results$model3_resample$aipw$bootstrap_samples),
  col = "green",
  lwd = 2
)
lines(
  density(results$model4_resample$aipw$bootstrap_samples),
  col = "purple",
  lwd = 2
)
legend(
  "topleft",
  legend = c("Model 1", "Model 2", "Model 3", "Model 4"),
  col = c("blue", "red", "green", "purple"),
  lwd = 2
)

# Comparison across models for AIPW (Reweight) - Increase vs Decrease
plot(
  density(results$model1_reweight$aipw$bootstrap_samples),
  main = "AIPW Estimates Across Models (Reweight) - Inc vs Dec",
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
  lwd = 2
)

# Comparison across models for DRS (Resample) - Increase vs Decrease
plot(
  density(results$model1_resample$drs$bootstrap_samples),
  main = "DRS Estimates Across Models (Resample) - Inc vs Dec",
  xlab = "ATT",
  col = "blue",
  lwd = 2,
  xlim = range(c(
    results$model1_resample$drs$bootstrap_samples,
    results$model2_resample$drs$bootstrap_samples,
    results$model3_resample$drs$bootstrap_samples,
    results$model4_resample$drs$bootstrap_samples
  ))
)
lines(
  density(results$model2_resample$drs$bootstrap_samples),
  col = "red",
  lwd = 2
)
lines(
  density(results$model3_resample$drs$bootstrap_samples),
  col = "green",
  lwd = 2
)
lines(
  density(results$model4_resample$drs$bootstrap_samples),
  col = "purple",
  lwd = 2
)
legend(
  "topleft",
  legend = c("Model 1", "Model 2", "Model 3", "Model 4"),
  col = c("blue", "red", "green", "purple"),
  lwd = 2
)

# Comparison across models for DRS (Reweight) - Increase vs Decrease
plot(
  density(results$model1_reweight$drs$bootstrap_samples),
  main = "DRS Estimates Across Models (Reweight) - Inc vs Dec",
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
  lwd = 2
)

# =============================================================================
# Save results
# =============================================================================

cat("\n\nSAVING RESULTS:\n")
cat("---------------\n")

# Save all individual results
saveRDS(results$model1_resample, "results/outcome/inc_dec/model1_resample.rds")
saveRDS(results$model1_reweight, "results/outcome/inc_dec/model1_reweight.rds")
saveRDS(results$model2_resample, "results/outcome/inc_dec/model2_resample.rds")
saveRDS(results$model2_reweight, "results/outcome/inc_dec/model2_reweight.rds")
saveRDS(results$model3_resample, "results/outcome/inc_dec/model3_resample.rds")
saveRDS(results$model3_reweight, "results/outcome/inc_dec/model3_reweight.rds")
saveRDS(results$model4_resample, "results/outcome/inc_dec/model4_resample.rds")
saveRDS(results$model4_reweight, "results/outcome/inc_dec/model4_reweight.rds")

# Save entire results list
saveRDS(results, "results/outcome/inc_dec/all_results.rds")

# Create comprehensive summary table with balance diagnostics
summary_table = data.frame(
  Model = rep(paste0("Model ", 1:4), each = 2),
  Bootstrap_Method = rep(c("Resample", "Reweight"), 4),
  AIPW_ATT = c(
    results$model1_resample$aipw$att,
    results$model1_reweight$aipw$att,
    results$model2_resample$aipw$att,
    results$model2_reweight$aipw$att,
    results$model3_resample$aipw$att,
    results$model3_reweight$aipw$att,
    results$model4_resample$aipw$att,
    results$model4_reweight$aipw$att
  ),
  AIPW_SE = c(
    results$model1_resample$aipw$se,
    results$model1_reweight$aipw$se,
    results$model2_resample$aipw$se,
    results$model2_reweight$aipw$se,
    results$model3_resample$aipw$se,
    results$model3_reweight$aipw$se,
    results$model4_resample$aipw$se,
    results$model4_reweight$aipw$se
  ),
  AIPW_CI_lower = c(
    results$model1_resample$aipw$ci_percentile[1],
    results$model1_reweight$aipw$ci_percentile[1],
    results$model2_resample$aipw$ci_percentile[1],
    results$model2_reweight$aipw$ci_percentile[1],
    results$model3_resample$aipw$ci_percentile[1],
    results$model3_reweight$aipw$ci_percentile[1],
    results$model4_resample$aipw$ci_percentile[1],
    results$model4_reweight$aipw$ci_percentile[1]
  ),
  AIPW_CI_upper = c(
    results$model1_resample$aipw$ci_percentile[2],
    results$model1_reweight$aipw$ci_percentile[2],
    results$model2_resample$aipw$ci_percentile[2],
    results$model2_reweight$aipw$ci_percentile[2],
    results$model3_resample$aipw$ci_percentile[2],
    results$model3_reweight$aipw$ci_percentile[2],
    results$model4_resample$aipw$ci_percentile[2],
    results$model4_reweight$aipw$ci_percentile[2]
  ),
  AIPW_pval = c(
    results$model1_resample$aipw$pval,
    results$model1_reweight$aipw$pval,
    results$model2_resample$aipw$pval,
    results$model2_reweight$aipw$pval,
    results$model3_resample$aipw$pval,
    results$model3_reweight$aipw$pval,
    results$model4_resample$aipw$pval,
    results$model4_reweight$aipw$pval
  ),
  DRS_ATT = c(
    results$model1_resample$drs$att,
    results$model1_reweight$drs$att,
    results$model2_resample$drs$att,
    results$model2_reweight$drs$att,
    results$model3_resample$drs$att,
    results$model3_reweight$drs$att,
    results$model4_resample$drs$att,
    results$model4_reweight$drs$att
  ),
  DRS_SE = c(
    results$model1_resample$drs$se,
    results$model1_reweight$drs$se,
    results$model2_resample$drs$se,
    results$model2_reweight$drs$se,
    results$model3_resample$drs$se,
    results$model3_reweight$drs$se,
    results$model4_resample$drs$se,
    results$model4_reweight$drs$se
  ),
  DRS_CI_lower = c(
    results$model1_resample$drs$ci_percentile[1],
    results$model1_reweight$drs$ci_percentile[1],
    results$model2_resample$drs$ci_percentile[1],
    results$model2_reweight$drs$ci_percentile[1],
    results$model3_resample$drs$ci_percentile[1],
    results$model3_reweight$drs$ci_percentile[1],
    results$model4_resample$drs$ci_percentile[1],
    results$model4_reweight$drs$ci_percentile[1]
  ),
  DRS_CI_upper = c(
    results$model1_resample$drs$ci_percentile[2],
    results$model1_reweight$drs$ci_percentile[2],
    results$model2_resample$drs$ci_percentile[2],
    results$model2_reweight$drs$ci_percentile[2],
    results$model3_resample$drs$ci_percentile[2],
    results$model3_reweight$drs$ci_percentile[2],
    results$model4_resample$drs$ci_percentile[2],
    results$model4_reweight$drs$ci_percentile[2]
  ),
  DRS_pval = c(
    results$model1_resample$drs$pval,
    results$model1_reweight$drs$pval,
    results$model2_resample$drs$pval,
    results$model2_reweight$drs$pval,
    results$model3_resample$drs$pval,
    results$model3_reweight$drs$pval,
    results$model4_resample$drs$pval,
    results$model4_reweight$drs$pval
  ),
  # Add balance diagnostics
  Avg_ASD = c(
    ifelse(
      !is.null(results$model1_resample$bootstrap_balance$avg_asd),
      mean(results$model1_resample$bootstrap_balance$avg_asd),
      NA
    ),
    ifelse(
      !is.null(results$model1_reweight$bootstrap_balance$avg_asd),
      mean(results$model1_reweight$bootstrap_balance$avg_asd),
      NA
    ),
    ifelse(
      !is.null(results$model2_resample$bootstrap_balance$avg_asd),
      mean(results$model2_resample$bootstrap_balance$avg_asd),
      NA
    ),
    ifelse(
      !is.null(results$model2_reweight$bootstrap_balance$avg_asd),
      mean(results$model2_reweight$bootstrap_balance$avg_asd),
      NA
    ),
    ifelse(
      !is.null(results$model3_resample$bootstrap_balance$avg_asd),
      mean(results$model3_resample$bootstrap_balance$avg_asd),
      NA
    ),
    ifelse(
      !is.null(results$model3_reweight$bootstrap_balance$avg_asd),
      mean(results$model3_reweight$bootstrap_balance$avg_asd),
      NA
    ),
    ifelse(
      !is.null(results$model4_resample$bootstrap_balance$avg_asd),
      mean(results$model4_resample$bootstrap_balance$avg_asd),
      NA
    ),
    ifelse(
      !is.null(results$model4_reweight$bootstrap_balance$avg_asd),
      mean(results$model4_reweight$bootstrap_balance$avg_asd),
      NA
    )
  ),
  Max_ASD = c(
    ifelse(
      !is.null(results$model1_resample$bootstrap_balance$max_asd),
      mean(results$model1_resample$bootstrap_balance$max_asd),
      NA
    ),
    ifelse(
      !is.null(results$model1_reweight$bootstrap_balance$max_asd),
      mean(results$model1_reweight$bootstrap_balance$max_asd),
      NA
    ),
    ifelse(
      !is.null(results$model2_resample$bootstrap_balance$max_asd),
      mean(results$model2_resample$bootstrap_balance$max_asd),
      NA
    ),
    ifelse(
      !is.null(results$model2_reweight$bootstrap_balance$max_asd),
      mean(results$model2_reweight$bootstrap_balance$max_asd),
      NA
    ),
    ifelse(
      !is.null(results$model3_resample$bootstrap_balance$max_asd),
      mean(results$model3_resample$bootstrap_balance$max_asd),
      NA
    ),
    ifelse(
      !is.null(results$model3_reweight$bootstrap_balance$max_asd),
      mean(results$model3_reweight$bootstrap_balance$max_asd),
      NA
    ),
    ifelse(
      !is.null(results$model4_resample$bootstrap_balance$max_asd),
      mean(results$model4_resample$bootstrap_balance$max_asd),
      NA
    ),
    ifelse(
      !is.null(results$model4_reweight$bootstrap_balance$max_asd),
      mean(results$model4_reweight$bootstrap_balance$max_asd),
      NA
    )
  ),
  Avg_ESS = c(
    ifelse(
      !is.null(results$model1_resample$bootstrap_balance$ess),
      mean(results$model1_resample$bootstrap_balance$ess),
      NA
    ),
    ifelse(
      !is.null(results$model1_reweight$bootstrap_balance$ess),
      mean(results$model1_reweight$bootstrap_balance$ess),
      NA
    ),
    ifelse(
      !is.null(results$model2_resample$bootstrap_balance$ess),
      mean(results$model2_resample$bootstrap_balance$ess),
      NA
    ),
    ifelse(
      !is.null(results$model2_reweight$bootstrap_balance$ess),
      mean(results$model2_reweight$bootstrap_balance$ess),
      NA
    ),
    ifelse(
      !is.null(results$model3_resample$bootstrap_balance$ess),
      mean(results$model3_resample$bootstrap_balance$ess),
      NA
    ),
    ifelse(
      !is.null(results$model3_reweight$bootstrap_balance$ess),
      mean(results$model3_reweight$bootstrap_balance$ess),
      NA
    ),
    ifelse(
      !is.null(results$model4_resample$bootstrap_balance$ess),
      mean(results$model4_resample$bootstrap_balance$ess),
      NA
    ),
    ifelse(
      !is.null(results$model4_reweight$bootstrap_balance$ess),
      mean(results$model4_reweight$bootstrap_balance$ess),
      NA
    )
  )
)

print(summary_table)

# Save as CSV
write.csv(
  summary_table,
  "results/outcome/inc_dec/bootstrap_DR_summary.csv",
  row.names = FALSE
)

cat("\nResults saved successfully!\n")

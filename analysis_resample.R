# Bootstrap Doubly Robust ATT Estimation Analysis - Resample Method
# Comparison: Increase vs Decrease
# This script runs DR_att() with 4 different outcome models using "resample" bootstrap method

# Setup
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
run = 123
source("utils.R")

# Load data
d = readRDS("data/data_ipwt2.rds")

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
      #change_res +
      #kinless +
      health_pre +
      chronic +
      death_due_covid +
      ppl_infected +
      income_loss +
      #job_loss +
      neighborhood +
      baseline_depr +
      baseline_lone
)

# Define 4 different outcome model formulas (RHS only, as character strings)
outcome_formula_1 = "1"

outcome_formula_2 = "baseline_lone "

outcome_formula_3 = "baseline_lone + baseline_depr "

outcome_formula_4 = "baseline_lone + baseline_depr + female + age_cat + edu +
                     emp_status + income + marital + coliving +
                     health_pre + chronic + death_due_covid +
                     ppl_infected + income_loss  + neighborhood"


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
  gbm_params = gbm_inc_dec,
  n_boot = 1000,
  seed = run,
  verbose = TRUE,
  parallel = T,
  n_cores = NULL,
  plot_diagnostics = TRUE,
  bootstrap_method = "resample"
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
  gbm_params = gbm_inc_dec,
  n_boot = 1000,
  seed = run,
  verbose = TRUE,
  parallel = TRUE,
  n_cores = NULL,
  plot_diagnostics = TRUE,
  bootstrap_method = "resample"
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
  gbm_params = gbm_inc_dec,
  n_boot = 1000,
  seed = run,
  verbose = TRUE,
  parallel = TRUE,
  n_cores = NULL,
  plot_diagnostics = TRUE,
  bootstrap_method = "resample"
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
  gbm_params = gbm_inc_dec,
  n_boot = 1000,
  seed = run,
  verbose = TRUE,
  parallel = TRUE,
  n_cores = NULL,
  plot_diagnostics = TRUE,
  bootstrap_method = "resample"
)

# =============================================================================
# Diagnostic plots - Increase vs Decrease (Resample Method)
# =============================================================================

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

# =============================================================================
# Save results
# =============================================================================

cat("\n\nSAVING RESULTS:\n")
cat("---------------\n")

# Save all individual results
saveRDS(results$model1_resample, "results/outcome/inc_dec/model1_resample.rds")
saveRDS(results$model2_resample, "results/outcome/inc_dec/model2_resample.rds")
saveRDS(results$model3_resample, "results/outcome/inc_dec/model3_resample.rds")
saveRDS(results$model4_resample, "results/outcome/inc_dec/model4_resample.rds")

# Save entire results list
saveRDS(results, "results/outcome/inc_dec/resample_results.rds")

# Create comprehensive summary table with balance diagnostics
summary_table = data.frame(
  Model = paste0("Model ", 1:4),
  Bootstrap_Method = "Resample",
  AIPW_ATT = c(
    results$model1_resample$aipw$att,
    results$model2_resample$aipw$att,
    results$model3_resample$aipw$att,
    results$model4_resample$aipw$att
  ),
  AIPW_SE = c(
    results$model1_resample$aipw$se,
    results$model2_resample$aipw$se,
    results$model3_resample$aipw$se,
    results$model4_resample$aipw$se
  ),
  AIPW_CI_lower = c(
    results$model1_resample$aipw$ci_percentile[1],
    results$model2_resample$aipw$ci_percentile[1],
    results$model3_resample$aipw$ci_percentile[1],
    results$model4_resample$aipw$ci_percentile[1]
  ),
  AIPW_CI_upper = c(
    results$model1_resample$aipw$ci_percentile[2],
    results$model2_resample$aipw$ci_percentile[2],
    results$model3_resample$aipw$ci_percentile[2],
    results$model4_resample$aipw$ci_percentile[2]
  ),
  AIPW_pval = c(
    results$model1_resample$aipw$pval,
    results$model2_resample$aipw$pval,
    results$model3_resample$aipw$pval,
    results$model4_resample$aipw$pval
  ),
  DRS_ATT = c(
    results$model1_resample$drs$att,
    results$model2_resample$drs$att,
    results$model3_resample$drs$att,
    results$model4_resample$drs$att
  ),
  DRS_SE = c(
    results$model1_resample$drs$se,
    results$model2_resample$drs$se,
    results$model3_resample$drs$se,
    results$model4_resample$drs$se
  ),
  DRS_CI_lower = c(
    results$model1_resample$drs$ci_percentile[1],
    results$model2_resample$drs$ci_percentile[1],
    results$model3_resample$drs$ci_percentile[1],
    results$model4_resample$drs$ci_percentile[1]
  ),
  DRS_CI_upper = c(
    results$model1_resample$drs$ci_percentile[2],
    results$model2_resample$drs$ci_percentile[2],
    results$model3_resample$drs$ci_percentile[2],
    results$model4_resample$drs$ci_percentile[2]
  ),
  DRS_pval = c(
    results$model1_resample$drs$pval,
    results$model2_resample$drs$pval,
    results$model3_resample$drs$pval,
    results$model4_resample$drs$pval
  ),
  # Add balance diagnostics
  Avg_ASD = c(
    ifelse(
      !is.null(results$model1_resample$bootstrap_balance$avg_asd),
      mean(results$model1_resample$bootstrap_balance$avg_asd),
      NA
    ),
    ifelse(
      !is.null(results$model2_resample$bootstrap_balance$avg_asd),
      mean(results$model2_resample$bootstrap_balance$avg_asd),
      NA
    ),
    ifelse(
      !is.null(results$model3_resample$bootstrap_balance$avg_asd),
      mean(results$model3_resample$bootstrap_balance$avg_asd),
      NA
    ),
    ifelse(
      !is.null(results$model4_resample$bootstrap_balance$avg_asd),
      mean(results$model4_resample$bootstrap_balance$avg_asd),
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
      !is.null(results$model2_resample$bootstrap_balance$max_asd),
      mean(results$model2_resample$bootstrap_balance$max_asd),
      NA
    ),
    ifelse(
      !is.null(results$model3_resample$bootstrap_balance$max_asd),
      mean(results$model3_resample$bootstrap_balance$max_asd),
      NA
    ),
    ifelse(
      !is.null(results$model4_resample$bootstrap_balance$max_asd),
      mean(results$model4_resample$bootstrap_balance$max_asd),
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
      !is.null(results$model2_resample$bootstrap_balance$ess),
      mean(results$model2_resample$bootstrap_balance$ess),
      NA
    ),
    ifelse(
      !is.null(results$model3_resample$bootstrap_balance$ess),
      mean(results$model3_resample$bootstrap_balance$ess),
      NA
    ),
    ifelse(
      !is.null(results$model4_resample$bootstrap_balance$ess),
      mean(results$model4_resample$bootstrap_balance$ess),
      NA
    )
  )
)

print(summary_table)

# Save as CSV
write.csv(
  summary_table,
  "results/outcome/inc_dec/resample_summary.csv",
  row.names = FALSE
)

cat("\nResults saved successfully!\n")

# SIMPLE/OR/DRS ATT ESTIMATION - COMPREHENSIVE ANALYSIS
#
# This script runs 12 analyses:
# - 2 pairwise comparisons: increase vs decrease, increase vs mix
# - 6 outcome models:
#   - Model 0: Simple difference in means (no covariates, no PS weights)
#   - Model 1: OR with full model (g-computation, no PS weights)
#   - Models 2-5: DRS with varying complexity (intercept only to full model)
# - No cross-fitting (k=1)
# - Stratified bootstrap with reweight method
#
# Total: 2 × 6 = 12 analyses

# Setup
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("utils.R")

# Load data
d = readRDS("data/data_iptw.rds")

# =============================================================================
# PARAMETERS
# =============================================================================

# Bootstrap parameters
n_boot = 5000
n_cores = NULL # Auto-detect
seed = 123
stratified = TRUE
sim = "balanced"
parallel = TRUE
ci_type = "all"
alpha = 0 # NO truncation

# Propensity score formula (same for both comparisons)
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

# Outcome model formulas (character strings for f.out parameter)
outcome_models = list(
  model0 = "1", # Placeholder for simple difference (f.out not used)
  model1 = "baseline_lone + baseline_depr + female + age_cat +
            edu + emp_status + income + marital + coliving +
            health_pre + chronic + death_due_covid +
            ppl_infected + income_loss + neighborhood", # Full model for OR (no PS weights)
  model2 = "1", # Intercept only for DRS
  model3 = "baseline_lone", # Baseline loneliness for DRS
  model4 = "baseline_lone + baseline_depr", # Baseline mental health for DRS
  model5 = "baseline_lone + baseline_depr + female + age_cat +
            edu + emp_status + income + marital + coliving +
            health_pre + chronic + death_due_covid +
            ppl_infected + income_loss + neighborhood" # Full model for DRS
)

# Outcome model descriptions (for reporting)
outcome_descriptions = c(
  "Simple: Difference in means",
  "OR: Full model (no PS weights)",
  "DRS: Intercept only",
  "DRS: Baseline loneliness",
  "DRS: Baseline loneliness + depression",
  "DRS: Full covariate adjustment"
)

# =============================================================================
# COMPARISON 1: INCREASE VS DECREASE
# =============================================================================

# Load tuned GBM parameters
gbm_params_inc_dec = readRDS("results/weighting/inc_dec_params.rds")

# Storage for all inc_dec results
results_inc_dec = list()

# Loop over models
for (model in 0:5) {
  cat(sprintf(
    "\n>>> Running Model %d: %s\n\n",
    model,
    outcome_descriptions[model + 1]
  ))

  # Build result path
  result_path = build_result_path(
    comparison = "inc_dec",
    model = model
  )

  # Determine estimator type and weight usage
  if (model == 0) {
    # Model 0: Simple difference in means
    estimator_type = "simple"
    use_weights = FALSE
  } else if (model == 1) {
    # Model 1: Outcome regression without PS weights
    estimator_type = "or"
    use_weights = FALSE
  } else {
    # Models 2-5: Doubly robust with PS weights
    estimator_type = "drs"
    use_weights = TRUE
  }

  # Run analysis
  result = DR_att(
    outcome = "loneliness",
    treatment = "remote_contact",
    treated_level = "increase",
    control_level = "decrease",
    f.ps = ps_formula,
    f.out = outcome_models[[model + 1]],
    data = d,
    gbm_params = gbm_params_inc_dec,
    bootstrap_method = "reweight",
    stratified = stratified,
    n_boot = n_boot,
    seed = seed,
    verbose = TRUE,
    parallel = parallel,
    n_cores = n_cores,
    plot_diagnostics = TRUE,
    ci_type = ci_type,
    sim = sim,
    save_to = result_path,
    alpha = alpha,
    use_weights = use_weights,
    estimator = estimator_type
  )

  # Store result with informative name
  result_name = sprintf("model%d", model)
  results_inc_dec[[result_name]] = result

  # Save individual result
  saveRDS(
    result,
    file.path(result_path, "result.rds")
  )

  cat(sprintf("\n✓ Completed: model=%d\n", model))
  cat(sprintf("  Saved to: %s\n", result_path))
}

# =============================================================================
# SUMMARY TABLES AND PLOTS - INC_DEC
# =============================================================================

# Create base directory for summaries
summary_dir_inc_dec = "results/outcome/DRS/inc_dec"
dir.create(summary_dir_inc_dec, showWarnings = FALSE, recursive = TRUE)

# Build master summary table
master_summary_inc_dec = build_master_summary(
  results_list = results_inc_dec,
  estimator_name = "DRS",
  comparison_name = "inc_dec"
)

# Save master summary
write.csv(
  master_summary_inc_dec,
  file.path(summary_dir_inc_dec, "DRS_master_summary_inc_dec.csv"),
  row.names = FALSE
)

cat(
  "\n✓ Master summary saved to:",
  file.path(summary_dir_inc_dec, "DRS_master_summary_inc_dec.csv"),
  "\n"
)

# Print summary
cat("\n=== MASTER SUMMARY TABLE (INC_DEC) ===\n\n")
print(master_summary_inc_dec, digits = 4)

# Create outcome model comparison boxplot
plot_outcome_model_boxplot(
  results_list = results_inc_dec,
  estimator_name = "drs",
  comparison_label = "Increase vs Decrease",
  save_path = file.path(summary_dir_inc_dec, "outcome_model_comparison.png")
)

cat(
  "✓ Boxplot saved to:",
  file.path(summary_dir_inc_dec, "outcome_model_comparison.png"),
  "\n"
)

# =============================================================================
# COMPARISON 2: INCREASE VS MIX
# =============================================================================

# Load tuned GBM parameters
gbm_params_inc_mix = readRDS("results/weighting/inc_mix_params.rds")

# Storage for all inc_mix results
results_inc_mix = list()

# Loop over models
for (model in 0:5) {
  cat(sprintf(
    "\n>>> Running Model %d: %s\n\n",
    model,
    outcome_descriptions[model + 1]
  ))

  # Build result path
  result_path = build_result_path(
    comparison = "inc_mix",
    model = model
  )

  # Determine estimator type and weight usage
  if (model == 0) {
    # Model 0: Simple difference in means
    estimator_type = "simple"
    use_weights = FALSE
  } else if (model == 1) {
    # Model 1: Outcome regression without PS weights
    estimator_type = "or"
    use_weights = FALSE
  } else {
    # Models 2-5: Doubly robust with PS weights
    estimator_type = "drs"
    use_weights = TRUE
  }

  # Run analysis
  result = DR_att(
    outcome = "loneliness",
    treatment = "remote_contact",
    treated_level = "increase",
    control_level = "mix",
    f.ps = ps_formula,
    f.out = outcome_models[[model + 1]],
    data = d,
    gbm_params = gbm_params_inc_mix,
    bootstrap_method = "reweight",
    stratified = stratified,
    n_boot = n_boot,
    seed = seed,
    verbose = TRUE,
    parallel = parallel,
    n_cores = n_cores,
    plot_diagnostics = TRUE,
    ci_type = ci_type,
    sim = sim,
    save_to = result_path,
    alpha = alpha,
    use_weights = use_weights,
    estimator = estimator_type
  )

  # Store result with informative name
  result_name = sprintf("model%d", model)
  results_inc_mix[[result_name]] = result

  # Save individual result
  saveRDS(
    result,
    file.path(result_path, "result.rds")
  )

  cat(sprintf("\n✓ Completed: model=%d\n", model))
  cat(sprintf("  Saved to: %s\n", result_path))
}

# =============================================================================
# SUMMARY TABLES AND PLOTS - INC_MIX
# =============================================================================

# Create base directory for summaries
summary_dir_inc_mix = "results/outcome/DRS/inc_mix"
dir.create(summary_dir_inc_mix, showWarnings = FALSE, recursive = TRUE)

# Build master summary table
master_summary_inc_mix = build_master_summary(
  results_list = results_inc_mix,
  estimator_name = "DRS",
  comparison_name = "inc_mix"
)

# Save master summary
write.csv(
  master_summary_inc_mix,
  file.path(summary_dir_inc_mix, "DRS_master_summary_inc_mix.csv"),
  row.names = FALSE
)

cat(
  "\n✓ Master summary saved to:",
  file.path(summary_dir_inc_mix, "DRS_master_summary_inc_mix.csv"),
  "\n"
)

# Print summary
cat("\n=== MASTER SUMMARY TABLE (INC_MIX) ===\n\n")
print(master_summary_inc_mix, digits = 4)

# Create outcome model comparison boxplot
plot_outcome_model_boxplot(
  results_list = results_inc_mix,
  estimator_name = "drs",
  comparison_label = "Increase vs Mix",
  save_path = file.path(summary_dir_inc_mix, "outcome_model_comparison.png")
)

cat(
  "✓ Boxplot saved to:",
  file.path(summary_dir_inc_mix, "outcome_model_comparison.png"),
  "\n"
)

# =============================================================================
# FINAL COMBINED SUMMARY
# =============================================================================

# Combine both comparisons
final_summary_drs = rbind(
  master_summary_inc_dec,
  master_summary_inc_mix
)

# Save final summary
dir.create("results/outcome/DRS", showWarnings = FALSE, recursive = TRUE)
write.csv(
  final_summary_drs,
  "results/outcome/DRS/DRS_master_summary_all_12.csv",
  row.names = FALSE
)

cat(
  "\n✓ Final Simple/OR/DRS master summary (12 analyses) saved to:",
  "results/outcome/DRS/DRS_master_summary_all_12.csv",
  "\n"
)

# Print final summary
cat("\n=== FINAL SIMPLE/OR/DRS SUMMARY (ALL 12 ANALYSES) ===\n\n")
print(final_summary_drs, digits = 4)

# Save complete workspace
saveRDS(
  list(
    inc_dec = results_inc_dec,
    inc_mix = results_inc_mix,
    master_summary_inc_dec = master_summary_inc_dec,
    master_summary_inc_mix = master_summary_inc_mix,
    final_summary = final_summary_drs
  ),
  "results/outcome/DRS/DRS_complete_workspace.rds"
)

cat(
  "\n✓ Complete DRS workspace saved to:",
  "results/outcome/DRS/DRS_complete_workspace.rds",
  "\n"
)

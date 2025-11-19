# DRS DOUBLY ROBUST ATT ESTIMATION - SIMPLIFIED ANALYSIS
#
# This script runs 8 DRS analyses:
# - 2 pairwise comparisons: increase vs decrease, increase vs mix
# - 4 outcome models: intercept only, baseline loneliness, + depression, full
# - No cross-fitting (k=1)
# - Stratified bootstrap with reweight method
#
# Total: 2 × 4 = 8 analyses

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
  model1 = "1", # Intercept only
  model2 = "baseline_lone", # Baseline loneliness
  model3 = "baseline_lone + baseline_depr", # Baseline mental health
  model4 = "baseline_lone + baseline_depr + female + age_cat +
            edu + emp_status + income + marital + coliving +
            health_pre + chronic + death_due_covid +
            ppl_infected + income_loss + neighborhood" # Full model
)

# Outcome model descriptions (for reporting)
outcome_descriptions = c(
  "Intercept only",
  "Baseline loneliness",
  "Baseline loneliness + depression",
  "Full covariate adjustment"
)

# =============================================================================
# COMPARISON 1: INCREASE VS DECREASE
# =============================================================================

cat("\n")
cat("=" %>% rep(80) %>% paste0(collapse = ""), "\n")
cat("COMPARISON 1: INCREASE VS DECREASE\n")
cat("=" %>% rep(80) %>% paste0(collapse = ""), "\n\n")

# Load tuned GBM parameters
gbm_params_inc_dec = readRDS("results/weighting/inc_dec_params.rds")

# Storage for all inc_dec results
results_inc_dec = list()

# Loop over models
for (model in 1:4) {
  cat(sprintf(
    "\n>>> Running Model %d: %s\n\n",
    model,
    outcome_descriptions[model]
  ))

  # Build result path
  result_path = build_result_path(
    comparison = "inc_dec",
    model = model
  )

  # Run analysis
  result = DR_att(
    outcome = "loneliness",
    treatment = "remote_contact",
    treated_level = "increase",
    control_level = "decrease",
    f.ps = ps_formula,
    f.out = outcome_models[[model]],
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
    alpha = alpha
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

cat("\n")
cat("=" %>% rep(80) %>% paste0(collapse = ""), "\n")
cat("CREATING SUMMARIES FOR INC_DEC\n")
cat("=" %>% rep(80) %>% paste0(collapse = ""), "\n\n")

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
cat("\n\nCreating outcome model comparison boxplot...\n")

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

cat("\n\n")
cat("=" %>% rep(80) %>% paste0(collapse = ""), "\n")
cat("COMPARISON 2: INCREASE VS MIX\n")
cat("=" %>% rep(80) %>% paste0(collapse = ""), "\n\n")

# Load tuned GBM parameters
gbm_params_inc_mix = readRDS("results/weighting/inc_mix_params.rds")

# Storage for all inc_mix results
results_inc_mix = list()

# Loop over models
for (model in 1:4) {
  cat(sprintf(
    "\n>>> Running Model %d: %s\n\n",
    model,
    outcome_descriptions[model]
  ))

  # Build result path
  result_path = build_result_path(
    comparison = "inc_mix",
    model = model
  )

  # Run analysis
  result = DR_att(
    outcome = "loneliness",
    treatment = "remote_contact",
    treated_level = "increase",
    control_level = "mix",
    f.ps = ps_formula,
    f.out = outcome_models[[model]],
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
    alpha = alpha
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

cat("\n")
cat("=" %>% rep(80) %>% paste0(collapse = ""), "\n")
cat("CREATING SUMMARIES FOR INC_MIX\n")
cat("=" %>% rep(80) %>% paste0(collapse = ""), "\n\n")

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
cat("\n\nCreating outcome model comparison boxplot...\n")

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

cat("\n\n")
cat("=" %>% rep(80) %>% paste0(collapse = ""), "\n")
cat("FINAL DRS SUMMARY (ALL 8 ANALYSES)\n")
cat("=" %>% rep(80) %>% paste0(collapse = ""), "\n\n")

# Combine both comparisons
final_summary_drs = rbind(
  master_summary_inc_dec,
  master_summary_inc_mix
)

# Save final summary
dir.create("results/outcome/DRS", showWarnings = FALSE, recursive = TRUE)
write.csv(
  final_summary_drs,
  "results/outcome/DRS/DRS_master_summary_all_8.csv",
  row.names = FALSE
)

cat(
  "\n✓ Final DRS master summary (8 analyses) saved to:",
  "results/outcome/DRS/DRS_master_summary_all_8.csv",
  "\n"
)

# Print final summary
cat("\n=== FINAL DRS SUMMARY (ALL 8 ANALYSES) ===\n\n")
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

cat("\n\n")
cat("=" %>% rep(80) %>% paste0(collapse = ""), "\n")
cat("DRS ANALYSIS COMPLETE!\n")
cat("=" %>% rep(80) %>% paste0(collapse = ""), "\n")
cat(sprintf("\nTotal analyses run: 8\n"))
cat(sprintf("Total bootstrap replications per analysis: %d\n", n_boot))
cat(sprintf("Results saved to: results/outcome/DRS/\n"))
cat("\n")

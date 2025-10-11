# DRS DOUBLY ROBUST ATT ESTIMATION - COMPREHENSIVE ANALYSIS
#
# This script runs 24 DRS analyses:
# - 2 pairwise comparisons: increase vs decrease, increase vs mix
# - 3 cross-fitting approaches: k=1 (no CF), k=2, k=3
# - 4 outcome models: intercept only, baseline loneliness, + depression, full
#
# Total: 2 × 3 × 4 = 24 analyses

# Setup
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("utils.R")


# Load data
d = readRDS("data/data_iptw.rds")

# =============================================================================
# PARAMETERS
# =============================================================================

# Bootstrap parameters
n_boot = 500
n_cores = NULL # Auto-detect
seed = 123
stratified = TRUE
sim = "balanced"
parallel = TRUE
ci_type = "all"
alpha = 0 # NO truncation

# Cross-fitting values
k_values = c(1, 2, 3) # k=1 means cross_fitting=FALSE

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

# Nested loops: k-fold × model
for (k in k_values) {
  cf_status = if (k == 1) {
    "No cross-fitting"
  } else {
    sprintf("%d-fold cross-fitting", k)
  }

  cat("\n")
  cat("-" %>% rep(80) %>% paste0(collapse = ""), "\n")
  cat(sprintf("K = %d (%s)\n", k, cf_status))
  cat("-" %>% rep(80) %>% paste0(collapse = ""), "\n\n")

  for (model in 1:4) {
    cat(sprintf(
      "\n>>> Running Model %d: %s\n\n",
      model,
      outcome_descriptions[model]
    ))

    # Build result path
    result_path = build_result_path(
      estimator = "drs",
      comparison = "inc_dec",
      k = k,
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
      estimator = "drs",
      bootstrap_method = "reweight",
      cross_fitting = (k > 1),
      k = k,
      stratify_folds = TRUE,
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
    result_name = sprintf("k%d_model%d", k, model)
    results_inc_dec[[result_name]] = result

    # Save individual result
    saveRDS(
      result,
      file.path(result_path, "result.rds")
    )

    cat(sprintf("\n✓ Completed: k=%d, model=%d\n", k, model))
    cat(sprintf("  Saved to: %s\n", result_path))
  }
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

# Create k-fold comparison plots for each k
cat("\n\nCreating k-fold comparison plot...\n")

# Extract results by k
results_k1_inc_dec = results_inc_dec[grep("^k1_", names(results_inc_dec))]
results_k2_inc_dec = results_inc_dec[grep("^k2_", names(results_inc_dec))]
results_k3_inc_dec = results_inc_dec[grep("^k3_", names(results_inc_dec))]

# Create comparison plot
plot_kfold_comparison(
  results_k1 = results_k1_inc_dec,
  results_k2 = results_k2_inc_dec,
  results_k3 = results_k3_inc_dec,
  estimator_name = "drs",
  comparison_label = "Increase vs Decrease",
  save_path = file.path(summary_dir_inc_dec, "kfold_comparison.png")
)

cat(
  "✓ K-fold comparison plot saved to:",
  file.path(summary_dir_inc_dec, "kfold_comparison.png"),
  "\n"
)

# Create per-k summary tables
for (k in k_values) {
  k_results = results_inc_dec[grep(sprintf("^k%d_", k), names(results_inc_dec))]

  # Build summary for this k
  k_summary = build_master_summary(
    results_list = k_results,
    estimator_name = "DRS",
    comparison_name = "inc_dec"
  )

  # Create directory for this k
  cf_suffix = if (k == 1) "nocf" else "cf"
  k_dir = file.path(summary_dir_inc_dec, sprintf("k%d_%s", k, cf_suffix))
  dir.create(k_dir, showWarnings = FALSE, recursive = TRUE)

  # Save k-specific summary
  write.csv(
    k_summary,
    file.path(k_dir, sprintf("summary_k%d.csv", k)),
    row.names = FALSE
  )

  cat(sprintf(
    "✓ K=%d summary saved to: %s\n",
    k,
    file.path(k_dir, sprintf("summary_k%d.csv", k))
  ))
}

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

# Nested loops: k-fold × model
for (k in k_values) {
  cf_status = if (k == 1) {
    "No cross-fitting"
  } else {
    sprintf("%d-fold cross-fitting", k)
  }

  cat("\n")
  cat("-" %>% rep(80) %>% paste0(collapse = ""), "\n")
  cat(sprintf("K = %d (%s)\n", k, cf_status))
  cat("-" %>% rep(80) %>% paste0(collapse = ""), "\n\n")

  for (model in 1:4) {
    cat(sprintf(
      "\n>>> Running Model %d: %s\n\n",
      model,
      outcome_descriptions[model]
    ))

    # Build result path
    result_path = build_result_path(
      estimator = "drs",
      comparison = "inc_mix",
      k = k,
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
      estimator = "drs",
      bootstrap_method = "reweight",
      cross_fitting = (k > 1),
      k = k,
      stratify_folds = TRUE,
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
    result_name = sprintf("k%d_model%d", k, model)
    results_inc_mix[[result_name]] = result

    # Save individual result
    saveRDS(
      result,
      file.path(result_path, "result.rds")
    )

    cat(sprintf("\n✓ Completed: k=%d, model=%d\n", k, model))
    cat(sprintf("  Saved to: %s\n", result_path))
  }
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

# Create k-fold comparison plots
cat("\n\nCreating k-fold comparison plot...\n")

# Extract results by k
results_k1_inc_mix = results_inc_mix[grep("^k1_", names(results_inc_mix))]
results_k2_inc_mix = results_inc_mix[grep("^k2_", names(results_inc_mix))]
results_k3_inc_mix = results_inc_mix[grep("^k3_", names(results_inc_mix))]

# Create comparison plot
plot_kfold_comparison(
  results_k1 = results_k1_inc_mix,
  results_k2 = results_k2_inc_mix,
  results_k3 = results_k3_inc_mix,
  estimator_name = "drs",
  comparison_label = "Increase vs Mix",
  save_path = file.path(summary_dir_inc_mix, "kfold_comparison.png")
)

cat(
  "✓ K-fold comparison plot saved to:",
  file.path(summary_dir_inc_mix, "kfold_comparison.png"),
  "\n"
)

# Create per-k summary tables
for (k in k_values) {
  k_results = results_inc_mix[grep(sprintf("^k%d_", k), names(results_inc_mix))]

  # Build summary for this k
  k_summary = build_master_summary(
    results_list = k_results,
    estimator_name = "DRS",
    comparison_name = "inc_mix"
  )

  # Create directory for this k
  cf_suffix = if (k == 1) "nocf" else "cf"
  k_dir = file.path(summary_dir_inc_mix, sprintf("k%d_%s", k, cf_suffix))
  dir.create(k_dir, showWarnings = FALSE, recursive = TRUE)

  # Save k-specific summary
  write.csv(
    k_summary,
    file.path(k_dir, sprintf("summary_k%d.csv", k)),
    row.names = FALSE
  )

  cat(sprintf(
    "✓ K=%d summary saved to: %s\n",
    k,
    file.path(k_dir, sprintf("summary_k%d.csv", k))
  ))
}

# =============================================================================
# FINAL COMBINED SUMMARY
# =============================================================================

cat("\n\n")
cat("=" %>% rep(80) %>% paste0(collapse = ""), "\n")
cat("FINAL DRS SUMMARY (ALL 24 ANALYSES)\n")
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
  "results/outcome/DRS/DRS_master_summary_all_24.csv",
  row.names = FALSE
)

cat(
  "\n✓ Final DRS master summary (24 analyses) saved to:",
  "results/outcome/DRS/DRS_master_summary_all_24.csv",
  "\n"
)

# Print final summary
cat("\n=== FINAL DRS SUMMARY (ALL 24 ANALYSES) ===\n\n")
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
cat(sprintf("\nTotal analyses run: 24\n"))
cat(sprintf("Total bootstrap replications per analysis: %d\n", n_boot))
cat(sprintf("Results saved to: results/outcome/DRS/\n"))
cat("\n")

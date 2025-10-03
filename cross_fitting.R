# ==============================================================================
# Cross-Fitting Comparison Analysis
# Compare standard AIPW vs Cross-fitted AIPW estimators
# Across 4 outcome model specifications
# ==============================================================================

# Setup -----------------------------------------------------------------------
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
set.seed(123)
source("utils.R")

# Load data with pre-computed IPTW weights
cat("Loading data...\n")
d = readRDS("data/data_iptw.rds")

# Define propensity score formula (for cross-fitting)
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

# Define 4 outcome model formulas (RHS only)
outcome_formulas = list(
  model1 = "1",
  model2 = "baseline_lone",
  model3 = "baseline_lone + baseline_depr",
  model4 = "baseline_lone + baseline_depr + female + age_cat + edu +
            emp_status + income + marital + coliving + change_res +
            kinless + health_pre + chronic + death_due_covid +
            ppl_infected + income_loss + job_loss + neighborhood"
)

# Load tuned GBM parameters
gbm_inc_dec = readRDS("results/weighting/inc_dec_params.rds")
gbm_inc_mix = readRDS("results/weighting/inc_mix_params.rds")

# Initialize results storage
results = list(
  inc_dec = list(standard = list(), cross_fitted = list()),
  inc_mix = list(standard = list(), cross_fitted = list())
)

# ==============================================================================
# PART 1: INCREASE VS DECREASE COMPARISON
# ==============================================================================

cat("\n=== INCREASE VS DECREASE COMPARISON ===\n\n")

# Prepare data for inc_dec comparison
inc_dec_analysis = d %>%
  filter(remote_contact %in% c("increase", "decrease")) %>%
  mutate(
    # Convert to binary for analysis functions
    remote_contact = ifelse(remote_contact == "increase", 1, 0),
    across(where(is.factor), droplevels)
  )

# Loop through each model specification
for (model_name in names(outcome_formulas)) {
  cat(sprintf("\n--- Processing %s (Inc vs Dec) ---\n", model_name))

  # 1. Standard AIPW (using pre-computed weights)
  cat("Running standard AIPW...\n")
  standard_result = aipw_att(
    outcome = "severe_loneliness",
    treatment = "remote_contact",
    f.out = outcome_formulas[[model_name]],
    wgt = "iptw",
    data = inc_dec_analysis,
    verbose = FALSE
  )
  results$inc_dec$standard[[model_name]] = standard_result
  cat(sprintf("Standard AIPW ATT: %.4f\n", standard_result$att))

  # 2. Cross-fitted AIPW
  cat("Running cross-fitted AIPW (5-fold)...\n")
  cf_result = cf_aipw_att(
    outcome = "severe_loneliness",
    treatment = "remote_contact",
    f.ps = ps_formula,
    f.out = outcome_formulas[[model_name]],
    data = inc_dec_analysis %>% select(-iptw), # Remove pre-computed weights for cross-fitting
    k = 5,
    stratify = TRUE,
    ps_params = gbm_inc_dec,
    seed = 123,
    verbose = FALSE
  )
  results$inc_dec$cross_fitted[[model_name]] = cf_result
  cat(sprintf("Cross-fitted AIPW ATT: %.4f\n", cf_result$att))

  # Calculate difference
  diff = cf_result$att - standard_result$att
  cat(sprintf("Difference (CF - Standard): %.4f\n", diff))
}

# ==============================================================================
# PART 2: INCREASE VS MIX COMPARISON
# ==============================================================================

cat("\n\n=== INCREASE VS MIX COMPARISON ===\n\n")

# Prepare data for inc_mix comparison
inc_mix_analysis = d %>%
  filter(remote_contact %in% c("increase", "mix")) %>%
  mutate(
    # Convert to binary for analysis functions
    remote_contact = ifelse(remote_contact == "increase", 1, 0),
    across(where(is.factor), droplevels)
  )

# Loop through each model specification
for (model_name in names(outcome_formulas)) {
  cat(sprintf("\n--- Processing %s (Inc vs Mix) ---\n", model_name))

  # 1. Standard AIPW (using pre-computed weights)
  cat("Running standard AIPW...\n")
  standard_result = aipw_att(
    outcome = "severe_loneliness",
    treatment = "remote_contact",
    f.out = outcome_formulas[[model_name]],
    wgt = "iptw",
    data = inc_mix_analysis,
    verbose = FALSE
  )
  results$inc_mix$standard[[model_name]] = standard_result
  cat(sprintf("Standard AIPW ATT: %.4f\n", standard_result$att))

  # 2. Cross-fitted AIPW
  cat("Running cross-fitted AIPW (5-fold)...\n")
  cf_result = cf_aipw_att(
    outcome = "severe_loneliness",
    treatment = "remote_contact",
    f.ps = ps_formula,
    f.out = outcome_formulas[[model_name]],
    data = inc_mix_analysis %>% select(-iptw), # Remove pre-computed weights for cross-fitting
    k = 5,
    stratify = TRUE,
    ps_params = gbm_inc_mix,
    seed = 123,
    verbose = FALSE
  )
  results$inc_mix$cross_fitted[[model_name]] = cf_result
  cat(sprintf("Cross-fitted AIPW ATT: %.4f\n", cf_result$att))

  # Calculate difference
  diff = cf_result$att - standard_result$att
  cat(sprintf("Difference (CF - Standard): %.4f\n", diff))
}

# ==============================================================================
# PART 3: CREATE SUMMARY TABLE
# ==============================================================================

cat("\n\n=== CREATING SUMMARY TABLE ===\n\n")

# Build comprehensive comparison table
comparison_table = data.frame()

# Process Inc vs Dec results
for (model_name in names(outcome_formulas)) {
  row_data = data.frame(
    Comparison = "Inc vs Dec",
    Model = gsub("model", "Model ", model_name),
    Standard_AIPW = results$inc_dec$standard[[model_name]]$att,
    CrossFitted_AIPW = results$inc_dec$cross_fitted[[model_name]]$att,
    Difference = results$inc_dec$cross_fitted[[model_name]]$att -
      results$inc_dec$standard[[model_name]]$att,
    Relative_Diff_Pct = 100 *
      (results$inc_dec$cross_fitted[[model_name]]$att -
        results$inc_dec$standard[[model_name]]$att) /
      abs(results$inc_dec$standard[[model_name]]$att)
  )
  comparison_table = rbind(comparison_table, row_data)
}

# Process Inc vs Mix results
for (model_name in names(outcome_formulas)) {
  row_data = data.frame(
    Comparison = "Inc vs Mix",
    Model = gsub("model", "Model ", model_name),
    Standard_AIPW = results$inc_mix$standard[[model_name]]$att,
    CrossFitted_AIPW = results$inc_mix$cross_fitted[[model_name]]$att,
    Difference = results$inc_mix$cross_fitted[[model_name]]$att -
      results$inc_mix$standard[[model_name]]$att,
    Relative_Diff_Pct = 100 *
      (results$inc_mix$cross_fitted[[model_name]]$att -
        results$inc_mix$standard[[model_name]]$att) /
      abs(results$inc_mix$standard[[model_name]]$att)
  )
  comparison_table = rbind(comparison_table, row_data)
}

# Round numeric columns
comparison_table[, 3:6] = round(comparison_table[, 3:6], 4)

# Display table
cat("Comparison Table:\n")
print(comparison_table, row.names = FALSE)

# ==============================================================================
# PART 4: VISUALIZATIONS
# ==============================================================================

cat("\n\n=== CREATING VISUALIZATIONS ===\n\n")

# Set up plotting parameters
par(mfrow = c(2, 2), mar = c(5, 4, 4, 2))

# Plot 1: Inc vs Dec comparison across models
plot_data_inc_dec = comparison_table[
  comparison_table$Comparison == "Inc vs Dec",
]
x_pos = 1:4
width = 0.35

plot(
  x = x_pos - width / 2,
  y = plot_data_inc_dec$Standard_AIPW,
  xlim = c(0.5, 4.5),
  ylim = range(c(
    plot_data_inc_dec$Standard_AIPW,
    plot_data_inc_dec$CrossFitted_AIPW
  )) *
    1.1,
  pch = 19,
  col = "blue",
  xlab = "Model",
  ylab = "ATT Estimate",
  main = "AIPW vs Cross-fitted AIPW: Inc vs Dec",
  xaxt = "n"
)
points(
  x = x_pos + width / 2,
  y = plot_data_inc_dec$CrossFitted_AIPW,
  pch = 19,
  col = "red"
)
axis(1, at = x_pos, labels = paste("Model", 1:4))
abline(h = 0, lty = 2, col = "gray")
legend(
  "topright",
  legend = c("Standard AIPW", "Cross-fitted AIPW"),
  col = c("blue", "red"),
  pch = 19,
  cex = 0.8
)

# Add connecting lines
for (i in 1:4) {
  segments(
    x_pos[i] - width / 2,
    plot_data_inc_dec$Standard_AIPW[i],
    x_pos[i] + width / 2,
    plot_data_inc_dec$CrossFitted_AIPW[i],
    col = "gray",
    lty = 2
  )
}

# Plot 2: Inc vs Mix comparison across models
plot_data_inc_mix = comparison_table[
  comparison_table$Comparison == "Inc vs Mix",
]

plot(
  x = x_pos - width / 2,
  y = plot_data_inc_mix$Standard_AIPW,
  xlim = c(0.5, 4.5),
  ylim = range(c(
    plot_data_inc_mix$Standard_AIPW,
    plot_data_inc_mix$CrossFitted_AIPW
  )) *
    1.1,
  pch = 19,
  col = "blue",
  xlab = "Model",
  ylab = "ATT Estimate",
  main = "AIPW vs Cross-fitted AIPW: Inc vs Mix",
  xaxt = "n"
)
points(
  x = x_pos + width / 2,
  y = plot_data_inc_mix$CrossFitted_AIPW,
  pch = 19,
  col = "red"
)
axis(1, at = x_pos, labels = paste("Model", 1:4))
abline(h = 0, lty = 2, col = "gray")
legend(
  "topright",
  legend = c("Standard AIPW", "Cross-fitted AIPW"),
  col = c("blue", "red"),
  pch = 19,
  cex = 0.8
)

# Add connecting lines
for (i in 1:4) {
  segments(
    x_pos[i] - width / 2,
    plot_data_inc_mix$Standard_AIPW[i],
    x_pos[i] + width / 2,
    plot_data_inc_mix$CrossFitted_AIPW[i],
    col = "gray",
    lty = 2
  )
}

# Plot 3: Differences across models
barplot(
  height = rbind(plot_data_inc_dec$Difference, plot_data_inc_mix$Difference),
  beside = TRUE,
  names.arg = paste("Model", 1:4),
  col = c("darkblue", "darkred"),
  main = "Difference: Cross-fitted - Standard AIPW",
  ylab = "Difference in ATT",
  ylim = c(
    min(comparison_table$Difference) * 1.2,
    max(comparison_table$Difference) * 1.2
  )
)
abline(h = 0, lty = 1)
legend(
  "topright",
  legend = c("Inc vs Dec", "Inc vs Mix"),
  fill = c("darkblue", "darkred"),
  cex = 0.8
)

# Plot 4: Relative differences (%)
barplot(
  height = rbind(
    plot_data_inc_dec$Relative_Diff_Pct,
    plot_data_inc_mix$Relative_Diff_Pct
  ),
  beside = TRUE,
  names.arg = paste("Model", 1:4),
  col = c("lightblue", "lightcoral"),
  main = "Relative Difference (%)",
  ylab = "Relative Difference (%)",
  ylim = c(
    min(comparison_table$Relative_Diff_Pct) * 1.2,
    max(comparison_table$Relative_Diff_Pct) * 1.2
  )
)
abline(h = 0, lty = 1)
legend(
  "topright",
  legend = c("Inc vs Dec", "Inc vs Mix"),
  fill = c("lightblue", "lightcoral"),
  cex = 0.8
)

# Reset plotting parameters
par(mfrow = c(1, 1))

# ==============================================================================
# PART 5: SAVE RESULTS
# ==============================================================================

cat("\n\n=== SAVING RESULTS ===\n\n")

# Create output directory if it doesn't exist
if (!dir.exists("results/outcome/cross_fitting")) {
  dir.create("results/outcome/cross_fitting", recursive = TRUE)
}

# Save results object
saveRDS(results, "results/outcome/cross_fitting/all_results.rds")
cat("Saved results object to results/outcome/cross_fitting/all_results.rds\n")

# Save comparison table
write.csv(
  comparison_table,
  "results/outcome/cross_fitting/comparison_table.csv",
  row.names = FALSE
)
cat(
  "Saved comparison table to results/outcome/cross_fitting/comparison_table.csv\n"
)

# Save plots to PDF
pdf(
  "results/outcome/cross_fitting/comparison_plots.pdf",
  width = 10,
  height = 10
)
par(mfrow = c(2, 2), mar = c(5, 4, 4, 2))

# Recreate all 4 plots for PDF
# Plot 1: Inc vs Dec
plot(
  x = x_pos - width / 2,
  y = plot_data_inc_dec$Standard_AIPW,
  xlim = c(0.5, 4.5),
  ylim = range(c(
    plot_data_inc_dec$Standard_AIPW,
    plot_data_inc_dec$CrossFitted_AIPW
  )) *
    1.1,
  pch = 19,
  col = "blue",
  xlab = "Model",
  ylab = "ATT Estimate",
  main = "AIPW vs Cross-fitted AIPW: Inc vs Dec",
  xaxt = "n"
)
points(
  x = x_pos + width / 2,
  y = plot_data_inc_dec$CrossFitted_AIPW,
  pch = 19,
  col = "red"
)
axis(1, at = x_pos, labels = paste("Model", 1:4))
abline(h = 0, lty = 2, col = "gray")
legend(
  "topright",
  legend = c("Standard AIPW", "Cross-fitted AIPW"),
  col = c("blue", "red"),
  pch = 19,
  cex = 0.8
)
for (i in 1:4) {
  segments(
    x_pos[i] - width / 2,
    plot_data_inc_dec$Standard_AIPW[i],
    x_pos[i] + width / 2,
    plot_data_inc_dec$CrossFitted_AIPW[i],
    col = "gray",
    lty = 2
  )
}

# Plot 2: Inc vs Mix
plot(
  x = x_pos - width / 2,
  y = plot_data_inc_mix$Standard_AIPW,
  xlim = c(0.5, 4.5),
  ylim = range(c(
    plot_data_inc_mix$Standard_AIPW,
    plot_data_inc_mix$CrossFitted_AIPW
  )) *
    1.1,
  pch = 19,
  col = "blue",
  xlab = "Model",
  ylab = "ATT Estimate",
  main = "AIPW vs Cross-fitted AIPW: Inc vs Mix",
  xaxt = "n"
)
points(
  x = x_pos + width / 2,
  y = plot_data_inc_mix$CrossFitted_AIPW,
  pch = 19,
  col = "red"
)
axis(1, at = x_pos, labels = paste("Model", 1:4))
abline(h = 0, lty = 2, col = "gray")
legend(
  "topright",
  legend = c("Standard AIPW", "Cross-fitted AIPW"),
  col = c("blue", "red"),
  pch = 19,
  cex = 0.8
)
for (i in 1:4) {
  segments(
    x_pos[i] - width / 2,
    plot_data_inc_mix$Standard_AIPW[i],
    x_pos[i] + width / 2,
    plot_data_inc_mix$CrossFitted_AIPW[i],
    col = "gray",
    lty = 2
  )
}

# Plot 3: Differences
barplot(
  height = rbind(plot_data_inc_dec$Difference, plot_data_inc_mix$Difference),
  beside = TRUE,
  names.arg = paste("Model", 1:4),
  col = c("darkblue", "darkred"),
  main = "Difference: Cross-fitted - Standard AIPW",
  ylab = "Difference in ATT",
  ylim = c(
    min(comparison_table$Difference) * 1.2,
    max(comparison_table$Difference) * 1.2
  )
)
abline(h = 0, lty = 1)
legend(
  "topright",
  legend = c("Inc vs Dec", "Inc vs Mix"),
  fill = c("darkblue", "darkred"),
  cex = 0.8
)

# Plot 4: Relative differences
barplot(
  height = rbind(
    plot_data_inc_dec$Relative_Diff_Pct,
    plot_data_inc_mix$Relative_Diff_Pct
  ),
  beside = TRUE,
  names.arg = paste("Model", 1:4),
  col = c("lightblue", "lightcoral"),
  main = "Relative Difference (%)",
  ylab = "Relative Difference (%)",
  ylim = c(
    min(comparison_table$Relative_Diff_Pct) * 1.2,
    max(comparison_table$Relative_Diff_Pct) * 1.2
  )
)
abline(h = 0, lty = 1)
legend(
  "topright",
  legend = c("Inc vs Dec", "Inc vs Mix"),
  fill = c("lightblue", "lightcoral"),
  cex = 0.8
)

dev.off()
cat("Saved plots to results/outcome/cross_fitting/comparison_plots.pdf\n")

# ==============================================================================
# SUMMARY
# ==============================================================================

cat("\n\n=== ANALYSIS COMPLETE ===\n\n")
cat("Key Findings:\n")
cat("--------------\n")

# Calculate average absolute difference
avg_abs_diff = mean(abs(comparison_table$Difference))
cat(sprintf("Average absolute difference: %.4f\n", avg_abs_diff))

# Find max difference
max_diff_row = comparison_table[which.max(abs(comparison_table$Difference)), ]
cat(sprintf(
  "Maximum difference: %.4f (%s, %s)\n",
  max_diff_row$Difference,
  max_diff_row$Comparison,
  max_diff_row$Model
))

# Report range of relative differences
cat(sprintf(
  "Relative difference range: %.2f%% to %.2f%%\n",
  min(comparison_table$Relative_Diff_Pct),
  max(comparison_table$Relative_Diff_Pct)
))

cat("\nAll results saved to results/outcome/cross_fitting/\n")

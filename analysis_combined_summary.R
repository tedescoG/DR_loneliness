# COMBINED SUMMARY: AIPW AND DRS ANALYSES
#
# This script loads all 48 analyses (24 AIPW + 24 DRS) and creates:
# - Master summary table with all 48 analyses
# - Comparison plots: AIPW vs DRS
# - Model complexity effects
# - Cross-fitting effects

# Setup
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("utils.R")

cat("\n")
cat("=" %>% rep(80) %>% paste0(collapse = ""), "\n")
cat("COMBINED SUMMARY: LOADING ALL RESULTS\n")
cat("=" %>% rep(80) %>% paste0(collapse = ""), "\n\n")

# =============================================================================
# LOAD ALL RESULTS
# =============================================================================

# Load AIPW workspace
cat("Loading AIPW results...\n")
aipw_workspace = readRDS("results/outcome/AIPW/AIPW_complete_workspace.rds")
aipw_summary = aipw_workspace$final_summary

# Load DRS workspace
cat("Loading DRS results...\n")
drs_workspace = readRDS("results/outcome/DRS/DRS_complete_workspace.rds")
drs_summary = drs_workspace$final_summary

cat("✓ All results loaded successfully\n\n")

# =============================================================================
# CREATE MASTER SUMMARY TABLE (ALL 48 ANALYSES)
# =============================================================================

cat("\n")
cat("=" %>% rep(80) %>% paste0(collapse = ""), "\n")
cat("CREATING MASTER SUMMARY TABLE (48 ANALYSES)\n")
cat("=" %>% rep(80) %>% paste0(collapse = ""), "\n\n")

# Combine AIPW and DRS summaries
master_summary_all = rbind(
  aipw_summary,
  drs_summary
)

# Sort by comparison, k, model, estimator
master_summary_all = master_summary_all[
  order(
    master_summary_all$Comparison,
    master_summary_all$K_Fold,
    master_summary_all$Model,
    master_summary_all$Estimator
  ),
]

# Create combined directory
dir.create("results/outcome/combined", showWarnings = FALSE, recursive = TRUE)

# Save master summary
write.csv(
  master_summary_all,
  "results/outcome/combined/master_summary_all_48.csv",
  row.names = FALSE
)

cat("✓ Master summary (48 analyses) saved to:",
    "results/outcome/combined/master_summary_all_48.csv\n\n")

# Print summary
cat("=== MASTER SUMMARY TABLE (ALL 48 ANALYSES) ===\n\n")
print(master_summary_all, digits = 4)

# =============================================================================
# COMPARISON PLOTS: AIPW vs DRS
# =============================================================================

cat("\n\n")
cat("=" %>% rep(80) %>% paste0(collapse = ""), "\n")
cat("CREATING COMPARISON PLOTS\n")
cat("=" %>% rep(80) %>% paste0(collapse = ""), "\n\n")

# Function to create AIPW vs DRS comparison plot
plot_estimator_comparison = function(
  aipw_results,
  drs_results,
  comparison_label,
  k_label,
  save_path,
  width = 3000,
  height = 2000,
  res = 300
) {
  #' Compare AIPW and DRS bootstrap distributions
  #'
  #' Creates a 2x4 plot: rows are estimators, columns are models

  # Extract bootstrap samples
  aipw_samples = lapply(aipw_results, function(x) x$aipw$bootstrap_samples)
  drs_samples = lapply(drs_results, function(x) x$drs$bootstrap_samples)

  # Extract point estimates
  aipw_estimates = sapply(aipw_results, function(x) x$aipw$att)
  drs_estimates = sapply(drs_results, function(x) x$drs$att)

  # Create plot function
  create_plot = function() {
    par(mfrow = c(2, 4), mar = c(4, 4, 3, 1))

    # Top row: AIPW
    for (model in 1:4) {
      dens = density(aipw_samples[[model]])

      hist(
        aipw_samples[[model]],
        breaks = 30,
        col = "lightblue",
        border = "white",
        main = sprintf("Model %d - AIPW", model),
        xlab = "ATT",
        ylab = "Frequency"
      )
      abline(v = aipw_estimates[model], col = "blue", lwd = 2, lty = 2)
    }

    # Bottom row: DRS
    for (model in 1:4) {
      dens = density(drs_samples[[model]])

      hist(
        drs_samples[[model]],
        breaks = 30,
        col = "lightgreen",
        border = "white",
        main = sprintf("Model %d - DRS", model),
        xlab = "ATT",
        ylab = "Frequency"
      )
      abline(v = drs_estimates[model], col = "darkgreen", lwd = 2, lty = 2)
    }

    # Add overall title
    mtext(
      sprintf("%s - %s", comparison_label, k_label),
      side = 3,
      line = -1.5,
      outer = TRUE,
      cex = 1.2
    )

    par(mfrow = c(1, 1))
  }

  # Display plot
  create_plot()

  # Save to file
  png(save_path, width = width, height = height, res = res)
  create_plot()
  dev.off()

  invisible(NULL)
}

# Create AIPW vs DRS plots for each comparison and k
comparisons = list(
  inc_dec = list(
    label = "Increase vs Decrease",
    aipw = aipw_workspace$inc_dec,
    drs = drs_workspace$inc_dec
  ),
  inc_mix = list(
    label = "Increase vs Mix",
    aipw = aipw_workspace$inc_mix,
    drs = drs_workspace$inc_mix
  )
)

for (comp_name in names(comparisons)) {
  comp = comparisons[[comp_name]]

  cat(sprintf("\nCreating plots for %s...\n", comp$label))

  for (k in c(1, 2, 3)) {
    # Extract results for this k
    k_pattern = sprintf("^k%d_", k)
    aipw_k = comp$aipw[grep(k_pattern, names(comp$aipw))]
    drs_k = comp$drs[grep(k_pattern, names(comp$drs))]

    # Create plot
    k_label = if (k == 1) "k=1 (No CF)" else sprintf("k=%d", k)
    save_path = sprintf(
      "results/outcome/combined/AIPW_vs_DRS_%s_k%d.png",
      comp_name,
      k
    )

    plot_estimator_comparison(
      aipw_results = aipw_k,
      drs_results = drs_k,
      comparison_label = comp$label,
      k_label = k_label,
      save_path = save_path
    )

    cat(sprintf("  ✓ Saved: %s\n", save_path))
  }
}

# =============================================================================
# MODEL COMPLEXITY EFFECTS
# =============================================================================

cat("\n\n")
cat("=" %>% rep(80) %>% paste0(collapse = ""), "\n")
cat("ANALYZING MODEL COMPLEXITY EFFECTS\n")
cat("=" %>% rep(80) %>% paste0(collapse = ""), "\n\n")

# Function to plot how ATT changes with model complexity
plot_model_complexity = function(
  summary_data,
  save_path,
  width = 3000,
  height = 2000,
  res = 300
) {
  #' Plot ATT estimates across model complexity
  #'
  #' Shows how estimates change from model 1 to 4

  # Create separate plots for each comparison and k
  create_plot = function() {
    par(mfrow = c(2, 3), mar = c(4, 4, 3, 2))

    # Loop through comparisons and k values
    for (comp in c("inc_dec", "inc_mix")) {
      comp_label = if (comp == "inc_dec") "Inc vs Dec" else "Inc vs Mix"

      for (k in c(1, 2, 3)) {
        # Subset data
        data_subset = summary_data[
          summary_data$Comparison == comp &
            summary_data$K_Fold == k,
        ]

        # Separate AIPW and DRS
        aipw_data = data_subset[data_subset$Estimator == "AIPW", ]
        drs_data = data_subset[data_subset$Estimator == "DRS", ]

        # Sort by model
        aipw_data = aipw_data[order(aipw_data$Model), ]
        drs_data = drs_data[order(drs_data$Model), ]

        # Plot
        plot(
          aipw_data$Model,
          aipw_data$ATT,
          type = "b",
          col = "blue",
          lwd = 2,
          pch = 16,
          ylim = range(c(
            aipw_data$CI_Perc_Lower,
            aipw_data$CI_Perc_Upper,
            drs_data$CI_Perc_Lower,
            drs_data$CI_Perc_Upper
          )),
          xlab = "Model",
          ylab = "ATT",
          main = sprintf("%s, k=%d", comp_label, k),
          xaxt = "n"
        )
        axis(1, at = 1:4)

        # Add confidence intervals for AIPW
        arrows(
          aipw_data$Model,
          aipw_data$CI_Perc_Lower,
          aipw_data$Model,
          aipw_data$CI_Perc_Upper,
          code = 3,
          angle = 90,
          length = 0.05,
          col = "blue"
        )

        # Add DRS line
        lines(
          drs_data$Model,
          drs_data$ATT,
          type = "b",
          col = "darkgreen",
          lwd = 2,
          pch = 17
        )

        # Add confidence intervals for DRS
        arrows(
          drs_data$Model,
          drs_data$CI_Perc_Lower,
          drs_data$Model,
          drs_data$CI_Perc_Upper,
          code = 3,
          angle = 90,
          length = 0.05,
          col = "darkgreen"
        )

        # Add horizontal line at 0
        abline(h = 0, lty = 3, col = "gray")

        # Add legend to first panel only
        if (comp == "inc_dec" && k == 1) {
          legend(
            "topright",
            legend = c("AIPW", "DRS"),
            col = c("blue", "darkgreen"),
            lwd = 2,
            pch = c(16, 17),
            cex = 0.8
          )
        }
      }
    }

    par(mfrow = c(1, 1))
  }

  # Display plot
  create_plot()

  # Save to file
  png(save_path, width = width, height = height, res = res)
  create_plot()
  dev.off()

  invisible(NULL)
}

cat("Creating model complexity plot...\n")
plot_model_complexity(
  summary_data = master_summary_all,
  save_path = "results/outcome/combined/model_complexity_effects.png"
)
cat("✓ Saved: results/outcome/combined/model_complexity_effects.png\n")

# =============================================================================
# CROSS-FITTING EFFECTS
# =============================================================================

cat("\n\n")
cat("=" %>% rep(80) %>% paste0(collapse = ""), "\n")
cat("ANALYZING CROSS-FITTING EFFECTS\n")
cat("=" %>% rep(80) %>% paste0(collapse = ""), "\n\n")

# Function to plot cross-fitting effects
plot_crossfitting_effects = function(
  summary_data,
  save_path,
  width = 3000,
  height = 2000,
  res = 300
) {
  #' Plot how estimates change with cross-fitting
  #'
  #' Shows k=1 vs k=2 vs k=3 for each model

  create_plot = function() {
    par(mfrow = c(2, 4), mar = c(4, 4, 3, 2))

    # Loop through comparisons
    for (comp in c("inc_dec", "inc_mix")) {
      comp_label = if (comp == "inc_dec") "Inc vs Dec" else "Inc vs Mix"

      # Loop through models
      for (model in 1:4) {
        # Subset data
        data_subset = summary_data[
          summary_data$Comparison == comp &
            summary_data$Model == model,
        ]

        # Separate AIPW and DRS
        aipw_data = data_subset[data_subset$Estimator == "AIPW", ]
        drs_data = data_subset[data_subset$Estimator == "DRS", ]

        # Sort by k
        aipw_data = aipw_data[order(aipw_data$K_Fold), ]
        drs_data = drs_data[order(drs_data$K_Fold), ]

        # Plot
        plot(
          aipw_data$K_Fold,
          aipw_data$ATT,
          type = "b",
          col = "blue",
          lwd = 2,
          pch = 16,
          ylim = range(c(
            aipw_data$CI_Perc_Lower,
            aipw_data$CI_Perc_Upper,
            drs_data$CI_Perc_Lower,
            drs_data$CI_Perc_Upper
          )),
          xlab = "K-fold",
          ylab = "ATT",
          main = sprintf("%s, Model %d", comp_label, model),
          xaxt = "n"
        )
        axis(1, at = c(1, 2, 3))

        # Add confidence intervals for AIPW
        arrows(
          aipw_data$K_Fold,
          aipw_data$CI_Perc_Lower,
          aipw_data$K_Fold,
          aipw_data$CI_Perc_Upper,
          code = 3,
          angle = 90,
          length = 0.05,
          col = "blue"
        )

        # Add DRS line
        lines(
          drs_data$K_Fold,
          drs_data$ATT,
          type = "b",
          col = "darkgreen",
          lwd = 2,
          pch = 17
        )

        # Add confidence intervals for DRS
        arrows(
          drs_data$K_Fold,
          drs_data$CI_Perc_Lower,
          drs_data$K_Fold,
          drs_data$CI_Perc_Upper,
          code = 3,
          angle = 90,
          length = 0.05,
          col = "darkgreen"
        )

        # Add horizontal line at 0
        abline(h = 0, lty = 3, col = "gray")

        # Add legend to first panel only
        if (comp == "inc_dec" && model == 1) {
          legend(
            "topright",
            legend = c("AIPW", "DRS"),
            col = c("blue", "darkgreen"),
            lwd = 2,
            pch = c(16, 17),
            cex = 0.8
          )
        }
      }
    }

    par(mfrow = c(1, 1))
  }

  # Display plot
  create_plot()

  # Save to file
  png(save_path, width = width, height = height, res = res)
  create_plot()
  dev.off()

  invisible(NULL)
}

cat("Creating cross-fitting effects plot...\n")
plot_crossfitting_effects(
  summary_data = master_summary_all,
  save_path = "results/outcome/combined/crossfitting_effects.png"
)
cat("✓ Saved: results/outcome/combined/crossfitting_effects.png\n")

# =============================================================================
# SUMMARY STATISTICS
# =============================================================================

cat("\n\n")
cat("=" %>% rep(80) %>% paste0(collapse = ""), "\n")
cat("SUMMARY STATISTICS\n")
cat("=" %>% rep(80) %>% paste0(collapse = ""), "\n\n")

# Overall summary by estimator
cat("=== Summary by Estimator ===\n\n")
aggregate_by_estimator = aggregate(
  ATT ~ Estimator,
  data = master_summary_all,
  FUN = function(x) {
    c(
      Mean = mean(x),
      SD = sd(x),
      Min = min(x),
      Max = max(x)
    )
  }
)
print(aggregate_by_estimator)

cat("\n\n=== Summary by Comparison ===\n\n")
aggregate_by_comparison = aggregate(
  ATT ~ Comparison + Estimator,
  data = master_summary_all,
  FUN = function(x) {
    c(
      Mean = mean(x),
      SD = sd(x),
      Min = min(x),
      Max = max(x)
    )
  }
)
print(aggregate_by_comparison)

cat("\n\n=== Summary by K-Fold ===\n\n")
aggregate_by_k = aggregate(
  ATT ~ K_Fold + Estimator,
  data = master_summary_all,
  FUN = function(x) {
    c(
      Mean = mean(x),
      SD = sd(x),
      Min = min(x),
      Max = max(x)
    )
  }
)
print(aggregate_by_k)

cat("\n\n=== Summary by Model ===\n\n")
aggregate_by_model = aggregate(
  ATT ~ Model + Estimator,
  data = master_summary_all,
  FUN = function(x) {
    c(
      Mean = mean(x),
      SD = sd(x),
      Min = min(x),
      Max = max(x)
    )
  }
)
print(aggregate_by_model)

# =============================================================================
# FINAL OUTPUT
# =============================================================================

cat("\n\n")
cat("=" %>% rep(80) %>% paste0(collapse = ""), "\n")
cat("COMBINED ANALYSIS COMPLETE!\n")
cat("=" %>% rep(80) %>% paste0(collapse = ""), "\n\n")

cat("Summary of outputs:\n")
cat("  - Master summary (48 analyses): results/outcome/combined/master_summary_all_48.csv\n")
cat("  - AIPW vs DRS plots: results/outcome/combined/AIPW_vs_DRS_*.png\n")
cat("  - Model complexity plot: results/outcome/combined/model_complexity_effects.png\n")
cat("  - Cross-fitting effects plot: results/outcome/combined/crossfitting_effects.png\n")
cat("\n")

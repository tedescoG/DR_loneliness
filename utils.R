# HELPER FUNCTIONS FOR ANALYSIS (REFACTORED)
# Unified bootstrap and doubly robust estimation functions

# Libraries
library(tidyverse)
library(gbm)
library(twang)
library(survey)
library(marginaleffects)
library(boot)
library(rsample)


# ASD ----------------------------------------------------------------------------------####
# Compute standardize absolute mean
std.diff = function(u, z, w) {
  #' Compute standardized absolute mean difference for covariate balance
  #'
  #' @param u numeric or factor - covariate to assess balance
  #' @param z numeric - binary treatment indicator (0/1)
  #' @param w numeric - propensity score weights
  #'
  #' @return numeric - absolute standardized difference (or vector for factors)

  # for variables other than unordered categorical variables
  # compute mean differences
  if (!is.factor(u)) {
    sd1 = sd(u[z == 1], na.rm = T) # calculates the standard deviation
    # for the treatment group
    if (sd1 > 0) {
      # calculate the absolute standardized mean difference
      result = abs(
        mean(u[z == 1], na.rm = TRUE) -
          weighted.mean(u[z == 0], w[z == 0], na.rm = TRUE)
      ) /
        sd1
    } else {
      result = 0
      warning("Covariate with standard deviation 0.")
    }
  } else {
    # for factors compute differences in percentages in each category
    result = NULL
    for (u.level in levels(u)) {
      # calculate the absolute standardized difference of the indicator variable
      result = c(result, std.diff(as.numeric(u == u.level), z, w))
    }
  }
  return(result)
}

# Average ASD across covariates
es.mean = function(i, gbm1, x, z, verbose) {
  #' Compute average absolute standardized difference across covariates
  #'
  #' @param i numeric - GBM iteration number
  #' @param gbm1 gbm.object - fitted GBM model
  #' @param x data.frame - covariate matrix
  #' @param z numeric - binary treatment indicator (0/1)
  #' @param verbose numeric - verbosity level (0/1/2)
  #'
  #' @return numeric - average standardized difference across all covariates

  # prints the iteration number
  i = floor(i) # makes sure that i is an integer

  # predict(gbm1, x, i) provides predicted values on the log-odds of
  # treatment for the gbm model with i iterations at the values of x
  # exp(predict(gbm1, x, i)) calculates the odds of treatment or the weight
  # from the predicted values
  w = exp(predict(gbm1, x, i))
  # assign treatment cases a weight of 1
  w[z == 1] = 1
  # sapply repeats calculation of std.diff for each variable (column) of x
  # unlist is an R function for managing data structures
  # mean(unlist(sapply(x, std.diff, z = z, w = w))) calculates the mean of the
  # standardized differences for all variables in x or ASAM
  if (verbose > 1) {
    cat(
      "\nASAM at iteration:",
      i,
      "=",
      mean(unlist(sapply(x, std.diff, z = z, w = w)))
    )
  }
  return(mean(unlist(sapply(x, std.diff, z = z, w = w))))
}

# Max ASD across covariates
es.max = function(i, gbm1, x, z, verbose) {
  #' Compute maximum absolute standardized difference across covariates
  #'
  #' @param i numeric - GBM iteration number
  #' @param gbm1 gbm.object - fitted GBM model
  #' @param x data.frame - covariate matrix
  #' @param z numeric - binary treatment indicator (0/1)
  #' @param verbose numeric - verbosity level (0/1/2)
  #'
  #' @return numeric - maximum standardized difference across all covariates

  # prints the iteration number
  i = floor(i) # makes sure that i is an integer

  # predict(gbm1, x, i) provides predicted values on the log-odds of
  # treatment for the gbm model with i iterations at the values of x
  # exp(predict(gbm1, x, i)) calculates the odds of treatment or the weight
  # from the predicted values
  w = exp(predict(gbm1, x, i))
  # assign treatment cases a weight of 1
  w[z == 1] = 1
  # sapply repeats calculation of std.diff for each variable (column) of x
  # unlist is an R function for managing data structures
  # mean(unlist(sapply(x, std.diff, z = z, w = w))) calculates the mean of the
  # standardized differences for all variables in x or ASAM
  if (verbose > 1) {
    cat(
      "\nASAM at iteration:",
      i,
      "=",
      mean(unlist(sapply(x, std.diff, z = z, w = w)))
    )
  }
  return(max(unlist(sapply(x, std.diff, z = z, w = w))))
}

# KS -----------------------------------------------------------------------------------####
# Compute KS
compute_KS = function(variable, treatment, weights) {
  #' Compute Kolmogorov-Smirnov statistic for covariate balance
  #'
  #' @param variable numeric - covariate values
  #' @param treatment numeric - binary treatment indicator (0/1)
  #' @param weights numeric - propensity score weights
  #'
  #' @return numeric - KS statistic (max difference between weighted CDFs)

  treatment_data = variable[treatment == 1]
  control_data = variable[treatment == 0]

  treatment_weights = weights[treatment == 1]
  control_weights = weights[treatment == 0]

  EDF1 = function(x) {
    sum(treatment_weights[treatment_data <= x]) / sum(treatment_weights)
  }
  EDF0 = function(x) {
    sum(control_weights[control_data <= x]) / sum(control_weights)
  }

  # Compute KS over a range of values in 'variable'
  ks_values = sapply(sort(unique(variable)), function(x) {
    abs(EDF1(x) - EDF0(x))
  })
  return(max(ks_values))
}

# Average KS across covariates
ks.mean = function(i, gbm1, x, z, verbose) {
  #' Compute average Kolmogorov-Smirnov statistic across covariates
  #'
  #' @param i numeric - GBM iteration number
  #' @param gbm1 gbm.object - fitted GBM model
  #' @param x data.frame - covariate matrix
  #' @param z numeric - binary treatment indicator (0/1)
  #' @param verbose numeric - verbosity level (0/1/2)
  #'
  #' @return numeric - average KS statistic across all covariates

  i = floor(i) # Ensure i is an integer

  # Compute propensity score weights
  w = exp(predict(gbm1, x, i))
  w[z == 1] = 1

  # Compute KS statistic for each covariate and take the mean
  # Return the max of the ks statistics for all variables
  ks_stats = sapply(x, function(covariate) compute_KS(covariate, z, w))
  # Take the average
  avg_ks = mean(ks_stats, na.rm = T)
  if (verbose > 1) {
    cat("Average KS at iteration:", i, "=", avg_ks, "\n")
  }
  return(avg_ks)
}

# Max KS across covariates
ks.max = function(i, gbm1, x, z, verbose) {
  #' Compute maximum Kolmogorov-Smirnov statistic across covariates
  #'
  #' @param i numeric - GBM iteration number
  #' @param gbm1 gbm.object - fitted GBM model
  #' @param x data.frame - covariate matrix
  #' @param z numeric - binary treatment indicator (0/1)
  #' @param verbose numeric - verbosity level (0/1/2)
  #'
  #' @return numeric - maximum KS statistic across all covariates

  i = floor(i) # Ensure i is an integer

  # Compute propensity score weights
  w = exp(predict(gbm1, x, i))
  w[z == 1] = 1

  # Compute KS statistic for each covariate and take the mean
  # Return the max of the ks statistics for all variables
  ks_stats = sapply(x, function(covariate) compute_KS(covariate, z, w))
  # Take the average
  max_ks = max(ks_stats, na.rm = T)
  if (verbose > 1) {
    cat("Average KS at iteration:", i, "=", max_ks, "\n")
  }
  return(max_ks)
}

# HYPER-PARAMETERS TUNING FUNC. --------------------------------------------------------####

# Helper function to process single combination of GBM's hyperparameters
process_tune_combination = function(
  i,
  param_combinations,
  x,
  y,
  method,
  seed
) {
  #' Process single hyperparameter combination for GBM tuning
  #'
  #' @param i integer - index of parameter combination to test
  #' @param param_combinations data.frame - grid of parameter combinations
  #' @param x data.frame - covariate matrix
  #' @param y numeric - binary treatment indicator (0/1)
  #' @param method character - balance metric ("es.mean", "es.max", "ks.mean", "ks.max")
  #'
  #' @return list with combination info, best balance value, and optimal iterations

  combination = param_combinations[i, ]
  gbm_params = list()

  # Build GBM parameters
  if ("distribution" %in% names(combination)) {
    gbm_params$distribution = as.character(combination$distribution)
  } else {
    gbm_params$distribution = "bernoulli" # Default for binary treatment
  }
  if ("n.trees" %in% names(combination)) {
    gbm_params$n.trees = combination$n.trees
    N_local = gbm_params$n.trees
  } else {
    N_local = 20000
  }
  if ("shrinkage" %in% names(combination)) {
    gbm_params$shrinkage = combination$shrinkage
  }
  if ("interaction.depth" %in% names(combination)) {
    gbm_params$interaction.depth = combination$interaction.depth
  }
  if ("bag.fraction" %in% names(combination)) {
    gbm_params$bag.fraction = combination$bag.fraction
  }
  if ("n.minobsinnode" %in% names(combination)) {
    gbm_params$n.minobsinnode = combination$n.minobsinnode
  }

  # Fit GBM model on full dataset
  model = do.call(
    "gbm.fit",
    c(list(x = x, y = y), gbm_params, verbose = FALSE)
  )

  # Select optimization function based on method
  method_matched = method
  optimization_function = ifelse(
    method_matched == "es.mean",
    es.mean,
    ifelse(
      method_matched == "es.max",
      es.max,
      ifelse(
        method_matched == "ks.mean",
        ks.mean,
        ks.max
      )
    )
  )

  # Find optimal number of iterations
  opt = optimize(
    optimization_function,
    interval = c(100, N_local),
    tol = 1,
    gbm1 = model,
    x = x,
    z = y,
    verbose = FALSE # Suppress verbose in parallel
  )

  return(list(
    combination_num = i,
    combination = combination,
    best_value = opt$objective,
    best_iter = opt$minimum
  ))
}

# Function to tune the GBM "chasing balance" (Griffin et al. 2017, Hirano et al. 2003)
tune.gbm = function(
  x,
  y,
  params.grid,
  method = c("es.mean", "es.max", "ks.mean", "ks.max"),
  verbose = 1,
  seed = 123,
  n.cores = NULL,
  parallel = TRUE
) {
  #' Tune GBM hyperparameters optimizing for covariate balance
  #'
  #' @param x data.frame - covariate matrix
  #' @param y numeric - binary treatment indicator (0/1)
  #' @param params.grid list - hyperparameters to test (n.trees, shrinkage, etc.)
  #' @param method character - balance metric to optimize
  #' @param verbose numeric - verbosity level (0/1/2)
  #' @param seed integer - random seed for reproducibility
  #' @param n.cores integer - number of cores for parallel processing (NULL = auto)
  #' @param parallel logical - enable parallel processing
  #'
  #' @return list with best parameters, performance, optimal iterations, and method
  # Load parallel library if requested
  if (parallel && !requireNamespace("parallel", quietly = TRUE)) {
    warning(
      "parallel package not available, falling back to sequential processing"
    )
    parallel = FALSE
  }

  # Core detection and setup
  if (parallel) {
    if (is.null(n.cores)) {
      n.cores = max(1, parallel::detectCores() - 1) # Leave one core free
    }
    n.cores = min(n.cores, parallel::detectCores())
    if (verbose >= 1) {
      cat("Using", n.cores, "cores for parallel processing\n")
    }
  }

  # Ensure y is numeric binary (0/1) and remove any names/attributes
  y = as.vector(as.numeric(y))
  if (max(y) > 1) {
    y = y - 1
  }

  # Generate all parameter combinations
  param_combinations = expand.grid(params.grid)

  if (verbose >= 1) {
    cat(
      "Starting hyperparameter tuning...\n",
      "Testing",
      nrow(param_combinations),
      "combinations\n",
      "Balance metric:",
      method,
      "\n\n"
    )
  }

  # Match method argument
  method = match.arg(method)
  RNGkind("L'Ecuyer-CMRG")
  set.seed(seed)
  # Execute parallel or sequential processing
  if (parallel && nrow(param_combinations) > 1) {
    # Parallel processing
    results = parallel::mclapply(
      1:nrow(param_combinations),
      function(i) {
        process_tune_combination(
          i,
          param_combinations,
          x,
          y,
          method
        )
      },
      mc.cores = n.cores,
      mc.set.seed = T
    )
  } else {
    # Sequential processing
    results = lapply(
      1:nrow(param_combinations),
      function(i) {
        process_tune_combination(
          i,
          param_combinations,
          x,
          y,
          method,
          seed = seed
        )
      }
    )
  }

  # Find best combination
  best_performance = Inf
  best_combination = NULL
  best_combination_num = NULL
  best_iteration = NULL

  for (result in results) {
    if (verbose >= 2) {
      cat(
        "Combination",
        result$combination_num,
        ":",
        method,
        "=",
        round(result$best_value, 4),
        "at iteration",
        result$best_iter,
        "\n"
      )
    }

    if (result$best_value < best_performance) {
      best_performance = result$best_value
      best_combination = result$combination
      best_combination_num = result$combination_num
      best_iteration = result$best_iter
    }
  }

  if (verbose >= 1) {
    cat(
      "\n=== BEST HYPERPARAMETERS ===\n",
      "Combination:",
      best_combination_num,
      "\n",
      method,
      ":",
      round(best_performance, 4),
      "\n",
      "Optimal iterations:",
      best_iteration,
      "\n"
    )
    cat("\nParameters:\n")
    print(best_combination)
  }

  # Return best parameters and performance
  return(list(
    params = as.list(best_combination),
    performance = best_performance,
    optimal_iter = best_iteration,
    method = method
  ))
}


# DRS Estimator -----------------------------------------------------------------------####
drs_att = function(outcome, treatment, f.out, wgt, data, verbose = T) {
  #' Compute ATT using g-computation (doubly robust standardization)
  #'
  #' @param outcome character string - name of outcome variable
  #' @param treatment character string - name of treatment variable (must be 0/1)
  #' @param f.out character string - formula for covariates (without outcome/treatment)
  #' @param wgt character string - name of weight variable (propensity score weights)
  #' @param data data.frame - dataset containing all variables
  #' @param verbose logical - whether to print results
  #'
  #' @return list with att estimate and predicted counterfactual outcomes

  # Extract variables
  Y = data %>% pull(all_of(outcome))
  Z = data %>% pull(all_of(treatment))

  # Create survey design with propensity score weights
  weighted_design = svydesign(
    ids = ~1,
    weights = as.formula(paste0("~", wgt)),
    data = data
  )

  # Build formula: outcome ~ treatment + covariates
  if (f.out != "1") {
    formula_str = paste(outcome, "~", treatment, "+", f.out)
  } else {
    formula_str = paste(outcome, "~", treatment)
  }

  # Fit weighted outcome model on full sample
  fit = svyglm(
    formula = as.formula(formula_str),
    family = 'quasibinomial',
    design = weighted_design
  )

  # Subset treated units
  treated_data = data %>% filter(!!sym(treatment) == 1)

  # Create counterfactual dataset: set treatment to 0 for treated units
  counterfactual_data = treated_data %>%
    mutate(!!sym(treatment) := 0)

  # Predict counterfactual outcomes mu_hat_0 for treated under control
  mu0 = predict(fit, newdata = counterfactual_data, type = "response")
  mu_hat_0 = mean(mu0)
  # Observed outcome mean for treated
  mu_hat_1 = mean(Y[Z == 1])

  # Compute ATT: difference between observed and counterfactual
  att_est = mu_hat_1 - mu_hat_0

  # Output
  if (verbose) {
    cat(sprintf("ATT estimate (g-computation): %.4f\n", att_est))
    cat(sprintf("Mean Y(0) [counterfactual]: %.4f\n", mu_hat_0))
  }

  # Return results
  return(list(
    att = att_est,
    mu_hat_0 = mu_hat_0,
    mu_hat_1 = mu_hat_1
  ))
}

# Outcome Regression ATT (without propensity score weights) ---------------------------####
or_att = function(outcome, treatment, f.out, data, verbose = T) {
  #' Compute ATT using outcome regression (g-computation without PS weights)
  #'
  #' @param outcome character string - name of outcome variable
  #' @param treatment character string - name of treatment variable (must be 0/1)
  #' @param f.out character string - formula for covariates (without outcome/treatment)
  #' @param data data.frame - dataset containing all variables
  #' @param verbose logical - whether to print results
  #'
  #' @return list with att estimate and predicted counterfactual outcomes
  #'
  #' @details This is the same as drs_att but without propensity score weighting.
  #' Used as Model 0 baseline to show the value of PS weighting.

  # Extract variables
  Y = data %>% pull(all_of(outcome))
  Z = data %>% pull(all_of(treatment))

  # Build formula: outcome ~ treatment + covariates
  if (f.out != "1") {
    formula_str = paste(outcome, "~", treatment, "+", f.out)
  } else {
    formula_str = paste(outcome, "~", treatment)
  }

  # Fit UNWEIGHTED outcome model on full sample
  fit = glm(
    formula = as.formula(formula_str),
    family = 'binomial',
    data = data
  )

  # Subset treated units
  treated_data = data %>% filter(!!sym(treatment) == 1)

  # Create counterfactual dataset: set treatment to 0 for treated units
  counterfactual_data = treated_data %>%
    mutate(!!sym(treatment) := 0)

  # Predict counterfactual outcomes mu_hat_0 for treated under control
  mu0 = predict(fit, newdata = counterfactual_data, type = "response")
  mu_hat_0 = mean(mu0)
  # Observed outcome mean for treated
  mu_hat_1 = mean(Y[Z == 1])

  # Compute ATT: difference between observed and counterfactual
  att_est = mu_hat_1 - mu_hat_0

  # Output
  if (verbose) {
    cat(sprintf("ATT estimate (outcome regression): %.4f\n", att_est))
    cat(sprintf("Mean Y(0) [counterfactual]: %.4f\n", mu_hat_0))
  }

  # Return results
  return(list(
    att = att_est,
    mu_hat_0 = mu_hat_0,
    mu_hat_1 = mu_hat_1
  ))
}

# Simple Difference in Means ATT (no covariates, no PS weights) -----------------------####
simple_att = function(outcome, treatment, data, verbose = T) {
  #' Compute ATT using simple difference in means (no covariates, no weights)
  #'
  #' @param outcome character string - name of outcome variable
  #' @param treatment character string - name of treatment variable (must be 0/1)
  #' @param data data.frame - dataset containing all variables
  #' @param verbose logical - whether to print results
  #'
  #' @return list with att estimate and group means
  #'
  #' @details Simple difference in means: E[Y|Z=1] - E[Y|Z=0]
  #' This is the most basic estimator with no covariate adjustment or weighting.

  # Extract variables
  Y = data %>% pull(all_of(outcome))
  Z = data %>% pull(all_of(treatment))

  # Compute group means
  mu_hat_1 = mean(Y[Z == 1]) # Treated mean
  mu_hat_0 = mean(Y[Z == 0]) # Control mean

  # Compute ATT: simple difference
  att_est = mu_hat_1 - mu_hat_0

  # Output
  if (verbose) {
    cat(sprintf("ATT estimate (simple difference): %.4f\n", att_est))
    cat(sprintf("Mean Y(1) [treated]: %.4f\n", mu_hat_1))
    cat(sprintf("Mean Y(0) [control]: %.4f\n", mu_hat_0))
  }

  # Return results (matching or_att and drs_att structure)
  return(list(
    att = att_est,
    mu_hat_0 = mu_hat_0,
    mu_hat_1 = mu_hat_1
  ))
}


# UNIFIED BOOTSTRAP FUNCTION -----------------------------------------------------------####

boot_iter = function(
  data,
  indices,
  outcome,
  treatment,
  treated_level,
  control_level,
  f.ps = NULL,
  f.out,
  gbm_params = NULL,
  bootstrap_method = c("resample", "reweight"),
  wgt = "iptw",
  seed_offset,
  alpha = 0,
  use_weights = TRUE,
  estimator = c("drs", "or", "simple")
) {
  #' Bootstrap iteration function for Simple/OR/DRS estimation
  #'
  #' Supports multiple estimators and bootstrap methods
  #'
  #' @param data data.frame - original dataset
  #' @param indices integer - bootstrap sample indices (generated by boot())
  #' @param outcome character - name of outcome variable
  #' @param treatment character - name of treatment variable
  #' @param treated_level string - level indicating treated group
  #' @param control_level string - level indicating control group
  #' @param f.ps formula - propensity score model formula (for reweight method)
  #' @param f.out character - outcome model formula (RHS only)
  #' @param gbm_params list - GBM parameters for ps() function (for reweight method)
  #' @param bootstrap_method character - "resample" (use existing weights) or "reweight" (re-estimate PS)
  #' @param wgt character - name of weight variable (for resample method)
  #' @param seed_offset integer - base seed for reproducibility
  #' @param alpha numeric - propensity score truncation level (0 = no truncation)
  #' @param use_weights logical - if FALSE, use or_att() instead of drs_att() (deprecated, use estimator instead)
  #' @param estimator character - "simple" (difference in means), "or" (outcome regression), or "drs" (doubly robust)
  #'
  #' @return numeric vector: c(att_estimate, avg_asd, max_asd, ess)

  # Wrap entire function in tryCatch for error handling
  tryCatch(
    {
      # Match arguments
      bootstrap_method = match.arg(bootstrap_method)

      # Resample data
      boot_data = data[indices, ]

      # Create binary treatment indicator
      boot_data[[treatment]] = as.numeric(
        boot_data[[treatment]] == treated_level
      )

      # Initialize weights vector
      weights_vec = NULL

      # Determine which estimator to use
      if (!missing(estimator)) {
        estimator = match.arg(estimator)
      } else {
        # Backward compatibility: use use_weights to determine estimator
        estimator = if (!use_weights) "or" else "drs"
      }

      # Route to appropriate estimator
      if (estimator == "simple") {
        # Simple difference in means (no weights, no covariates)
        simple_result = simple_att(
          outcome = outcome,
          treatment = treatment,
          data = boot_data,
          verbose = FALSE
        )

        drs_est = simple_result$att
        weights_vec = NULL
      } else if (estimator == "or") {
        # Outcome regression without weights
        or_result = or_att(
          outcome = outcome,
          treatment = treatment,
          f.out = f.out,
          data = boot_data,
          verbose = FALSE
        )

        drs_est = or_result$att
        weights_vec = NULL
      } else if (estimator == "drs" && bootstrap_method == "resample") {
        # Use existing weights from data
        drs_result = drs_att(
          outcome = outcome,
          treatment = treatment,
          f.out = f.out,
          wgt = wgt,
          data = boot_data,
          verbose = FALSE
        )

        drs_est = drs_result$att
        weights_vec = boot_data[[wgt]]
      } else if (estimator == "drs" && bootstrap_method == "reweight") {
        # Reweight: re-estimate propensity scores
        set.seed(seed_offset)
        ps_fit = do.call(
          "ps",
          c(
            list(
              formula = f.ps,
              data = boot_data,
              estimand = "ATT",
              stop.method = "es.mean",
              verbose = FALSE
            ),
            gbm_params
          )
        )

        # Extract weights
        boot_data$ps_wgt = get.weights(ps_fit, stop.method = "es.mean")

        # Compute DRS estimate
        drs_result = drs_att(
          outcome = outcome,
          treatment = treatment,
          f.out = f.out,
          wgt = "ps_wgt",
          data = boot_data,
          verbose = FALSE
        )

        drs_est = drs_result$att
        weights_vec = boot_data$ps_wgt
      }

      # ========== COMPUTE BALANCE STATISTICS ==========
      avg_asd = NA_real_
      max_asd = NA_real_
      ess = NA_real_

      tryCatch(
        {
          # Extract all covariates (exclude treatment, outcome, and weight columns)
          exclude_cols = c(treatment, outcome)
          if (bootstrap_method == "resample") {
            exclude_cols = c(exclude_cols, wgt)
          }
          if (bootstrap_method == "reweight") {
            exclude_cols = c(exclude_cols, "ps_wgt")
          }

          covariate_cols = setdiff(names(boot_data), exclude_cols)

          if (length(covariate_cols) > 0 && !is.null(weights_vec)) {
            covariates = boot_data[, covariate_cols, drop = FALSE]

            # Ensure treatment is numeric 0/1
            treat_vec = as.numeric(boot_data[[treatment]])

            # Ensure weights are numeric
            weights_vec = as.numeric(weights_vec)

            # Check for valid weights
            if (all(is.finite(weights_vec)) && all(weights_vec > 0)) {
              # Compute standardized differences
              std_diffs = unlist(sapply(
                covariates,
                std.diff,
                z = treat_vec,
                w = weights_vec
              ))
              std_diffs = std_diffs[!is.na(std_diffs)]

              if (length(std_diffs) > 0) {
                avg_asd = mean(std_diffs)
                max_asd = max(std_diffs)
              }

              # Compute effective sample size for control group
              # ESS = (sum of weights)^2 / sum of weights^2
              control_weights = weights_vec[treat_vec == 0]
              ess = sum(control_weights)^2 / sum(control_weights^2)
            }
          }
        },
        error = function(e) {
          # If balance calculation fails, leave as NA
        }
      )

      # Return results: DRS estimate + balance metrics
      return(c(
        drs = drs_est,
        avg_asd = avg_asd,
        max_asd = max_asd,
        ess = ess
      ))
    },
    error = function(e) {
      # Return NAs if anything fails during bootstrap iteration
      return(c(
        drs = NA_real_,
        avg_asd = NA_real_,
        max_asd = NA_real_,
        ess = NA_real_
      ))
    }
  )
}


# DR ATT ANALYSIS FUNCTION -------------------------------------------------------------####

DR_att = function(
  outcome,
  treatment,
  treated_level,
  control_level,
  f.ps,
  f.out,
  data,
  gbm_params = list(
    n.trees = 10000,
    interaction.depth = 2,
    shrinkage = 0.01,
    bag.fraction = 0.5,
    n.minobsinnode = 10
  ),
  bootstrap_method = c("resample", "reweight"),
  stratified = TRUE,
  wgt = "iptw",
  n_boot = 1000,
  seed = 123,
  verbose = TRUE,
  parallel = TRUE,
  n_cores = NULL,
  plot_diagnostics = TRUE,
  ci_type = c("all", "norm", "basic", "perc"),
  sim = c("ordinary", "balanced"),
  save_to = NULL,
  alpha = 0,
  use_weights = TRUE,
  estimator = c("drs", "or", "simple")
) {
  #' Bootstrap ATT Estimation for Simple/OR/DRS estimators
  #'
  #' Supports multiple estimators with bootstrap inference:
  #' - simple: simple difference in means (no covariates, no weights)
  #' - or: outcome regression with covariates (no propensity scores)
  #' - drs: doubly robust with propensity scores and outcome modeling
  #'
  #' @param outcome character - name of outcome variable
  #' @param treatment character - name of treatment variable
  #' @param treated_level string - level of treatment variable indicating treated
  #' @param control_level string - level of treatment variable indicating control
  #' @param f.ps formula - propensity score model formula (required for reweight method)
  #' @param f.out character - outcome model formula (RHS only)
  #' @param data data.frame - dataset (must contain wgt if bootstrap_method = "resample")
  #' @param gbm_params list - GBM parameters for ps() function
  #' @param bootstrap_method character - "resample" (faster) or "reweight" (more conservative)
  #' @param stratified logical - use stratified bootstrap to maintain treatment/control proportions
  #' @param wgt character - name of propensity score weight column in data (for resample method)
  #' @param n_boot integer - number of bootstrap replications
  #' @param seed integer - random seed
  #' @param verbose logical - print progress and results
  #' @param parallel logical - use parallel processing
  #' @param n_cores integer - number of cores (NULL = auto-detect)
  #' @param plot_diagnostics logical - generate diagnostic plots
  #' @param ci_type character - confidence interval types: "all", "norm", "basic", "perc"
  #' @param sim character - bootstrap simulation: "ordinary" or "balanced"
  #' @param save_to character - directory path to save diagnostic plots (NULL = don't save)
  #' @param alpha numeric - propensity score truncation level (0 = no truncation)
  #' @param use_weights logical - if FALSE, use outcome regression without PS weights (deprecated, use estimator instead)
  #' @param estimator character - "simple" (difference in means), "or" (outcome regression), or "drs" (doubly robust)
  #'
  #' @return list with ATT estimates, SEs, CIs, p-values, balance diagnostics, and bootstrap samples

  # Match and validate arguments
  bootstrap_method = match.arg(bootstrap_method)
  ci_type = match.arg(ci_type)
  sim = match.arg(sim)
  estimator = match.arg(estimator)

  # Validate inputs based on estimator and bootstrap method
  if (estimator == "drs") {
    if (bootstrap_method == "resample") {
      if (!wgt %in% names(data)) {
        stop(sprintf(
          "Weight variable '%s' not found in data. Available columns: %s",
          wgt,
          paste(names(data), collapse = ", ")
        ))
      }
    }

    if (bootstrap_method == "reweight" && is.null(f.ps)) {
      stop(
        "f.ps formula required when bootstrap_method = 'reweight' and estimator = 'drs'"
      )
    }
  }

  if (verbose) {
    estimator_name = switch(
      estimator,
      simple = "Simple",
      or = "OR",
      drs = "DRS"
    )
    cat(sprintf("=== Bootstrap %s ATT Estimation ===\n", estimator_name))
    cat("Estimator:", estimator, "\n")
    cat("Outcome:", outcome, "\n")
    cat("Treatment:", treatment, "\n")
    cat("Treated level:", treated_level, "\n")
    cat("Control level:", control_level, "\n")
    if (estimator == "drs") {
      cat("Bootstrap method:", bootstrap_method, "\n")
      if (bootstrap_method == "resample") {
        cat("PS weight variable:", wgt, "\n")
      }
    }
    cat("Bootstrap replications:", n_boot, "\n")
    cat("Bootstrap simulation type:", sim, "\n")
    cat("Stratified bootstrap:", stratified, "\n")
    cat("Parallel processing:", parallel, "\n\n")
  }

  # Subset data to treated and control groups
  analysis_data = data %>%
    filter(!!sym(treatment) %in% c(treated_level, control_level)) %>%
    mutate(across(where(is.factor), droplevels))

  n = nrow(analysis_data)

  # Set up parallel processing
  parallel_type = "no"
  if (parallel) {
    if (!requireNamespace("parallel", quietly = TRUE)) {
      warning("parallel package not available, using sequential processing")
    } else {
      if (is.null(n_cores)) {
        n_cores = max(1, parallel::detectCores() - 1)
      }
      n_cores = min(n_cores, parallel::detectCores())
      parallel_type = "multicore"
      if (verbose) {
        cat("Using", n_cores, "cores\n\n")
      }
    }
  }

  # Prepare strata for stratified bootstrap
  strata_vec = NULL
  if (stratified) {
    strata_vec = as.numeric(analysis_data[[treatment]] == treated_level)
    if (verbose) {
      n_treated = sum(strata_vec)
      n_control = length(strata_vec) - n_treated
      cat(sprintf(
        "Stratification: %d treated, %d control (%.1f%% treated)\n\n",
        n_treated,
        n_control,
        100 * n_treated / n
      ))
    }
  }

  # Bootstrap
  if (verbose) {
    cat("Running bootstrap using boot package...\n\n")
  }

  # Set seed for reproducibility
  set.seed(seed)

  # Run boot() with boot_iter function
  boot_result = boot::boot(
    data = analysis_data,
    statistic = boot_iter,
    R = n_boot,
    sim = sim,
    strata = strata_vec,
    parallel = parallel_type,
    ncpus = if (parallel) n_cores else 1,
    # Additional arguments passed to boot_iter
    outcome = outcome,
    treatment = treatment,
    treated_level = treated_level,
    control_level = control_level,
    f.ps = f.ps,
    f.out = f.out,
    gbm_params = gbm_params,
    bootstrap_method = bootstrap_method,
    wgt = wgt,
    seed_offset = seed,
    alpha = alpha,
    use_weights = use_weights,
    estimator = estimator
  )

  # Extract point estimates and bootstrap samples
  # boot_result has 4 columns: drs, avg_asd, max_asd, ess
  drs_point_att = boot_result$t0[1]

  if (verbose) {
    estimator_name = switch(
      estimator,
      simple = "Simple",
      or = "OR",
      drs = "DRS"
    )
    cat("\nPoint Estimate:\n")
    cat(sprintf("  %s ATT: %.4f\n", estimator_name, drs_point_att))
    cat("\n")
  }

  # Extract bootstrap samples
  drs_boots = boot_result$t[, 1]

  # Extract balance metrics if available
  if (ncol(boot_result$t) >= 4) {
    avg_asd_boots = boot_result$t[, 2]
    max_asd_boots = boot_result$t[, 3]
    ess_boots = boot_result$t[, 4]
  } else {
    avg_asd_boots = NULL
    max_asd_boots = NULL
    ess_boots = NULL
  }

  # Remove failed iterations
  valid_idx = !is.na(drs_boots)
  n_failed = sum(!valid_idx)

  if (n_failed > 0 && verbose) {
    cat(sprintf("\nWarning: %d bootstrap iterations failed\n", n_failed))
  }

  # Filter valid bootstrap samples
  drs_boots = drs_boots[valid_idx]

  # Also filter balance metrics if they exist
  if (!is.null(avg_asd_boots)) {
    avg_asd_boots = avg_asd_boots[valid_idx]
    max_asd_boots = max_asd_boots[valid_idx]
    ess_boots = ess_boots[valid_idx]
  }

  # Compute standard error
  drs_se = sd(drs_boots)

  # Perform Shapiro-Wilk normality test for DRS
  drs_shapiro = tryCatch(
    shapiro.test(drs_boots),
    error = function(e) {
      warning(sprintf("Shapiro-Wilk test failed for DRS: %s", e$message))
      list(
        statistic = NA,
        p.value = NA,
        method = "Shapiro-Wilk normality test"
      )
    }
  )

  # Compute confidence intervals using boot.ci()
  # We need to filter the boot object to remove NA rows
  boot_result_filtered = boot_result
  boot_result_filtered$t = boot_result$t[valid_idx, ]
  boot_result_filtered$R = sum(valid_idx)

  # Determine which CI types to compute
  if (ci_type == "all") {
    ci_types_to_compute = c("norm", "basic", "perc")
  } else {
    ci_types_to_compute = ci_type
  }

  # Compute CIs for DRS (index = 1)
  drs_ci = tryCatch(
    boot::boot.ci(
      boot_result_filtered,
      type = ci_types_to_compute,
      index = 1
    ),
    error = function(e) {
      warning("boot.ci() failed for DRS, falling back to manual computation")
      NULL
    }
  )

  # Extract CIs or fall back to manual computation
  if (!is.null(drs_ci)) {
    drs_ci_normal = if ("normal" %in% names(drs_ci)) {
      drs_ci$normal[2:3]
    } else {
      c(drs_point_att - 1.96 * drs_se, drs_point_att + 1.96 * drs_se)
    }
    drs_ci_basic = if ("basic" %in% names(drs_ci)) {
      drs_ci$basic[4:5]
    } else {
      NULL
    }
    drs_ci_percentile = if ("percent" %in% names(drs_ci)) {
      drs_ci$percent[4:5]
    } else {
      quantile(drs_boots, probs = c(0.025, 0.975))
    }
  } else {
    drs_ci_normal = c(
      drs_point_att - 1.96 * drs_se,
      drs_point_att + 1.96 * drs_se
    )
    drs_ci_basic = NULL
    drs_ci_percentile = quantile(drs_boots, probs = c(0.025, 0.975))
  }

  # Compute p-value using normal approximation
  drs_z = drs_point_att / drs_se
  drs_pval = 2 * pnorm(-abs(drs_z))

  # Print results
  if (verbose) {
    estimator_name = switch(
      estimator,
      simple = "Simple",
      or = "OR",
      drs = "DRS"
    )
    cat(sprintf("\n=== %s Results ===\n", estimator_name))
    cat(sprintf("ATT estimate:     %.4f\n", drs_point_att))
    cat(sprintf("Standard error:   %.4f\n", drs_se))
    cat(sprintf("p-value:          %.4f\n", drs_pval))
    cat("\nConfidence Intervals (95%):\n")
    if (!is.null(drs_ci_normal)) {
      cat(sprintf(
        "  Normal:      [%.4f, %.4f]\n",
        drs_ci_normal[1],
        drs_ci_normal[2]
      ))
    }
    if (!is.null(drs_ci_basic)) {
      cat(sprintf(
        "  Basic:       [%.4f, %.4f]\n",
        drs_ci_basic[1],
        drs_ci_basic[2]
      ))
    }
    if (!is.null(drs_ci_percentile)) {
      cat(sprintf(
        "  Percentile:  [%.4f, %.4f]\n",
        drs_ci_percentile[1],
        drs_ci_percentile[2]
      ))
    }

    # Display Shapiro-Wilk test results for DRS
    cat("\nNormality Test (Shapiro-Wilk):\n")
    if (!is.na(drs_shapiro$p.value)) {
      cat(sprintf(
        "  W statistic: %.4f\n",
        drs_shapiro$statistic
      ))
      cat(sprintf(
        "  p-value:     %.4f %s\n",
        drs_shapiro$p.value,
        ifelse(
          drs_shapiro$p.value < 0.05,
          "(reject normality at 0.05 level)",
          "(fail to reject normality)"
        )
      ))
    } else {
      cat("  Test failed (likely due to sample size or other issues)\n")
    }

    # Display balance diagnostics if available
    if (!is.null(avg_asd_boots)) {
      balance_summary = list(
        avg_asd = list(
          mean = mean(avg_asd_boots, na.rm = TRUE),
          sd = sd(avg_asd_boots, na.rm = TRUE),
          min = min(avg_asd_boots, na.rm = TRUE),
          max = max(avg_asd_boots, na.rm = TRUE)
        ),
        max_asd = list(
          mean = mean(max_asd_boots, na.rm = TRUE),
          sd = sd(max_asd_boots, na.rm = TRUE),
          min = min(max_asd_boots, na.rm = TRUE),
          max = max(max_asd_boots, na.rm = TRUE)
        ),
        ess = list(
          mean = mean(ess_boots, na.rm = TRUE),
          sd = sd(ess_boots, na.rm = TRUE),
          min = min(ess_boots, na.rm = TRUE),
          max = max(ess_boots, na.rm = TRUE)
        )
      )

      cat("\n=== Bootstrap Balance Diagnostics ===\n")
      cat(sprintf(
        "Average ASD across bootstrap samples: %.4f (SD: %.4f)\n",
        balance_summary$avg_asd$mean,
        balance_summary$avg_asd$sd
      ))
      cat(sprintf(
        "  Range: [%.4f, %.4f]\n",
        balance_summary$avg_asd$min,
        balance_summary$avg_asd$max
      ))
      cat(sprintf(
        "Maximum ASD across bootstrap samples: %.4f (SD: %.4f)\n",
        balance_summary$max_asd$mean,
        balance_summary$max_asd$sd
      ))
      cat(sprintf(
        "  Range: [%.4f, %.4f]\n",
        balance_summary$max_asd$min,
        balance_summary$max_asd$max
      ))
      cat(sprintf(
        "Effective Sample Size: %.1f (SD: %.1f)\n",
        balance_summary$ess$mean,
        balance_summary$ess$sd
      ))
      cat(sprintf(
        "  Range: [%.1f, %.1f]\n",
        balance_summary$ess$min,
        balance_summary$ess$max
      ))
    }
  }

  # Generate diagnostic plots with color coding
  if (plot_diagnostics) {
    # Set colors based on estimator
    hist_col = switch(
      estimator,
      simple = "lightgrey",
      or = "lightblue",
      drs = "lightgreen"
    )

    plot_title_prefix = switch(
      estimator,
      simple = "Simple",
      or = "OR",
      drs = "DRS"
    )

    par(mfrow = c(1, 2), mar = c(4, 4, 2, 1))

    # Histogram with color based on estimator
    hist(
      drs_boots,
      main = paste(plot_title_prefix, "Bootstrap Distribution"),
      xlab = "ATT Estimate",
      col = hist_col,
      border = "white",
      breaks = 30
    )
    abline(v = drs_point_att, col = "red", lwd = 2, lty = 2)

    # Q-Q plot
    qqnorm(drs_boots, main = paste(plot_title_prefix, "Q-Q Plot"))
    qqline(drs_boots, col = "red", lwd = 2)

    par(mfrow = c(1, 1))

    # Save to file if save_to is specified
    if (!is.null(save_to)) {
      # Create directory if it doesn't exist
      if (!dir.exists(save_to)) {
        dir.create(save_to, recursive = TRUE)
      }

      # Construct filename
      filename = sprintf(
        "DR_bootstrap_%s_vs_%s_%s%s_n%d.png",
        treated_level,
        control_level,
        sim,
        if (alpha > 0) paste0("_alpha", alpha) else "",
        n_boot
      )
      filepath = file.path(save_to, filename)

      # Open high-resolution PNG device for academic paper quality
      png(filepath, width = 3000, height = 3000, res = 300)

      par(mfrow = c(1, 2), mar = c(4, 4, 2, 1))

      # Histogram with color based on estimator
      hist(
        drs_boots,
        main = paste(plot_title_prefix, "Bootstrap Distribution"),
        xlab = "ATT Estimate",
        col = hist_col,
        border = "white",
        breaks = 30
      )
      abline(v = drs_point_att, col = "red", lwd = 2, lty = 2)

      # Q-Q plot
      qqnorm(drs_boots, main = paste(plot_title_prefix, "Q-Q Plot"))
      qqline(drs_boots, col = "red", lwd = 2)

      par(mfrow = c(1, 1))

      dev.off()
      if (verbose) {
        cat(sprintf("\nDiagnostic plot saved to: %s\n", filepath))
      }
    }
  }

  # Create balance_diagnostics if metrics are available
  balance_diagnostics = NULL
  bootstrap_balance = NULL

  if (!is.null(avg_asd_boots) && length(avg_asd_boots) > 0) {
    balance_diagnostics = list(
      avg_asd = list(
        mean = mean(avg_asd_boots, na.rm = TRUE),
        sd = sd(avg_asd_boots, na.rm = TRUE),
        min = min(avg_asd_boots, na.rm = TRUE),
        max = max(avg_asd_boots, na.rm = TRUE)
      ),
      max_asd = list(
        mean = mean(max_asd_boots, na.rm = TRUE),
        sd = sd(max_asd_boots, na.rm = TRUE),
        min = min(max_asd_boots, na.rm = TRUE),
        max = max(max_asd_boots, na.rm = TRUE)
      ),
      ess = list(
        mean = mean(ess_boots, na.rm = TRUE),
        sd = sd(ess_boots, na.rm = TRUE),
        min = min(ess_boots, na.rm = TRUE),
        max = max(ess_boots, na.rm = TRUE)
      )
    )

    bootstrap_balance = data.frame(
      iteration = seq_len(length(avg_asd_boots)),
      avg_asd = avg_asd_boots,
      max_asd = max_asd_boots,
      ess = ess_boots
    )
  }

  # Return results
  results = list(
    drs = list(
      att = drs_point_att,
      se = drs_se,
      ci_normal = drs_ci_normal,
      ci_basic = drs_ci_basic,
      ci_percentile = drs_ci_percentile,
      pval = drs_pval,
      bootstrap_samples = drs_boots,
      boot_ci_object = drs_ci,
      shapiro_test = drs_shapiro
    ),
    n_boot = length(drs_boots),
    n_failed = n_failed,
    balance_diagnostics = balance_diagnostics,
    bootstrap_balance = bootstrap_balance,
    boot_object = boot_result
  )

  return(results)
}


# NEW HELPER FUNCTIONS FOR DEFINITIVE ANALYSIS -----------------------------------------####

build_result_path = function(
  comparison,
  model,
  base = "results/outcome/DRS"
) {
  #' Build consistent directory path for DRS analysis results
  #'
  #' @param comparison character - "inc_dec" or "inc_mix"
  #' @param model integer - model number (1-4)
  #' @param base character - base directory path (defaults to DRS folder)
  #'
  #' @return character - full path to results directory (created if doesn't exist)

  path = file.path(
    base,
    comparison,
    sprintf("model%d", model)
  )
  dir.create(path, showWarnings = FALSE, recursive = TRUE)
  return(path)
}

build_master_summary = function(results_list, estimator_name, comparison_name) {
  #' Build comprehensive summary table from multiple DRS analysis results
  #'
  #' @param results_list named list - analysis results, names should indicate model
  #' @param estimator_name character - "DRS" (kept for consistency)
  #' @param comparison_name character - "inc_dec" or "inc_mix"
  #'
  #' @return data.frame - comprehensive summary with all analyses

  summary_rows = list()

  for (result_name in names(results_list)) {
    result = results_list[[result_name]]

    # Extract model from name (e.g., "model1", "model2", etc.)
    model_num = as.integer(sub("model", "", result_name))

    # Extract DRS results
    est_lower = tolower(estimator_name)
    if (!est_lower %in% names(result)) {
      warning(sprintf(
        "Estimator '%s' not found in result '%s'",
        est_lower,
        result_name
      ))
      next
    }

    est = result[[est_lower]]

    # Build row from DRS results (no nested structure handling)
    row = data.frame(
      Estimator = estimator_name,
      Comparison = comparison_name,
      Model = model_num,
      ATT = est$att,
      SE = est$se,
      CI_Normal_Lower = est$ci_normal[1],
      CI_Normal_Upper = est$ci_normal[2],
      CI_Basic_Lower = ifelse(
        !is.null(est$ci_basic),
        est$ci_basic[1],
        NA
      ),
      CI_Basic_Upper = ifelse(
        !is.null(est$ci_basic),
        est$ci_basic[2],
        NA
      ),
      CI_Perc_Lower = est$ci_percentile[1],
      CI_Perc_Upper = est$ci_percentile[2],
      p_value = est$pval,
      Shapiro_W = ifelse(
        !is.null(est$shapiro_test),
        est$shapiro_test$statistic,
        NA
      ),
      Shapiro_p = ifelse(
        !is.null(est$shapiro_test),
        est$shapiro_test$p.value,
        NA
      ),
      n_boot = result$n_boot,
      n_failed = result$n_failed,
      stringsAsFactors = FALSE
    )

    # Add balance diagnostics if available
    if (!is.null(result$balance_diagnostics)) {
      row$Avg_ASD_Mean = result$balance_diagnostics$avg_asd$mean
      row$Avg_ASD_SD = result$balance_diagnostics$avg_asd$sd
      row$Max_ASD_Mean = result$balance_diagnostics$max_asd$mean
      row$Max_ASD_SD = result$balance_diagnostics$max_asd$sd
      row$ESS_Mean = result$balance_diagnostics$ess$mean
      row$ESS_SD = result$balance_diagnostics$ess$sd
    }

    summary_rows[[length(summary_rows) + 1]] = row
  }

  if (length(summary_rows) == 0) {
    stop("No valid results found in results_list")
  }

  summary_df = do.call(rbind, summary_rows)
  rownames(summary_df) = NULL

  return(summary_df)
}


# OUTCOME MODEL COMPARISON PLOT --------------------------------------------------------####
plot_outcome_model_boxplot = function(
  results_list,
  estimator_name,
  comparison_label,
  save_path,
  width = 2400,
  height = 1800,
  res = 300
) {
  #' Create boxplot comparison of Simple/OR/DRS ATT estimates across 6 outcome models
  #'
  #' @param results_list named list - results for model0, model1, model2, model3, model4, model5
  #' @param estimator_name character - "DRS" (kept for consistency)
  #' @param comparison_label character - label for plot title
  #' @param save_path character - path to save PNG file
  #' @param width numeric - plot width in pixels
  #' @param height numeric - plot height in pixels
  #' @param res numeric - resolution in DPI
  #'
  #' @return NULL (displays and saves plot)

  est = tolower(estimator_name)

  # Helper function to extract bootstrap samples (no nested structure)
  extract_boots = function(x) {
    x[[est]]$bootstrap_samples
  }

  # Helper function to extract point estimates (no nested structure)
  extract_att = function(x) {
    x[[est]]$att
  }

  # Extract data for all 6 models (0-5)
  model_names = paste0("model", 0:5)
  bootstrap_data = lapply(model_names, function(m) {
    extract_boots(results_list[[m]])
  })
  point_estimates = sapply(model_names, function(m) {
    extract_att(results_list[[m]])
  })

  # Model labels - removed "Model *:" prefix
  model_labels = c(
    "Avg Diff",
    "Full",
    "Intercept",
    "Baseline",
    "Base+Depression",
    "Full"
  )

  # Color vectors for 6 models
  box_colors = c("lightgrey", "lightblue", rep("lightgreen", 4))
  border_colors = c("darkgrey", "darkblue", rep("darkgreen", 4))

  # Create plot function
  create_plot = function() {
    par(mar = c(6, 5, 4, 2), font.lab = 2) # Adjusted bottom margin

    # First create the boxplot without data (to get axis limits)
    temp_plot = boxplot(
      bootstrap_data,
      names = model_labels,
      main = sprintf("%s - ATT Comparison", comparison_label),
      xlab = "",
      ylab = "ATT Estimate",
      plot = FALSE
    )

    # Now create empty plot with axis
    plot(
      1,
      type = "n",
      xlim = c(0.5, 6.5),
      ylim = range(temp_plot$stats),
      xaxt = "n",
      main = sprintf("%s - ATT Comparison", comparison_label),
      xlab = "",
      ylab = "ATT Estimate",
      cex.lab = 1.2
    ) # Bigger axis label

    # Add grid in the background
    grid(nx = NA, ny = NULL, col = rgb(0, 0, 0, alpha = 0.1), lty = 1)

    # Add x-axis labels with bigger size
    axis(
      1,
      at = 1:6,
      labels = model_labels,
      cex.axis = 1.1,
      padj = 0.5,
      mgp = c(3, 1, 0)
    )

    # Add x-axis title
    mtext("Outcome Model", side = 1, line = 4, font = 2, cex = 1.2)

    # Now add the boxplot on top
    boxplot(
      bootstrap_data,
      names = model_labels,
      col = box_colors,
      border = border_colors,
      outline = TRUE,
      las = 1,
      add = TRUE,
      xaxt = "n",
      yaxt = "n",
      cex.axis = 1.1 # Bigger y-axis labels
    )

    # Add horizontal line at zero
    abline(h = 0, col = "black", lty = 2, lwd = 1.5)

    # Overlay point estimates
    points(1:6, point_estimates, col = "red", pch = 19, cex = 1.5)

    # Add legend with tight-fitting box
    legend(
      "topright",
      legend = c("Simple", "OR", "DRS", "Point Estimate"),
      col = c("lightgrey", "lightblue", "lightgreen", "red"),
      pch = c(15, 15, 15, 19),
      lty = c(NA, NA, NA, NA),
      lwd = c(NA, NA, NA, NA),
      cex = 1.0,
      pt.cex = 1.2, # Bigger colored squares
      bg = "white", # White background
      box.lty = 1, # Solid box border
      box.lwd = 1, # Box line width
      x.intersp = 0.5, # Horizontal spacing
      y.intersp = 0.8, # Vertical spacing
      inset = 0 # Small inset from corner
    )
  }

  # Display plot
  create_plot()

  # Save to file
  png(save_path, width = width, height = height, res = res)
  create_plot()
  dev.off()

  invisible(NULL)
}

# CF VS NO-CF HISTOGRAM COMPARISON -----------------------------------------------------####

# HELPER X RESULTS ---------------------------------------------------------------------####
extract_results = function(
  result,
  model_num,
  method,
  estimator
) {
  #' Extract results from DRS result structure
  #'
  #' @param result list - analysis result object
  #' @param model_num integer - model number
  #' @param method character - method name
  #' @param estimator character - estimator name ("DRS")
  #'
  #' @return data.frame - extracted results

  est_data = result[[tolower(estimator)]]

  data.frame(
    Model = model_num,
    Method = method,
    Estimator = estimator,
    ATT = est_data$att,
    SE = est_data$se,
    CI_Normal_Lower = est_data$ci_normal[1],
    CI_Normal_Upper = est_data$ci_normal[2],
    CI_Basic_Lower = ifelse(
      !is.null(est_data$ci_basic),
      est_data$ci_basic[1],
      NA
    ),
    CI_Basic_Upper = ifelse(
      !is.null(est_data$ci_basic),
      est_data$ci_basic[2],
      NA
    ),
    CI_Perc_Lower = est_data$ci_percentile[1],
    CI_Perc_Upper = est_data$ci_percentile[2],
    p_value = est_data$pval,
    Shapiro_W = ifelse(
      !is.null(est_data$shapiro_test),
      est_data$shapiro_test$statistic,
      NA
    ),
    Shapiro_p = ifelse(
      !is.null(est_data$shapiro_test),
      est_data$shapiro_test$p.value,
      NA
    ),
    n_boot = result$n_boot,
    n_failed = result$n_failed
  )
}

# Add balance diagnostics
extract_balance = function(result, model_num, method) {
  if (!is.null(result$balance_diagnostics)) {
    data.frame(
      Model = model_num,
      Method = method,
      Avg_ASD_Mean = result$balance_diagnostics$avg_asd$mean,
      Avg_ASD_SD = result$balance_diagnostics$avg_asd$sd,
      Max_ASD_Mean = result$balance_diagnostics$max_asd$mean,
      Max_ASD_SD = result$balance_diagnostics$max_asd$sd,
      ESS_Mean = result$balance_diagnostics$ess$mean,
      ESS_SD = result$balance_diagnostics$ess$sd
    )
  } else {
    NULL
  }
}

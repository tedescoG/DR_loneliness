# HELPER FUNCTIONS FOR ANALYSIS

# Libraries
library(tidyverse)
library(gbm)
library(twang)
library(survey)
library(marginaleffects)
library(boot)

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
  #' @param seed integer - random seed for reproducibility
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
  set.seed(seed)
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
          method,
          seed = seed
        )
      },
      mc.cores = n.cores,
      mc.set.seed = TRUE
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

# AIPW Estimator -----------------------------------------------------------------------####

aipw_att = function(outcome, treatment, f.out, wgt, data, verbose = T) {
  #' Compute AIPW estimator for Average Treatment Effect on Treated (ATT)
  #'
  #' @param outcome character - name of outcome variable
  #' @param treatment character - name of treatment variable (must be 0/1)
  #' @param f.out character - outcome model formula (RHS only)
  #' @param wgt character - name of weight variable (propensity scores)
  #' @param data data.frame - dataset containing all variables
  #' @param verbose logical - whether to print results
  #'
  #' @return list with ATT estimate and counterfactual means
  # Extract variables
  Y = data %>% pull(all_of(outcome))
  Z = data %>% pull(all_of(treatment))
  weights = data %>% pull(all_of(wgt))
  N1 = sum(Z)

  # Fit outcome model mu0 on control units only
  control_data = data %>% filter(Z == 0)
  mu0 = glm(
    formula = as.formula(paste(outcome, "~", f.out)),
    data = control_data,
    family = "quasibinomial",
    control = glm.control(maxit = 50) # avoid convergence problems
  )

  # Predict mu0 for all observations
  mu_hat_0 = predict(mu0, newdata = data, type = "response")

  # Moodie (2018) AIPW ATT estimator
  att_est = (1 / N1) *
    sum(
      (Z - (1 - Z) * weights) * (Y - mu_hat_0)
    )

  # output
  if (verbose) {
    cat(sprintf("ATT estimate (AIPW): %.4f\n", att_est))
    cat(sprintf("Mean Y(0) [counterfactual]: %.4f\n", mu_hat_0))
  }

  # Return results
  return(list(
    att = att_est,
    mu_hat_0 = mu_hat_0,
    mu_hat_1 = mean(Y[Z == 1])
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
  if (nzchar(f.out)) {
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


# BOOTSTRAP DR ESTIMATORS ------------------------------------------------------------------####

# Helper function for single bootstrap iteration with weight resampling (no PS re-estimation)
boot_iter_resample = function(
  indices,
  data,
  outcome,
  treatment,
  treated_level,
  control_level,
  f.out,
  wgt
) {
  #' Bootstrap iteration using resampled data with existing weights
  #'
  #' @param indices integer - bootstrap sample indices
  #' @param data data.frame - full dataset
  #' @param outcome character - name of outcome variable
  #' @param treatment character - name of treatment variable
  #' @param treated_level string - level indicating treated group
  #' @param control_level string - level indicating control group
  #' @param f.out character - outcome model formula (RHS only)
  #' @param wgt character - name of weight variable
  #'
  #' @return numeric vector with AIPW, DRS estimates and balance metrics
  # Resample data
  boot_data = data[indices, ]

  # Create binary treatment indicator
  boot_data[[treatment]] = as.numeric(
    boot_data[[treatment]] == treated_level
  )

  # Compute AIPW estimate with resampled weights (these work fine)
  aipw_result = aipw_att(
    outcome = outcome,
    treatment = treatment,
    f.out = f.out,
    wgt = wgt,
    data = boot_data,
    verbose = FALSE
  )

  # Compute DRS estimate with resampled weights (these work fine)
  drs_result = drs_att(
    outcome = outcome,
    treatment = treatment,
    f.out = f.out,
    wgt = wgt,
    data = boot_data,
    verbose = FALSE
  )

  # Compute balance diagnostics using our std.diff function
  avg_asd = NA
  max_asd = NA
  ess = NA

  tryCatch(
    {
      # Extract all covariates (exclude treatment, outcome, and weight)
      covariate_cols = setdiff(names(boot_data), c(treatment, outcome, wgt))

      if (length(covariate_cols) > 0) {
        covariates = boot_data[, covariate_cols, drop = FALSE]

        # Ensure treatment is numeric 0/1
        treat_vec = as.numeric(boot_data[[treatment]])

        # Ensure weights are numeric
        weights_vec = as.numeric(boot_data[[wgt]])

        # Check for valid weights
        if (all(is.finite(weights_vec)) && all(weights_vec > 0)) {
          # Compute standardized differences using our function
          std_diffs = unlist(sapply(covariates, std.diff, z = treat_vec, w = weights_vec))
          std_diffs = std_diffs[!is.na(std_diffs)]

          if (length(std_diffs) > 0) {
            avg_asd = mean(std_diffs)
            max_asd = max(std_diffs)
          }

          # Compute effective sample size
          # ESS = (sum of weights)^2 / sum of weights^2
          ess = sum(weights_vec)^2 / sum(weights_vec^2)
        }
      }
    },
    error = function(e) {
      # If calculation fails, leave as NA
    }
  )

  return(c(
    aipw = aipw_result$att,
    drs = drs_result$att,
    avg_asd = avg_asd,
    max_asd = max_asd,
    ess = ess
  ))
}

# Helper function for single bootstrap iteration with PS re-estimation
boot_iter_reweight = function(
  indices,
  data,
  outcome,
  treatment,
  treated_level,
  control_level,
  f.ps,
  f.out,
  ps_params,
  seed_offset
) {
  #' Bootstrap iteration with propensity score re-estimation
  #'
  #' @param indices integer - bootstrap sample indices
  #' @param data data.frame - full dataset
  #' @param outcome character - name of outcome variable
  #' @param treatment character - name of treatment variable
  #' @param treated_level string - level indicating treated group
  #' @param control_level string - level indicating control group
  #' @param f.ps formula - propensity score model formula
  #' @param f.out character - outcome model formula (RHS only)
  #' @param ps_params list - GBM parameters for ps() function
  #' @param seed_offset integer - offset for random seed
  #'
  #' @return numeric vector with AIPW, DRS estimates and balance metrics
  # Resample data
  boot_data = data[indices, ]

  # Create binary treatment indicator
  boot_data[[treatment]] = as.numeric(
    boot_data[[treatment]] == treated_level
  )

  # Fit propensity score model (this works fine)
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
      ps_params
    )
  )

  # Extract weights
  boot_data$ps_wgt = get.weights(ps_fit, stop.method = "es.mean")

  # Compute AIPW estimate (this works fine)
  aipw_result = aipw_att(
    outcome = outcome,
    treatment = treatment,
    f.out = f.out,
    wgt = "ps_wgt",
    data = boot_data,
    verbose = FALSE
  )

  # Compute DRS estimate (this works fine)
  drs_result = drs_att(
    outcome = outcome,
    treatment = treatment,
    f.out = f.out,
    wgt = "ps_wgt",
    data = boot_data,
    verbose = FALSE
  )

  # Compute balance diagnostics using our std.diff function
  avg_asd = NA
  max_asd = NA
  ess = NA

  tryCatch(
    {
      # Extract all covariates (exclude treatment, outcome, and weight)
      covariate_cols = setdiff(names(boot_data), c(treatment, outcome, "ps_wgt"))

      if (length(covariate_cols) > 0) {
        covariates = boot_data[, covariate_cols, drop = FALSE]

        # Ensure treatment is numeric 0/1
        treat_vec = as.numeric(boot_data[[treatment]])

        # Ensure weights are numeric
        weights_vec = as.numeric(boot_data$ps_wgt)

        # Check for valid weights
        if (all(is.finite(weights_vec)) && all(weights_vec > 0)) {
          # Compute standardized differences using our function
          std_diffs = unlist(sapply(covariates, std.diff, z = treat_vec, w = weights_vec))
          std_diffs = std_diffs[!is.na(std_diffs)]

          if (length(std_diffs) > 0) {
            avg_asd = mean(std_diffs)
            max_asd = max(std_diffs)
          }

          # Compute effective sample size
          # ESS = (sum of weights)^2 / sum of weights^2
          ess = sum(weights_vec)^2 / sum(weights_vec^2)
        }
      }
    },
    error = function(e) {
      # If calculation fails, leave as NA
    }
  )

  return(c(
    aipw = aipw_result$att,
    drs = drs_result$att,
    avg_asd = avg_asd,
    max_asd = max_asd,
    ess = ess
  ))
}

# Main bootstrap DR function
DR_att = function(
  outcome,
  treatment,
  treated_level,
  control_level,
  f.ps,
  f.out,
  data,
  ps_params = list(
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
  plot_diagnostics = TRUE
) {
  #' Bootstrap Doubly Robust ATT Estimation
  #'
  #' @param outcome character - name of outcome variable
  #' @param treatment character - name of treatment variable
  #' @param treated_level string - level of treatment variable indicating treated
  #' @param control_level string - level of treatment variable indicating control
  #' @param f.ps formula - propensity score model formula (only used if bootstrap_method = "reweight")
  #' @param f.out character - outcome model formula (RHS only)
  #' @param data data.frame - dataset (must contain wgt if bootstrap_method = "resample")
  #' @param ps_params list - GBM parameters for ps() function (only used if bootstrap_method = "reweight")
  #' @param bootstrap_method character - bootstrap method: "resample" (default, faster) or "reweight" (more conservative)
  #' @param stratified logical - use stratified bootstrap to maintain treatment/control proportions (default TRUE)
  #' @param wgt character - name of propensity score weight column in data
  #' @param n_boot integer - number of bootstrap replications
  #' @param seed integer - random seed
  #' @param verbose logical - print progress
  #' @param parallel logical - use parallel processing
  #' @param n_cores integer - number of cores (NULL = auto-detect)
  #' @param plot_diagnostics logical - generate diagnostic plots
  #'
  #' @return list with estimates, SEs, CIs, p-values, and bootstrap samples

  # Match and validate bootstrap method
  bootstrap_method = match.arg(bootstrap_method)

  # Validate inputs based on bootstrap method
  if (bootstrap_method == "resample") {
    if (!wgt %in% names(data)) {
      stop(sprintf(
        "Weight variable '%s' not found in data. Available columns: %s",
        wgt,
        paste(names(data), collapse = ", ")
      ))
    }
  }

  if (verbose) {
    cat("=== Bootstrap Doubly Robust ATT Estimation ===\n")
    cat("Outcome:", outcome, "\n")
    cat("Treatment:", treatment, "\n")
    cat("Treated level:", treated_level, "\n")
    cat("Control level:", control_level, "\n")
    cat("Bootstrap method:", bootstrap_method, "\n")
    if (bootstrap_method == "resample") {
      cat("PS weight variable:", wgt, "\n")
    }
    cat("Bootstrap replications:", n_boot, "\n")
    cat("Stratified bootstrap:", stratified, "\n")
    cat("Parallel processing:", parallel, "\n\n")
  }

  # Subset data to treated and control groups
  analysis_data = data %>%
    filter(!!sym(treatment) %in% c(treated_level, control_level)) %>%
    mutate(across(where(is.factor), droplevels))

  n = nrow(analysis_data)

  # Set up parallel processing
  if (parallel) {
    if (!requireNamespace("parallel", quietly = TRUE)) {
      warning("parallel package not available, using sequential processing")
      parallel = FALSE
    } else {
      if (is.null(n_cores)) {
        n_cores = max(1, parallel::detectCores() - 1)
      }
      n_cores = min(n_cores, parallel::detectCores())
      if (verbose) {
        cat("Using", n_cores, "cores\n\n")
      }
    }
  }

  # Point estimates on full sample
  if (verbose) {
    cat("Computing point estimates on full sample...\n")
  }

  # Create binary treatment for full sample
  full_data = analysis_data
  full_data[[treatment]] = as.numeric(
    full_data[[treatment]] == treated_level
  )

  # Point estimates
  aipw_point = aipw_att(
    outcome = outcome,
    treatment = treatment,
    f.out = f.out,
    wgt = wgt,
    data = full_data,
    verbose = FALSE
  )

  drs_point = drs_att(
    outcome = outcome,
    treatment = treatment,
    f.out = f.out,
    wgt = wgt,
    data = full_data,
    verbose = FALSE
  )

  if (verbose) {
    cat("\nPoint Estimates:\n")
    cat(sprintf("  AIPW ATT: %.4f\n", aipw_point$att))
    cat(sprintf("  DRS ATT:  %.4f\n\n", drs_point$att))
  }

  # Bootstrap
  if (verbose) {
    cat("Running bootstrap...\n")
  }

  # Set up RNG for reproducibility
  RNGkind("L'Ecuyer-CMRG")
  set.seed(seed)

  # Pre-compute treatment group indices for stratified bootstrap
  if (stratified) {
    treated_indices = which(analysis_data[[treatment]] == treated_level)
    control_indices = which(analysis_data[[treatment]] == control_level)
    n_treated = length(treated_indices)
    n_control = length(control_indices)

    if (verbose) {
      cat(sprintf(
        "Stratification: %d treated, %d control (%.1f%% treated)\n\n",
        n_treated,
        n_control,
        100 * n_treated / n
      ))
    }
  }

  # Create progress bar for sequential processing
  if (!parallel && verbose) {
    pb = txtProgressBar(min = 0, max = n_boot, style = 3)
  }

  # Bootstrap function wrapper
  boot_fn = function(b) {
    # Generate bootstrap indices (stratified or unstratified)
    if (stratified) {
      # Stratified bootstrap: sample separately from treated and control
      boot_treated = sample(treated_indices, n_treated, replace = TRUE)
      boot_control = sample(control_indices, n_control, replace = TRUE)
      indices = c(boot_treated, boot_control)
    } else {
      # Unstratified bootstrap: simple random sampling
      indices = sample(1:n, n, replace = TRUE)
    }

    # Select bootstrap method
    if (bootstrap_method == "resample") {
      result = boot_iter_resample(
        indices = indices,
        data = analysis_data,
        outcome = outcome,
        treatment = treatment,
        treated_level = treated_level,
        control_level = control_level,
        f.out = f.out,
        wgt = wgt
      )
    } else {
      # reestimate_ps
      result = boot_iter_reweight(
        indices = indices,
        data = analysis_data,
        outcome = outcome,
        treatment = treatment,
        treated_level = treated_level,
        control_level = control_level,
        f.ps = f.ps,
        f.out = f.out,
        ps_params = ps_params,
        seed_offset = b
      )
    }

    if (!parallel && verbose) {
      setTxtProgressBar(pb, b)
    }
    return(result)
  }

  # Run bootstrap
  if (parallel) {
    boot_results = parallel::mclapply(
      1:n_boot,
      boot_fn,
      mc.cores = n_cores,
      mc.set.seed = TRUE
    )
  } else {
    boot_results = lapply(1:n_boot, boot_fn)
  }

  if (!parallel && verbose) {
    close(pb)
  }

  # Extract results
  boot_matrix = do.call(rbind, boot_results)
  aipw_boots = boot_matrix[, "aipw"]
  drs_boots = boot_matrix[, "drs"]

  # Extract balance metrics if available (for resample and reweight methods)
  if (ncol(boot_matrix) > 2) {
    avg_asd_boots = boot_matrix[, "avg_asd"]
    max_asd_boots = boot_matrix[, "max_asd"]
    ess_boots = boot_matrix[, "ess"]
  } else {
    avg_asd_boots = NULL
    max_asd_boots = NULL
    ess_boots = NULL
  }

  # Remove failed iterations
  n_failed = sum(is.na(aipw_boots))
  if (n_failed > 0 && verbose) {
    cat(sprintf("\nWarning: %d bootstrap iterations failed\n", n_failed))
  }

  # Keep track of valid indices for all metrics
  valid_idx = !is.na(aipw_boots)
  aipw_boots = aipw_boots[valid_idx]
  drs_boots = drs_boots[valid_idx]

  # Also filter balance metrics if they exist
  if (!is.null(avg_asd_boots)) {
    avg_asd_boots = avg_asd_boots[valid_idx]
    max_asd_boots = max_asd_boots[valid_idx]
    ess_boots = ess_boots[valid_idx]
  }

  # Compute statistics
  # AIPW
  aipw_se = sd(aipw_boots)
  aipw_ci_percentile = quantile(aipw_boots, probs = c(0.025, 0.975))
  aipw_ci_normal = c(
    aipw_point$att - 1.96 * aipw_se,
    aipw_point$att + 1.96 * aipw_se
  )
  aipw_z = aipw_point$att / aipw_se
  aipw_pval = 2 * pnorm(-abs(aipw_z))

  # DRS
  drs_se = sd(drs_boots)
  drs_ci_percentile = quantile(drs_boots, probs = c(0.025, 0.975))
  drs_ci_normal = c(
    drs_point$att - 1.96 * drs_se,
    drs_point$att + 1.96 * drs_se
  )
  drs_z = drs_point$att / drs_se
  drs_pval = 2 * pnorm(-abs(drs_z))

  # Print results
  if (verbose) {
    cat("\n=== AIPW Results ===\n")
    cat(sprintf("ATT estimate:     %.4f\n", aipw_point$att))
    cat(sprintf("Standard error:   %.4f\n", aipw_se))
    cat(sprintf(
      "95%% CI (percentile): [%.4f, %.4f]\n",
      aipw_ci_percentile[1],
      aipw_ci_percentile[2]
    ))
    cat(sprintf(
      "95%% CI (normal):     [%.4f, %.4f]\n",
      aipw_ci_normal[1],
      aipw_ci_normal[2]
    ))
    cat(sprintf("p-value:          %.4f\n", aipw_pval))

    cat("\n=== DRS Results ===\n")
    cat(sprintf("ATT estimate:     %.4f\n", drs_point$att))
    cat(sprintf("Standard error:   %.4f\n", drs_se))
    cat(sprintf(
      "95%% CI (percentile): [%.4f, %.4f]\n",
      drs_ci_percentile[1],
      drs_ci_percentile[2]
    ))
    cat(sprintf(
      "95%% CI (normal):     [%.4f, %.4f]\n",
      drs_ci_normal[1],
      drs_ci_normal[2]
    ))
    cat(sprintf("p-value:          %.4f\n", drs_pval))

    # Display balance diagnostics if available
    if (!is.null(avg_asd_boots)) {
      # Compute balance summary statistics
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
      cat(sprintf("Average ASD across bootstrap samples: %.4f (SD: %.4f)\n",
                  balance_summary$avg_asd$mean, balance_summary$avg_asd$sd))
      cat(sprintf("  Range: [%.4f, %.4f]\n",
                  balance_summary$avg_asd$min, balance_summary$avg_asd$max))
      cat(sprintf("Maximum ASD across bootstrap samples: %.4f (SD: %.4f)\n",
                  balance_summary$max_asd$mean, balance_summary$max_asd$sd))
      cat(sprintf("  Range: [%.4f, %.4f]\n",
                  balance_summary$max_asd$min, balance_summary$max_asd$max))
      cat(sprintf("Effective Sample Size: %.1f (SD: %.1f)\n",
                  balance_summary$ess$mean, balance_summary$ess$sd))
      cat(sprintf("  Range: [%.1f, %.1f]\n",
                  balance_summary$ess$min, balance_summary$ess$max))
    }
  }

  # Generate diagnostic plots
  if (plot_diagnostics) {
    par(mfrow = c(2, 2), mar = c(4, 4, 2, 1))

    # AIPW histogram
    hist(
      aipw_boots,
      main = sprintf("AIPW Bootstrap Distribution (%s)", bootstrap_method),
      xlab = "ATT Estimate",
      col = "lightblue",
      border = "white",
      breaks = 30
    )
    abline(v = aipw_point$att, col = "red", lwd = 2, lty = 2)

    # AIPW QQ-plot
    qqnorm(aipw_boots, main = "AIPW Q-Q Plot")
    qqline(aipw_boots, col = "red", lwd = 2)

    # DRS histogram
    hist(
      drs_boots,
      main = sprintf("DRS Bootstrap Distribution (%s)", bootstrap_method),
      xlab = "ATT Estimate",
      col = "lightgreen",
      border = "white",
      breaks = 30
    )
    abline(v = drs_point$att, col = "red", lwd = 2, lty = 2)

    # DRS QQ-plot
    qqnorm(drs_boots, main = "DRS Q-Q Plot")
    qqline(drs_boots, col = "red", lwd = 2)

    par(mfrow = c(1, 1))
  }

  # Create balance_diagnostics if metrics are available
  balance_diagnostics = NULL
  bootstrap_balance = NULL

  if (!is.null(avg_asd_boots)) {
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
      iteration = 1:length(avg_asd_boots),
      avg_asd = avg_asd_boots,
      max_asd = max_asd_boots,
      ess = ess_boots
    )
  }

  # Return results
  return(list(
    aipw = list(
      att = aipw_point$att,
      se = aipw_se,
      ci_percentile = aipw_ci_percentile,
      ci_normal = aipw_ci_normal,
      pval = aipw_pval,
      bootstrap_samples = aipw_boots
    ),
    drs = list(
      att = drs_point$att,
      se = drs_se,
      ci_percentile = drs_ci_percentile,
      ci_normal = drs_ci_normal,
      pval = drs_pval,
      bootstrap_samples = drs_boots
    ),
    n_boot = length(aipw_boots),
    n_failed = n_failed,
    balance_diagnostics = balance_diagnostics,
    bootstrap_balance = bootstrap_balance
  ))
}

# diagnostics for ps
# ps_inc_dec$desc$es.mean.ATT$ess.ctrl
# ps_inc_dec$desc$es.mean.ATT$mean.es
# hist(ps_inc_dec$ps$es.mean.ATT)
# hist(ps_inc_dec$w$es.mean.ATT)

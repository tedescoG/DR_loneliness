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
      mc.set.seed = F
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

## Cross-fitted AIPW -------------------------------------------------------------------####
cf_aipw_att = function(
  outcome,
  treatment,
  f.ps,
  f.out,
  data,
  k = 5,
  stratify = TRUE,
  gbm_params = list(
    n.trees = 10000,
    interaction.depth = 2,
    shrinkage = 0.01,
    bag.fraction = 0.5,
    n.minobsinnode = 10
  ),
  seed = 123,
  verbose = TRUE,
  alpha = 0
) {
  #' Cross-fitted AIPW estimator for ATT
  #'
  #' @param outcome character - outcome variable name
  #' @param treatment character - treatment variable name (binary 0/1)
  #' @param f.ps formula - propensity score formula
  #' @param f.out character - outcome model formula (RHS only)
  #' @param data data.frame - dataset
  #' @param k integer - number of folds
  #' @param stratify logical - stratify folds by treatment
  #' @param gbm_params list - GBM hyperparameters (no internal tuning)
  #' @param seed integer - random seed
  #' @param verbose logical - print progress
  #' @param alpha numeric - propensity score truncation level (0 = no truncation)
  #'
  #' @return list with ATT estimate and counterfactual means

  # Extract variables
  Y = data[[outcome]]
  Z = data[[treatment]]
  n = nrow(data)
  N1 = sum(Z)

  # Create folds
  stratify_var = if (stratify) treatment else NULL
  set.seed(seed)
  folds = vfold_cv(data, v = k, strata = all_of(stratify_var))

  if (verbose) {
    cat(sprintf("Cross-fitted AIPW with %d folds\n", k))
    if (stratify) cat("Using stratified folding by treatment\n")
  }

  # Initialize storage for fold-specific predictions and ATT estimates
  mu_hat_0 = numeric(n)
  weights = numeric(n)
  att_estimates = numeric(k)

  # Cross-fitting loop
  for (fold in 1:k) {
    if (verbose) {
      cat(sprintf("Processing fold %d/%d...\n", fold, k))
    }

    # Split data: predict on fold k, train on all other folds
    test_idx = folds$splits[[fold]]$in_id
    train_idx = setdiff(1:n, test_idx)

    train_data = data[train_idx, ]
    test_data = data[test_idx, ]

    # Fit propensity score model on K-1 folds
    set.seed(seed + fold) # Ensure reproducibility across folds
    ps_fit = do.call(
      "ps",
      c(
        list(
          formula = f.ps,
          data = train_data,
          estimand = "ATT",
          stop.method = "es.mean",
          verbose = FALSE
        ),
        gbm_params
      )
    )

    # Predict propensity scores for held-out fold k
    ps_pred = predict(
      ps_fit$gbm.obj,
      newdata = test_data,
      n.trees = ps_fit$desc$es.mean$n.trees,
      type = "response"
    )

    # Truncate propensity scores if alpha > 0
    if (alpha > 0) {
      ps_pred = pmax(alpha, pmin(ps_pred, 1 - alpha))
    }

    # Compute IPT weights for held-out fold (ATT weights)
    weights[test_idx] = ifelse(
      Z[test_idx] == 1,
      1,
      ps_pred / (1 - ps_pred)
    )

    # Fit outcome model on control units from K-1 folds
    control_train = train_data[train_data[[treatment]] == 0, ]
    mu0_fit = glm(
      formula = as.formula(paste(outcome, "~", f.out)),
      data = control_train,
      family = "quasibinomial",
      control = glm.control(maxit = 50)
    )

    # Predict counterfactual outcomes for held-out fold k
    mu_hat_0[test_idx] = predict(
      mu0_fit,
      newdata = test_data,
      type = "response"
    )

    # Compute fold-specific AIPW estimate using only fold k data
    Y_k = Y[test_idx]
    Z_k = Z[test_idx]
    weights_k = weights[test_idx]
    mu_hat_0_k = mu_hat_0[test_idx]
    N1_k = sum(Z_k)

    att_estimates[fold] = (1 / N1_k) *
      sum((Z_k - (1 - Z_k) * weights_k) * (Y_k - mu_hat_0_k))
  }

  # Average fold-specific ATT estimates for final cross-fitted estimate
  att_est = mean(att_estimates)

  # Output
  if (verbose) {
    cat(sprintf("\nCross-fitted AIPW ATT: %.4f\n", att_est))
    cat(sprintf("Mean Y(0) [counterfactual]: %.4f\n", mean(mu_hat_0[Z == 0])))
    cat(sprintf("Mean Y(1) [observed]: %.4f\n", mean(Y[Z == 1])))
  }

  # Return results (same structure as aipw_att for compatibility)
  return(list(
    att = att_est,
    mu_hat_0 = mean(mu_hat_0[Z == 0]),
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

## Cross-fitted DRS -------------------------------------------------------------------####
cf_drs_att = function(
  outcome,
  treatment,
  f.ps,
  f.out,
  data,
  k = 5,
  stratify = TRUE,
  gbm_params = list(
    n.trees = 10000,
    interaction.depth = 2,
    shrinkage = 0.01,
    bag.fraction = 0.5,
    n.minobsinnode = 10
  ),
  seed = 123,
  verbose = TRUE,
  alpha = 0
) {
  #' Cross-fitted DRS (Doubly Robust Standardization) ATT estimator
  #'
  #' Implements K-fold cross-fitting for DRS estimation to reduce overfitting bias
  #' when using flexible machine learning methods (GBM) for propensity scores.
  #'
  #' @param outcome character - name of outcome variable
  #' @param treatment character - name of treatment variable (must be 0/1 numeric)
  #' @param f.ps formula - propensity score model formula
  #' @param f.out character - outcome model formula (RHS only, excluding outcome and treatment)
  #' @param data data.frame - dataset (treatment should already be 0/1 numeric)
  #' @param k integer - number of cross-fitting folds
  #' @param stratify logical - stratify folds by treatment status
  #' @param gbm_params list - GBM parameters for ps() function
  #' @param seed integer - random seed for fold creation
  #' @param verbose logical - print progress and results
  #' @param alpha numeric - propensity score truncation level (0 = no truncation)
  #'
  #' @return list with att estimate and fold-specific estimates

  if (verbose) {
    cat("=== Cross-Fitted DRS ATT Estimation ===\n")
    cat(sprintf("K-fold cross-fitting: k=%d\n", k))
    cat(sprintf("Stratified folds: %s\n\n", stratify))
  }

  n = nrow(data)

  # Extract outcome and treatment
  Y = data[[outcome]]
  Z = data[[treatment]]

  # Verify treatment is 0/1
  if (!all(Z %in% c(0, 1))) {
    stop("Treatment variable must be 0/1 numeric")
  }

  # Create K folds
  stratify_var = if (stratify) treatment else NULL
  set.seed(seed)
  folds = vfold_cv(data, v = k, strata = all_of(stratify_var))

  # Initialize storage
  weights = numeric(n)
  drs_fold_estimates = numeric(k)

  if (verbose) {
    cat("Processing folds:\n")
  }

  # Cross-fitting loop
  for (fold in 1:k) {
    if (verbose) {
      cat(sprintf("  Fold %d/%d...\n", fold, k))
    }

    # Split data
    test_idx = folds$splits[[fold]]$in_id
    train_idx = setdiff(1:n, test_idx)

    train_data = data[train_idx, ]
    test_data = data[test_idx, ]

    # Fit propensity score model on K-1 folds
    set.seed(seed)
    ps_fit = do.call(
      "ps",
      c(
        list(
          formula = f.ps,
          data = train_data,
          estimand = "ATT",
          stop.method = "es.mean",
          verbose = FALSE
        ),
        gbm_params
      )
    )

    # Predict propensity scores for held-out fold
    ps_pred = predict(
      ps_fit$gbm.obj,
      newdata = test_data,
      n.trees = ps_fit$desc$es.mean$n.trees,
      type = "response"
    )

    # Truncate propensity scores if alpha > 0
    if (alpha > 0) {
      ps_pred = pmax(alpha, pmin(ps_pred, 1 - alpha))
    }

    # Compute ATT weights for held-out fold
    weights[test_idx] = ifelse(
      Z[test_idx] == 1,
      1,
      ps_pred / (1 - ps_pred)
    )

    # DRS estimation on held-out fold
    test_data_wgt = test_data
    test_data_wgt$ps_wgt = weights[test_idx]

    weighted_design = svydesign(
      ids = ~1,
      weights = ~ps_wgt,
      data = test_data_wgt
    )

    # Build outcome model formula
    if (f.out != "1") {
      formula_str = paste(outcome, "~", treatment, "+", f.out)
    } else {
      formula_str = paste(outcome, "~", treatment)
    }

    # Fit weighted outcome model on held-out fold
    fit = svyglm(
      formula = as.formula(formula_str),
      family = 'quasibinomial',
      design = weighted_design
    )

    # Extract fold-specific outcomes
    Y_k = Y[test_idx]
    Z_k = Z[test_idx]

    # Predict counterfactuals for treated units in fold
    treated_k = test_data_wgt[test_data_wgt[[treatment]] == 1, ]
    counterfactual_k = treated_k
    counterfactual_k[[treatment]] = 0

    mu0_k = predict(fit, newdata = counterfactual_k, type = "response")
    mu_hat_0_k = mean(mu0_k)
    mu_hat_1_k = mean(Y_k[Z_k == 1])

    # Fold-specific ATT
    drs_fold_estimates[fold] = mu_hat_1_k - mu_hat_0_k

    if (verbose) {
      cat(sprintf("    Fold %d ATT: %.4f\n", fold, drs_fold_estimates[fold]))
    }
  }

  # Average across folds for final estimate
  att_est = mean(drs_fold_estimates)

  if (verbose) {
    cat(sprintf("\nCross-fitted DRS ATT: %.4f\n", att_est))
  }

  # Return results (compatible with drs_att structure)
  return(list(
    att = att_est
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
  cross_fitting = FALSE,
  k = 5,
  stratify_folds = TRUE,
  wgt = "iptw",
  seed_offset,
  alpha = 0
) {
  #' Unified bootstrap iteration function supporting all estimation methods
  #'
  #' Supports: resample/reweight × non-cross-fitted/cross-fitted (4 combinations)
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
  #' @param cross_fitting logical - use cross-fitting or not
  #' @param k integer - number of cross-fitting folds
  #' @param stratify_folds logical - stratify folds by treatment (for cross-fitting)
  #' @param wgt character - name of weight variable (for resample method)
  #' @param seed_offset integer - base seed for reproducibility
  #' @param alpha numeric - propensity score truncation level (0 = no truncation)
  #'
  #' @return numeric vector: c(aipw, drs, avg_asd, max_asd, ess)

  # Wrap entire function in tryCatch for error handling
  tryCatch(
    {
      # Match bootstrap method
      bootstrap_method = match.arg(bootstrap_method)

      # Resample data
      boot_data = data[indices, ]

      # Create binary treatment indicator
      boot_data[[treatment]] = as.numeric(
        boot_data[[treatment]] == treated_level
      )

      # Initialize estimates and weights
      aipw_est = NA_real_
      drs_est = NA_real_
      weights_vec = NULL

      if (!cross_fitting) {
        # ========== NON-CROSS-FITTED ESTIMATION ==========

        if (bootstrap_method == "resample") {
          # Use existing weights from data
          aipw_result = aipw_att(
            outcome = outcome,
            treatment = treatment,
            f.out = f.out,
            wgt = wgt,
            data = boot_data,
            verbose = FALSE
          )

          drs_result = drs_att(
            outcome = outcome,
            treatment = treatment,
            f.out = f.out,
            wgt = wgt,
            data = boot_data,
            verbose = FALSE
          )

          weights_vec = boot_data[[wgt]]
        } else {
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

          # Compute estimates
          aipw_result = aipw_att(
            outcome = outcome,
            treatment = treatment,
            f.out = f.out,
            wgt = "ps_wgt",
            data = boot_data,
            verbose = FALSE
          )

          drs_result = drs_att(
            outcome = outcome,
            treatment = treatment,
            f.out = f.out,
            wgt = "ps_wgt",
            data = boot_data,
            verbose = FALSE
          )

          weights_vec = boot_data$ps_wgt
        }

        aipw_est = aipw_result$att
        drs_est = drs_result$att
      } else {
        # ========== CROSS-FITTED ESTIMATION ==========

        n = nrow(boot_data)
        Y = boot_data[[outcome]]
        Z = boot_data[[treatment]]
        N1 = sum(Z)

        # Create K folds with stratification
        # Use indices to create unique seed for each bootstrap iteration
        stratify_var = if (stratify_folds) treatment else NULL
        set.seed(seed_offset + min(indices))
        folds = vfold_cv(
          boot_data,
          v = k,
          strata = all_of(stratify_var)
        )

        # Initialize storage for fold-specific estimates
        mu_hat_0 = numeric(n)
        weights = numeric(n)
        aipw_fold_estimates = numeric(k)
        drs_fold_estimates = numeric(k)

        # Cross-fitting loop
        for (fold in 1:k) {
          # Split data: predict on fold k, train on K-1 folds
          test_idx = folds$splits[[fold]]$in_id
          train_idx = setdiff(1:n, test_idx)

          train_data = boot_data[train_idx, ]
          test_data = boot_data[test_idx, ]

          # Fit propensity score model on K-1 folds
          set.seed(seed_offset) # Fixed seed for GBM within bootstrap
          ps_fit = do.call(
            "ps",
            c(
              list(
                formula = f.ps,
                data = train_data,
                estimand = "ATT",
                stop.method = "es.mean",
                verbose = FALSE
              ),
              gbm_params
            )
          )

          # Predict propensity scores for held-out fold k
          ps_pred = predict(
            ps_fit$gbm.obj,
            newdata = test_data,
            n.trees = ps_fit$desc$es.mean$n.trees,
            type = "response"
          )

          # Truncate propensity scores if alpha > 0
          if (alpha > 0) {
            ps_pred = pmax(alpha, pmin(ps_pred, 1 - alpha))
          }

          # Compute IPT weights for held-out fold (ATT weights)
          weights[test_idx] = ifelse(
            Z[test_idx] == 1,
            1,
            ps_pred / (1 - ps_pred)
          )

          # Fit outcome model on control units from K-1 folds
          control_train = train_data[train_data[[treatment]] == 0, ]
          mu0_fit = glm(
            formula = as.formula(paste(outcome, "~", f.out)),
            data = control_train,
            family = "quasibinomial",
            control = glm.control(maxit = 50)
          )

          # Predict counterfactual outcomes for held-out fold k
          mu_hat_0[test_idx] = predict(
            mu0_fit,
            newdata = test_data,
            type = "response"
          )

          # Compute fold-specific AIPW estimate
          Y_k = Y[test_idx]
          Z_k = Z[test_idx]
          weights_k = weights[test_idx]
          mu_hat_0_k = mu_hat_0[test_idx]
          N1_k = sum(Z_k)

          aipw_fold_estimates[fold] = (1 / N1_k) *
            sum((Z_k - (1 - Z_k) * weights_k) * (Y_k - mu_hat_0_k))

          # Compute fold-specific DRS estimate
          # Create weighted survey design for fold k
          test_data_wgt = test_data
          test_data_wgt$ps_wgt = weights[test_idx]

          weighted_design = svydesign(
            ids = ~1,
            weights = ~ps_wgt,
            data = test_data_wgt
          )

          # Build outcome model formula
          if (f.out != "1") {
            formula_str = paste(outcome, "~", treatment, "+", f.out)
          } else {
            formula_str = paste(outcome, "~", treatment)
          }

          # Fit weighted outcome model on fold k
          fit = svyglm(
            formula = as.formula(formula_str),
            family = 'quasibinomial',
            design = weighted_design
          )

          # Predict counterfactuals for treated units in fold k
          treated_k = test_data_wgt[test_data_wgt[[treatment]] == 1, ]
          counterfactual_k = treated_k
          counterfactual_k[[treatment]] = 0

          mu0_drs_k = predict(
            fit,
            newdata = counterfactual_k,
            type = "response"
          )
          mu_hat_0_drs = mean(mu0_drs_k)
          mu_hat_1_drs = mean(Y_k[Z_k == 1])

          drs_fold_estimates[fold] = mu_hat_1_drs - mu_hat_0_drs
        }

        # Average across folds for final bootstrap estimates
        aipw_est = mean(aipw_fold_estimates)
        drs_est = mean(drs_fold_estimates)

        # Use full sample weights for balance diagnostics
        weights_vec = weights
      }

      # ========== COMPUTE BALANCE STATISTICS ==========
      avg_asd = NA_real_
      max_asd = NA_real_
      ess = NA_real_

      tryCatch(
        {
          # Extract all covariates (exclude treatment, outcome, and weight columns)
          exclude_cols = c(treatment, outcome)
          if (!cross_fitting && bootstrap_method == "resample") {
            exclude_cols = c(exclude_cols, wgt)
          }
          if (!cross_fitting && bootstrap_method == "reweight") {
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

      # Return results
      return(c(
        aipw = aipw_est,
        drs = drs_est,
        avg_asd = avg_asd,
        max_asd = max_asd,
        ess = ess
      ))
    },
    error = function(e) {
      # Return NAs if anything fails during bootstrap iteration
      return(c(
        aipw = NA_real_,
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
  cross_fitting = FALSE,
  k = 5,
  stratify_folds = TRUE,
  n_boot = 1000,
  seed = 123,
  verbose = TRUE,
  parallel = TRUE,
  n_cores = NULL,
  plot_diagnostics = TRUE,
  ci_type = c("all", "norm", "basic", "perc"),
  sim = c("ordinary", "balanced"),
  save_to = NULL,
  alpha = 0
) {
  #' Unified Bootstrap Doubly Robust ATT Estimation
  #'
  #' Supports all combinations of:
  #' - Bootstrap methods: resample (use existing weights) or reweight (re-estimate PS)
  #' - Simulation schemes: ordinary or balanced bootstrap
  #' - Cross-fitting: TRUE or FALSE
  #' - Confidence intervals: normal, basic, percentile (no BCa)
  #'
  #' @param outcome character - name of outcome variable
  #' @param treatment character - name of treatment variable
  #' @param treated_level string - level of treatment variable indicating treated
  #' @param control_level string - level of treatment variable indicating control
  #' @param f.ps formula - propensity score model formula (required for reweight and cross-fitting)
  #' @param f.out character - outcome model formula (RHS only)
  #' @param data data.frame - dataset (must contain wgt if bootstrap_method = "resample")
  #' @param gbm_params list - GBM parameters for ps() function
  #' @param bootstrap_method character - "resample" (faster) or "reweight" (more conservative)
  #' @param stratified logical - use stratified bootstrap to maintain treatment/control proportions
  #' @param wgt character - name of propensity score weight column in data (for resample method)
  #' @param cross_fitting logical - use cross-fitting to reduce overfitting bias
  #' @param k integer - number of cross-fitting folds
  #' @param stratify_folds logical - stratify folds by treatment (within cross-fitting)
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
  #'
  #' @return list with estimates, SEs, CIs, p-values, balance diagnostics, and bootstrap samples

  # Match and validate arguments
  bootstrap_method = match.arg(bootstrap_method)
  ci_type = match.arg(ci_type)
  sim = match.arg(sim)

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

  if (cross_fitting && is.null(f.ps)) {
    stop("f.ps formula required when cross_fitting = TRUE")
  }

  if (bootstrap_method == "reweight" && is.null(f.ps)) {
    stop("f.ps formula required when bootstrap_method = 'reweight'")
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
    cat("Cross-fitting:", cross_fitting, "\n")
    if (cross_fitting) {
      cat("  Number of folds:", k, "\n")
      cat("  Stratify folds:", stratify_folds, "\n")
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
    cat("Running bootstrap using boot package...\n")
    if (cross_fitting) {
      cat("Note: Each bootstrap iteration performs", k, "cross-fitting folds\n")
      cat("Total model fits:", n_boot, "×", k, "=", n_boot * k, "\n\n")
    }
  }

  # Set seed for reproducibility
  set.seed(seed)

  # Run boot() with unified boot_iter function
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
    cross_fitting = cross_fitting,
    k = k,
    stratify_folds = stratify_folds,
    wgt = wgt,
    seed_offset = seed,
    alpha = alpha
  )

  # Extract point estimates from boot object (t0)
  aipw_point_att = boot_result$t0[1]
  drs_point_att = boot_result$t0[2]

  if (verbose) {
    cat("\nPoint Estimates:\n")
    cat(sprintf("  AIPW ATT: %.4f\n", aipw_point_att))
    cat(sprintf("  DRS ATT:  %.4f\n\n", drs_point_att))
  }

  # Extract bootstrap samples
  aipw_boots = boot_result$t[, 1]
  drs_boots = boot_result$t[, 2]

  # Extract balance metrics if available
  if (ncol(boot_result$t) >= 5) {
    avg_asd_boots = boot_result$t[, 3]
    max_asd_boots = boot_result$t[, 4]
    ess_boots = boot_result$t[, 5]
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

  # Compute standard errors
  aipw_se = sd(aipw_boots)
  drs_se = sd(drs_boots)

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

  # Compute CIs for AIPW (index = 1)
  aipw_ci = tryCatch(
    {
      boot::boot.ci(boot_result_filtered, type = ci_types_to_compute, index = 1)
    },
    error = function(e) {
      warning("boot.ci() failed for AIPW, falling back to manual computation")
      NULL
    }
  )

  # Compute CIs for DRS (index = 2)
  drs_ci = tryCatch(
    {
      boot::boot.ci(boot_result_filtered, type = ci_types_to_compute, index = 2)
    },
    error = function(e) {
      warning("boot.ci() failed for DRS, falling back to manual computation")
      NULL
    }
  )

  # Extract CIs or fall back to manual computation
  if (!is.null(aipw_ci)) {
    aipw_ci_normal = if ("normal" %in% names(aipw_ci)) {
      aipw_ci$normal[2:3]
    } else {
      c(aipw_point_att - 1.96 * aipw_se, aipw_point_att + 1.96 * aipw_se)
    }
    aipw_ci_basic = if ("basic" %in% names(aipw_ci)) {
      aipw_ci$basic[4:5]
    } else {
      NULL
    }
    aipw_ci_percentile = if ("percent" %in% names(aipw_ci)) {
      aipw_ci$percent[4:5]
    } else {
      quantile(aipw_boots, probs = c(0.025, 0.975))
    }
  } else {
    aipw_ci_normal = c(
      aipw_point_att - 1.96 * aipw_se,
      aipw_point_att + 1.96 * aipw_se
    )
    aipw_ci_basic = NULL
    aipw_ci_percentile = quantile(aipw_boots, probs = c(0.025, 0.975))
  }

  if (!is.null(drs_ci)) {
    drs_ci_normal = if ("normal" %in% names(drs_ci)) {
      drs_ci$normal[2:3]
    } else {
      c(drs_point_att - 1.96 * drs_se, drs_point_att + 1.96 * drs_se)
    }
    drs_ci_basic = if ("basic" %in% names(drs_ci)) drs_ci$basic[4:5] else NULL
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

  # Compute p-values using normal approximation
  aipw_z = aipw_point_att / aipw_se
  aipw_pval = 2 * pnorm(-abs(aipw_z))

  drs_z = drs_point_att / drs_se
  drs_pval = 2 * pnorm(-abs(drs_z))

  # Print results
  if (verbose) {
    cat("\n=== AIPW Results ===\n")
    cat(sprintf("ATT estimate:     %.4f\n", aipw_point_att))
    cat(sprintf("Standard error:   %.4f\n", aipw_se))
    cat(sprintf("p-value:          %.4f\n", aipw_pval))
    cat("\nConfidence Intervals (95%):\n")
    if (!is.null(aipw_ci_normal)) {
      cat(sprintf(
        "  Normal:      [%.4f, %.4f]\n",
        aipw_ci_normal[1],
        aipw_ci_normal[2]
      ))
    }
    if (!is.null(aipw_ci_basic)) {
      cat(sprintf(
        "  Basic:       [%.4f, %.4f]\n",
        aipw_ci_basic[1],
        aipw_ci_basic[2]
      ))
    }
    if (!is.null(aipw_ci_percentile)) {
      cat(sprintf(
        "  Percentile:  [%.4f, %.4f]\n",
        aipw_ci_percentile[1],
        aipw_ci_percentile[2]
      ))
    }

    cat("\n=== DRS Results ===\n")
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

  # Generate diagnostic plots
  if (plot_diagnostics) {
    # Display plots first
    par(mfrow = c(2, 2), mar = c(4, 4, 2, 1))

    # AIPW histogram
    hist(
      aipw_boots,
      main = sprintf(
        "AIPW Bootstrap Distribution (%s%s)",
        bootstrap_method,
        if (cross_fitting) paste0(", CF k=", k) else ""
      ),
      xlab = "ATT Estimate",
      col = "lightblue",
      border = "white",
      breaks = 30
    )
    abline(v = aipw_point_att, col = "red", lwd = 2, lty = 2)

    # AIPW QQ-plot
    qqnorm(aipw_boots, main = "AIPW Q-Q Plot")
    qqline(aipw_boots, col = "red", lwd = 2)

    # DRS histogram
    hist(
      drs_boots,
      main = sprintf(
        "DRS Bootstrap Distribution (%s%s)",
        bootstrap_method,
        if (cross_fitting) paste0(", CF k=", k) else ""
      ),
      xlab = "ATT Estimate",
      col = "lightgreen",
      border = "white",
      breaks = 30
    )
    abline(v = drs_point_att, col = "red", lwd = 2, lty = 2)

    # DRS QQ-plot
    qqnorm(drs_boots, main = "DRS Q-Q Plot")
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
        "DR_bootstrap_%s_vs_%s_%s_%s%s_n%d.png",
        treated_level,
        control_level,
        bootstrap_method,
        sim,
        if (cross_fitting) paste0("_CF", k) else "",
        n_boot
      )
      filepath = file.path(save_to, filename)

      # Open high-resolution PNG device for academic paper quality
      png(filepath, width = 3000, height = 3000, res = 300)

      par(mfrow = c(2, 2), mar = c(4, 4, 2, 1))

      # AIPW histogram
      hist(
        aipw_boots,
        main = sprintf(
          "AIPW Bootstrap Distribution (%s%s)",
          bootstrap_method,
          if (cross_fitting) paste0(", CF k=", k) else ""
        ),
        xlab = "ATT Estimate",
        col = "lightblue",
        border = "white",
        breaks = 30
      )
      abline(v = aipw_point_att, col = "red", lwd = 2, lty = 2)

      # AIPW QQ-plot
      qqnorm(aipw_boots, main = "AIPW Q-Q Plot")
      qqline(aipw_boots, col = "red", lwd = 2)

      # DRS histogram
      hist(
        drs_boots,
        main = sprintf(
          "DRS Bootstrap Distribution (%s%s)",
          bootstrap_method,
          if (cross_fitting) paste0(", CF k=", k) else ""
        ),
        xlab = "ATT Estimate",
        col = "lightgreen",
        border = "white",
        breaks = 30
      )
      abline(v = drs_point_att, col = "red", lwd = 2, lty = 2)

      # DRS QQ-plot
      qqnorm(drs_boots, main = "DRS Q-Q Plot")
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
      att = aipw_point_att,
      se = aipw_se,
      ci_normal = aipw_ci_normal,
      ci_basic = aipw_ci_basic,
      ci_percentile = aipw_ci_percentile,
      pval = aipw_pval,
      bootstrap_samples = aipw_boots,
      boot_ci_object = aipw_ci
    ),
    drs = list(
      att = drs_point_att,
      se = drs_se,
      ci_normal = drs_ci_normal,
      ci_basic = drs_ci_basic,
      ci_percentile = drs_ci_percentile,
      pval = drs_pval,
      bootstrap_samples = drs_boots,
      boot_ci_object = drs_ci
    ),
    n_boot = length(aipw_boots),
    n_failed = n_failed,
    balance_diagnostics = balance_diagnostics,
    bootstrap_balance = bootstrap_balance,
    boot_object = boot_result # Include boot object for additional analyses
  ))
}

# HELPER X RESULTS ---------------------------------------------------------------------####
extract_results = function(result, model_num, method, estimator) {
  est = result[[tolower(estimator)]]

  data.frame(
    Model = model_num,
    Method = method,
    Estimator = estimator,
    ATT = est$att,
    SE = est$se,
    CI_Normal_Lower = est$ci_normal[1],
    CI_Normal_Upper = est$ci_normal[2],
    CI_Basic_Lower = ifelse(!is.null(est$ci_basic), est$ci_basic[1], NA),
    CI_Basic_Upper = ifelse(!is.null(est$ci_basic), est$ci_basic[2], NA),
    CI_Perc_Lower = est$ci_percentile[1],
    CI_Perc_Upper = est$ci_percentile[2],
    p_value = est$pval,
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

plot_model_comparison <- function(
  results,
  save_path,
  width = 3000,
  height = 2000,
  res = 300,
  colors = c("blue", "red", "green", "purple"),
  ylim_padding = 1.05,
  aipw_title = "AIPW: Reweight Method",
  drs_title = "DRS: Reweight Method"
) {
  #' Plot Model Comparison Distributions
  #'
  #' Creates a comparison plot of bootstrap distributions for multiple models,
  #' showing both AIPW and DRS methods side by side. Displays the plot and saves to file.
  #'
  #' @param results A list of model results, each containing `aipw` and `drs`
  #'   components with `bootstrap_samples`
  #' @param save_path Path where the plot should be saved (PNG format)
  #' @param width Plot width in pixels (default: 3000)
  #' @param height Plot height in pixels (default: 2000)
  #' @param res Resolution in DPI (default: 300)
  #' @param colors Vector of colors for each model (default: blue, red, green, purple)
  #' @param ylim_padding Padding factor for y-axis limit (default: 1.05 for 5% padding)
  #' @param aipw_title Title for AIPW panel (default: "AIPW: Reweight Method")
  #' @param drs_title Title for DRS panel (default: "DRS: Reweight Method")
  #'
  #' @return NULL (displays plot and saves to file)

  # Extract model names
  model_names <- names(results)
  n_models <- length(model_names)

  # Ensure we have enough colors
  if (length(colors) < n_models) {
    colors <- rep(colors, length.out = n_models)
  }

  # === AIPW Panel ===

  # Extract all AIPW bootstrap samples
  aipw_samples <- lapply(results, function(x) {
    x$aipw$bootstrap_samples
  })

  # Calculate densities
  aipw_densities <- lapply(aipw_samples, density)

  # Calculate y-axis limit
  aipw_ylim <- max(sapply(aipw_densities, function(d) max(d$y)))

  # Calculate x-axis limit
  aipw_xlim <- range(unlist(aipw_samples))

  # === DRS Panel ===

  # Extract all DRS bootstrap samples
  drs_samples <- lapply(results, function(x) {
    x$drs$bootstrap_samples
  })

  # Calculate densities
  drs_densities <- lapply(drs_samples, density)

  # Calculate y-axis limit
  drs_ylim <- max(sapply(drs_densities, function(d) max(d$y)))

  # Calculate x-axis limit
  drs_xlim <- range(unlist(drs_samples))

  # Create plot function (used for both display and save)
  create_plot <- function() {
    par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))

    # AIPW plot
    plot(
      aipw_densities[[1]],
      main = aipw_title,
      xlab = "ATT",
      col = colors[1],
      lwd = 2,
      xlim = aipw_xlim,
      ylim = c(0, aipw_ylim * ylim_padding)
    )

    if (n_models > 1) {
      for (i in 2:n_models) {
        lines(aipw_densities[[i]], col = colors[i], lwd = 2)
      }
    }

    legend(
      "topleft",
      legend = paste("Model", 1:n_models),
      col = colors[1:n_models],
      lwd = 2,
      cex = 0.8
    )

    # DRS plot
    plot(
      drs_densities[[1]],
      main = drs_title,
      xlab = "ATT",
      col = colors[1],
      lwd = 2,
      xlim = drs_xlim,
      ylim = c(0, drs_ylim * ylim_padding)
    )

    if (n_models > 1) {
      for (i in 2:n_models) {
        lines(drs_densities[[i]], col = colors[i], lwd = 2)
      }
    }

    legend(
      "topleft",
      legend = paste("Model", 1:n_models),
      col = colors[1:n_models],
      lwd = 2,
      cex = 0.8
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

# MULTINOMIAL TREATMENT x SEVERE LONELINESS (BINARY)

# Setup
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("utils.R")

# Import data
d = readRDS(
  "/Users/gaetanotedesco/Desktop/Research/Thesis/DR_loneliness/data/data_iptw.rds"
)
d$remote_contact = relevel(d$remote_contact, ref = "increase")
treated = subset(d, remote_contact == "increase")

# outcome model
weighted_design = svydesign(ids = ~1, weights = ~iptw, data = d)


### Model 1: Outcome ~ Treatment ####
# Fitting weighted regression
fit1 = svyglm(
  severe_loneliness ~ remote_contact,
  family = 'quasibinomial',
  design = weighted_design
)

# invariant-decrease is significant increasing the prob. of loneliness
# Estimating potential outcomes
avg_comparisons(
  fit1,
  variables = "remote_contact",
  newdata = treated,
  wts = "iptw"
)


### Model 2: Outcome ~ Treatment + Baseline ####
# Fitting weighted regression
fit2 = svyglm(
  severe_loneliness ~ remote_contact + baseline_lone,
  family = 'quasibinomial',
  design = weighted_design
)

# invariant-decrease is significant increasing the prob. of loneliness
# while increase-decrease not.

# Estimating potential outcomes
avg_comparisons(
  fit2,
  variables = "remote_contact",
  newdata = treated,
  wts = "iptw"
)


### Model 3: Outcome ~ Treatment + Baseline + strong predictor ####
# Fitting weighted regression
fit3 = svyglm(
  severe_loneliness ~ remote_contact + baseline_lone + baseline_depr,
  family = 'quasibinomial',
  design = weighted_design
)

# invariant-decrease is significant increasing the prob. of loneliness
# while increase-decrease not.

# Estimating potential outcomes
avg_comparisons(
  fit3,
  variables = "remote_contact",
  newdata = treated,
  wts = "iptw"
)


### Model 4: Outcome ~ Treatment + all covariates ####
# Fitting weighted regression
fit4 = svyglm(
  severe_loneliness ~
    remote_contact +
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
      baseline_depr +
      baseline_lone +
      neighborhood,
  family = 'quasibinomial',
  design = weighted_design
)

# invariant-decrease is significant increasing the prob. of loneliness
# while increase-decrease not.

# Estimating potential outcomes
avg_comparisons(
  fit4,
  variables = "remote_contact",
  newdata = treated,
  wts = "iptw"
)

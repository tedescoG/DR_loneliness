# IPTW x MULTINOMIAL TREATMENT

# Setup
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("utils.R")
run = 123

# Parameter for ps function
params_ps = list(
  estimand = "ATT",
  stop.method = c("es.mean"),
  verbose = F,
  n.grid = 50
)

# Grid of hyperparameters for GBM
gbm.grid = list(
  distribution = c("bernoulli"),
  n.trees = c(5000, 10000),
  interaction.depth = c(1),
  shrinkage = seq(0.01, 0.05, by = 0.005),
  bag.fraction = c(.3, .5, 1, .75),
  n.minobsinnode = c(10, 15)
)

# Import data
d = readRDS("~/Desktop/Research/Thesis/DR_loneliness/data/data_clean.rds")
# Define ps function
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


# INCREASE VS DECREASE ----------------------------------------------------------------------------####

inc_dec = subset(d, remote_contact %in% c("increase", "decrease")) %>%
  mutate(across(where(is.factor), droplevels)) %>%
  mutate(remote_contact = ifelse(remote_contact == "increase", 1, 0))


# gbm optimization
system.time({
  tuned_gbm = tune.gbm(
    x = inc_dec %>%
      select(
        -c(
          "severe_loneliness",
          "remote_contact"
        )
      ),
    y = inc_dec[, "remote_contact"],
    params.grid = gbm.grid,
    method = "es.mean",
    parallel = T
  )
})


# Save hyperparams combination
saveRDS(
  tuned_gbm$params,
  "/Users/gaetanotedesco/Desktop/Research/Thesis/DR_loneliness/results/weighting/inc_dec_params.rds"
)

set.seed(run)
ps_inc_dec = do.call(
  "ps",
  c(list(formula = ps_formula, data = inc_dec), params_ps, tuned_gbm$params)
)

# Exploring/ saving results
(weight_summary_inc_dec = summary(ps_inc_dec))
saveRDS(
  weight_summary_inc_dec,
  "/Users/gaetanotedesco/Desktop/Research/Thesis/DR_loneliness/results/weighting/inc_dec/weighting_summary.rds"
)

# Global graphical assessment
# Save plots with high resolution for academic publication
# Convergence plot
png(
  "/Users/gaetanotedesco/Desktop/Research/Thesis/DR_loneliness/results/weighting/inc_dec/convergence_plot.png",
  width = 8,
  height = 6,
  units = "in",
  res = 300
)
plot(ps_inc_dec, plots = 1)
dev.off()

# Overlap plot
png(
  "/Users/gaetanotedesco/Desktop/Research/Thesis/DR_loneliness/results/weighting/inc_dec/overlap_plot.png",
  width = 8,
  height = 6,
  units = "in",
  res = 300
)
plot(ps_inc_dec, plots = 2)
dev.off()

# Balance plot
png(
  "/Users/gaetanotedesco/Desktop/Research/Thesis/DR_loneliness/results/weighting/inc_dec/balance_plot.png",
  width = 8,
  height = 6,
  units = "in",
  res = 300
)
plot(ps_inc_dec, plots = 3)
dev.off()

# Weights plot
png(
  "/Users/gaetanotedesco/Desktop/Research/Thesis/DR_loneliness/results/weighting/inc_dec/weights_plot.png",
  width = 8,
  height = 6,
  units = "in",
  res = 300
)
plot(ps_inc_dec, plots = 6)
dev.off()

# Also display in R session
plot(ps_inc_dec, plots = 1) #convergence
plot(ps_inc_dec, plots = 2) #overlap
plot(ps_inc_dec, plots = 3) #balance
plot(ps_inc_dec, plots = 6) #weights


# Balance table
(bal_tab = bal.table(ps_inc_dec))
saveRDS(
  bal_tab,
  "/Users/gaetanotedesco/Desktop/Research/Thesis/DR_loneliness/results/weighting/inc_dec/balance_table.rds"
)

# Saving weights in df
inc_dec$ipsw = get.weights(ps_inc_dec, stop.method = "es.mean")


# INCREASE VS MIX ---------------------------------------------------------------------------------####
inc_mix = subset(d, remote_contact %in% c("increase", "mix")) %>%
  mutate(across(where(is.factor), droplevels)) %>%
  mutate(remote_contact = ifelse(remote_contact == "increase", 1, 0))


# gbm optimization
system.time({
  tuned_gbm = tune.gbm(
    x = inc_mix %>%
      select(
        -c(
          "severe_loneliness",
          "remote_contact",
        )
      ),
    y = inc_mix[, "remote_contact"],
    params.grid = gbm.grid,
    method = "es.mean",
    parallel = T
  )
})


# Save hyperparams combination
saveRDS(
  tuned_gbm$params,
  "/Users/gaetanotedesco/Desktop/Research/Thesis/DR_loneliness/results/weighting/inc_mix_params.rds"
)

set.seed(run)
ps_inc_mix = do.call(
  "ps",
  c(list(formula = ps_formula, data = inc_mix), params_ps, tuned_gbm$params)
)

# Exploring/ saving results
(weight_summary_inc_mix = summary(ps_inc_mix))
saveRDS(
  weight_summary_inc_mix,
  "/Users/gaetanotedesco/Desktop/Research/Thesis/DR_loneliness/results/weighting/inc_mix/weighting_summury.rds"
)

# Global graphical assessment
# Save plots with high resolution for academic publication
# Convergence plot
png(
  "/Users/gaetanotedesco/Desktop/Research/Thesis/DR_loneliness/results/weighting/inc_mix/convergence_plot.png",
  width = 8,
  height = 6,
  units = "in",
  res = 300
)
plot(ps_inc_mix, plots = 1)
dev.off()

# Overlap plot
png(
  "/Users/gaetanotedesco/Desktop/Research/Thesis/DR_loneliness/results/weighting/inc_mix/overlap_plot.png",
  width = 8,
  height = 6,
  units = "in",
  res = 300
)
plot(ps_inc_mix, plots = 2)
dev.off()

# Balance plot
png(
  "/Users/gaetanotedesco/Desktop/Research/Thesis/DR_loneliness/results/weighting/inc_mix/balance_plot.png",
  width = 8,
  height = 6,
  units = "in",
  res = 300
)
plot(ps_inc_mix, plots = 3)
dev.off()

# Weights plot
png(
  "/Users/gaetanotedesco/Desktop/Research/Thesis/DR_loneliness/results/weighting/inc_mix/weights_plot.png",
  width = 8,
  height = 6,
  units = "in",
  res = 300
)
plot(ps_inc_mix, plots = 6)
dev.off()

plot(ps_inc_mix, plots = 1) #convergence
plot(ps_inc_mix, plots = 2) #overlap
plot(ps_inc_mix, plots = 3) #balance
plot(ps_inc_mix, plots = 6) #weights


(bal_tab = bal.table(ps_inc_mix))
saveRDS(
  bal_tab,
  "/Users/gaetanotedesco/Desktop/Research/Thesis/DR_loneliness/results/weighting/inc_mix/balance_table.rds"
)

# Saving weights in df
inc_mix$ipsw = get.weights(ps_inc_mix, stop.method = "es.mean")


# Save data with weights
d$iptw = ifelse(
  d$remote_contact == "increase",
  1,
  ifelse(d$remote_contact == "decrease", inc_dec$ipsw, inc_mix$ipsw)
)

saveRDS(
  d,
  "/Users/gaetanotedesco/Desktop/Research/Thesis/DR_loneliness/data/data_iptw.rds"
)

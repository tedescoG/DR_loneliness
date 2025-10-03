# Setup
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("utils.R")

# Load the data
d <- readRDS("data/data_iptw.rds")
d1 = readRDS("data/data_clean.rds")

d = d %>%
  mutate(emp_status = d1$emp_status)

# INCREASE vs DECREASE -----------------------------------------------####
inc_dec <- d %>%
  filter(remote_contact %in% c("increase", "decrease")) %>%
  mutate(across(where(is.factor), droplevels)) %>%
  mutate(remote_contact = as.numeric(remote_contact == "increase"))

m1 = "1"
m2 = "baseline_lone"
m3 = "baseline_lone + baseline_depr"
m4 = "female +
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
      baseline_lone"


att_inc_dec1 = aipw_att2(
  outcome = "severe_loneliness",
  treatment = "remote_contact",
  f.out = m1,
  wgt = "iptw",
  data = inc_dec
)

att_inc_dec2 = aipw_att2(
  outcome = "severe_loneliness",
  treatment = "remote_contact",
  f.out = m2,
  wgt = "iptw",
  data = inc_dec
)


att_inc_dec3 = aipw_att2(
  outcome = "severe_loneliness",
  treatment = "remote_contact",
  f.out = m3,
  wgt = "iptw",
  data = inc_dec
)


att_inc_dec4 = aipw_att2(
  outcome = "severe_loneliness",
  treatment = "remote_contact",
  f.out = m4,
  wgt = "iptw",
  data = inc_dec
)


# INCREASE vs MIX
inc_mix <- d %>%
  filter(remote_contact %in% c("increase", "mix")) %>%
  mutate(across(where(is.factor), droplevels)) %>%
  mutate(remote_contact = as.numeric(remote_contact == "increase"))

att_inc_mix1 = aipw_att(
  outcome = "severe_loneliness",
  treatment = "remote_contact",
  f.out = m1,
  wgt = "iptw",
  data = inc_mix
)


att_inc_mix2 = aipw_att(
  outcome = "severe_loneliness",
  treatment = "remote_contact",
  f.out = m2,
  wgt = "iptw",
  data = inc_mix
)


att_inc_mix3 = aipw_att(
  outcome = "severe_loneliness",
  treatment = "remote_contact",
  f.out = m3,
  wgt = "iptw",
  data = inc_mix
)


att_inc_mix4 = aipw_att(
  outcome = "severe_loneliness",
  treatment = "remote_contact",
  f.out = m4,
  wgt = "iptw",
  data = inc_mix
)

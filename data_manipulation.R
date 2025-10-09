# DATA MANIPULATION
# Variable construction from questionnaire's items

# Setup
library(haven)
library(tidyverse)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Data Import
df = read_dta("data/pooled_ita_w.dta") %>%
  select(-c(1:6)) %>%
  as.data.frame()

# Covariate definitions ----------------------------------------------------------------------------------####

# Sex
table(df$F2, useNA = "ifany")
female = ifelse(df$F2 == 1, 1, 0)

# Age
sum(is.na(df$F1))
age_num = 2020 - as.numeric(df$F1)
age_cat = cut(
  age_num,
  breaks = c(18, 25, 35, 49, 59, 69, Inf),
  labels = c("18-25", "26-35", "36-49", "50-59", "60-69", "70+"),
  right = F,
  include.lowest = T
)
table(age_cat)

# Education
table(df$F7, useNA = "ifany")
edu = cut(
  as.numeric(df$F7),
  breaks = c(0, 3, 4, Inf),
  labels = c("low", "medium", "high")
)

# Employment Status
table(df$F8, useNA = "ifany")
emp_status = factor(
  df$F8,
  levels = 1:7,
  labels = c(
    "employed",
    "student",
    'unemployed',
    'other',
    'retired',
    'home care',
    'other'
  )
)
table(emp_status)

# Income class
table(df$F9, useNA = "ifany")
income = factor(
  df$F9,
  levels = 1:4,
  labels = c('comfortable', 'sufficient', 'struggling', 'desperate')
)
table(income)

# Marital status
table(df$F6, useNA = "ifany")
marital = factor(
  df$F6,
  levels = 1:6,
  labels = c(
    "single",
    'married',
    'domestic partnership',
    'not cohabiting with partner',
    'separated divorced',
    'widowed'
  )
)
table(marital)

# Kin available
table(df$A1_10, useNA = "ifany")
kinless = ifelse(df$A1_10 == 1, 1, 0)
table(kinless)

# Co-living
# check whether respondent lived alone before and during lockdown
coliving = ifelse((df$A2_10 == 1) & (df$A3_1 == 1), 0, 1)
table(coliving)

# Self-reported health condition pre-covid
table(df$F10, useNA = "ifany")
# aggregation of `bad` and `very bad`
df$F10[df$F10 == 5] = 4
health_pre = factor(
  df$F10,
  levels = 1:4,
  labels = c('very good', 'good', 'acceptable', 'poor')
)
table(health_pre)

# Chronic illness
table(df$F11, useNA = "ifany")
chronic = ifelse(df$F11 == 2, 0, 1)
table(chronic)

# Death of relative/friend due to COVID
table(df$E1_7, useNA = "ifany")
death_due_covid = ifelse(replace_na(df$E1_7, 0) == 1, 1, 0)

# Relative/friend infected
table(df$E1_8, useNA = "ifany")
ppl_infected = ifelse(replace_na(df$E1_8, 0) == 1, 1, 0)

# Income loss during COVID
table(df$E1_4, useNA = "ifany")
income_loss = ifelse(replace_na(df$E1_4, 0) == 1, 1, 0)

# Job loss during COVID
table(df$E1_5, useNA = "ifany")
job_loss = ifelse(replace_na(df$E1_5, 0) == 1, 1, 0)
table(job_loss, useNA = "ifany")

# Received help/support during COVID
table(df$E2_7, useNA = "ifany")
support = ifelse(df$E2_7 == 1, 0, 1)
table(support)

# Self-reported depression pre-COVID (Baseline)
table(df$C2, useNA = "ifany")
baseline_depr = factor(
  df$C2,
  levels = 1:4,
  labels = c("more often", "as usual", "less often", "never")
)
table(baseline_depr)

# Self-reported loneliness pre_COVID (Baseline)
table(df$C4, useNA = "ifany")
baseline_lone = factor(
  df$C4,
  levels = 1:4,
  labels = c("more often", "as usual", "less often", "never")
)
table(baseline_lone)

# Self-reported neighborhood description
table(df$F5, useNA = "ifany")
neighborhood = factor(
  df$F5,
  levels = 1:5,
  labels = c("big city", "suburb", "town", "village", "countryside")
)

# OUTCOME ------------------------------------------------------------------------------------------------####
table(df$C3, useNA = "ifany")
loneliness = ifelse(df$C3 %in% c(4, 3), 1, 0)
table(loneliness)


# MULTINOMIAL TREATMENT ----------------------------------------------------------------------------------####
# Questionnaire items about increase
B5 = df[, grep("^B5", names(df))]
B5 = data.frame(lapply(B5, function(x) as.numeric(x)))

# Questionnaire items about decrease
B6 = df[, grep("^B6", names(df))]
B6 = data.frame(lapply(B6, function(x) as.numeric(x)))

# function for generating the treatment
multinomial_treatment = function(df1, df2) {
  n = max(nrow(df1), nrow(df2))
  treat = rep(0, n)
  for (i in 1:n) {
    if (
      (rowSums(df1[i, 1:6], na.rm = T) >= 1) &
        ((rowSums(df2[i, 1:6], na.rm = T) == 0) | (df2[i, 7] == 1))
    ) {
      treat[i] = 3 # only increased contacts
    } else if (
      (rowSums(df1[i, 1:6], na.rm = T) >= 1) &
        ((rowSums(df2[i, 1:6], na.rm = T) >= 1))
    ) {
      treat[i] = 2 # mix pattern: increase + decrease
    } else if ((df1[i, 7] == 1) & ((rowSums(df2[i, 1:6], na.rm = T) >= 0))) {
      treat[i] = 1 # decrease: unchanged or decreased contacts
    }
  }
  return(treat)
}

remote_contact = factor(
  multinomial_treatment(B5, B6),
  levels = c(3, 2, 1),
  labels = c("increase", "mix", "decrease")
)


# Saving the final data ----------------------------------------------------------------------------------####
data_clean = data.frame(
  female,
  age_cat,
  edu,
  emp_status,
  income,
  marital,
  coliving,
  health_pre,
  kinless,
  chronic,
  death_due_covid,
  ppl_infected,
  income_loss,
  job_loss,
  support,
  neighborhood,
  baseline_depr,
  baseline_lone,
  loneliness,
  remote_contact
) %>%
  filter(age_cat %in% c("50-59", "60-69", "70+")) %>%
  filter(emp_status != "student") %>%
  filter(!rowSums(is.na(.))) %>%
  mutate(across(where(is.factor), droplevels))


# Save data in .rds
saveRDS(
  data_clean,
  "data/data_clean.rds"
)

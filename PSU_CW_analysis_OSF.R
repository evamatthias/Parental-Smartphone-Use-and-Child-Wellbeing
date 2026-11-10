################################################################################
# Analysis script: Parental smartphone use and Child Wellbeing

# Copyright (c) 2026 Eva Matthias
# Licensed under CC BY 4.0. See https://creativecommons.org/licenses/by/4.0/
# This file is part of Parental smartphone use and Child Wellbeing. 
# If you redistribute this file, you must retain this copyright notice and the license information.

# Data structure: dyadic long format (one row per person and day).
#   Suffix .a = actor , .p = partner
#   _wc = within-person centered, _bc /_bc = between-person centered
################################################################################

# Setup ####

## Set working directory
# If using RStudio:
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Otherwise, set the path manually:
# setwd("path/to/folder")

## Load required packages
library(tidyverse)   
library(psych)      
library(gtsummary)   
library(MissMech)    
library(effectsize)  
library(gt)         
library(nlme)       
library(sjPlot)     
library(emmeans)     

## Import data####

data1 <- read.csv2("PSU_CW_data_OSF.csv")

demographics <- data1 %>%
  distinct(couplecode, id_person.a, .keep_all = TRUE)%>%
  filter(!is.na(id_person.a))


# Data preparation ####

demographics <- demographics %>%
  mutate(
    # Mean daily working hours extrapolated to a working week (4.5 days)
    eworkh_week = eworkh_mean_person * 4.5
  )

# Person-level data with mothers as actors (.a = mother, .p = father)
demo_female <- demographics %>%
  filter(gender.a == 1)

# Daily data split by gender
data_female <- data1 %>% filter(gender_f.a == 1)
data_male   <- data1 %>% filter(gender_m.a == 1)


# Descriptive statistics ####

## Overall sample
describe(demographics[, c("kid_age", "subjSES.a",
                          "eworkh_mean_person", "eworkh_week")])
describe(demo_female$kid_age)   # target child age (one value per family)

table(demographics$kidtgender.a)


## Daily variables
describe(data1[, c("ephubc.A_bc", "ephubc.a_wc", "ephubc_needthreat_mean.a",
                   "pos_emo_child.a", "neg_emo_child.a", "neg_emo.a")])

## Continuous variables by gender (1 = mothers, 2 = fathers)
summarise_by_gender <- function(data, var) {
  data %>%
    group_by(gender.a) %>%
    summarise(
      Variable = var,
      n        = sum(!is.na(.data[[var]])),
      M        = mean(.data[[var]], na.rm = TRUE),
      SD       = sd(.data[[var]], na.rm = TRUE),
      Median   = median(.data[[var]], na.rm = TRUE),
      Min      = min(.data[[var]], na.rm = TRUE),
      Max      = max(.data[[var]], na.rm = TRUE),
      .groups  = "drop"
    )
}

# Working hours only for employed participants (occupation codes 5-7 excluded)
demographics_employed <- demographics %>% filter(!occu.a %in% c(5, 6, 7))

desc_by_gender <- bind_rows(
  map_dfr(c( "kid_age", "subjSES.a"),
          ~ summarise_by_gender(demographics, .x)),
  map_dfr(c("eworkh_mean_person", "eworkh_week"),
          ~ summarise_by_gender(demographics_employed, .x))
)
desc_by_gender



# Missing data analysis ####

vars_ema <- c("ephubc.a", "neg_emo.a", "neg_emo_child.a",
              "pos_emo_child.a")

## Compliance per person (% of answered prompts)
compliance_person <- data1 %>%
  group_by(id_person.a) %>%
  summarise(across(all_of(vars_ema), ~ mean(!is.na(.x)) * 100))

compliance_person %>%
  summarise(across(all_of(vars_ema),
                   list(Mean = mean, SD = sd, Min = min, Max = max),
                   .names = "{.col}_{.fn}"))

## Overall compliance and missingness
missing_summary <- tibble(
  Variable        = vars_ema,
  Valid_n         = colSums(!is.na(data1[, vars_ema])),
  Missing_n       = colSums(is.na(data1[, vars_ema])),
  Total_n         = nrow(data1),
  Compliance_pct  = round(colMeans(!is.na(data1[, vars_ema])) * 100, 2),
  Missing_pct     = round(colMeans(is.na(data1[, vars_ema])) * 100, 2)
)
missing_summary

## Need threat: missingness only on days with smartphone use (ephubc.a != 1)
data1 %>%
  filter(ephubc.a != 1) %>%
  summarise(
    n_available     = n(),
    missing_n       = sum(is.na(ephubc_needthreat_mean.a)),
    missing_percent = round(mean(is.na(ephubc_needthreat_mean.a)) * 100, 2)
  )

## Mobile sensing data (participants with sensing data only)
vars_sensing   <- c("duration_screen_h.a", "frequency_screen.a")
data1_movisens <- data1 %>% filter(!is.na(id_movisens.a))

tibble(
  Variable    = vars_sensing,
  Missing_n   = colSums(is.na(data1_movisens[, vars_sensing])),
  Total_n     = nrow(data1_movisens),
  Missing_pct = round(colMeans(is.na(data1_movisens[, vars_sensing])) * 100, 2)
)

## Test of MCAR 
TestMCARNormality(data1[, c("ephubc.a", "ephubc.p", "neg_emo.a")])


# Paired t-tests: mothers vs. fathers #####

# Variable pairs: mothers (.a) vs. fathers (.p), person level
var_pairs <- list(
  "Self reported Smartphone Use in the presence of the child" = c("ephubc.A_mean_person", "ephubc.P_mean_person"),
  "Screen duration"                = c("duration_screen.a_mean_person", "duration_screen.p_mean_person"),
  "Screen frequency"               = c("frequency_screen.a_mean_person","frequency_screen.p_mean_person"),
  "Need threat"                    = c("ephubc_needthreat_mean.a",      "ephubc_needthreat_mean.p"),
  "Positive emotions toward child" = c("pos_emo_child.a",               "pos_emo_child.p"),
  "Negative emotions toward child" = c("neg_emo_child.a",               "neg_emo_child.p")
)

# Paired t-test with Cohen's d (complete pairs only)
run_ttest <- function(data, var1, var2, label) {
  d <- data[complete.cases(data[, c(var1, var2)]), c(var1, var2)]
  t_res <- t.test(d[[var1]], d[[var2]], paired = TRUE)
  eff   <- cohens_d(d[[var1]], d[[var2]], paired = TRUE)

  tibble(
    Variable  = label,
    M_Mother  = mean(d[[var1]]),
    SD_Mother = sd(d[[var1]]),
    M_Father  = mean(d[[var2]]),
    SD_Father = sd(d[[var2]]),
    t         = unname(t_res$statistic),
    df        = unname(t_res$parameter),
    p         = ifelse(t_res$p.value < .001, "< .001", sprintf("%.3f", t_res$p.value)),
    d         = eff$Cohens_d,
    N         = nrow(d)
  )
}

ttest_table <- imap_dfr(var_pairs, ~ run_ttest(demo_female, .x[1], .x[2], .y))
ttest_table

# Correlation matrix ####
# Note. Fathers' correlations are shown above the diagonal, mothers' correlations
# below the diagonal, and mother-father intercorrelations on the diagonal.

vars_mothers <- c("ephubc.A_bc", "duration_screen.a_mean_person",
                  "frequency_screen.a_mean_person", "ephubc_needthreat_mean.a",
                  "pos_emo_child.a", "neg_emo_child.a")
vars_fathers <- c("ephubc.P_bc", "duration_screen.p_mean_person",
                  "frequency_screen.p_mean_person", "ephubc_needthreat_mean.p",
                  "pos_emo_child.p", "neg_emo_child.p")
construct_names <- c("Self-reported Smartphone Use", "Screen duration", "Screen frequency",
                     "Need threat", "Positive emotion", "Negative emotion")

k <- length(vars_mothers)
as_n_matrix <- function(n) if (is.matrix(n)) n else matrix(n, k, k)

cor_mothers <- corr.test(demo_female[, vars_mothers], use = "pairwise", adjust = "none")
cor_fathers <- corr.test(demo_female[, vars_fathers], use = "pairwise", adjust = "none")

# Fathers: upper triangle
R <- cor_fathers$r
P <- cor_fathers$p
N <- as_n_matrix(cor_fathers$n)

# Mothers: lower triangle
lower <- lower.tri(R)
R[lower] <- cor_mothers$r[lower]
P[lower] <- cor_mothers$p[lower]
N[lower] <- as_n_matrix(cor_mothers$n)[lower]

# Mother-father correlations: diagonal
dyad_cor <- map(seq_len(k), ~ cor.test(demo_female[[vars_mothers[.x]]],
                                       demo_female[[vars_fathers[.x]]]))
diag(R) <- map_dbl(dyad_cor, ~ .x$estimate)
diag(P) <- map_dbl(dyad_cor, ~ .x$p.value)
diag(N) <- map_dbl(dyad_cor, ~ .x$parameter + 2)

# Table with r, significance stars and n
stars <- ifelse(P < .001, "***", ifelse(P < .01, "**", ifelse(P < .05, "*", "")))
cor_table <- matrix(sprintf("%.2f%s (n=%d)", R, stars, as.integer(N)),
                    nrow = k, dimnames = list(construct_names, construct_names))
cor_table

# Exact p-values (same layout: fathers above, mothers below, dyad on diagonal)
dimnames(P) <- list(construct_names, construct_names)
round(P, 3)

# Table with r, exact p and n in each cell
cor_table_p <- matrix(
  sprintf("%.2f (p %s, n=%d)", R,
          ifelse(P < .001, "< .001", paste("=", sub("^0", "", sprintf("%.3f", P)))),
          as.integer(N)),
  nrow = k, dimnames = list(construct_names, construct_names))
cor_table_p


## Multilevel models with random intercept and random slopes ####

# create seperate female and male datasets####
data_male <- data1 %>%
  filter(gender_m.a==1)

data_female <- data1 %>%
  filter(gender_f.a==1)


# Models including child age as a moderator ####
# NOTE: Deviation from the preregistration. The preregistered analyses did not
# include child age. Child age was added as a moderator during the peer-review
# process.
# To reproduce the preregistered models, remove "kid_age*" from the fixed
# effects (e.g., "+ kid_age*ephubc.a_wc + kid_age*ephubc.A_bc" becomes
# "+ ephubc.a_wc + ephubc.A_bc"). All other model specifications remain unchanged.


need_threat_female <- lme(fixed = ephubc_needthreat_mean.a ~ 
                            + time 
                          + weekend
                          + eworkh.a_wc + eworkh.A_bc 
                          + kid_age*ephubc.a_wc + kid_age*ephubc.A_bc,
                          random = ~ 1  
                          #+ ephubc.a_wc 
                          | id_person.a, 
                          correlation = corAR1(form = ~ time | id_person.a), 
                          data=data_female,
                          control = lmeControl(msMaxIter = 5000, msMaxEval = 5000),
                          na.action = na.omit) 
summary(need_threat_female)
tab_model(need_threat_female, dv.labels = "Need threat female")


need_threat_male <- lme(fixed = ephubc_needthreat_mean.a ~ 
                          + time 
                        +weekend
                        + eworkh.a_wc + eworkh.A_bc 
                        + kid_age*ephubc.a_wc +  kid_age*ephubc.A_bc ,
                        random = ~ 1  
                        #+ ephubc.a_wc 
                        | id_person.a, 
                        correlation = corAR1(form = ~ time | id_person.a), 
                        data=data_male,
                        control = lmeControl(msMaxIter = 5000, msMaxEval = 5000),
                        na.action = na.omit) 
summary(need_threat_male)
tab_model(need_threat_male, dv.labels = "Need threat male")

need_threat_female_duration<- lme(fixed = ephubc_needthreat_mean.a ~ 
                                    + time 
                                  +weekend
                                  + eworkh.a_wc + eworkh.A_bc 
                                  +  kid_age*duration_screen_h.a_wc +  kid_age*duration_screen.a_bc, 
                                  random = ~ 1  
                                  + duration_screen_h.a_wc 
                                  | id_person.a, 
                                  correlation = corAR1(form = ~ time | id_person.a), 
                                  data=data_female,
                                  control = lmeControl(msMaxIter = 5000, msMaxEval = 5000),
                                  na.action = na.omit) 
summary(need_threat_female_duration)
tab_model(need_threat_female_duration, dv.labels = "Need threat female Duration")

need_threat_male_duration <- lme(fixed = ephubc_needthreat_mean.a ~ 
                                   + time 
                                 + weekend
                                 + eworkh.a_wc + eworkh.A_bc 
                                 +  kid_age*duration_screen_h.a_wc +  kid_age*duration_screen.a_bc , 
                                 random = ~ 1  
                                 # + duration_screen_h.a_wc 
                                 | id_person.a, 
                                 correlation = corAR1(form = ~ time | id_person.a), 
                                 data=data_male,
                                 control = lmeControl(msMaxIter = 5000, msMaxEval = 5000),
                                 na.action = na.omit) 
summary(need_threat_male_duration)
tab_model(need_threat_male_duration, dv.labels = "Need threat male Duration")

need_threat_female_frequency<- lme(fixed = ephubc_needthreat_mean.a ~ 
                                     + time 
                                   + weekend
                                   + eworkh.a_wc + eworkh.A_bc 
                                   +  kid_age*frequency_screen.a_wc +  kid_age*frequency_screen.a_bc , 
                                   random = ~ 1  
                                   + frequency_screen.a_wc 
                                   | id_person.a, 
                                   correlation = corAR1(form = ~ time | id_person.a), 
                                   data=data_female,
                                   control = lmeControl(msMaxIter = 5000, msMaxEval = 5000),
                                   na.action = na.omit) 
summary(need_threat_female_frequency)
tab_model(need_threat_female_frequency, dv.labels = "Need threat female frequency")

need_threat_male_frequency <- lme(fixed = ephubc_needthreat_mean.a ~ 
                                    + time
                                  + weekend
                                  + eworkh.a_wc + eworkh.A_bc 
                                  +  kid_age*frequency_screen.a_wc +  kid_age*frequency_screen.a_bc , 
                                  random = ~ 1  
                                  # + frequency_screen.a_wc 
                                  | id_person.a, 
                                  correlation = corAR1(form = ~ time | id_person.a), 
                                  data=data_male,
                                  control = lmeControl(msMaxIter = 5000, msMaxEval = 5000),
                                  na.action = na.omit) 
summary(need_threat_male_frequency)
tab_model(need_threat_male_frequency, dv.labels = "Need threat male frequency")

pos_emo_female <- lme(fixed = pos_emo_child.a ~ 
                        + time 
                      + weekend
                      + eworkh.a_wc + eworkh.A_bc 
                      +  kid_age*ephubc.a_wc +  kid_age*ephubc.A_bc ,
                      random = ~ 1  
                      #  + ephubc.a_wc
                      | id_person.a, 
                      correlation = corAR1(form = ~ time | id_person.a), 
                      data=data_female,
                      control = lmeControl(msMaxIter = 5000, msMaxEval = 5000),
                      na.action = na.omit) 
summary(pos_emo_female)
tab_model(pos_emo_female, dv.labels = "Pos Emotions Female")

pos_emo_male <- lme(fixed = pos_emo_child.a ~ 
                      + time 
                    + weekend
                    + eworkh.a_wc + eworkh.A_bc 
                    +  kid_age*ephubc.a_wc +  kid_age*ephubc.A_bc , 
                    random = ~ 1  
                    # + ephubc.a_wc
                    | id_person.a, 
                    correlation = corAR1(form = ~ time | id_person.a), 
                    data=data_male,
                    control = lmeControl(msMaxIter = 5000, msMaxEval = 5000),
                    na.action = na.omit) 
summary(pos_emo_male)
tab_model(pos_emo_male, dv.labels = "Pos Emotions male")

pos_emo_female_duration <- lme(fixed = pos_emo_child.a ~ 
                                 + time 
                               + weekend
                               #+ eworkh.a_wc + eworkh.A_bc 
                               +  kid_age*duration_screen_h.a_wc +   kid_age*duration_screen.a_bc , 
                               random = ~ 1  
                               #+duration_screen_h.a_wc 
                               | id_person.a, 
                               correlation = corAR1(form = ~ time | id_person.a), 
                               data=data_female,
                               control = lmeControl(msMaxIter = 5000, msMaxEval = 5000),
                               na.action = na.omit) 
summary(pos_emo_female_duration)
tab_model(pos_emo_female_duration, dv.labels = "Pos Emotions female Duration")

pos_emo_male_duration <- lme(fixed = pos_emo_child.a ~ 
                               + time 
                             + weekend
                             + eworkh.a_wc + eworkh.A_bc 
                             +  kid_age*duration_screen_h.a_wc +   kid_age*duration_screen.a_bc, 
                             random = ~ 1  
                             #+duration_screen_h.a_wc 
                             | id_person.a, 
                             correlation = corAR1(form = ~ time | id_person.a), 
                             data=data_male,
                             control = lmeControl(msMaxIter = 5000, msMaxEval = 5000),
                             na.action = na.omit) 
summary(pos_emo_male_duration)
tab_model(pos_emo_male_duration, dv.labels = "Pos Emotions male Duration")

pos_emo_female_frequency <- lme(fixed = pos_emo_child.a ~ 
                                  + time 
                                + weekend
                                + eworkh.a_wc + eworkh.A_bc 
                                +  kid_age*frequency_screen.a_wc +   kid_age*frequency_screen.a_bc , 
                                random = ~ 1  
                                +frequency_screen.a_wc 
                                | id_person.a, 
                                correlation = corAR1(form = ~ time | id_person.a), 
                                data=data_female,
                                control = lmeControl(msMaxIter = 5000, msMaxEval = 5000),
                                na.action = na.omit) 
summary(pos_emo_female_frequency)
tab_model(pos_emo_female_frequency, dv.labels = "Pos Emotions female frequency")

pos_emo_male_frequency <- lme(fixed = pos_emo_child.a ~ 
                                + time
                              + weekend
                              + eworkh.a_wc + eworkh.A_bc 
                              +  kid_age*frequency_screen.a_wc +   kid_age*frequency_screen.a_bc, 
                              random = ~ 1  
                              # +frequency_screen.a_wc 
                              | id_person.a, 
                              correlation = corAR1(form = ~ time | id_person.a), 
                              data=data_male,
                              control = lmeControl(msMaxIter = 5000, msMaxEval = 5000),
                              na.action = na.omit) 
summary(pos_emo_male_frequency)
tab_model(pos_emo_male_frequency, dv.labels = "Pos Emotions male frequency")

neg_emo_female <- lme(fixed = neg_emo_child.a ~ 
                        + time 
                      + weekend
                      + eworkh.a_wc + eworkh.A_bc 
                      +  kid_age*ephubc.a_wc +  kid_age*ephubc.A_bc, 
                      random = ~ 1  
                      + ephubc.a_wc
                      | id_person.a, 
                      correlation = corAR1(form = ~ time | id_person.a), 
                      data=data_female,
                      control = lmeControl(msMaxIter = 5000, msMaxEval = 5000),
                      na.action = na.omit) 
summary(neg_emo_female)
tab_model(neg_emo_female, dv.labels = "Neg Emotions Female")

neg_emo_male <- lme(fixed = neg_emo_child.a ~ 
                      + time 
                    + weekend
                    + eworkh.a_wc + eworkh.A_bc 
                    +  kid_age*ephubc.a_wc +  kid_age*ephubc.A_bc, 
                    random = ~ 1 
                    + ephubc.a_wc 
                    | id_person.a, 
                    correlation = corAR1(form = ~ time | id_person.a), 
                    data=data_male,
                    control = lmeControl(msMaxIter = 5000, msMaxEval = 5000),
                    na.action = na.omit) 
summary(neg_emo_male)
tab_model(neg_emo_male, dv.labels = "Neg Emotions male")

neg_emo_female_duration <- lme(fixed = neg_emo_child.a ~ 
                                 + time 
                               + weekend
                               + eworkh.a_wc + eworkh.A_bc 
                               +  kid_age*duration_screen_h.a_wc +   kid_age*duration_screen.a_bc, 
                               random = ~ 1 
                               # + duration_screen_h.a_wc
                               | id_person.a, 
                               correlation = corAR1(form = ~ time | id_person.a), 
                               data=data_female,
                               control = lmeControl(msMaxIter = 5000, msMaxEval = 5000),
                               na.action = na.omit) 
summary(neg_emo_female_duration)
tab_model(neg_emo_female_duration, dv.labels = "Neg Emotions Female Duration")

neg_emo_male_duration <- lme(fixed = neg_emo_child.a ~ 
                               + time 
                             + weekend
                             + eworkh.a_wc + eworkh.A_bc 
                             +  kid_age*duration_screen_h.a_wc +   kid_age*duration_screen.a_bc, 
                             random = ~ 1  
                             + duration_screen_h.a_wc
                             | id_person.a, 
                             correlation = corAR1(form = ~ time | id_person.a), 
                             data=data_male,
                             control = lmeControl(msMaxIter = 5000, msMaxEval = 5000),
                             na.action = na.omit) 
summary(neg_emo_male_duration)
tab_model(neg_emo_male_duration, dv.labels = "Neg Emotions male Duration")

neg_emo_female_frequency <- lme(fixed = neg_emo_child.a ~ 
                                  + time 
                                + weekend
                                + eworkh.a_wc + eworkh.A_bc 
                                +  kid_age*frequency_screen.a_wc +  kid_age*frequency_screen.a_bc , 
                                random = ~ 1  
                                #+frequency_screen.a_wc 
                                | id_person.a, 
                                correlation = corAR1(form = ~ time | id_person.a), 
                                data=data_female,
                                control = lmeControl(msMaxIter = 5000, msMaxEval = 5000),
                                na.action = na.omit) 
summary(neg_emo_female_frequency)
tab_model(neg_emo_female_frequency, dv.labels = "Neg Emotions Female frequency")

neg_emo_male_frequency <- lme(fixed = neg_emo_child.a ~ 
                                + time 
                              + weekend
                              + eworkh.a_wc + eworkh.A_bc 
                              +  kid_age*frequency_screen.a_wc +   kid_age*frequency_screen.a_bc, 
                              random = ~ 1  
                              # +frequency_screen.a_wc 
                              | id_person.a, 
                              correlation = corAR1(form = ~ time | id_person.a), 
                              data=data_male,
                              control = lmeControl(msMaxIter = 5000, msMaxEval = 5000),
                              na.action = na.omit) 
summary(neg_emo_male_frequency)
tab_model(neg_emo_male_frequency, dv.labels = "Neg Emotions male frequency")

#Save output models including child age moderator####

tab_model(
  need_threat_female,
  need_threat_female_frequency,
  need_threat_female_duration,
  pos_emo_female,
  pos_emo_female_frequency,
  pos_emo_female_duration,
  neg_emo_female,
  neg_emo_female_frequency,
  neg_emo_female_duration,
  dv.labels = c(
    "Need Threat SR",
    "Need Threat Frequency",
    "Need Threat Duration",
    "Positive Emotions SR",
    "Positive Emotions Frequency",
    "Positive Emotions Duration",
    "Negative Emotions SR",
    "Negative Emotions Frequency",
    "Negative Emotions Duration"
  ),
  file = "Mothers_interaction age.doc"
)

tab_model(
  need_threat_male,
  need_threat_male_frequency,
  need_threat_male_duration,
  pos_emo_male,
  pos_emo_male_frequency,
  pos_emo_male_duration,
  neg_emo_male,
  neg_emo_male_frequency,
  neg_emo_male_duration,
  dv.labels = c(
    "Need Threat SR",
    "Need Threat Frequency",
    "Need Threat Duration",
    "Positive Emotions SR",
    "Positive Emotions Frequency",
    "Positive Emotions Duration",
    "Negative Emotions SR",
    "Negative Emotions Frequency",
    "Negative Emotions Duration"
  ),
  file = "Fathers_interaction age.doc"
)

#Simple slope analysis####
# note if the CI crosses 0 than it is not significant!

m <- mean(data_female$kid_age, na.rm = TRUE)
s <- sd(data_female$kid_age, na.rm = TRUE)

#needthreat female

# Simple slopes bei MW ± 1 SD
slopes <- emtrends(
  need_threat_female,
  var = "ephubc.A_bc",
  at = list(kid_age = c(m - s, m, m + s))
)

summary(slopes)
confint(slopes)

# Simple slopes bei MW ± 1 SD
slopes <- emtrends(
  need_threat_female,
  var = "ephubc.a_wc",
  at = list(kid_age = c(m - s, m, m + s))
)

summary(slopes)
confint(slopes)

# Simple slopes bei MW ± 1 SD
slopes <- emtrends(
  need_threat_female_frequency,
  var = "frequency_screen.a_bc",
  at = list(kid_age = c(m - s, m, m + s))
)

summary(slopes)
confint(slopes)

# Simple slopes bei MW ± 1 SD
slopes <- emtrends(
  need_threat_female_frequency,
  var = "frequency_screen.a_wc",
  at = list(kid_age = c(m - s, m, m + s))
)

summary(slopes)
confint(slopes)


# Simple slopes bei MW ± 1 SD
slopes <- emtrends(
  need_threat_female_duration,
  var = "duration_screen.a_bc",
  at = list(kid_age = c(m - s, m, m + s))
)

summary(slopes)
confint(slopes)

# Simple slopes bei MW ± 1 SD
slopes <- emtrends(
  need_threat_female_duration,
  var = "duration_screen_h.a_wc",
  at = list(kid_age = c(m - s, m, m + s))
)

summary(slopes)
confint(slopes)

#pos_emo_female
# Simple slopes bei MW ± 1 SD
slopes <- emtrends(
  pos_emo_female,
  var = "ephubc.A_bc",
  at = list(kid_age = c(m - s, m, m + s))
)

summary(slopes)
confint(slopes)

# Simple slopes bei MW ± 1 SD
slopes <- emtrends(
  pos_emo_female,
  var = "ephubc.a_wc",
  at = list(kid_age = c(m - s, m, m + s))
)

summary(slopes)
confint(slopes)

# Simple slopes bei MW ± 1 SD
slopes <- emtrends(
  pos_emo_female_frequency,
  var = "frequency_screen.a_bc",
  at = list(kid_age = c(m - s, m, m + s))
)

summary(slopes)
confint(slopes)

# Simple slopes bei MW ± 1 SD
slopes <- emtrends(
  pos_emo_female_frequency,
  var = "frequency_screen.a_wc",
  at = list(kid_age = c(m - s, m, m + s))
)

summary(slopes)
confint(slopes)


# Simple slopes bei MW ± 1 SD
slopes <- emtrends(
  pos_emo_female_duration,
  var = "duration_screen.a_bc",
  at = list(kid_age = c(m - s, m, m + s))
)

summary(slopes)
confint(slopes)

# Simple slopes bei MW ± 1 SD
slopes <- emtrends(
  pos_emo_female_duration,
  var = "duration_screen_h.a_wc",
  at = list(kid_age = c(m - s, m, m + s))
)

summary(slopes)
confint(slopes)

#neg_emo_female
# Simple slopes bei MW ± 1 SD
slopes <- emtrends(
  neg_emo_female,
  var = "ephubc.A_bc",
  at = list(kid_age = c(m - s, m, m + s))
)

summary(slopes)
confint(slopes)

# Simple slopes bei MW ± 1 SD
slopes <- emtrends(
  neg_emo_female,
  var = "ephubc.a_wc",
  at = list(kid_age = c(m - s, m, m + s))
)

summary(slopes)
confint(slopes)

# Simple slopes bei MW ± 1 SD
slopes <- emtrends(
  neg_emo_female_frequency,
  var = "frequency_screen.a_bc",
  at = list(kid_age = c(m - s, m, m + s))
)

summary(slopes)
confint(slopes)

# Simple slopes bei MW ± 1 SD
slopes <- emtrends(
  neg_emo_female_frequency,
  var = "frequency_screen.a_wc",
  at = list(kid_age = c(m - s, m, m + s))
)

summary(slopes)
confint(slopes)


# Simple slopes bei MW ± 1 SD
slopes <- emtrends(
  neg_emo_female_duration,
  var = "duration_screen.a_bc",
  at = list(kid_age = c(m - s, m, m + s))
)

summary(slopes)
confint(slopes)

# Simple slopes bei MW ± 1 SD
slopes <- emtrends(
  neg_emo_female_duration,
  var = "duration_screen_h.a_wc",
  at = list(kid_age = c(m - s, m, m + s))
)

summary(slopes)
confint(slopes)

#slope analysis male####
#need threat male
m <- mean(data_male$kid_age, na.rm = TRUE)
s <- sd(data_male$kid_age, na.rm = TRUE)

# Simple slopes bei MW ± 1 SD
slopes <- emtrends(
  need_threat_male,
  var = "ephubc.A_bc",
  at = list(kid_age = c(m - s, m, m + s))
)

summary(slopes)
confint(slopes)

# Simple slopes bei MW ± 1 SD
slopes <- emtrends(
  need_threat_male,
  var = "ephubc.a_wc",
  at = list(kid_age = c(m - s, m, m + s))
)

summary(slopes)
confint(slopes)


# Simple slopes bei MW ± 1 SD
slopes <- emtrends(
  need_threat_male_frequency,
  var = "frequency_screen.a_bc",
  at = list(kid_age = c(m - s, m, m + s))
)

summary(slopes)
confint(slopes)

# Simple slopes bei MW ± 1 SD
slopes <- emtrends(
  need_threat_male_frequency,
  var = "frequency_screen.a_wc",
  at = list(kid_age = c(m - s, m, m + s))
)

summary(slopes)
confint(slopes)


# Simple slopes bei MW ± 1 SD
slopes <- emtrends(
  need_threat_male_duration,
  var = "duration_screen.a_bc",
  at = list(kid_age = c(m - s, m, m + s))
)

summary(slopes)
confint(slopes)

# Simple slopes bei MW ± 1 SD
slopes <- emtrends(
  need_threat_male_duration,
  var = "duration_screen_h.a_wc",
  at = list(kid_age = c(m - s, m, m + s))
)

summary(slopes)
confint(slopes)


#Simple slopes bei MW ± 1 SD
slopes <- emtrends(
  neg_emo_male,
  var = "ephubc.A_bc",
  at = list(kid_age = c(m - s, m, m + s))
)

summary(slopes)
confint(slopes)

# Simple slopes bei MW ± 1 SD
slopes <- emtrends(
  neg_emo_male,
  var = "ephubc.a_wc",
  at = list(kid_age = c(m - s, m, m + s))
)

summary(slopes)
confint(slopes)

# Simple slopes bei MW ± 1 SD
slopes <- emtrends(
  neg_emo_male_frequency,
  var = "frequency_screen.a_bc",
  at = list(kid_age = c(m - s, m, m + s))
)

summary(slopes)
confint(slopes)

# Simple slopes bei MW ± 1 SD
slopes <- emtrends(
  neg_emo_male_frequency,
  var = "frequency_screen.a_wc",
  at = list(kid_age = c(m - s, m, m + s))
)

summary(slopes)
confint(slopes)


# Simple slopes bei MW ± 1 SD
slopes <- emtrends(
  neg_emo_male_duration,
  var = "duration_screen.a_bc",
  at = list(kid_age = c(m - s, m, m + s))
)

summary(slopes)
confint(slopes)

# Simple slopes bei MW ± 1 SD
slopes <- emtrends(
  neg_emo_male_duration,
  var = "duration_screen_h.a_wc",
  at = list(kid_age = c(m - s, m, m + s))
)
summary(slopes)
confint(slopes)

#pos_emo_male
# Simple slopes bei MW ± 1 SD
slopes <- emtrends(
  pos_emo_male,
  var = "ephubc.A_bc",
  at = list(kid_age = c(m - s, m, m + s))
)

summary(slopes)
confint(slopes)

# Simple slopes bei MW ± 1 SD
slopes <- emtrends(
  pos_emo_male,
  var = "ephubc.a_wc",
  at = list(kid_age = c(m - s, m, m + s))
)

summary(slopes)
confint(slopes)

# Simple slopes bei MW ± 1 SD
slopes <- emtrends(
  pos_emo_male_frequency,
  var = "frequency_screen.a_bc",
  at = list(kid_age = c(m - s, m, m + s))
)

summary(slopes)
confint(slopes)

# Simple slopes bei MW ± 1 SD
slopes <- emtrends(
  pos_emo_male_frequency,
  var = "frequency_screen.a_wc",
  at = list(kid_age = c(m - s, m, m + s))
)

summary(slopes)
confint(slopes)


# Simple slopes bei MW ± 1 SD
slopes <- emtrends(
  pos_emo_male_duration,
  var = "duration_screen.a_bc",
  at = list(kid_age = c(m - s, m, m + s))
)

summary(slopes)
confint(slopes)

# Simple slopes bei MW ± 1 SD
slopes <- emtrends(
  pos_emo_male_duration,
  var = "duration_screen_h.a_wc",
  at = list(kid_age = c(m - s, m, m + s))
)

summary(slopes)
confint(slopes)

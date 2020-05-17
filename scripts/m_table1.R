#' ########################################
#' make Table 1 to display baseline (day 0) characteristics
#' ########################################
#' 
#' depends: micu_cohort
#' 

library(tidyverse)
library(gt)
library(gtsummary)
set.seed(16)

potential_cohort <- read_csv("./tabs/potential_subject_summary.csv")
micu_cohort <- read_csv("./tabs/micu_subject_summary.csv")



#' table 1 for pilot cohort
micu_cohort %>%
  filter(grepl("pilot",subject_id)) %>%
  #filter(!is.na(Sa_dom)) %>% # filter subjects with no enrollment nares Sa data
  select(AN_duration,
         Sa_detect,
         Sa_reads,
         Sa_prop,
         Sa_dom,
         r_Sa_30after_enrollment,
         age,
         gender,
         race,
         contains('admit'),
         copd,
         asthma,
         ild,
         lymphleuk,
         dm,
         chf,
         cirrhosis,
         contains('vanco_before0_7'),
         contains('metro_before0_7'),
         contains('linez_before0_7'),
         contains('dapto_before0_7'),
         contains('cefaz_before0_7'),
         contains('piptaz_before0_7'),
         contains('cefep_before0_7'),
         contains('mero_before0_7')) %>%
  tbl_summary(
    data = .,
    by = Sa_detect,
    type = list("AN_duration" ~ "continuous", "Sa_reads" ~ "continuous", "Sa_prop" ~ "continuous", "admit_fio2" ~ "continuous", "admit_peep" ~ "continuous")
  ) %>%
  add_n() %>%
  add_overall() %>%
  add_p(pvalue_fun = function(x) style_pvalue(x, digits = 2)) %>%
  gtsummary::as_gt() %>%
  gt::tab_header(title = "Table 1: Pilot Cohort", subtitle = "Primary Exposure: Detectable Anterior Nares Staphylococcus aureus") -> t1_Sa_pilot

t1_Sa_pilot

# t1_Sa_pilot %>%
#   gt::as_raw_html() %>%
#   write_lines(path = "./tabs/t1_Sa_pilot.html")



#' table 1 for validation cohort
micu_cohort %>%
  filter(grepl("validation",subject_id)) %>%
  #filter(!is.na(Sa_dom)) %>% # filter subjects with no enrollment nares Sa data
  select(AN_duration,
         Sa_detect,
         Sa_reads,
         Sa_prop,
         Sa_dom,
         r_Sa_30after_enrollment,
         age,
         gender,
         race,
         contains('admit'),
         copd,
         asthma,
         ild,
         lymphleuk,
         dm,
         chf,
         cirrhosis,
         contains('vanco_before0_7'),
         contains('metro_before0_7'),
         contains('linez_before0_7'),
         contains('dapto_before0_7'),
         contains('cefaz_before0_7'),
         contains('piptaz_before0_7'),
         contains('cefep_before0_7'),
         contains('mero_before0_7')) %>%
  tbl_summary(
    data = .,
    by = Sa_detect,
    type = list("AN_duration" ~ "continuous", "Sa_reads" ~ "continuous", "Sa_prop" ~ "continuous", "admit_fio2" ~ "continuous", "admit_peep" ~ "continuous")
  ) %>%
  add_n() %>%
  add_overall() %>%
  add_p(pvalue_fun = function(x) style_pvalue(x, digits = 2)) %>%
  gtsummary::as_gt() %>%
  gt::tab_header(title = "Table 1: Validation Cohort", subtitle = "Primary Exposure: Detectable Anterior Nares Staphylococcus aureus") -> t1_Sa_validation

t1_Sa_validation

# t1_Sa_validation %>%
#   gt::as_raw_html() %>%
#   write_lines(path = "./tabs/t1_Sa_validation.html")


#' table 1 for combined cohorts
micu_cohort %>%
  #filter(!is.na(Sa_dom)) %>% # filter subjects with no enrollment nares Sa data
  select(AN_duration,
         Sa_detect,
         Sa_reads,
         Sa_prop,
         Sa_dom,
         r_Sa_30after_enrollment,
         age,
         gender,
         race,
         contains('admit'),
         copd,
         asthma,
         ild,
         lymphleuk,
         dm,
         chf,
         cirrhosis,
         contains('vanco_before0_7'),
         contains('metro_before0_7'),
         contains('linez_before0_7'),
         contains('dapto_before0_7'),
         contains('cefaz_before0_7'),
         contains('piptaz_before0_7'),
         contains('cefep_before0_7'),
         contains('mero_before0_7')) %>%
  tbl_summary(
    data = .,
    by = Sa_detect,
    type = list("AN_duration" ~ "continuous", "Sa_reads" ~ "continuous", "Sa_prop" ~ "continuous", "admit_fio2" ~ "continuous", "admit_peep" ~ "continuous")
  ) %>%
  add_n() %>%
  add_overall() %>%
  add_p(pvalue_fun = function(x) style_pvalue(x, digits = 2)) %>%
  gtsummary::as_gt() %>%
  gt::tab_header(title = "Table 1: Combined Cohorts", subtitle = "Primary Exposure: Detectable Anterior Nares Staphylococcus aureus") -> t1_Sa_combined

t1_Sa_combined

# t1_Sa_combined %>%
#   gt::as_raw_html() %>%
#   write_lines(path = "./tabs/t1_Sa_combined_cohorts.html")






#
###
#####
###
#








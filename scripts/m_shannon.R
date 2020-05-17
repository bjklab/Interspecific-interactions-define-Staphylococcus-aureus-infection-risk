#' ########################################
#' comparison of Sa across ET and AN sites
#' ########################################
#' 
#' depends: asv_micu_complete & micu_cohort
#' 

library(tidyverse)
library(tidybayes)
library(bayesplot)
library(ggnewscale)
library(ggsci)
library(cowplot)
library(gtsummary)
set.seed(16)


micu_cohort <- read_csv("./tabs/micu_subject_summary.csv")
micu_cohort


list(read_csv("./tabs/AN_asv_summary.csv.gz",
              col_types = cols(`respiratory_cx_pathogen` = col_character(),
                               `blood_cx_pathogen` = col_character())),
     read_csv("./tabs/ET_asv_summary.csv.gz",
              col_types = cols(`respiratory_cx_pathogen` = col_character(),
                               `blood_cx_pathogen` = col_character()))) %>%
  bind_rows() -> asv_micu_complete
asv_micu_complete


asv_micu_complete %>%
  select(seqvar_id, blast_family, blast_genus, best) %>%
  distinct() %>%
  mutate(asv_id = paste0("asv_",seqvar_id)) -> asv_key
asv_key





asv_micu_complete %>%
  filter(subject_id %in% micu_cohort$subject_id) %>% # included subjects only
  mutate(fever_yn = temp_f_max >= 100.4,
         #AN_duration = as.numeric(last_AN_date - first_AN_date) + 1,
         Sa_rcx = grepl("staphylococcus_aureus", respiratory_cx_pathogen),
         Sa_rcx = replace(Sa_rcx,is.na(Sa_rcx),FALSE),
         Sa_bcx = grepl("staphylococcus_aureus", blood_cx_pathogen),
         Sa_bcx = replace(Sa_bcx,is.na(Sa_bcx),FALSE),
         Sa_fever_rcx = (Sa_rcx & fever_yn) | (Sa_rcx & lag(fever_yn,1)) | (Sa_rcx & lag(fever_yn,2)) | (Sa_rcx & lead(fever_yn,1)) | (Sa_rcx & lead(fever_yn,2)),
         Sa_fever_bcx = (Sa_bcx & fever_yn) | (Sa_bcx & lag(fever_yn,1)) | (Sa_bcx & lag(fever_yn,2)) | (Sa_bcx & lead(fever_yn,1)) | (Sa_bcx & lead(fever_yn,2)),
         Sa_infection = Sa_rcx | Sa_bcx,
         Sa_fever_infection = Sa_fever_rcx | Sa_fever_bcx) %>%
  group_by(subject_id) %>%
  mutate(Sa_prop = unique(sumstaph[subject_day == 0 & specimen_type == "AN"] / specimen_read_total[subject_day == 0 & specimen_type == "AN"]),
         Sa_reads = unique(sumstaph[subject_day == 0 & specimen_type == "AN"]),
         #AN_duration = unique(AN_duration),
         Sa_dom = Sa_prop > 1/3,
         Sa_detect = Sa_prop > 0) %>%
  ungroup() %>%
  filter(subject_day == 0 & specimen_type == "AN") %>%
  select(subject_id,specimen_shannon,Sa_detect) %>%
  distinct() %>%
  left_join(mutate(select(micu_cohort,subject_id,r_Sa_30after_enrollment),subject_id = as.character(subject_id)), by = "subject_id") -> Sa_AN0_shannon

Sa_AN0_shannon


Sa_AN0_shannon %>%
  mutate(study = stringr::str_extract(string = subject_id, pattern = "pilot|validation")) %>%
  filter(study == "pilot") %>%
  select(Sa_detect,specimen_shannon) %>%
  tbl_summary(
    data = .,
    by = Sa_detect,
    type = list("specimen_shannon" ~ "continuous")) %>%
  add_n() %>%
  add_overall() %>%
  add_p(pvalue_fun = function(x) style_pvalue(x, digits = 2)) %>%
  gtsummary::as_gt() %>%
  gt::tab_header(title = "Anterior Nares Shannon Diversity at Day 0", subtitle = "Pilot Cohort: detectable Staphylococcus aureus at enrollment?") -> shannon_pilot_exposure
shannon_pilot_exposure


Sa_AN0_shannon %>%
  mutate(study = stringr::str_extract(string = subject_id, pattern = "pilot|validation")) %>%
  filter(study == "validation") %>%
  select(Sa_detect,specimen_shannon) %>%
  tbl_summary(
    data = .,
    by = Sa_detect,
    type = list("specimen_shannon" ~ "continuous")) %>%
  add_n() %>%
  add_overall() %>%
  add_p(pvalue_fun = function(x) style_pvalue(x, digits = 2)) %>%
  gtsummary::as_gt() %>%
  gt::tab_header(title = "Anterior Nares Shannon Diversity at Day 0", subtitle = "Validation Cohort: detectable Staphylococcus aureus at enrollment?") -> shannon_validation_exposure
shannon_validation_exposure


Sa_AN0_shannon %>%
  mutate(study = stringr::str_extract(string = subject_id, pattern = "pilot|validation")) %>%
  filter(study == "pilot") %>%
  select(r_Sa_30after_enrollment,specimen_shannon) %>%
  tbl_summary(
    data = .,
    by = r_Sa_30after_enrollment,
    type = list("specimen_shannon" ~ "continuous")) %>%
  add_n() %>%
  add_overall() %>%
  add_p(pvalue_fun = function(x) style_pvalue(x, digits = 2)) %>%
  gtsummary::as_gt() %>%
  gt::tab_header(title = "Anterior Nares Shannon Diversity at Day 0", subtitle = "Pilot Cohort: detectable Staphylococcus aureus at enrollment?") -> shannon_pilot_outcome
shannon_pilot_outcome


Sa_AN0_shannon %>%
  mutate(study = stringr::str_extract(string = subject_id, pattern = "pilot|validation")) %>%
  filter(study == "validation") %>%
  select(r_Sa_30after_enrollment,specimen_shannon) %>%
  tbl_summary(
    data = .,
    by = r_Sa_30after_enrollment,
    type = list("specimen_shannon" ~ "continuous")) %>%
  add_n() %>%
  add_overall() %>%
  add_p(pvalue_fun = function(x) style_pvalue(x, digits = 2)) %>%
  gtsummary::as_gt() %>%
  gt::tab_header(title = "Anterior Nares Shannon Diversity at Day 0", subtitle = "Validation Cohort: detectable Staphylococcus aureus at enrollment?") -> shannon_validation_outcome
shannon_validation_outcome




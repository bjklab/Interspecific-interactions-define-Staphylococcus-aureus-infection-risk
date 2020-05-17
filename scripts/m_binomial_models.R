#' ########################################
#' binomial models relating Sa LRTI ~ AN microbiome (day 0)
#' ########################################
#' 
#' depends: asv_micu_complete & Sa_cohort & sa_linear_cohorts
#' 

library(tidyverse)
library(rstan)
library(tidybayes)
library(bayesplot)


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





micu_cohort %>%
  mutate(subject_id = as.character(subject_id)) %>%
  select(subject_id, r_Sa_30after_enrollment, contains("before0_7")) %>%
  rename(outcome = r_Sa_30after_enrollment) %>%
  right_join(filter(asv_micu_complete, subject_id %in% .$subject_id), by = "subject_id") %>%
  filter(specimen_type == "AN" & subject_day == 0) %>%
  #count(subject_id) #check: 90 subjects = 14 pilot + 76 validation
  distinct() -> an0_data
an0_data



#' ########################################
#' ########################################
#' SIMPLE MODELS
#' ########################################
#' ########################################




#' ########################################
#' SIMPLE binomial: Sa LRTI ~ AN Sa ASVs (day 0) - POOLED
#' ########################################


an0_data %>%
  #filter() %>% #select cohort
  filter(grepl("aureus",best)) %>%
  select(subject_id,subject_day,outcome,seqvar_id,read_count,specimen_read_total) %>%
  group_by(seqvar_id) %>%
  filter(sum(read_count, na.rm = TRUE) > 0) %>%
  ungroup() %>%
  distinct() %>%
  select(-subject_id,-subject_day) %>%
  mutate_at(.vars = vars(contains('read')), .funs = ~ as.vector(scale(.x, center=TRUE, scale = TRUE))) %>% # scale read count and specimen read total
  identity() -> d_tib
d_tib



d_tib %>%
  tidybayes::compose_data() -> d_list
str(d_list)



d_tib %>%
  pull(seqvar_id) %>%
  unique() %>%
  sort() %>%
  enframe() %>%
  rename(seqvar_num = name, seqvar_id = value) %>%
  left_join(select(asv_micu_complete,seqvar_id,best,blast_genus), by = "seqvar_id") %>%
  distinct() %>%
  mutate(best = ifelse(best == "unidentified" | best == "uncultured bacterium", blast_genus, best)) -> seqvar_key
seqvar_key



m <- stan(file = './models/binomial/staph_lrti_model-AN_Sa_ASVs_only-reghs_scaled_reads.stan',
          data = d_list,
          iter = 1000,
          chains = 8,
          cores = 8,
          verbose = TRUE,
          seed = 16,
          control=list(adapt_delta=0.95, max_treedepth=10))

rstan::check_hmc_diagnostics(m)
#[1] "n_eff / iter looks reasonable for all parameters"
#[1] "Rhat looks reasonable for all parameters"
#[1] "0 of 500 iterations ended with a divergence (0%)"
#[1] "0 of 500 iterations saturated the maximum tree depth of 10 (0%)"
#[1] "E-BFMI indicated no pathological behavior"



summary(m, probs = c(0.025, 0.975, 0.05, 0.95, 0.1, 0.9, 0.25, 0.75, 0.2, 0.8, 0.4, 0.6, 0.5))$summary %>%
  as_tibble(.,rownames = "id") %>%
  rename(parameter = id) %>%
  mutate(typeS_2.5_pos = `2.5%` > 0 & `97.5%` > 0,
         typeS_2.5_neg = `2.5%` < 0 & `97.5%` < 0,
         typeS_5_pos = `5%` > 0 & `95%` > 0,
         typeS_5_neg = `5%` < 0 & `95%` < 0,
         typeS_10_pos = `10%` > 0 & `90%` > 0,
         typeS_10_neg = `10%` < 0 & `90%` < 0,
         typeS_25_pos = `25%` > 0 & `75%` > 0,
         typeS_25_neg = `25%` < 0 & `75%` < 0) %>%
  select(parameter, n_eff, Rhat, mean, se_mean, sd, contains("typeS"), everything()) %>%
  arrange(desc(parameter))



summary(m, probs = c(0.025, 0.975, 0.05, 0.95, 0.1, 0.9, 0.25, 0.75, 0.2, 0.8, 0.4, 0.6, 0.5))$summary %>%
  as_tibble(.,rownames = "id") %>%
  rename(parameter = id) %>%
  mutate(typeS_2.5_pos = `2.5%` > 0 & `97.5%` > 0,
         typeS_2.5_neg = `2.5%` < 0 & `97.5%` < 0,
         typeS_5_pos = `5%` > 0 & `95%` > 0,
         typeS_5_neg = `5%` < 0 & `95%` < 0,
         typeS_10_pos = `10%` > 0 & `90%` > 0,
         typeS_10_neg = `10%` < 0 & `90%` < 0,
         typeS_25_pos = `25%` > 0 & `75%` > 0,
         typeS_25_neg = `25%` < 0 & `75%` < 0) %>%
  select(parameter, n_eff, Rhat, mean, se_mean, sd, contains("typeS"), everything()) %>%
  arrange(desc(parameter)) %>%
  write_csv(path = "./models/binomial/staph_lrti_model-AN_Sa_ASVs_only-reghs_scaled_reads_POOLED.csv")



summary(m, probs = c(0.025, 0.975, 0.05, 0.95, 0.1, 0.9, 0.25, 0.75, 0.2, 0.8, 0.4, 0.6, 0.5))$summary %>%
  as_tibble(.,rownames = "id") %>%
  rename(parameter = id) %>%
  mutate(typeS_2.5_pos = `2.5%` > 0 & `97.5%` > 0,
         typeS_2.5_neg = `2.5%` < 0 & `97.5%` < 0,
         typeS_5_pos = `5%` > 0 & `95%` > 0,
         typeS_5_neg = `5%` < 0 & `95%` < 0,
         typeS_10_pos = `10%` > 0 & `90%` > 0,
         typeS_10_neg = `10%` < 0 & `90%` < 0,
         typeS_25_pos = `25%` > 0 & `75%` > 0,
         typeS_25_neg = `25%` < 0 & `75%` < 0) %>%
  select(parameter, n_eff, Rhat, mean, se_mean, sd, contains("typeS"), everything()) %>%
  arrange(desc(parameter)) %>%
  filter(grepl("beta",parameter)) %>%
  ggplot(data = .) +
  geom_segment(aes(x = `2.5%`, xend = `97.5%`, y = parameter, yend = parameter)) +
  geom_point(aes(x = mean, y = parameter))



summary(m, probs = c(0.025, 0.975, 0.05, 0.95, 0.1, 0.9, 0.25, 0.75, 0.2, 0.8, 0.4, 0.6, 0.5))$summary %>%
  as_tibble(.,rownames = "id") %>%
  rename(parameter = id) %>%
  mutate(typeS_2.5_pos = `2.5%` > 0 & `97.5%` > 0,
         typeS_2.5_neg = `2.5%` < 0 & `97.5%` < 0,
         typeS_5_pos = `5%` > 0 & `95%` > 0,
         typeS_5_neg = `5%` < 0 & `95%` < 0,
         typeS_10_pos = `10%` > 0 & `90%` > 0,
         typeS_10_neg = `10%` < 0 & `90%` < 0,
         typeS_25_pos = `25%` > 0 & `75%` > 0,
         typeS_25_neg = `25%` < 0 & `75%` < 0) %>%
  select(parameter, n_eff, Rhat, mean, se_mean, sd, contains("typeS"), everything()) %>%
  arrange(desc(parameter)) %>%
  filter(grepl("beta",parameter)) %>%
  mutate(seqvar_num = readr::parse_number(parameter)) %>%
  left_join(seqvar_key, by = 'seqvar_num') %>%
  mutate(param_lab = paste(parameter,seqvar_id,sep = '\n')) %T>%
  write_csv(path = "./models/binomial/staph_lrti_model-AN_Sa_ASVs_only-reghs_scaled_reads_POOLED_seqvar.csv") %>%
  ggplot(data = .) +
  geom_segment(aes(x = `2.5%`, xend = `97.5%`, y = param_lab, yend = param_lab, colour = typeS_10_pos)) +
  geom_point(aes(x = mean, y = param_lab, colour = typeS_10_pos))




#' ########################################
#' SIMPLE binomial: Sa LRTI ~ AN Sa ASVs (day 0) - PILOT
#' ########################################


an0_data %>%
  filter(grepl("pilot",subject_id)) %>% #select cohort
  filter(grepl("aureus",best)) %>%
  select(subject_id,subject_day,outcome,seqvar_id,read_count,specimen_read_total) %>%
  group_by(seqvar_id) %>%
  filter(sum(read_count, na.rm = TRUE) > 0) %>%
  ungroup() %>%
  distinct() %>%
  select(-subject_id,-subject_day) %>%
  mutate_at(.vars = vars(contains('read')), .funs = ~ as.vector(scale(.x, center=TRUE, scale = TRUE))) %>% # scale read count and specimen read total
  identity() -> d_tib
d_tib



d_tib %>%
  tidybayes::compose_data() -> d_list
str(d_list)



d_tib %>%
  pull(seqvar_id) %>%
  unique() %>%
  sort() %>%
  enframe() %>%
  rename(seqvar_num = name, seqvar_id = value) %>%
  left_join(select(asv_micu_complete,seqvar_id,best,blast_genus), by = "seqvar_id") %>%
  distinct() %>%
  mutate(best = ifelse(best == "unidentified" | best == "uncultured bacterium", blast_genus, best)) -> seqvar_key
seqvar_key



m <- stan(file = './models/binomial/staph_lrti_model-AN_Sa_ASVs_only-reghs_scaled_reads.stan',
          data = d_list,
          iter = 1000,
          chains = 8,
          cores = 8,
          verbose = TRUE,
          seed = 16,
          control=list(adapt_delta=0.95, max_treedepth=10))

rstan::check_hmc_diagnostics(m)
#[1] "n_eff / iter looks reasonable for all parameters"
#[1] "Rhat looks reasonable for all parameters"
#[1] "0 of 500 iterations ended with a divergence (0%)"
#[1] "0 of 500 iterations saturated the maximum tree depth of 10 (0%)"
#[1] "E-BFMI indicated no pathological behavior"



summary(m, probs = c(0.025, 0.975, 0.05, 0.95, 0.1, 0.9, 0.25, 0.75, 0.2, 0.8, 0.4, 0.6, 0.5))$summary %>%
  as_tibble(.,rownames = "id") %>%
  rename(parameter = id) %>%
  mutate(typeS_2.5_pos = `2.5%` > 0 & `97.5%` > 0,
         typeS_2.5_neg = `2.5%` < 0 & `97.5%` < 0,
         typeS_5_pos = `5%` > 0 & `95%` > 0,
         typeS_5_neg = `5%` < 0 & `95%` < 0,
         typeS_10_pos = `10%` > 0 & `90%` > 0,
         typeS_10_neg = `10%` < 0 & `90%` < 0,
         typeS_25_pos = `25%` > 0 & `75%` > 0,
         typeS_25_neg = `25%` < 0 & `75%` < 0) %>%
  select(parameter, n_eff, Rhat, mean, se_mean, sd, contains("typeS"), everything()) %>%
  arrange(desc(parameter))



summary(m, probs = c(0.025, 0.975, 0.05, 0.95, 0.1, 0.9, 0.25, 0.75, 0.2, 0.8, 0.4, 0.6, 0.5))$summary %>%
  as_tibble(.,rownames = "id") %>%
  rename(parameter = id) %>%
  mutate(typeS_2.5_pos = `2.5%` > 0 & `97.5%` > 0,
         typeS_2.5_neg = `2.5%` < 0 & `97.5%` < 0,
         typeS_5_pos = `5%` > 0 & `95%` > 0,
         typeS_5_neg = `5%` < 0 & `95%` < 0,
         typeS_10_pos = `10%` > 0 & `90%` > 0,
         typeS_10_neg = `10%` < 0 & `90%` < 0,
         typeS_25_pos = `25%` > 0 & `75%` > 0,
         typeS_25_neg = `25%` < 0 & `75%` < 0) %>%
  select(parameter, n_eff, Rhat, mean, se_mean, sd, contains("typeS"), everything()) %>%
  arrange(desc(parameter)) %>%
  write_csv(path = "./models/binomial/staph_lrti_model-AN_Sa_ASVs_only-reghs_scaled_reads_PILOT.csv")



summary(m, probs = c(0.025, 0.975, 0.05, 0.95, 0.1, 0.9, 0.25, 0.75, 0.2, 0.8, 0.4, 0.6, 0.5))$summary %>%
  as_tibble(.,rownames = "id") %>%
  rename(parameter = id) %>%
  mutate(typeS_2.5_pos = `2.5%` > 0 & `97.5%` > 0,
         typeS_2.5_neg = `2.5%` < 0 & `97.5%` < 0,
         typeS_5_pos = `5%` > 0 & `95%` > 0,
         typeS_5_neg = `5%` < 0 & `95%` < 0,
         typeS_10_pos = `10%` > 0 & `90%` > 0,
         typeS_10_neg = `10%` < 0 & `90%` < 0,
         typeS_25_pos = `25%` > 0 & `75%` > 0,
         typeS_25_neg = `25%` < 0 & `75%` < 0) %>%
  select(parameter, n_eff, Rhat, mean, se_mean, sd, contains("typeS"), everything()) %>%
  arrange(desc(parameter)) %>%
  filter(grepl("beta",parameter)) %>%
  ggplot(data = .) +
  geom_segment(aes(x = `2.5%`, xend = `97.5%`, y = parameter, yend = parameter)) +
  geom_point(aes(x = mean, y = parameter))



summary(m, probs = c(0.025, 0.975, 0.05, 0.95, 0.1, 0.9, 0.25, 0.75, 0.2, 0.8, 0.4, 0.6, 0.5))$summary %>%
  as_tibble(.,rownames = "id") %>%
  rename(parameter = id) %>%
  mutate(typeS_2.5_pos = `2.5%` > 0 & `97.5%` > 0,
         typeS_2.5_neg = `2.5%` < 0 & `97.5%` < 0,
         typeS_5_pos = `5%` > 0 & `95%` > 0,
         typeS_5_neg = `5%` < 0 & `95%` < 0,
         typeS_10_pos = `10%` > 0 & `90%` > 0,
         typeS_10_neg = `10%` < 0 & `90%` < 0,
         typeS_25_pos = `25%` > 0 & `75%` > 0,
         typeS_25_neg = `25%` < 0 & `75%` < 0) %>%
  select(parameter, n_eff, Rhat, mean, se_mean, sd, contains("typeS"), everything()) %>%
  arrange(desc(parameter)) %>%
  filter(grepl("beta",parameter)) %>%
  mutate(seqvar_num = readr::parse_number(parameter)) %>%
  left_join(seqvar_key, by = 'seqvar_num') %>%
  mutate(param_lab = paste(parameter,seqvar_id,sep = '\n')) %T>%
  write_csv(path = "./models/binomial/staph_lrti_model-AN_Sa_ASVs_only-reghs_scaled_reads_PILOT_seqvar.csv") %>%
  ggplot(data = .) +
  geom_segment(aes(x = `2.5%`, xend = `97.5%`, y = param_lab, yend = param_lab, colour = typeS_10_pos)) +
  geom_point(aes(x = mean, y = param_lab, colour = typeS_10_pos))




#' ########################################
#' SIMPLE binomial: Sa LRTI ~ AN Sa ASVs (day 0) - VALIDATION
#' ########################################


an0_data %>%
  filter(grepl("validation",subject_id)) %>% #select cohort
  filter(grepl("aureus",best)) %>%
  select(subject_id,subject_day,outcome,seqvar_id,read_count,specimen_read_total) %>%
  group_by(seqvar_id) %>%
  filter(sum(read_count, na.rm = TRUE) > 0) %>%
  ungroup() %>%
  distinct() %>%
  select(-subject_id,-subject_day) %>%
  mutate_at(.vars = vars(contains('read')), .funs = ~ as.vector(scale(.x, center=TRUE, scale = TRUE))) %>% # scale read count and specimen read total
  identity() -> d_tib
d_tib



d_tib %>%
  tidybayes::compose_data() -> d_list
str(d_list)



d_tib %>%
  pull(seqvar_id) %>%
  unique() %>%
  sort() %>%
  enframe() %>%
  rename(seqvar_num = name, seqvar_id = value) %>%
  left_join(select(asv_micu_complete,seqvar_id,best,blast_genus), by = "seqvar_id") %>%
  distinct() %>%
  mutate(best = ifelse(best == "unidentified" | best == "uncultured bacterium", blast_genus, best)) -> seqvar_key
seqvar_key



m <- stan(file = './models/binomial/staph_lrti_model-AN_Sa_ASVs_only-reghs_scaled_reads.stan',
          data = d_list,
          iter = 1000,
          chains = 8,
          cores = 8,
          verbose = TRUE,
          seed = 16,
          control=list(adapt_delta=0.95, max_treedepth=10))

rstan::check_hmc_diagnostics(m)
#[1] "n_eff / iter looks reasonable for all parameters"
#[1] "Rhat looks reasonable for all parameters"
#[1] "0 of 500 iterations ended with a divergence (0%)"
#[1] "0 of 500 iterations saturated the maximum tree depth of 10 (0%)"
#[1] "E-BFMI indicated no pathological behavior"



summary(m, probs = c(0.025, 0.975, 0.05, 0.95, 0.1, 0.9, 0.25, 0.75, 0.2, 0.8, 0.4, 0.6, 0.5))$summary %>%
  as_tibble(.,rownames = "id") %>%
  rename(parameter = id) %>%
  mutate(typeS_2.5_pos = `2.5%` > 0 & `97.5%` > 0,
         typeS_2.5_neg = `2.5%` < 0 & `97.5%` < 0,
         typeS_5_pos = `5%` > 0 & `95%` > 0,
         typeS_5_neg = `5%` < 0 & `95%` < 0,
         typeS_10_pos = `10%` > 0 & `90%` > 0,
         typeS_10_neg = `10%` < 0 & `90%` < 0,
         typeS_25_pos = `25%` > 0 & `75%` > 0,
         typeS_25_neg = `25%` < 0 & `75%` < 0) %>%
  select(parameter, n_eff, Rhat, mean, se_mean, sd, contains("typeS"), everything()) %>%
  arrange(desc(parameter))



summary(m, probs = c(0.025, 0.975, 0.05, 0.95, 0.1, 0.9, 0.25, 0.75, 0.2, 0.8, 0.4, 0.6, 0.5))$summary %>%
  as_tibble(.,rownames = "id") %>%
  rename(parameter = id) %>%
  mutate(typeS_2.5_pos = `2.5%` > 0 & `97.5%` > 0,
         typeS_2.5_neg = `2.5%` < 0 & `97.5%` < 0,
         typeS_5_pos = `5%` > 0 & `95%` > 0,
         typeS_5_neg = `5%` < 0 & `95%` < 0,
         typeS_10_pos = `10%` > 0 & `90%` > 0,
         typeS_10_neg = `10%` < 0 & `90%` < 0,
         typeS_25_pos = `25%` > 0 & `75%` > 0,
         typeS_25_neg = `25%` < 0 & `75%` < 0) %>%
  select(parameter, n_eff, Rhat, mean, se_mean, sd, contains("typeS"), everything()) %>%
  arrange(desc(parameter)) %>%
  write_csv(path = "./models/binomial/staph_lrti_model-AN_Sa_ASVs_only-reghs_scaled_reads_VALIDATION.csv")



summary(m, probs = c(0.025, 0.975, 0.05, 0.95, 0.1, 0.9, 0.25, 0.75, 0.2, 0.8, 0.4, 0.6, 0.5))$summary %>%
  as_tibble(.,rownames = "id") %>%
  rename(parameter = id) %>%
  mutate(typeS_2.5_pos = `2.5%` > 0 & `97.5%` > 0,
         typeS_2.5_neg = `2.5%` < 0 & `97.5%` < 0,
         typeS_5_pos = `5%` > 0 & `95%` > 0,
         typeS_5_neg = `5%` < 0 & `95%` < 0,
         typeS_10_pos = `10%` > 0 & `90%` > 0,
         typeS_10_neg = `10%` < 0 & `90%` < 0,
         typeS_25_pos = `25%` > 0 & `75%` > 0,
         typeS_25_neg = `25%` < 0 & `75%` < 0) %>%
  select(parameter, n_eff, Rhat, mean, se_mean, sd, contains("typeS"), everything()) %>%
  arrange(desc(parameter)) %>%
  filter(grepl("beta",parameter)) %>%
  ggplot(data = .) +
  geom_segment(aes(x = `2.5%`, xend = `97.5%`, y = parameter, yend = parameter)) +
  geom_point(aes(x = mean, y = parameter))



summary(m, probs = c(0.025, 0.975, 0.05, 0.95, 0.1, 0.9, 0.25, 0.75, 0.2, 0.8, 0.4, 0.6, 0.5))$summary %>%
  as_tibble(.,rownames = "id") %>%
  rename(parameter = id) %>%
  mutate(typeS_2.5_pos = `2.5%` > 0 & `97.5%` > 0,
         typeS_2.5_neg = `2.5%` < 0 & `97.5%` < 0,
         typeS_5_pos = `5%` > 0 & `95%` > 0,
         typeS_5_neg = `5%` < 0 & `95%` < 0,
         typeS_10_pos = `10%` > 0 & `90%` > 0,
         typeS_10_neg = `10%` < 0 & `90%` < 0,
         typeS_25_pos = `25%` > 0 & `75%` > 0,
         typeS_25_neg = `25%` < 0 & `75%` < 0) %>%
  select(parameter, n_eff, Rhat, mean, se_mean, sd, contains("typeS"), everything()) %>%
  arrange(desc(parameter)) %>%
  filter(grepl("beta",parameter)) %>%
  mutate(seqvar_num = readr::parse_number(parameter)) %>%
  left_join(seqvar_key, by = 'seqvar_num') %>%
  mutate(param_lab = paste(parameter,seqvar_id,sep = '\n')) %T>%
  write_csv(path = "./models/binomial/staph_lrti_model-AN_Sa_ASVs_only-reghs_scaled_reads_VALIDATION_seqvar.csv") %>%
  ggplot(data = .) +
  geom_segment(aes(x = `2.5%`, xend = `97.5%`, y = param_lab, yend = param_lab, colour = typeS_10_pos)) +
  geom_point(aes(x = mean, y = param_lab, colour = typeS_10_pos))






#' ########################################
#' ########################################
#' INTEGRATIVE MODELS (fixed Sa + random ASV)
#' ########################################
#' ########################################



#' ########################################
#' INTEGRATIVE: LRTI ~ day 0 AN Sa (fixed) + AN microbiome (random) - POOLED
#' ########################################

lrti_binomial_Sa_cohorts %>%
  mutate(cohort = factor(cohort, levels = c('Pilot','Validation','Pooled')),
         typeS10 = typeS_10_neg | typeS_10_pos,
         typeS5 = typeS_5_neg | typeS_5_pos) %>%
  group_by(seqvar_id) %>%
  mutate(checktypeS = sum(typeS5) > 0) %>%
  ungroup() %>%
  filter(checktypeS == TRUE & is.na(seqvar_id) == FALSE) %>%
  pull(seqvar_id) %>%
  unique() -> sa_seqvars
sa_seqvars



an0_data %>%
  filter(seqvar_id %in% sa_seqvars) %>%
  mutate(seqvar_id = paste0("asv_",seqvar_id)) %>%
  select(subject_id,subject_day,seqvar_id,read_count) %>%
  distinct() %>%
  spread(key = seqvar_id, value = read_count) %>%
  right_join(filter(an0_data, !seqvar_id %in% sa_seqvars), by = c("subject_id", "subject_day")) %>%
  distinct() -> an0_data_sa_asv
an0_data_sa_asv



an0_data_sa_asv %>%
  #filter(grepl("",subject_id)) %>% #select cohort
  select(subject_id,subject_day,outcome,seqvar_id,read_count,specimen_read_total,contains("asv_")) %>%
  group_by(seqvar_id) %>%
  filter(sum(read_count, na.rm = TRUE) > 0) %>% # only include ASVs detectable in the selected specimens
  ungroup() %>%
  distinct() %>%
  select(-subject_id,-subject_day) %>%
  mutate_at(.vars = vars(contains('read'),contains("asv_")), .funs = ~ as.vector(scale(.x, center=TRUE, scale = TRUE))) %>% # scale read count and specimen read total
  identity() -> d_tib
d_tib



d_tib %>%
  tidybayes::compose_data() -> d_list
str(d_list)



d_tib %>%
  pull(seqvar_id) %>%
  unique() %>%
  sort() %>%
  enframe() %>%
  rename(seqvar_num = name, seqvar_id = value) %>%
  left_join(select(asv_micu_complete,seqvar_id,best,blast_genus), by = "seqvar_id") %>%
  distinct() %>%
  mutate(best = ifelse(best == "unidentified" | best == "uncultured bacterium", blast_genus, best)) -> seqvar_key
seqvar_key


  
m <- stan(file = './models/binomial/staph_lrti_model-AN_Sa_ASVs_microbiome-reghs_scaled_reads.stan',
          data = d_list,
          iter = 1000,
          chains = 8,
          cores = 8,
          verbose = TRUE,
          seed = 16,
          control=list(adapt_delta=0.95, max_treedepth=10))

rstan::check_hmc_diagnostics(m)
#[1] "n_eff / iter looks reasonable for all parameters"
#[1] "Rhat looks reasonable for all parameters"
#[1] "0 of 500 iterations ended with a divergence (0%)"
#[1] "0 of 500 iterations saturated the maximum tree depth of 10 (0%)"
#[1] "E-BFMI indicated no pathological behavior"



summary(m, probs = c(0.025, 0.975, 0.05, 0.95, 0.1, 0.9, 0.25, 0.75, 0.2, 0.8, 0.4, 0.6, 0.5))$summary %>%
  as_tibble(.,rownames = "id") %>%
  rename(parameter = id) %>%
  mutate(typeS_2.5_pos = `2.5%` > 0 & `97.5%` > 0,
         typeS_2.5_neg = `2.5%` < 0 & `97.5%` < 0,
         typeS_5_pos = `5%` > 0 & `95%` > 0,
         typeS_5_neg = `5%` < 0 & `95%` < 0,
         typeS_10_pos = `10%` > 0 & `90%` > 0,
         typeS_10_neg = `10%` < 0 & `90%` < 0,
         typeS_25_pos = `25%` > 0 & `75%` > 0,
         typeS_25_neg = `25%` < 0 & `75%` < 0) %>%
  select(parameter, n_eff, Rhat, mean, se_mean, sd, contains("typeS"), everything()) %>%
  arrange(desc(parameter))



summary(m, probs = c(0.025, 0.975, 0.05, 0.95, 0.1, 0.9, 0.25, 0.75, 0.2, 0.8, 0.4, 0.6, 0.5))$summary %>%
  as_tibble(.,rownames = "id") %>%
  rename(parameter = id) %>%
  #read_csv(file = "./models/binomial/staph_lrti_model-AN_Sa_ASVs_microbiome-reghs_scaled_reads_POOLED.csv") %>%
  mutate(typeS_2.5_pos = `2.5%` > 0 & `97.5%` > 0,
         typeS_2.5_neg = `2.5%` < 0 & `97.5%` < 0,
         typeS_5_pos = `5%` > 0 & `95%` > 0,
         typeS_5_neg = `5%` < 0 & `95%` < 0,
         typeS_10_pos = `10%` > 0 & `90%` > 0,
         typeS_10_neg = `10%` < 0 & `90%` < 0,
         typeS_25_pos = `25%` > 0 & `75%` > 0,
         typeS_25_neg = `25%` < 0 & `75%` < 0) %>%
  select(parameter, n_eff, Rhat, mean, se_mean, sd, contains("typeS"), everything()) %>%
  arrange(desc(parameter)) %>%
  write_csv(path = "./models/binomial/staph_lrti_model-AN_Sa_ASVs_microbiome-reghs_scaled_reads_POOLED.csv")



summary(m, probs = c(0.025, 0.975, 0.05, 0.95, 0.1, 0.9, 0.25, 0.75, 0.2, 0.8, 0.4, 0.6, 0.5))$summary %>%
  as_tibble(.,rownames = "id") %>%
  rename(parameter = id) %>%
  mutate(typeS_2.5_pos = `2.5%` > 0 & `97.5%` > 0,
         typeS_2.5_neg = `2.5%` < 0 & `97.5%` < 0,
         typeS_5_pos = `5%` > 0 & `95%` > 0,
         typeS_5_neg = `5%` < 0 & `95%` < 0,
         typeS_10_pos = `10%` > 0 & `90%` > 0,
         typeS_10_neg = `10%` < 0 & `90%` < 0,
         typeS_25_pos = `25%` > 0 & `75%` > 0,
         typeS_25_neg = `25%` < 0 & `75%` < 0) %>%
  select(parameter, n_eff, Rhat, mean, se_mean, sd, contains("typeS"), everything()) %>%
  arrange(desc(parameter)) %>%
  filter(grepl("beta",parameter)) %>%
  ggplot(data = .) +
  geom_segment(aes(x = `2.5%`, xend = `97.5%`, y = parameter, yend = parameter)) +
  geom_point(aes(x = mean, y = parameter))



summary(m, probs = c(0.025, 0.975, 0.05, 0.95, 0.1, 0.9, 0.25, 0.75, 0.2, 0.8, 0.4, 0.6, 0.5))$summary %>%
  as_tibble(.,rownames = "id") %>%
  rename(parameter = id) %>%
  #read_csv(file = "./models/binomial/staph_lrti_model-AN_Sa_ASVs_microbiome-reghs_scaled_reads_POOLED.csv") %>%
  mutate(typeS_2.5_pos = `2.5%` > 0 & `97.5%` > 0,
         typeS_2.5_neg = `2.5%` < 0 & `97.5%` < 0,
         typeS_5_pos = `5%` > 0 & `95%` > 0,
         typeS_5_neg = `5%` < 0 & `95%` < 0,
         typeS_10_pos = `10%` > 0 & `90%` > 0,
         typeS_10_neg = `10%` < 0 & `90%` < 0,
         typeS_25_pos = `25%` > 0 & `75%` > 0,
         typeS_25_neg = `25%` < 0 & `75%` < 0) %>%
  select(parameter, n_eff, Rhat, mean, se_mean, sd, contains("typeS"), everything()) %>%
  arrange(desc(parameter)) %>%
  filter(grepl("beta",parameter)) %>%
  mutate(seqvar_num = readr::parse_number(parameter)) %>%
  left_join(seqvar_key, by = 'seqvar_num') %>%
  mutate(param_lab = paste(parameter,seqvar_id,sep = '\n')) %T>%
  write_csv(path = "./models/binomial/staph_lrti_model-AN_Sa_ASVs_microbiome-reghs_scaled_reads_POOLED_seqvar.csv") %>%
  ggplot(data = .) +
  geom_segment(aes(x = `2.5%`, xend = `97.5%`, y = param_lab, yend = param_lab, colour = typeS_10_pos)) +
  geom_point(aes(x = mean, y = param_lab, colour = typeS_10_pos))



summary(m, probs = c(0.025, 0.975, 0.05, 0.95, 0.1, 0.9, 0.25, 0.75, 0.2, 0.8, 0.4, 0.6, 0.5))$summary %>%
  as_tibble(.,rownames = "id") %>%
  rename(parameter = id) %>%
  mutate(typeS_2.5_pos = `2.5%` > 0 & `97.5%` > 0,
         typeS_2.5_neg = `2.5%` < 0 & `97.5%` < 0,
         typeS_5_pos = `5%` > 0 & `95%` > 0,
         typeS_5_neg = `5%` < 0 & `95%` < 0,
         typeS_10_pos = `10%` > 0 & `90%` > 0,
         typeS_10_neg = `10%` < 0 & `90%` < 0,
         typeS_25_pos = `25%` > 0 & `75%` > 0,
         typeS_25_neg = `25%` < 0 & `75%` < 0) %>%
  select(parameter, n_eff, Rhat, mean, se_mean, sd, contains("typeS"), everything()) %>%
  arrange(desc(parameter)) %>%
  filter(grepl("beta",parameter)) %>%
  mutate(seqvar_num = readr::parse_number(parameter)) %>%
  left_join(seqvar_key, by = 'seqvar_num') %>%
  mutate(param_lab = paste(parameter,seqvar_id,sep = '\n')) %>%
  filter(typeS_10_neg | typeS_10_pos) %>%
  select(seqvar_id, blast_genus, `50%`, best) %>%
  distinct()



#' ########################################
#' INTEGRATIVE: LRTI ~ day 0 AN Sa (fixed) + AN microbiome (random) - PILOT
#' ########################################

lrti_binomial_Sa_cohorts %>%
  mutate(cohort = factor(cohort, levels = c('Pilot','Validation','Pooled')),
         typeS10 = typeS_10_neg | typeS_10_pos,
         typeS5 = typeS_5_neg | typeS_5_pos) %>%
  group_by(seqvar_id) %>%
  mutate(checktypeS = sum(typeS5) > 0) %>%
  ungroup() %>%
  filter(checktypeS == TRUE & is.na(seqvar_id) == FALSE) %>%
  pull(seqvar_id) %>%
  unique() -> sa_seqvars
sa_seqvars



an0_data %>%
  filter(seqvar_id %in% sa_seqvars) %>%
  mutate(seqvar_id = paste0("asv_",seqvar_id)) %>%
  select(subject_id,subject_day,seqvar_id,read_count) %>%
  distinct() %>%
  spread(key = seqvar_id, value = read_count) %>%
  right_join(filter(an0_data, !seqvar_id %in% sa_seqvars), by = c("subject_id", "subject_day")) %>%
  distinct() -> an0_data_sa_asv
an0_data_sa_asv



an0_data_sa_asv %>%
  filter(grepl("pilot",subject_id)) %>% #select cohort: PILOT
  select(subject_id,subject_day,outcome,seqvar_id,read_count,specimen_read_total,contains("asv_")) %>%
  group_by(seqvar_id) %>%
  filter(sum(read_count, na.rm = TRUE) > 0) %>% # only include ASVs detectable in the selected specimens
  ungroup() %>%
  distinct() %>%
  select(-subject_id,-subject_day) %>%
  mutate_at(.vars = vars(contains('read'),contains("asv_")), .funs = ~ as.vector(scale(.x, center=TRUE, scale = TRUE))) %>% # scale read count and specimen read total
  identity() -> d_tib
d_tib



d_tib %>%
  tidybayes::compose_data() -> d_list
str(d_list)



d_tib %>%
  pull(seqvar_id) %>%
  unique() %>%
  sort() %>%
  enframe() %>%
  rename(seqvar_num = name, seqvar_id = value) %>%
  left_join(select(asv_micu_complete,seqvar_id,best,blast_genus), by = "seqvar_id") %>%
  distinct() %>%
  mutate(best = ifelse(best == "unidentified" | best == "uncultured bacterium", blast_genus, best)) -> seqvar_key
seqvar_key



m <- stan(file = './models/binomial/staph_lrti_model-AN_Sa_ASVs_microbiome-reghs_scaled_reads.stan',
          data = d_list,
          iter = 1000,
          chains = 8,
          cores = 8,
          verbose = TRUE,
          seed = 16,
          control=list(adapt_delta=0.95, max_treedepth=10))

rstan::check_hmc_diagnostics(m)



summary(m, probs = c(0.025, 0.975, 0.05, 0.95, 0.1, 0.9, 0.25, 0.75, 0.2, 0.8, 0.4, 0.6, 0.5))$summary %>%
  as_tibble(.,rownames = "id") %>%
  rename(parameter = id) %>%
  mutate(typeS_2.5_pos = `2.5%` > 0 & `97.5%` > 0,
         typeS_2.5_neg = `2.5%` < 0 & `97.5%` < 0,
         typeS_5_pos = `5%` > 0 & `95%` > 0,
         typeS_5_neg = `5%` < 0 & `95%` < 0,
         typeS_10_pos = `10%` > 0 & `90%` > 0,
         typeS_10_neg = `10%` < 0 & `90%` < 0,
         typeS_25_pos = `25%` > 0 & `75%` > 0,
         typeS_25_neg = `25%` < 0 & `75%` < 0) %>%
  select(parameter, n_eff, Rhat, mean, se_mean, sd, contains("typeS"), everything()) %>%
  arrange(desc(parameter))



summary(m, probs = c(0.025, 0.975, 0.05, 0.95, 0.1, 0.9, 0.25, 0.75, 0.2, 0.8, 0.4, 0.6, 0.5))$summary %>%
  as_tibble(.,rownames = "id") %>%
  rename(parameter = id) %>%
  #read_csv(file = "./models/binomial/staph_lrti_model-AN_Sa_ASVs_microbiome-reghs_scaled_reads_PILOT.csv") %>%
  mutate(typeS_2.5_pos = `2.5%` > 0 & `97.5%` > 0,
         typeS_2.5_neg = `2.5%` < 0 & `97.5%` < 0,
         typeS_5_pos = `5%` > 0 & `95%` > 0,
         typeS_5_neg = `5%` < 0 & `95%` < 0,
         typeS_10_pos = `10%` > 0 & `90%` > 0,
         typeS_10_neg = `10%` < 0 & `90%` < 0,
         typeS_25_pos = `25%` > 0 & `75%` > 0,
         typeS_25_neg = `25%` < 0 & `75%` < 0) %>%
  select(parameter, n_eff, Rhat, mean, se_mean, sd, contains("typeS"), everything()) %>%
  arrange(desc(parameter)) %>%
  write_csv(path = "./models/binomial/staph_lrti_model-AN_Sa_ASVs_microbiome-reghs_scaled_reads_PILOT.csv")



summary(m, probs = c(0.025, 0.975, 0.05, 0.95, 0.1, 0.9, 0.25, 0.75, 0.2, 0.8, 0.4, 0.6, 0.5))$summary %>%
  as_tibble(.,rownames = "id") %>%
  rename(parameter = id) %>%
  mutate(typeS_2.5_pos = `2.5%` > 0 & `97.5%` > 0,
         typeS_2.5_neg = `2.5%` < 0 & `97.5%` < 0,
         typeS_5_pos = `5%` > 0 & `95%` > 0,
         typeS_5_neg = `5%` < 0 & `95%` < 0,
         typeS_10_pos = `10%` > 0 & `90%` > 0,
         typeS_10_neg = `10%` < 0 & `90%` < 0,
         typeS_25_pos = `25%` > 0 & `75%` > 0,
         typeS_25_neg = `25%` < 0 & `75%` < 0) %>%
  select(parameter, n_eff, Rhat, mean, se_mean, sd, contains("typeS"), everything()) %>%
  arrange(desc(parameter)) %>%
  filter(grepl("beta",parameter)) %>%
  ggplot(data = .) +
  geom_segment(aes(x = `2.5%`, xend = `97.5%`, y = parameter, yend = parameter)) +
  geom_point(aes(x = mean, y = parameter))



summary(m, probs = c(0.025, 0.975, 0.05, 0.95, 0.1, 0.9, 0.25, 0.75, 0.2, 0.8, 0.4, 0.6, 0.5))$summary %>%
  as_tibble(.,rownames = "id") %>%
  rename(parameter = id) %>%
  #read_csv(file = "./models/binomial/staph_lrti_model-AN_Sa_ASVs_microbiome-reghs_scaled_reads_PILOT.csv") %>%
  mutate(typeS_2.5_pos = `2.5%` > 0 & `97.5%` > 0,
         typeS_2.5_neg = `2.5%` < 0 & `97.5%` < 0,
         typeS_5_pos = `5%` > 0 & `95%` > 0,
         typeS_5_neg = `5%` < 0 & `95%` < 0,
         typeS_10_pos = `10%` > 0 & `90%` > 0,
         typeS_10_neg = `10%` < 0 & `90%` < 0,
         typeS_25_pos = `25%` > 0 & `75%` > 0,
         typeS_25_neg = `25%` < 0 & `75%` < 0) %>%
  select(parameter, n_eff, Rhat, mean, se_mean, sd, contains("typeS"), everything()) %>%
  arrange(desc(parameter)) %>%
  filter(grepl("beta",parameter)) %>%
  mutate(seqvar_num = readr::parse_number(parameter)) %>%
  left_join(seqvar_key, by = 'seqvar_num') %>%
  mutate(param_lab = paste(parameter,seqvar_id,sep = '\n')) %T>%
  write_csv(path = "./models/binomial/staph_lrti_model-AN_Sa_ASVs_microbiome-reghs_scaled_reads_PILOT_seqvar.csv") %>%
  ggplot(data = .) +
  geom_segment(aes(x = `2.5%`, xend = `97.5%`, y = param_lab, yend = param_lab, colour = typeS_10_pos)) +
  geom_point(aes(x = mean, y = param_lab, colour = typeS_10_pos))



summary(m, probs = c(0.025, 0.975, 0.05, 0.95, 0.1, 0.9, 0.25, 0.75, 0.2, 0.8, 0.4, 0.6, 0.5))$summary %>%
  as_tibble(.,rownames = "id") %>%
  rename(parameter = id) %>%
  mutate(typeS_2.5_pos = `2.5%` > 0 & `97.5%` > 0,
         typeS_2.5_neg = `2.5%` < 0 & `97.5%` < 0,
         typeS_5_pos = `5%` > 0 & `95%` > 0,
         typeS_5_neg = `5%` < 0 & `95%` < 0,
         typeS_10_pos = `10%` > 0 & `90%` > 0,
         typeS_10_neg = `10%` < 0 & `90%` < 0,
         typeS_25_pos = `25%` > 0 & `75%` > 0,
         typeS_25_neg = `25%` < 0 & `75%` < 0) %>%
  select(parameter, n_eff, Rhat, mean, se_mean, sd, contains("typeS"), everything()) %>%
  arrange(desc(parameter)) %>%
  filter(grepl("beta",parameter)) %>%
  mutate(seqvar_num = readr::parse_number(parameter)) %>%
  left_join(seqvar_key, by = 'seqvar_num') %>%
  mutate(param_lab = paste(parameter,seqvar_id,sep = '\n')) %>%
  filter(typeS_10_neg | typeS_10_pos) %>%
  select(seqvar_id, blast_genus, `50%`, best) %>%
  distinct()





#' ########################################
#' INTEGRATIVE: LRTI ~ day 0 AN Sa (fixed) + AN microbiome (random) - VALIDATION
#' ########################################

lrti_binomial_Sa_cohorts %>%
  mutate(cohort = factor(cohort, levels = c('Pilot','Validation','Pooled')),
         typeS10 = typeS_10_neg | typeS_10_pos,
         typeS5 = typeS_5_neg | typeS_5_pos) %>%
  group_by(seqvar_id) %>%
  mutate(checktypeS = sum(typeS5) > 0) %>%
  ungroup() %>%
  filter(checktypeS == TRUE & is.na(seqvar_id) == FALSE) %>%
  pull(seqvar_id) %>%
  unique() -> sa_seqvars
sa_seqvars



an0_data %>%
  filter(seqvar_id %in% sa_seqvars) %>%
  mutate(seqvar_id = paste0("asv_",seqvar_id)) %>%
  select(subject_id,subject_day,seqvar_id,read_count) %>%
  distinct() %>%
  spread(key = seqvar_id, value = read_count) %>%
  right_join(filter(an0_data, !seqvar_id %in% sa_seqvars), by = c("subject_id", "subject_day")) %>%
  distinct() -> an0_data_sa_asv
an0_data_sa_asv



an0_data_sa_asv %>%
  filter(grepl("validation",subject_id)) %>% #select cohort: VALIDATION
  select(subject_id,subject_day,outcome,seqvar_id,read_count,specimen_read_total,contains("asv_")) %>%
  group_by(seqvar_id) %>%
  filter(sum(read_count, na.rm = TRUE) > 0) %>% # only include ASVs detectable in the selected specimens
  ungroup() %>%
  distinct() %>%
  select(-subject_id,-subject_day) %>%
  mutate_at(.vars = vars(contains('read'),contains("asv_")), .funs = ~ as.vector(scale(.x, center=TRUE, scale = TRUE))) %>% # scale read count and specimen read total
  identity() -> d_tib
d_tib



d_tib %>%
  select(-asv_241b2b5f1ea0f1065ca7027534068198) %>% # validation cohort does not contain asv_241b2b5f1ea0f1065ca7027534068198
  distinct() -> d_tib
d_tib


d_tib %>%
  tidybayes::compose_data() -> d_list
str(d_list)



d_tib %>%
  pull(seqvar_id) %>%
  unique() %>%
  sort() %>%
  enframe() %>%
  rename(seqvar_num = name, seqvar_id = value) %>%
  left_join(select(asv_micu_complete,seqvar_id,best,blast_genus), by = "seqvar_id") %>%
  distinct() %>%
  mutate(best = ifelse(best == "unidentified" | best == "uncultured bacterium", blast_genus, best)) -> seqvar_key
seqvar_key



m <- stan(file = './models/binomial/staph_lrti_model-AN_Sa_ASVs_microbiome-reghs_scaled_reads_validation.stan', # note: model code altered given missing Sa ASV
          data = d_list,
          iter = 1000,
          chains = 8,
          cores = 8,
          verbose = TRUE,
          seed = 16,
          control=list(adapt_delta=0.95, max_treedepth=10))

rstan::check_hmc_diagnostics(m)



summary(m, probs = c(0.025, 0.975, 0.05, 0.95, 0.1, 0.9, 0.25, 0.75, 0.2, 0.8, 0.4, 0.6, 0.5))$summary %>%
  as_tibble(.,rownames = "id") %>%
  rename(parameter = id) %>%
  mutate(typeS_2.5_pos = `2.5%` > 0 & `97.5%` > 0,
         typeS_2.5_neg = `2.5%` < 0 & `97.5%` < 0,
         typeS_5_pos = `5%` > 0 & `95%` > 0,
         typeS_5_neg = `5%` < 0 & `95%` < 0,
         typeS_10_pos = `10%` > 0 & `90%` > 0,
         typeS_10_neg = `10%` < 0 & `90%` < 0,
         typeS_25_pos = `25%` > 0 & `75%` > 0,
         typeS_25_neg = `25%` < 0 & `75%` < 0) %>%
  select(parameter, n_eff, Rhat, mean, se_mean, sd, contains("typeS"), everything()) %>%
  arrange(desc(parameter))



summary(m, probs = c(0.025, 0.975, 0.05, 0.95, 0.1, 0.9, 0.25, 0.75, 0.2, 0.8, 0.4, 0.6, 0.5))$summary %>%
  as_tibble(.,rownames = "id") %>%
  rename(parameter = id) %>%
  #read_csv(file = "./models/binomial/staph_lrti_model-AN_Sa_ASVs_microbiome-reghs_scaled_reads_VALIDATION.csv") %>%
  mutate(typeS_2.5_pos = `2.5%` > 0 & `97.5%` > 0,
         typeS_2.5_neg = `2.5%` < 0 & `97.5%` < 0,
         typeS_5_pos = `5%` > 0 & `95%` > 0,
         typeS_5_neg = `5%` < 0 & `95%` < 0,
         typeS_10_pos = `10%` > 0 & `90%` > 0,
         typeS_10_neg = `10%` < 0 & `90%` < 0,
         typeS_25_pos = `25%` > 0 & `75%` > 0,
         typeS_25_neg = `25%` < 0 & `75%` < 0) %>%
  select(parameter, n_eff, Rhat, mean, se_mean, sd, contains("typeS"), everything()) %>%
  arrange(desc(parameter)) %>%
  write_csv(path = "./models/binomial/staph_lrti_model-AN_Sa_ASVs_microbiome-reghs_scaled_reads_VALIDATION.csv")



summary(m, probs = c(0.025, 0.975, 0.05, 0.95, 0.1, 0.9, 0.25, 0.75, 0.2, 0.8, 0.4, 0.6, 0.5))$summary %>%
  as_tibble(.,rownames = "id") %>%
  rename(parameter = id) %>%
  mutate(typeS_2.5_pos = `2.5%` > 0 & `97.5%` > 0,
         typeS_2.5_neg = `2.5%` < 0 & `97.5%` < 0,
         typeS_5_pos = `5%` > 0 & `95%` > 0,
         typeS_5_neg = `5%` < 0 & `95%` < 0,
         typeS_10_pos = `10%` > 0 & `90%` > 0,
         typeS_10_neg = `10%` < 0 & `90%` < 0,
         typeS_25_pos = `25%` > 0 & `75%` > 0,
         typeS_25_neg = `25%` < 0 & `75%` < 0) %>%
  select(parameter, n_eff, Rhat, mean, se_mean, sd, contains("typeS"), everything()) %>%
  arrange(desc(parameter)) %>%
  filter(grepl("beta",parameter)) %>%
  ggplot(data = .) +
  geom_segment(aes(x = `2.5%`, xend = `97.5%`, y = parameter, yend = parameter)) +
  geom_point(aes(x = mean, y = parameter))



summary(m, probs = c(0.025, 0.975, 0.05, 0.95, 0.1, 0.9, 0.25, 0.75, 0.2, 0.8, 0.4, 0.6, 0.5))$summary %>%
  as_tibble(.,rownames = "id") %>%
  rename(parameter = id) %>%
  #read_csv(file = "./models/binomial/staph_lrti_model-AN_Sa_ASVs_microbiome-reghs_scaled_reads_VALIDATION.csv") %>%
  mutate(typeS_2.5_pos = `2.5%` > 0 & `97.5%` > 0,
         typeS_2.5_neg = `2.5%` < 0 & `97.5%` < 0,
         typeS_5_pos = `5%` > 0 & `95%` > 0,
         typeS_5_neg = `5%` < 0 & `95%` < 0,
         typeS_10_pos = `10%` > 0 & `90%` > 0,
         typeS_10_neg = `10%` < 0 & `90%` < 0,
         typeS_25_pos = `25%` > 0 & `75%` > 0,
         typeS_25_neg = `25%` < 0 & `75%` < 0) %>%
  select(parameter, n_eff, Rhat, mean, se_mean, sd, contains("typeS"), everything()) %>%
  arrange(desc(parameter)) %>%
  filter(grepl("beta",parameter)) %>%
  mutate(seqvar_num = readr::parse_number(parameter)) %>%
  left_join(seqvar_key, by = 'seqvar_num') %>%
  mutate(param_lab = paste(parameter,seqvar_id,sep = '\n')) %T>%
  write_csv(path = "./models/binomial/staph_lrti_model-AN_Sa_ASVs_microbiome-reghs_scaled_reads_VALIDATION_seqvar.csv") %>%
  ggplot(data = .) +
  geom_segment(aes(x = `2.5%`, xend = `97.5%`, y = param_lab, yend = param_lab, colour = typeS_10_pos)) +
  geom_point(aes(x = mean, y = param_lab, colour = typeS_10_pos))



summary(m, probs = c(0.025, 0.975, 0.05, 0.95, 0.1, 0.9, 0.25, 0.75, 0.2, 0.8, 0.4, 0.6, 0.5))$summary %>%
  as_tibble(.,rownames = "id") %>%
  rename(parameter = id) %>%
  mutate(typeS_2.5_pos = `2.5%` > 0 & `97.5%` > 0,
         typeS_2.5_neg = `2.5%` < 0 & `97.5%` < 0,
         typeS_5_pos = `5%` > 0 & `95%` > 0,
         typeS_5_neg = `5%` < 0 & `95%` < 0,
         typeS_10_pos = `10%` > 0 & `90%` > 0,
         typeS_10_neg = `10%` < 0 & `90%` < 0,
         typeS_25_pos = `25%` > 0 & `75%` > 0,
         typeS_25_neg = `25%` < 0 & `75%` < 0) %>%
  select(parameter, n_eff, Rhat, mean, se_mean, sd, contains("typeS"), everything()) %>%
  arrange(desc(parameter)) %>%
  filter(grepl("beta",parameter)) %>%
  mutate(seqvar_num = readr::parse_number(parameter)) %>%
  left_join(seqvar_key, by = 'seqvar_num') %>%
  mutate(param_lab = paste(parameter,seqvar_id,sep = '\n')) %>%
  filter(typeS_10_neg | typeS_10_pos) %>%
  select(seqvar_id, blast_genus, `50%`, best) %>%
  distinct()













#' ########################################
#' ########################################
#' COMPLETE MODELS (fixed Sa + random ASV)
#' ########################################
#' ########################################



#' ########################################
#' COMPLETE: LRTI ~ day 0 AN Sa (fixed) + AN microbiome (random) - POOLED
#' ########################################

lrti_binomial_Sa_cohorts %>%
  mutate(cohort = factor(cohort, levels = c('Pilot','Validation','Pooled')),
         typeS10 = typeS_10_neg | typeS_10_pos,
         typeS5 = typeS_5_neg | typeS_5_pos) %>%
  group_by(seqvar_id) %>%
  mutate(checktypeS = sum(typeS5) > 0) %>%
  ungroup() %>%
  filter(checktypeS == TRUE & is.na(seqvar_id) == FALSE) %>%
  pull(seqvar_id) %>%
  unique() -> sa_seqvars
sa_seqvars



an0_data %>%
  filter(seqvar_id %in% sa_seqvars) %>%
  mutate(seqvar_id = paste0("asv_",seqvar_id)) %>%
  select(subject_id,subject_day,seqvar_id,read_count) %>%
  distinct() %>%
  spread(key = seqvar_id, value = read_count) %>%
  right_join(filter(an0_data, !seqvar_id %in% sa_seqvars), by = c("subject_id", "subject_day")) %>%
  distinct() -> an0_data_sa_asv
an0_data_sa_asv



an0_data_sa_asv %>%
  #filter(grepl("",subject_id)) %>% #select cohort
  select(subject_id,subject_day,outcome,seqvar_id,read_count,specimen_read_total,contains("asv_"),contains("before0_7")) %>%
  rename_all(.funs = ~ gsub("_before0_7","",.x)) %>%
  group_by(seqvar_id) %>%
  filter(sum(read_count, na.rm = TRUE) > 0) %>% # only include ASVs detectable in the selected specimens
  ungroup() %>%
  distinct() %>%
  select(-subject_id,-subject_day) %>%
  mutate_at(.vars = vars(contains('read'),contains("asv_")), .funs = ~ as.vector(scale(.x, center=TRUE, scale = TRUE))) %>% # scale read count and specimen read total
  identity() -> d_tib
d_tib



d_tib %>%
  tidybayes::compose_data() -> d_list
str(d_list)



d_tib %>%
  pull(seqvar_id) %>%
  unique() %>%
  sort() %>%
  enframe() %>%
  rename(seqvar_num = name, seqvar_id = value) %>%
  left_join(select(asv_micu_complete,seqvar_id,best,blast_genus), by = "seqvar_id") %>%
  distinct() %>%
  mutate(best = ifelse(best == "unidentified" | best == "uncultured bacterium", blast_genus, best)) -> seqvar_key
seqvar_key



m <- stan(file = './models/binomial/staph_lrti_model-AN_Sa_ASVs_microbiome_abx-reghs_scaled_reads.stan',
          data = d_list,
          iter = 1000,
          chains = 8,
          cores = 8,
          verbose = TRUE,
          seed = 16,
          control=list(adapt_delta=0.95, max_treedepth=10))

rstan::check_hmc_diagnostics(m)



summary(m, probs = c(0.025, 0.975, 0.05, 0.95, 0.1, 0.9, 0.25, 0.75, 0.2, 0.8, 0.4, 0.6, 0.5))$summary %>%
  as_tibble(.,rownames = "id") %>%
  rename(parameter = id) %>%
  mutate(typeS_2.5_pos = `2.5%` > 0 & `97.5%` > 0,
         typeS_2.5_neg = `2.5%` < 0 & `97.5%` < 0,
         typeS_5_pos = `5%` > 0 & `95%` > 0,
         typeS_5_neg = `5%` < 0 & `95%` < 0,
         typeS_10_pos = `10%` > 0 & `90%` > 0,
         typeS_10_neg = `10%` < 0 & `90%` < 0) %>%
  select(parameter, n_eff, Rhat, mean, se_mean, sd, contains("typeS"), everything()) %>%
  arrange(desc(parameter))



summary(m, probs = c(0.025, 0.975, 0.05, 0.95, 0.1, 0.9, 0.25, 0.75, 0.2, 0.8, 0.4, 0.6, 0.5))$summary %>%
  as_tibble(.,rownames = "id") %>%
  rename(parameter = id) %>%
  #read_csv(file = "./models/binomial/staph_lrti_model-AN_Sa_ASVs_microbiome_abx-reghs_scaled_reads_POOLED.csv") %>%
  mutate(typeS_2.5_pos = `2.5%` > 0 & `97.5%` > 0,
         typeS_2.5_neg = `2.5%` < 0 & `97.5%` < 0,
         typeS_5_pos = `5%` > 0 & `95%` > 0,
         typeS_5_neg = `5%` < 0 & `95%` < 0,
         typeS_10_pos = `10%` > 0 & `90%` > 0,
         typeS_10_neg = `10%` < 0 & `90%` < 0,
         typeS_25_pos = `25%` > 0 & `75%` > 0,
         typeS_25_neg = `25%` < 0 & `75%` < 0) %>%
  select(parameter, n_eff, Rhat, mean, se_mean, sd, contains("typeS"), everything()) %>%
  arrange(desc(parameter)) %>%
  write_csv(path = "./models/binomial/staph_lrti_model-AN_Sa_ASVs_microbiome_abx-reghs_scaled_reads_POOLED.csv")



summary(m, probs = c(0.025, 0.975, 0.05, 0.95, 0.1, 0.9, 0.25, 0.75, 0.2, 0.8, 0.4, 0.6, 0.5))$summary %>%
  as_tibble(.,rownames = "id") %>%
  rename(parameter = id) %>%
  mutate(typeS_2.5_pos = `2.5%` > 0 & `97.5%` > 0,
         typeS_2.5_neg = `2.5%` < 0 & `97.5%` < 0,
         typeS_5_pos = `5%` > 0 & `95%` > 0,
         typeS_5_neg = `5%` < 0 & `95%` < 0,
         typeS_10_pos = `10%` > 0 & `90%` > 0,
         typeS_10_neg = `10%` < 0 & `90%` < 0,
         typeS_25_pos = `25%` > 0 & `75%` > 0,
         typeS_25_neg = `25%` < 0 & `75%` < 0) %>%
  select(parameter, n_eff, Rhat, mean, se_mean, sd, contains("typeS"), everything()) %>%
  arrange(desc(parameter)) %>%
  filter(grepl("beta",parameter)) %>%
  ggplot(data = .) +
  geom_segment(aes(x = `2.5%`, xend = `97.5%`, y = parameter, yend = parameter)) +
  geom_point(aes(x = mean, y = parameter))



summary(m, probs = c(0.025, 0.975, 0.05, 0.95, 0.1, 0.9, 0.25, 0.75, 0.2, 0.8, 0.4, 0.6, 0.5))$summary %>%
  as_tibble(.,rownames = "id") %>%
  rename(parameter = id) %>%
  #read_csv(file = "./models/binomial/staph_lrti_model-AN_Sa_ASVs_microbiome_abx-reghs_scaled_reads_POOLED.csv") %>%
  mutate(typeS_2.5_pos = `2.5%` > 0 & `97.5%` > 0,
         typeS_2.5_neg = `2.5%` < 0 & `97.5%` < 0,
         typeS_5_pos = `5%` > 0 & `95%` > 0,
         typeS_5_neg = `5%` < 0 & `95%` < 0,
         typeS_10_pos = `10%` > 0 & `90%` > 0,
         typeS_10_neg = `10%` < 0 & `90%` < 0,
         typeS_25_pos = `25%` > 0 & `75%` > 0,
         typeS_25_neg = `25%` < 0 & `75%` < 0) %>%
  select(parameter, n_eff, Rhat, mean, se_mean, sd, contains("typeS"), everything()) %>%
  arrange(desc(parameter)) %>%
  filter(grepl("beta",parameter)) %>%
  mutate(seqvar_num = readr::parse_number(parameter)) %>%
  left_join(seqvar_key, by = 'seqvar_num') %>%
  mutate(param_lab = paste(parameter,seqvar_id,sep = '\n')) %T>%
  write_csv(path = "./models/binomial/staph_lrti_model-AN_Sa_ASVs_microbiome_abx-reghs_scaled_reads_POOLED_seqvar.csv") %>%
  ggplot(data = .) +
  geom_segment(aes(x = `2.5%`, xend = `97.5%`, y = param_lab, yend = param_lab, colour = typeS_10_pos)) +
  geom_point(aes(x = mean, y = param_lab, colour = typeS_10_pos))



summary(m, probs = c(0.025, 0.975, 0.05, 0.95, 0.1, 0.9, 0.25, 0.75, 0.2, 0.8, 0.4, 0.6, 0.5))$summary %>%
  as_tibble(.,rownames = "id") %>%
  rename(parameter = id) %>%
  mutate(typeS_2.5_pos = `2.5%` > 0 & `97.5%` > 0,
         typeS_2.5_neg = `2.5%` < 0 & `97.5%` < 0,
         typeS_5_pos = `5%` > 0 & `95%` > 0,
         typeS_5_neg = `5%` < 0 & `95%` < 0,
         typeS_10_pos = `10%` > 0 & `90%` > 0,
         typeS_10_neg = `10%` < 0 & `90%` < 0,
         typeS_25_pos = `25%` > 0 & `75%` > 0,
         typeS_25_neg = `25%` < 0 & `75%` < 0) %>%
  select(parameter, n_eff, Rhat, mean, se_mean, sd, contains("typeS"), everything()) %>%
  arrange(desc(parameter)) %>%
  filter(grepl("beta",parameter)) %>%
  mutate(seqvar_num = readr::parse_number(parameter)) %>%
  left_join(seqvar_key, by = 'seqvar_num') %>%
  mutate(param_lab = paste(parameter,seqvar_id,sep = '\n')) %>%
  filter(typeS_10_neg | typeS_10_pos) %>%
  select(seqvar_id, blast_genus, `50%`, best) %>%
  distinct()



#' ########################################
#' COMPLETE: LRTI ~ day 0 AN Sa (fixed) + AN microbiome (random) - PILOT
#' ########################################

lrti_binomial_Sa_cohorts %>%
  mutate(cohort = factor(cohort, levels = c('Pilot','Validation','Pooled')),
         typeS10 = typeS_10_neg | typeS_10_pos,
         typeS5 = typeS_5_neg | typeS_5_pos) %>%
  group_by(seqvar_id) %>%
  mutate(checktypeS = sum(typeS5) > 0) %>%
  ungroup() %>%
  filter(checktypeS == TRUE & is.na(seqvar_id) == FALSE) %>%
  pull(seqvar_id) %>%
  unique() -> sa_seqvars
sa_seqvars



an0_data %>%
  filter(seqvar_id %in% sa_seqvars) %>%
  mutate(seqvar_id = paste0("asv_",seqvar_id)) %>%
  select(subject_id,subject_day,seqvar_id,read_count) %>%
  distinct() %>%
  spread(key = seqvar_id, value = read_count) %>%
  right_join(filter(an0_data, !seqvar_id %in% sa_seqvars), by = c("subject_id", "subject_day")) %>%
  distinct() -> an0_data_sa_asv
an0_data_sa_asv



an0_data_sa_asv %>%
  filter(grepl("pilot",subject_id)) %>% #select cohort: PILOT
  select(subject_id,subject_day,outcome,seqvar_id,read_count,specimen_read_total,contains("asv_"),contains("before0_7")) %>%
  rename_all(.funs = ~ gsub("_before0_7","",.x)) %>%
  group_by(seqvar_id) %>%
  filter(sum(read_count, na.rm = TRUE) > 0) %>% # only include ASVs detectable in the selected specimens
  ungroup() %>%
  distinct() %>%
  select(-subject_id,-subject_day) %>%
  mutate_at(.vars = vars(contains('read'),contains("asv_")), .funs = ~ as.vector(scale(.x, center=TRUE, scale = TRUE))) %>% # scale read count and specimen read total
  identity() -> d_tib
d_tib



d_tib %>%
  tidybayes::compose_data() -> d_list
str(d_list)



d_tib %>%
  pull(seqvar_id) %>%
  unique() %>%
  sort() %>%
  enframe() %>%
  rename(seqvar_num = name, seqvar_id = value) %>%
  left_join(select(asv_micu_complete,seqvar_id,best,blast_genus), by = "seqvar_id") %>%
  distinct() %>%
  mutate(best = ifelse(best == "unidentified" | best == "uncultured bacterium", blast_genus, best)) -> seqvar_key
seqvar_key



m <- stan(file = './models/binomial/staph_lrti_model-AN_Sa_ASVs_microbiome_abx-reghs_scaled_reads.stan',
          data = d_list,
          iter = 1000,
          chains = 8,
          cores = 8,
          verbose = TRUE,
          seed = 16,
          control=list(adapt_delta=0.95, max_treedepth=10))

rstan::check_hmc_diagnostics(m)



summary(m, probs = c(0.025, 0.975, 0.05, 0.95, 0.1, 0.9, 0.25, 0.75, 0.2, 0.8, 0.4, 0.6, 0.5))$summary %>%
  as_tibble(.,rownames = "id") %>%
  rename(parameter = id) %>%
  mutate(typeS_2.5_pos = `2.5%` > 0 & `97.5%` > 0,
         typeS_2.5_neg = `2.5%` < 0 & `97.5%` < 0,
         typeS_5_pos = `5%` > 0 & `95%` > 0,
         typeS_5_neg = `5%` < 0 & `95%` < 0,
         typeS_10_pos = `10%` > 0 & `90%` > 0,
         typeS_10_neg = `10%` < 0 & `90%` < 0,
         typeS_25_pos = `25%` > 0 & `75%` > 0,
         typeS_25_neg = `25%` < 0 & `75%` < 0) %>%
  select(parameter, n_eff, Rhat, mean, se_mean, sd, contains("typeS"), everything()) %>%
  arrange(desc(parameter))



summary(m, probs = c(0.025, 0.975, 0.05, 0.95, 0.1, 0.9, 0.25, 0.75, 0.2, 0.8, 0.4, 0.6, 0.5))$summary %>%
  as_tibble(.,rownames = "id") %>%
  rename(parameter = id) %>%
  #read_csv(file = "./models/binomial/staph_lrti_model-AN_Sa_ASVs_microbiome_abx-reghs_scaled_reads_PILOT.csv") %>%
  mutate(typeS_2.5_pos = `2.5%` > 0 & `97.5%` > 0,
         typeS_2.5_neg = `2.5%` < 0 & `97.5%` < 0,
         typeS_5_pos = `5%` > 0 & `95%` > 0,
         typeS_5_neg = `5%` < 0 & `95%` < 0,
         typeS_10_pos = `10%` > 0 & `90%` > 0,
         typeS_10_neg = `10%` < 0 & `90%` < 0,
         typeS_25_pos = `25%` > 0 & `75%` > 0,
         typeS_25_neg = `25%` < 0 & `75%` < 0) %>%
  select(parameter, n_eff, Rhat, mean, se_mean, sd, contains("typeS"), everything()) %>%
  arrange(desc(parameter)) %>%
  write_csv(path = "./models/binomial/staph_lrti_model-AN_Sa_ASVs_microbiome_abx-reghs_scaled_reads_PILOT.csv")



summary(m, probs = c(0.025, 0.975, 0.05, 0.95, 0.1, 0.9, 0.25, 0.75, 0.2, 0.8, 0.4, 0.6, 0.5))$summary %>%
  as_tibble(.,rownames = "id") %>%
  rename(parameter = id) %>%
  mutate(typeS_2.5_pos = `2.5%` > 0 & `97.5%` > 0,
         typeS_2.5_neg = `2.5%` < 0 & `97.5%` < 0,
         typeS_5_pos = `5%` > 0 & `95%` > 0,
         typeS_5_neg = `5%` < 0 & `95%` < 0,
         typeS_10_pos = `10%` > 0 & `90%` > 0,
         typeS_10_neg = `10%` < 0 & `90%` < 0,
         typeS_25_pos = `25%` > 0 & `75%` > 0,
         typeS_25_neg = `25%` < 0 & `75%` < 0) %>%
  select(parameter, n_eff, Rhat, mean, se_mean, sd, contains("typeS"), everything()) %>%
  arrange(desc(parameter)) %>%
  filter(grepl("beta",parameter)) %>%
  ggplot(data = .) +
  geom_segment(aes(x = `2.5%`, xend = `97.5%`, y = parameter, yend = parameter)) +
  geom_point(aes(x = mean, y = parameter))



summary(m, probs = c(0.025, 0.975, 0.05, 0.95, 0.1, 0.9, 0.25, 0.75, 0.2, 0.8, 0.4, 0.6, 0.5))$summary %>%
  as_tibble(.,rownames = "id") %>%
  rename(parameter = id) %>%
  #read_csv(file = "./models/binomial/staph_lrti_model-AN_Sa_ASVs_microbiome_abx-reghs_scaled_reads_PILOT.csv") %>%
  mutate(typeS_2.5_pos = `2.5%` > 0 & `97.5%` > 0,
         typeS_2.5_neg = `2.5%` < 0 & `97.5%` < 0,
         typeS_5_pos = `5%` > 0 & `95%` > 0,
         typeS_5_neg = `5%` < 0 & `95%` < 0,
         typeS_10_pos = `10%` > 0 & `90%` > 0,
         typeS_10_neg = `10%` < 0 & `90%` < 0,
         typeS_25_pos = `25%` > 0 & `75%` > 0,
         typeS_25_neg = `25%` < 0 & `75%` < 0) %>%
  select(parameter, n_eff, Rhat, mean, se_mean, sd, contains("typeS"), everything()) %>%
  arrange(desc(parameter)) %>%
  filter(grepl("beta",parameter)) %>%
  mutate(seqvar_num = readr::parse_number(parameter)) %>%
  left_join(seqvar_key, by = 'seqvar_num') %>%
  mutate(param_lab = paste(parameter,seqvar_id,sep = '\n')) %T>%
  write_csv(path = "./models/binomial/staph_lrti_model-AN_Sa_ASVs_microbiome_abx-reghs_scaled_reads_PILOT_seqvar.csv") %>%
  ggplot(data = .) +
  geom_segment(aes(x = `2.5%`, xend = `97.5%`, y = param_lab, yend = param_lab, colour = typeS_10_pos)) +
  geom_point(aes(x = mean, y = param_lab, colour = typeS_10_pos))



summary(m, probs = c(0.025, 0.975, 0.05, 0.95, 0.1, 0.9, 0.25, 0.75, 0.2, 0.8, 0.4, 0.6, 0.5))$summary %>%
  as_tibble(.,rownames = "id") %>%
  rename(parameter = id) %>%
  mutate(typeS_2.5_pos = `2.5%` > 0 & `97.5%` > 0,
         typeS_2.5_neg = `2.5%` < 0 & `97.5%` < 0,
         typeS_5_pos = `5%` > 0 & `95%` > 0,
         typeS_5_neg = `5%` < 0 & `95%` < 0,
         typeS_10_pos = `10%` > 0 & `90%` > 0,
         typeS_10_neg = `10%` < 0 & `90%` < 0,
         typeS_25_pos = `25%` > 0 & `75%` > 0,
         typeS_25_neg = `25%` < 0 & `75%` < 0) %>%
  select(parameter, n_eff, Rhat, mean, se_mean, sd, contains("typeS"), everything()) %>%
  arrange(desc(parameter)) %>%
  filter(grepl("beta",parameter)) %>%
  mutate(seqvar_num = readr::parse_number(parameter)) %>%
  left_join(seqvar_key, by = 'seqvar_num') %>%
  mutate(param_lab = paste(parameter,seqvar_id,sep = '\n')) %>%
  filter(typeS_10_neg | typeS_10_pos) %>%
  select(seqvar_id, blast_genus, `50%`, best) %>%
  distinct()





#' ########################################
#' COMPLETE: LRTI ~ day 0 AN Sa (fixed) + AN microbiome (random) - VALIDATION
#' ########################################

lrti_binomial_Sa_cohorts %>%
  mutate(cohort = factor(cohort, levels = c('Pilot','Validation','Pooled')),
         typeS10 = typeS_10_neg | typeS_10_pos,
         typeS5 = typeS_5_neg | typeS_5_pos) %>%
  group_by(seqvar_id) %>%
  mutate(checktypeS = sum(typeS5) > 0) %>%
  ungroup() %>%
  filter(checktypeS == TRUE & is.na(seqvar_id) == FALSE) %>%
  pull(seqvar_id) %>%
  unique() -> sa_seqvars
sa_seqvars



an0_data %>%
  filter(seqvar_id %in% sa_seqvars) %>%
  mutate(seqvar_id = paste0("asv_",seqvar_id)) %>%
  select(subject_id,subject_day,seqvar_id,read_count) %>%
  distinct() %>%
  spread(key = seqvar_id, value = read_count) %>%
  right_join(filter(an0_data, !seqvar_id %in% sa_seqvars), by = c("subject_id", "subject_day")) %>%
  distinct() -> an0_data_sa_asv
an0_data_sa_asv



an0_data_sa_asv %>%
  filter(grepl("validation",subject_id)) %>% #select cohort: VALIDATION
  select(subject_id,subject_day,outcome,seqvar_id,read_count,specimen_read_total,contains("asv_"),contains("before0_7")) %>%
  rename_all(.funs = ~ gsub("_before0_7","",.x)) %>%
  group_by(seqvar_id) %>%
  filter(sum(read_count, na.rm = TRUE) > 0) %>% # only include ASVs detectable in the selected specimens
  ungroup() %>%
  distinct() %>%
  select(-subject_id,-subject_day) %>%
  mutate_at(.vars = vars(contains('read'),contains("asv_")), .funs = ~ as.vector(scale(.x, center=TRUE, scale = TRUE))) %>% # scale read count and specimen read total
  identity() -> d_tib
d_tib



d_tib %>%
  select(-asv_241b2b5f1ea0f1065ca7027534068198) %>% # validation cohort does not contain asv_241b2b5f1ea0f1065ca7027534068198
  distinct() -> d_tib
d_tib


d_tib %>%
  tidybayes::compose_data() -> d_list
str(d_list)



d_tib %>%
  pull(seqvar_id) %>%
  unique() %>%
  sort() %>%
  enframe() %>%
  rename(seqvar_num = name, seqvar_id = value) %>%
  left_join(select(asv_micu_complete,seqvar_id,best,blast_genus), by = "seqvar_id") %>%
  distinct() %>%
  mutate(best = ifelse(best == "unidentified" | best == "uncultured bacterium", blast_genus, best)) -> seqvar_key
seqvar_key



m <- stan(file = './models/binomial/staph_lrti_model-AN_Sa_ASVs_microbiome_abx-reghs_scaled_reads_validation.stan', # note: model code altered given missing Sa ASV
          data = d_list,
          iter = 1000,
          chains = 8,
          cores = 8,
          verbose = TRUE,
          seed = 16,
          control=list(adapt_delta=0.95, max_treedepth=10))

rstan::check_hmc_diagnostics(m)



summary(m, probs = c(0.025, 0.975, 0.05, 0.95, 0.1, 0.9, 0.25, 0.75, 0.2, 0.8, 0.4, 0.6, 0.5))$summary %>%
  as_tibble(.,rownames = "id") %>%
  rename(parameter = id) %>%
  mutate(typeS_2.5_pos = `2.5%` > 0 & `97.5%` > 0,
         typeS_2.5_neg = `2.5%` < 0 & `97.5%` < 0,
         typeS_5_pos = `5%` > 0 & `95%` > 0,
         typeS_5_neg = `5%` < 0 & `95%` < 0,
         typeS_10_pos = `10%` > 0 & `90%` > 0,
         typeS_10_neg = `10%` < 0 & `90%` < 0,
         typeS_25_pos = `25%` > 0 & `75%` > 0,
         typeS_25_neg = `25%` < 0 & `75%` < 0) %>%
  select(parameter, n_eff, Rhat, mean, se_mean, sd, contains("typeS"), everything()) %>%
  arrange(desc(parameter))



summary(m, probs = c(0.025, 0.975, 0.05, 0.95, 0.1, 0.9, 0.25, 0.75, 0.2, 0.8, 0.4, 0.6, 0.5))$summary %>%
  as_tibble(.,rownames = "id") %>%
  rename(parameter = id) %>%
  #read_csv(file = "./models/binomial/staph_lrti_model-AN_Sa_ASVs_microbiome_abx-reghs_scaled_reads_VALIDATION.csv") %>%
  mutate(typeS_2.5_pos = `2.5%` > 0 & `97.5%` > 0,
         typeS_2.5_neg = `2.5%` < 0 & `97.5%` < 0,
         typeS_5_pos = `5%` > 0 & `95%` > 0,
         typeS_5_neg = `5%` < 0 & `95%` < 0,
         typeS_10_pos = `10%` > 0 & `90%` > 0,
         typeS_10_neg = `10%` < 0 & `90%` < 0,
         typeS_25_pos = `25%` > 0 & `75%` > 0,
         typeS_25_neg = `25%` < 0 & `75%` < 0) %>%
  select(parameter, n_eff, Rhat, mean, se_mean, sd, contains("typeS"), everything()) %>%
  arrange(desc(parameter)) %>%
  write_csv(path = "./models/binomial/staph_lrti_model-AN_Sa_ASVs_microbiome_abx-reghs_scaled_reads_VALIDATION.csv")



summary(m, probs = c(0.025, 0.975, 0.05, 0.95, 0.1, 0.9, 0.25, 0.75, 0.2, 0.8, 0.4, 0.6, 0.5))$summary %>%
  as_tibble(.,rownames = "id") %>%
  rename(parameter = id) %>%
  mutate(typeS_2.5_pos = `2.5%` > 0 & `97.5%` > 0,
         typeS_2.5_neg = `2.5%` < 0 & `97.5%` < 0,
         typeS_5_pos = `5%` > 0 & `95%` > 0,
         typeS_5_neg = `5%` < 0 & `95%` < 0,
         typeS_10_pos = `10%` > 0 & `90%` > 0,
         typeS_10_neg = `10%` < 0 & `90%` < 0,
         typeS_25_pos = `25%` > 0 & `75%` > 0,
         typeS_25_neg = `25%` < 0 & `75%` < 0) %>%
  select(parameter, n_eff, Rhat, mean, se_mean, sd, contains("typeS"), everything()) %>%
  arrange(desc(parameter)) %>%
  filter(grepl("beta",parameter)) %>%
  ggplot(data = .) +
  geom_segment(aes(x = `2.5%`, xend = `97.5%`, y = parameter, yend = parameter)) +
  geom_point(aes(x = mean, y = parameter))



summary(m, probs = c(0.025, 0.975, 0.05, 0.95, 0.1, 0.9, 0.25, 0.75, 0.2, 0.8, 0.4, 0.6, 0.5))$summary %>%
  as_tibble(.,rownames = "id") %>%
  rename(parameter = id) %>%
  #read_csv(file = "./models/binomial/staph_lrti_model-AN_Sa_ASVs_microbiome_abx-reghs_scaled_reads_VALIDATION.csv") %>%
  mutate(typeS_2.5_pos = `2.5%` > 0 & `97.5%` > 0,
         typeS_2.5_neg = `2.5%` < 0 & `97.5%` < 0,
         typeS_5_pos = `5%` > 0 & `95%` > 0,
         typeS_5_neg = `5%` < 0 & `95%` < 0,
         typeS_10_pos = `10%` > 0 & `90%` > 0,
         typeS_10_neg = `10%` < 0 & `90%` < 0,
         typeS_25_pos = `25%` > 0 & `75%` > 0,
         typeS_25_neg = `25%` < 0 & `75%` < 0) %>%
  select(parameter, n_eff, Rhat, mean, se_mean, sd, contains("typeS"), everything()) %>%
  arrange(desc(parameter)) %>%
  filter(grepl("beta",parameter)) %>%
  mutate(seqvar_num = readr::parse_number(parameter)) %>%
  left_join(seqvar_key, by = 'seqvar_num') %>%
  mutate(param_lab = paste(parameter,seqvar_id,sep = '\n')) %T>%
  write_csv(path = "./models/binomial/staph_lrti_model-AN_Sa_ASVs_microbiome_abx-reghs_scaled_reads_VALIDATION_seqvar.csv") %>%
  ggplot(data = .) +
  geom_segment(aes(x = `2.5%`, xend = `97.5%`, y = param_lab, yend = param_lab, colour = typeS_10_pos)) +
  geom_point(aes(x = mean, y = param_lab, colour = typeS_10_pos))



summary(m, probs = c(0.025, 0.975, 0.05, 0.95, 0.1, 0.9, 0.25, 0.75, 0.2, 0.8, 0.4, 0.6, 0.5))$summary %>%
  as_tibble(.,rownames = "id") %>%
  rename(parameter = id) %>%
  mutate(typeS_2.5_pos = `2.5%` > 0 & `97.5%` > 0,
         typeS_2.5_neg = `2.5%` < 0 & `97.5%` < 0,
         typeS_5_pos = `5%` > 0 & `95%` > 0,
         typeS_5_neg = `5%` < 0 & `95%` < 0,
         typeS_10_pos = `10%` > 0 & `90%` > 0,
         typeS_10_neg = `10%` < 0 & `90%` < 0,
         typeS_25_pos = `25%` > 0 & `75%` > 0,
         typeS_25_neg = `25%` < 0 & `75%` < 0) %>%
  select(parameter, n_eff, Rhat, mean, se_mean, sd, contains("typeS"), everything()) %>%
  arrange(desc(parameter)) %>%
  filter(grepl("beta",parameter)) %>%
  mutate(seqvar_num = readr::parse_number(parameter)) %>%
  left_join(seqvar_key, by = 'seqvar_num') %>%
  mutate(param_lab = paste(parameter,seqvar_id,sep = '\n')) %>%
  filter(typeS_10_neg | typeS_10_pos) %>%
  select(seqvar_id, blast_genus, `50%`, best) %>%
  distinct()






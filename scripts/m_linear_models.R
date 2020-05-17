#' ########################################
#' linear models relating ET ~ AN Staph aureus
#' ########################################
#' 
#' depends: asv_micu_complete & micu_complete
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




#' ########################################
#' POOLED COHORTS
#' ########################################

asv_micu_complete %>%
  filter(subject_id %in% micu_cohort$subject_id) %>% # included subjects only
  mutate(Sa_prop = sumstaph / specimen_read_total) %>%
  select(subject_id,subject_day,specimen_type,Sa_prop) %>%
  distinct() %>%
  filter(specimen_type %in% c("AN","ET")) -> et_an_comparison

et_an_comparison %>%
  group_by(subject_id,subject_day) %>%
  filter(paste(sort(specimen_type), collapse = " ") == "AN ET") %>% # only include same-day same-subject pairs
  spread(key = specimen_type, value = Sa_prop) -> et_an_same_day
et_an_same_day

et_an_same_day %>%
  rename(AN_Sa = AN,
         ET_Sa = ET) -> et_an_same_day_Sa
et_an_same_day_Sa



asv_micu_complete %>%
  filter(subject_id %in% micu_cohort$subject_id) %>% # included subjects only
  filter(grepl("aureus",best)) %>% # S.aureus ASVs only
  filter(specimen_type %in% c("AN","ET")) %>% # AN and ET specimens only
  select(subject_id,subject_day,specimen_type,seqvar_id,read_count,specimen_read_total) %>%
  distinct() %>%
  group_by(subject_id,subject_day,specimen_type,seqvar_id) %>%
  summarise(Sa_ASV_prop = read_count / specimen_read_total) %>%
  ungroup() %>%
  distinct() %>%
  group_by(seqvar_id) %>%
  filter(sum(Sa_ASV_prop, na.rm = TRUE) > 0) %>% # filter out any Staph aureus ASVs not present in this specimen set
  ungroup() %>%
  distinct() -> et_an_asv_comparison
et_an_asv_comparison

et_an_asv_comparison %>%
  group_by(subject_id,subject_day,seqvar_id) %>%
  filter(paste(sort(specimen_type), collapse = " ") == "AN ET") %>% # only include same-day same-subject pairs
  spread(key = specimen_type, value = Sa_ASV_prop) %>%
  ungroup() %>%
  distinct() -> et_an_asv_same_day
et_an_asv_same_day

et_an_asv_same_day %>%
  rename(AN_asv = AN,
         ET_asv = ET) -> et_an_asv_same_day_Sa
et_an_asv_same_day_Sa


et_an_asv_same_day_Sa %>%
  left_join(et_an_same_day_Sa, by = c("subject_id", "subject_day")) %>%
  #mutate_at(.vars = vars(contains("AN")), .funs = ~ log(.x)) %>%
  #mutate_at(.vars = vars(contains("ET")), .funs = ~ log(.x)) %>%
  select(seqvar_id,ET_Sa,AN_asv) %>%
  na.omit() %>%
  #filter(abs(ET_Sa) != Inf & abs(AN_asv) != Inf) %>%
  group_by(seqvar_id) %>%
  filter(sum(AN_asv, na.rm = TRUE) > 0) %>% # remove ASVs not found in this specimen set
  ungroup() %>%
  tidybayes::compose_data() -> d_list
str(d_list)


et_an_asv_same_day_Sa %>%
  left_join(et_an_same_day_Sa, by = c("subject_id", "subject_day")) %>%
  #mutate_at(.vars = vars(contains("AN")), .funs = ~ log(.x)) %>%
  #mutate_at(.vars = vars(contains("ET")), .funs = ~ log(.x)) %>%
  select(seqvar_id,ET_Sa,AN_asv) %>%
  na.omit() %>%
  #filter(abs(ET_Sa) != Inf & abs(AN_asv) != Inf) %>%
  group_by(seqvar_id) %>%
  filter(sum(AN_asv, na.rm = TRUE) > 0) %>% # remove ASVs not found in this specimen set
  ungroup() %>%
  pull(seqvar_id) %>%
  unique() %>%
  sort() %>%
  enframe() %>%
  rename(seqvar_num = name, seqvar_id = value) %>%
  left_join(select(asv_micu_complete,seqvar_id,best,blast_genus), by = "seqvar_id") %>%
  distinct() %>%
  mutate(best = ifelse(best == "unidentified" | best == "uncultured bacterium", blast_genus, best)) -> seqvar_key
seqvar_key


m <- stan(file = './models/linear/staph_model_prop_ET_v_ANsv_reghs.stan',
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
  #write_csv(path = "./models/linear/staph_model_prop_ET_v_ANsv_reghs_POOLED.csv") %>%
  identity()

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
  mutate(param_lab = paste(parameter,seqvar_id,sep = '\n')) %>%
  #write_csv(path = "./models/linear/staph_model_prop_ET_v_ANsv_reghs_POOLED_seqvar.csv") %>%
  ggplot(data = .) +
  geom_segment(aes(x = `2.5%`, xend = `97.5%`, y = param_lab, yend = param_lab, colour = typeS_10_pos)) +
  geom_point(aes(x = mean, y = param_lab, colour = typeS_10_pos))






#' ########################################
#' PILOT COHORT
#' ########################################

asv_micu_complete %>%
  filter(subject_id %in% micu_cohort$subject_id) %>% # included subjects only
  mutate(Sa_prop = sumstaph / specimen_read_total) %>%
  select(subject_id,subject_day,specimen_type,Sa_prop) %>%
  distinct() %>%
  filter(specimen_type %in% c("AN","ET")) -> et_an_comparison

et_an_comparison %>%
  group_by(subject_id,subject_day) %>%
  filter(paste(sort(specimen_type), collapse = " ") == "AN ET") %>% # only include same-day same-subject pairs
  spread(key = specimen_type, value = Sa_prop) -> et_an_same_day
et_an_same_day

et_an_same_day %>%
  rename(AN_Sa = AN,
         ET_Sa = ET) -> et_an_same_day_Sa
et_an_same_day_Sa



asv_micu_complete %>%
  filter(subject_id %in% micu_cohort$subject_id) %>% # included subjects only
  filter(grepl("aureus",best)) %>% # S.aureus ASVs only
  filter(specimen_type %in% c("AN","ET")) %>% # AN and ET specimens only
  select(subject_id,subject_day,specimen_type,seqvar_id,read_count,specimen_read_total) %>%
  distinct() %>%
  group_by(subject_id,subject_day,specimen_type,seqvar_id) %>%
  summarise(Sa_ASV_prop = read_count / specimen_read_total) %>%
  ungroup() %>%
  distinct() %>%
  group_by(seqvar_id) %>%
  filter(sum(Sa_ASV_prop, na.rm = TRUE) > 0) %>% # filter out any Staph aureus ASVs not present in this specimen set
  ungroup() %>%
  distinct() -> et_an_asv_comparison
et_an_asv_comparison

et_an_asv_comparison %>%
  group_by(subject_id,subject_day,seqvar_id) %>%
  filter(paste(sort(specimen_type), collapse = " ") == "AN ET") %>% # only include same-day same-subject pairs
  spread(key = specimen_type, value = Sa_ASV_prop) %>%
  ungroup() %>%
  distinct() -> et_an_asv_same_day
et_an_asv_same_day

et_an_asv_same_day %>%
  rename(AN_asv = AN,
         ET_asv = ET) -> et_an_asv_same_day_Sa
et_an_asv_same_day_Sa


et_an_asv_same_day_Sa %>%
  left_join(et_an_same_day_Sa, by = c("subject_id", "subject_day")) %>%
  filter(grepl("pilot", subject_id)) %>%
  #mutate_at(.vars = vars(contains("AN")), .funs = ~ log(.x)) %>%
  #mutate_at(.vars = vars(contains("ET")), .funs = ~ log(.x)) %>%
  select(seqvar_id,ET_Sa,AN_asv) %>%
  na.omit() %>%
  #filter(abs(ET_Sa) != Inf & abs(AN_asv) != Inf) %>%
  group_by(seqvar_id) %>%
  filter(sum(AN_asv, na.rm = TRUE) > 0) %>% # remove ASVs not found in this specimen set
  ungroup() %>%
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



m <- stan(file = './models/linear/staph_model_prop_ET_v_ANsv_reghs.stan',
          data = d_list,
          iter = 1000,
          chains = 8,
          cores = 8,
          verbose = TRUE,
          seed = 16,
          control=list(adapt_delta=0.95, max_treedepth=10))

rstan::check_hmc_diagnostics(m)
# save(m, file = "./models/linear/staph_model_prop_ET_v_ANsv_reghs_PILOT.rdata",compress = TRUE)

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
  #write_csv(path = "./models/linear/staph_model_prop_ET_v_ANsv_reghs_PILOT.csv") %>%
  identity()

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
  mutate(param_lab = paste(parameter,seqvar_id,sep = '\n')) %>%
  #write_csv(path = "./models/linear/staph_model_prop_ET_v_ANsv_reghs_PILOT_seqvar.csv") %>%
  ggplot(data = .) +
  geom_segment(aes(x = `2.5%`, xend = `97.5%`, y = param_lab, yend = param_lab, colour = typeS_10_pos)) +
  geom_point(aes(x = mean, y = param_lab, colour = typeS_10_pos))



#' ########################################
#' VALIDATION COHORT
#' ########################################

asv_micu_complete %>%
  filter(subject_id %in% micu_cohort$subject_id) %>% # included subjects only
  mutate(Sa_prop = sumstaph / specimen_read_total) %>%
  select(subject_id,subject_day,specimen_type,Sa_prop) %>%
  distinct() %>%
  filter(specimen_type %in% c("AN","ET")) -> et_an_comparison

et_an_comparison %>%
  group_by(subject_id,subject_day) %>%
  filter(paste(sort(specimen_type), collapse = " ") == "AN ET") %>% # only include same-day same-subject pairs
  spread(key = specimen_type, value = Sa_prop) -> et_an_same_day
et_an_same_day

et_an_same_day %>%
  rename(AN_Sa = AN,
         ET_Sa = ET) -> et_an_same_day_Sa
et_an_same_day_Sa



asv_micu_complete %>%
  filter(subject_id %in% micu_cohort$subject_id) %>% # included subjects only
  filter(grepl("aureus",best)) %>% # S.aureus ASVs only
  filter(specimen_type %in% c("AN","ET")) %>% # AN and ET specimens only
  select(subject_id,subject_day,specimen_type,seqvar_id,read_count,specimen_read_total) %>%
  distinct() %>%
  group_by(subject_id,subject_day,specimen_type,seqvar_id) %>%
  summarise(Sa_ASV_prop = read_count / specimen_read_total) %>%
  ungroup() %>%
  distinct() %>%
  group_by(seqvar_id) %>%
  filter(sum(Sa_ASV_prop, na.rm = TRUE) > 0) %>% # filter out any Staph aureus ASVs not present in this specimen set
  ungroup() %>%
  distinct() -> et_an_asv_comparison
et_an_asv_comparison

et_an_asv_comparison %>%
  group_by(subject_id,subject_day,seqvar_id) %>%
  filter(paste(sort(specimen_type), collapse = " ") == "AN ET") %>% # only include same-day same-subject pairs
  spread(key = specimen_type, value = Sa_ASV_prop) %>%
  ungroup() %>%
  distinct() -> et_an_asv_same_day
et_an_asv_same_day

et_an_asv_same_day %>%
  rename(AN_asv = AN,
         ET_asv = ET) -> et_an_asv_same_day_Sa
et_an_asv_same_day_Sa


et_an_asv_same_day_Sa %>%
  left_join(et_an_same_day_Sa, by = c("subject_id", "subject_day")) %>%
  filter(grepl("validation", subject_id)) %>%
  #mutate_at(.vars = vars(contains("AN")), .funs = ~ log(.x)) %>%
  #mutate_at(.vars = vars(contains("ET")), .funs = ~ log(.x)) %>%
  select(seqvar_id,ET_Sa,AN_asv) %>%
  na.omit() %>%
  #filter(abs(ET_Sa) != Inf & abs(AN_asv) != Inf) %>%
  group_by(seqvar_id) %>%
  filter(sum(AN_asv, na.rm = TRUE) > 0) %>% # remove ASVs not found in this specimen set
  ungroup() %>%
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



m <- stan(file = './models/linear/staph_model_prop_ET_v_ANsv_reghs.stan',
          data = d_list,
          iter = 1000,
          chains = 8,
          verbose = TRUE,
          seed = 16,
          control=list(adapt_delta=0.95, max_treedepth=10))

rstan::check_hmc_diagnostics(m)
# save(m, file = "./models/linear/staph_model_prop_ET_v_ANsv_reghs_VALIDATION.rdata",compress = TRUE)

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
  #write_csv(path = "./models/linear/staph_model_prop_ET_v_ANsv_reghs_VALIDATION.csv") %>%
  identity()

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
  mutate(param_lab = paste(parameter,seqvar_id,sep = '\n')) %>%
  #write_csv(path = "./models/linear/staph_model_prop_ET_v_ANsv_reghs_VALIDATION_seqvar.csv") %>%
  ggplot(data = .) +
  geom_segment(aes(x = `2.5%`, xend = `97.5%`, y = param_lab, yend = param_lab, colour = typeS_10_pos)) +
  geom_point(aes(x = mean, y = param_lab, colour = typeS_10_pos))



#' ########################################
#' COMBINE COHORTS
#' ########################################

sa_linear_cohorts <- bind_rows(mutate(read_csv("./models/linear/staph_model_prop_ET_v_ANsv_reghs_POOLED_seqvar.csv"), cohort = 'Pooled'),
                               mutate(read_csv("./models/linear/staph_model_prop_ET_v_ANsv_reghs_PILOT_seqvar.csv"), cohort = 'Pilot'),
                               mutate(read_csv("./models/linear/staph_model_prop_ET_v_ANsv_reghs_VALIDATION_seqvar.csv"), cohort = 'Validation')) %>%
  filter(grepl('beta',parameter))
sa_linear_cohorts


sa_linear_cohorts %>%
  #mutate(cohort = factor(cohort, levels = c('Pilot','Validation','Pooled')),
  #       typeS10 = typeS_10_neg | typeS_10_pos,
  #       typeS5 = typeS_5_neg | typeS_5_pos) %>%
  mutate(cohort = factor(cohort, levels = c('Pilot','Validation','Pooled')),
         typeS10 = typeS_10_neg | typeS_10_pos,
         typeS5 = typeS_5_neg | typeS_5_pos,
         typeS25 = typeS_25_neg | typeS_25_pos,
         typeScol = ifelse(typeS5, "< 0.05", ifelse(typeS10, "< 0.1", ifelse(typeS25, "< 0.25", ">= 0.25")))) %>% #count(typeScol)
  group_by(seqvar_id) %>%
  mutate(checktypeS = sum(typeS10)) %>%
  ungroup() %>%
  #filter(checktypeS > 0) %>%
  mutate(blast_genus = gsub(" 1","",blast_genus),
         seqvar_lab = paste0(blast_genus,"\n(",seqvar_id,")"),
         typeScol = factor(typeScol, levels = c("< 0.05", "< 0.1", "< 0.25", ">= 0.25"))) %>% #count(typeScol)
  ggplot(data = .) +
  geom_vline(xintercept = 0, linetype = 2) +
  geom_segment(aes(x = `25%`, xend = `75%`, y = seqvar_id, yend = seqvar_id, colour = typeScol)) +
  geom_point(aes(x = `50%`, y = seqvar_id, colour = typeScol)) +
  facet_wrap(~cohort, nrow = 1) +
  #scale_colour_manual(values = c('TRUE' = 'red', 'FALSE' = 'black'), guide = FALSE) +
  scale_colour_manual(values = c('< 0.05' = '#cb181d', '< 0.1' = '#fb6a4a', '< 0.25' = '#fcae91', '>= 0.25' = 'darkgrey'), drop = FALSE) +
  scale_x_continuous(breaks = seq(from = -1, to = 5, by = 1)) +
  theme_bw() +
  theme(strip.background = element_blank(),
        axis.text.x = element_text(colour = 'black'),
        axis.text.y = element_text(colour = 'black', size = 6),
        plot.title = element_text(hjust = 0),
        plot.subtitle = element_text(hjust = 0, size = 10),
        legend.position = 'bottom') +
  labs(x = expression(paste("ET ", italic("S. aureus"), " Attributable to AN ASV (", beta, " posterior median and 50% CI)")),
       y = expression(paste(italic("S. aureus"), " Anterior Nares ASVs")),
       colour = "Type S Error") -> p_sa_linear_cohorts

p_sa_linear_cohorts

# ggsave(plot = p_sa_linear_cohorts, filename = "./figs/supp/p_sa_linear_cohorts.png", height = 4, width = 8, units = "in", dpi = 600)
# ggsave(plot = p_sa_linear_cohorts, filename = "./figs/supp/p_sa_linear_cohorts.svg", height = 4, width = 8, units = "in")
# ggsave(plot = p_sa_linear_cohorts, filename = "./figs/supp/p_sa_linear_cohorts.pdf", height = 4, width = 8, units = "in")




# test version with viridis plasma color scale:
p_sa_linear_cohorts +
  scale_colour_manual(values = c('< 0.05' = '#7301A8FF', '< 0.1' = '#BD3786FF', '< 0.25' = '#ED7953FF', '>= 0.25' = 'darkgrey'), drop = FALSE)





#
###
#####
###
#








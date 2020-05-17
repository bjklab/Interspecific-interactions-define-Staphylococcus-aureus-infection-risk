#' ########################################
#' joint species distribution model of family-level & ASV-level (all ASVs)
#' ########################################
#' 
#' depends: model output from jsdm script
#' 

library(tidyverse)
library(patchwork)
set.seed(16)



#' #########################################
#' 
#' read posterior betas from pilot and validation models
#' 
#' #########################################

list.files(path = "./models/jsdm/", pattern = "pilot_posterior_betas", full.names = TRUE) %>%
  map(.f = ~ read_rds(.x)) %>%
  bind_rows() %>%
  mutate(cohort = "pilot") -> pilot_posterior_betas


list.files(path = "./models/jsdm/", pattern = "valid_posterior_betas", full.names = TRUE) %>%
  map(.f = ~ read_rds(.x)) %>%
  bind_rows() %>%
  mutate(cohort = "validation") -> valid_posterior_betas


bind_rows(pilot_posterior_betas, valid_posterior_betas) %>%
  identity() -> posterior_betas

rm(pilot_posterior_betas, valid_posterior_betas)




#' #########################################
#' 
#' read posterior betas from pilot and validation models
#' 
#' #########################################

posterior_betas %>%
  gather(key = "beta", value = "posterior_draws", -asv_number, -asv_id, -seqvar_id, -blast_family, -blast_genus, -best, -chain, -cohort) %>%
  group_by(cohort, asv_number, asv_id, seqvar_id, blast_family, blast_genus, best, beta) %>%
  summarise_at(.vars = vars(posterior_draws), .funs = list("hpdi_low" = ~ rethinking::HPDI(.x, prob = 0.95)[1],
                                                           "hpdi_high" = ~ rethinking::HPDI(.x, prob = 0.95)[2],
                                                           "pi_low" = ~ rethinking::PI(.x, prob = 0.95)[1],
                                                           "pi_high" = ~ rethinking::PI(.x, prob = 0.95)[2],
                                                           "median" = ~ quantile(.x, 0.5),
                                                           "type_S" = ~ min(sum(.x > 0, na.rm = TRUE)/length(.x), sum(.x < 0, na.rm = TRUE)/length(.x), na.rm = TRUE))) %>%
  ungroup() %>%
  mutate(posterior_certainty = type_S < 0.05) %>%
  identity() -> posterior_asv_summary

posterior_asv_summary



posterior_asv_summary %>%
  filter(grepl("aureus",best)) %>%
  ggplot(data = .) +
  geom_segment(aes(x = hpdi_low, xend = hpdi_high, y = seqvar_id, yend = seqvar_id, color = posterior_certainty)) +
  geom_point(aes(x = median, y = seqvar_id, color = posterior_certainty)) +
  facet_wrap(facets = ~ beta + cohort, scales = "free_x") +
  labs(title = "Type S HPDI Low")










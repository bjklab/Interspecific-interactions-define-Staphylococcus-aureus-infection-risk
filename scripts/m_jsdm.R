#' ########################################
#' joint species distribution model of family-level & ASV-level (all ASVs)
#' ########################################
#' 
#' depends: micu_cohort, asv_micu_complete
#' 

library(tidyverse)
library(jSDM)
library(parallel)
#library(patchwork)
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






#' ########################################
#' AN ADMISSION JSDM
#' ########################################

#' outcome matrix for probit regression
asv_micu_complete %>%
  filter(!is.na(seqvar_id)) %>%
  filter(specimen_type == "AN") %>%
  filter(grepl("pilot",subject_id)) %>%
  filter(subject_day == 0) %>%
  select(subject_id, seqvar_id, read_count) %>%
  distinct() %>%
  #qplot(data = ., x = read_count) + scale_x_log10() + geom_vline(xintercept = 10, color = "red", linetype = 2)
  mutate(seqvar_id = paste0("asv_",seqvar_id),
         read_count = read_count > 10, # for probit regression, presence > 10 reads
         read_count = as.numeric(read_count)) %>%
  group_by(seqvar_id) %>%
  filter(sum(read_count) > 0) %>% # filter out unobserved ASVs
  ungroup() %>%
  #pivot_wider(id_cols = "subject_id", names_from = "seqvar_id", values_from = "read_count") %>%
  reshape2::acast(subject_id ~ seqvar_id, value.var = "read_count") -> y_pilot
str(y_pilot)




asv_micu_complete %>%
  filter(!is.na(seqvar_id)) %>%
  filter(specimen_type == "AN") %>%
  filter(grepl("validation",subject_id)) %>%
  filter(subject_day == 0) %>%
  select(subject_id, seqvar_id, read_count) %>%
  distinct() %>%
  #qplot(data = ., x = read_count) + scale_x_log10() + geom_vline(xintercept = 10, color = "red", linetype = 2)
  mutate(seqvar_id = paste0("asv_",seqvar_id),
         read_count = read_count > 10, # for probit regression, presence > 10 reads
         read_count = as.numeric(read_count)) %>%
  group_by(seqvar_id) %>%
  filter(sum(read_count) > 0) %>% # filter out unobserved ASVs
  ungroup() %>%
  #pivot_wider(id_cols = "subject_id", names_from = "seqvar_id", values_from = "read_count") %>%
  reshape2::acast(subject_id ~ seqvar_id, value.var = "read_count") -> y_valid
str(y_valid)



asv_micu_complete %>%
  filter(!is.na(seqvar_id)) %>%
  filter(specimen_type == "AN") %>%
  #filter(grepl("validation",subject_id)) %>%
  filter(subject_day == 0) %>%
  select(subject_id, seqvar_id, read_count) %>%
  distinct() %>%
  #qplot(data = ., x = read_count) + scale_x_log10() + geom_vline(xintercept = 10, color = "red", linetype = 2)
  mutate(seqvar_id = paste0("asv_",seqvar_id),
         read_count = read_count > 10, # for probit regression, presence > 10 reads
         read_count = as.numeric(read_count)) %>%
  group_by(seqvar_id) %>%
  filter(sum(read_count) > 0) %>% # filter out unobserved ASVs
  ungroup() %>%
  #pivot_wider(id_cols = "subject_id", names_from = "seqvar_id", values_from = "read_count") %>%
  reshape2::acast(subject_id ~ seqvar_id, value.var = "read_count") -> y_pooled
str(y_pooled)




#' exposure matrix for probit regression
micu_cohort %>%
  filter(grepl("pilot",subject_id)) %>%
  select(age, gender, race, copd, asthma, ild, lymphleuk, dm, cirrhosis, chf, contains("before0_7")) %>%
  mutate_if(.predicate = ~ is.character(.x), .funs = ~ factor(.x)) %>%
  identity() -> x_pilot
x_pilot

micu_cohort %>%
  filter(grepl("validation",subject_id)) %>%
  select(age, gender, race, copd, asthma, ild, lymphleuk, dm, cirrhosis, chf, contains("before0_7")) %>%
  mutate_if(.predicate = ~ is.character(.x), .funs = ~ factor(.x)) %>%
  identity() -> x_valid
x_valid  


micu_cohort %>%
  #filter(grepl("validation",subject_id)) %>%
  select(age, gender, race, copd, asthma, ild, lymphleuk, dm, cirrhosis, chf, contains("before0_7")) %>%
  mutate_if(.predicate = ~ is.character(.x), .funs = ~ factor(.x)) %>%
  identity() -> x_pooled
x_pooled  




#' #############################################
#' 
#' RUN MODELS: jSDM package
#' 
#' #############################################

#' jSDM package demonstration dataset
# data(frogs, package="jSDM")
# head(frogs)
# 
# # data.obs
# PA_frogs <- frogs[,4:12]
# 
# # Normalized continuous variables
# Env_frogs <- cbind(scale(frogs[,1]),frogs[,2],scale(frogs[,3]))
# colnames(Env_frogs) <- colnames(frogs[,1:3])
# 
# mod_jSDM_block_frogs <- jSDM_probit_block (
#   # Response variable 
#   presence_site_sp = as.matrix(PA_frogs), 
#   # Explanatory variables 
#   site_suitability = ~.,   
#   site_data = as.data.frame(Env_frogs), n_latent=2,
#   # Chains
#   burnin=20000, mcmc=5000, thin=5,
#   # Starting values
#   alpha_start=0, beta_start=0,
#   lambda_start=0, W_start=0,
#   V_alpha_start=1, 
#   # Priors
#   shape=0.5, rate=0.0005,
#   mu_beta=0, V_beta=1.0E6,
#   mu_lambda=0, V_lambda=10,
#   # Various 
#   seed=1234, verbose=1)
# 
# #plot(coda::as.mcmc(mod_jSDM_block_frogs$mcmc.alpha[,1:2]))
# 
# plot_residual_cor(mod_jSDM_block_frogs)


#' ######################################
#' 
#' jSDM pilot data - test
#' 
#' ######################################

# x_pilot %>%
#   select(age,
#          copd,
#          asthma,
#          ild,
#          lymphleuk,
#          dm,
#          cirrhosis,
#          chf,
#          contains("before0_7")) %>%
#   #mutate(age = age = scale(age)[,1]) %>%
#   mutate_if(.predicate = ~ is.numeric(.x), .funs = ~ scale(.x)[,1]) %>%
#   mutate_if(.predicate = ~ is.logical(.x), .funs = ~ as.numeric(.x)) %>%
#   identity() -> dx_pilot
# dx_pilot
# 
# mod_jSDM_block_pilot <- jSDM_probit_block (
#   # Response variable 
#   presence_site_sp = as.matrix(y_pilot), 
#   # Explanatory variables 
#   site_suitability = ~.,   
#   site_data = as.data.frame(dx_pilot), n_latent=2,
#   # Chains
#   burnin=20000, mcmc=10000, thin=5,
#   # Starting values
#   alpha_start=0, beta_start=0,
#   lambda_start=0, W_start=0,
#   V_alpha_start=1, 
#   # Priors
#   shape=0.5, rate=0.0005,
#   mu_beta=0, V_beta=1.0E6,
#   mu_lambda=0, V_lambda=10,
#   # Various 
#   seed=16, verbose=1)
# 
# mod_jSDM_block_pilot %>%
#   write_rds(path = "./models/jsdm/mod_jSDM_block_pilot.rds.gz", compress = "gz")
# 
# mod_jSDM_block_pilot <- read_rds(path = "./models/jsdm/mod_jSDM_block_pilot.rds.gz")
# 
# str(mod_jSDM_block_pilot)
# 
# str(mod_jSDM_block_pilot$mcmc.alpha)
# 
# # plot(coda::as.mcmc(mod_jSDM_block_pilot$mcmc.alpha[,1:14]))
# 
# plot(mod_jSDM_block_pilot$mcmc.Deviance,main = "Deviance")
# 
# #plot_residual_cor(mod_jSDM_block_pilot)
# 
# get_residual_cor(mod_jSDM_block_pilot) -> residcor_jSDM_block_pilot
# 
# residcor_jSDM_block_pilot$cor.median %>%
#   as_tibble() %>%
#   mutate(asv_x = colnames(y_pilot)) %>%
#   gather(key = "asv_y", value = "cor", -asv_x) %>%
#   mutate(asv_y = structure(.Data = colnames(y_pilot), .Names = paste0("V",seq(ncol(y_pilot))))[asv_y]) -> resid_cor_pilot
# resid_cor_pilot
# 
# resid_cor_pilot %>%
#   qplot(data = ., x = asv_x, y = asv_y, geom = "tile", fill = cor) +
#   scale_fill_viridis_c(option = "magma") +
#   coord_equal() +
#   theme_classic() +
#   theme(axis.text.x = element_text(angle = 270, vjust = 0.5, hjust = 0))
# 
# resid_cor_pilot %>%
#   filter(asv_x %in% pull(filter(asv_key, grepl("aureus", best)),asv_id)) %>%
#   mutate(genus_y = structure(.Data = asv_key$blast_genus, .Names = asv_key$asv_id)[asv_y],
#          asv_genus = paste0(genus_y,"\n(",gsub("asv_","",asv_y),")")) %>%
#   arrange(cor) %>%
#   print(n=100)
# 
# resid_cor_pilot %>%
#   filter(asv_x %in% pull(filter(asv_key, grepl("aureus", best)),asv_id)) %>%
#   mutate(genus_y = structure(.Data = asv_key$blast_genus, .Names = asv_key$asv_id)[asv_y],
#          asv_genus = paste0(genus_y,"\n(",gsub("asv_","",asv_y),")")) %>%
#   arrange(desc(cor)) %>%
#   print(n=100)



#' ######################################
#' 
#' jSDM pilot data - parallel chains
#' 
#' ######################################

x_pilot %>%
  select(#age,
         #copd,
         #asthma,
         #ild,
         #lymphleuk,
         #dm,
         #cirrhosis,
         #chf,
         contains("before0_7")) %>%
  #mutate(age = age = scale(age)[,1]) %>%
  mutate_if(.predicate = ~ is.numeric(.x), .funs = ~ scale(.x)[,1]) %>%
  mutate_if(.predicate = ~ is.logical(.x), .funs = ~ as.numeric(.x)) %>%
  identity() -> dx_pilot
dx_pilot


run_jsdm <- function(x,y,seed_number) {
  my_jsdm <- jSDM_probit_block (
  # Response variable 
  presence_site_sp = as.matrix(y), 
  # Explanatory variables 
  site_suitability = ~.,   
  site_data = as.data.frame(x), n_latent=2,
  # Chains
  burnin=20000, mcmc=5000, thin=5,
  # Starting values
  alpha_start=0, beta_start=0,
  lambda_start=0, W_start=0,
  V_alpha_start=1, 
  # Priors
  shape=0.5, rate=0.0005,
  mu_beta=0, V_beta=1.0E6,
  mu_lambda=0, V_lambda=10,
  # Various 
  seed=seed_number,
  verbose=1)
  return(my_jsdm)
}


parallel::mcMap(f = function(x,y,z) run_jsdm(x,y,z),
                y = list(y_pilot, y_pilot, y_pilot, y_pilot),
                x = list(dx_pilot, dx_pilot, dx_pilot, dx_pilot),
                z = as.list(seq(4)*16),
                mc.cores = 4
                ) -> pmod_jSDM_block_pilot

names(pmod_jSDM_block_pilot) <- paste0("chain_",seq_along(pmod_jSDM_block_pilot))

str(pmod_jSDM_block_pilot[1])

pmod_jSDM_block_pilot %>%
  map2(.x = .,
       .y = as.list(paste0("./models/jsdm/pmod_jSDM_block_pilot_",names(.),".rds.gz")),
       .f = ~ write_rds(x = .x, path = .y, compress = "gz")
       )

as.list(structure(.Data = c("chain_1","chain_2","chain_3","chain_4"), .Names = c("chain_1","chain_2","chain_3","chain_4"))) %>%
  map(.f = ~ read_rds(path = paste0("./models/jsdm/pmod_jSDM_block_pilot_",.x,".rds.gz"))) -> pmod_jSDM_block_pilot

str(pmod_jSDM_block_pilot[[1]]$mcmc.sp)



  
#' POSTERIOR BETAS & DIAGNOSTICS

# map2(.f = ~ rbind(.x,.y),
#      .x = list(matrix(seq(5,8),nrow=2),matrix(seq(5,8),nrow=2)),
#      .y = list(matrix(seq(4),nrow=2),matrix(seq(4),nrow=2)))
# 
# pmap(.l = list(list(matrix(seq(5,8),nrow=2),matrix(seq(5,8),nrow=2)),
#                list(matrix(seq(4),nrow=2),matrix(seq(4),nrow=2)),
#                list(matrix(seq(9,12),nrow=2),matrix(seq(9,12),nrow=2))),
#      .f = rbind)

map(.x = pmod_jSDM_block_pilot, .f = ~ pluck(.x, "mcmc.sp")) %>%
  pmap(.l = ., .f = rbind) %>%
  map(.x = ., .f = ~ as_tibble(.x)) %>%
  bind_rows(.id = "asv_number") %>%
  mutate(asv_id = structure(.Names = paste0("sp_",seq(ncol(y_pilot))), .Data = colnames(y_pilot))[asv_number]) %>%
  group_by(asv_number) %>%
  mutate(chain = rep(seq(4), each = 1000)) %>%
  ungroup() %>%
  left_join(asv_key, by = "asv_id") %>%
  select(contains("asv"), contains("seqvar"), contains("blast"), best, chain, everything()) %>%
  identity() -> pilot_posterior_betas

pilot_posterior_betas

map(.x = as.list(structure(.Names = paste0("chain_",seq(4)), .Data = seq(4))),
    .f = ~ filter(pilot_posterior_betas, chain == .x)) %>%
  map2(.x = ., .y = names(.), .f = ~ write_rds(x = .x, path = paste0("./models/jsdm/pilot_posterior_betas_",.y,".rds.gz"), compress = "gz"))


#' convergence diagnostics for posterior betas

pilot_posterior_betas %>%
  gather(key = "beta", value = "posterior_sample", -asv_number, -asv_id, -seqvar_id, -blast_family, -blast_genus, -best, -chain) %>%
  select(chain, beta, asv_number, posterior_sample) %>%
  mutate(chain = paste0("chain_",chain)) %>%
  group_by(chain, beta, asv_number) %>%
  mutate(iter = seq_along(posterior_sample)) %>%
  ungroup() %>%
  spread(key = chain, value = posterior_sample) %>%
  select(-iter) %>%
  group_by(beta, asv_number) %>%
  nest() %>%
  summarise(rhat = map_dbl(.x = data, .f = ~ rstan::Rhat(as.matrix(.x)))) %>%
  identity() -> pilot_beta_rhat

pilot_beta_rhat

pilot_beta_rhat %>%
  qplot(data = ., x = asv_number, y = rhat, color = beta, facets = ~ beta)



#' examine posterior betas

pilot_posterior_betas %>%
  gather(key = "beta", value = "posterior_draws", -asv_number, -asv_id, -seqvar_id, -blast_family, -blast_genus, -best, -chain) %>%
  group_by(asv_number, asv_id, seqvar_id, blast_family, blast_genus, best, beta) %>%
  summarise_at(.vars = vars(posterior_draws), .funs = list("hpdi_low" = ~ rethinking::HPDI(.x, prob = 0.95)[1],
                                                           "hpdi_high" = ~ rethinking::HPDI(.x, prob = 0.95)[2],
                                                           "pi_low" = ~ rethinking::PI(.x, prob = 0.95)[1],
                                                           "pi_high" = ~ rethinking::PI(.x, prob = 0.95)[2],
                                                           "median" = ~ quantile(.x, 0.5))) %>%
  ungroup() %>%
  identity() -> pilot_posterior_asv_summary

pilot_posterior_asv_summary

pilot_posterior_asv_summary %>%
  filter(grepl("Corynebacterium|Staphylococcus",blast_genus)) %>%
  ggplot(data = .) +
  geom_segment(aes(x = hpdi_low, xend = hpdi_high, y = seqvar_id, yend = seqvar_id)) +
  geom_point(aes(x = median, y = seqvar_id)) +
  facet_wrap(facets = ~ beta, scales = "free_x")



#' POSTERIOR RESIDUAL CORRELATION

pmod_jSDM_block_pilot %>%
  map(get_residual_cor) -> residcor_jSDM_block_pilot

residcor_jSDM_block_pilot %>%
  map(.x = ., .f = ~ pluck(.x, "cor.median") %>%
        as_tibble() %>%
        mutate(asv_x = colnames(y_pilot)) %>%
        gather(key = "asv_y", value = "cor", -asv_x) %>%
        mutate(asv_y = structure(.Data = colnames(y_pilot), .Names = paste0("V",seq(ncol(y_pilot))))[asv_y])
      ) %>%
  bind_rows(.id = "chain") -> resid_cor_pilot
resid_cor_pilot


resid_cor_pilot %>%
  write_csv("./models/jsdm/pilot_residual_cor.csv.gz")


resid_cor_pilot %>%
  group_by(asv_x,asv_y) %>%
  summarise(interchain_cor_variance = var(cor)) %>%
  arrange(desc(interchain_cor_variance))


resid_cor_pilot %>%
  qplot(data = ., x = asv_x, y = asv_y, geom = "tile", fill = cor) +
  facet_wrap(~ chain) +
  scale_fill_viridis_c(option = "magma", limits = c(-1,1)) +
  coord_equal() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 270, vjust = 0.5, hjust = 0))


resid_cor_pilot %>%
  filter(asv_x %in% pull(filter(asv_key, grepl("aureus", best)),asv_id)) %>%
  mutate(genus_y = structure(.Data = asv_key$blast_genus, .Names = asv_key$asv_id)[asv_y],
         asv_genus = paste0(genus_y,"\n(",gsub("asv_","",asv_y),")")) %>%
  arrange(cor) %>%
  print(n=100)







#' ######################################
#' 
#' jSDM valid data - parallel chains
#' 
#' ######################################

x_valid %>%
  select(#age,
         #copd,
         #asthma,
         #ild,
         #lymphleuk,
         #dm,
         #cirrhosis,
         #chf,
         contains("before0_7")) %>%
  #mutate(age = age = scale(age)[,1]) %>%
  mutate_if(.predicate = ~ is.numeric(.x), .funs = ~ scale(.x)[,1]) %>%
  mutate_if(.predicate = ~ is.logical(.x), .funs = ~ as.numeric(.x)) %>%
  identity() -> dx_valid
dx_valid


run_jsdm <- function(x,y,seed_number) {
  my_jsdm <- jSDM_probit_block (
    # Response variable 
    presence_site_sp = as.matrix(y), 
    # Explanatory variables 
    site_suitability = ~.,   
    site_data = as.data.frame(x), n_latent=2,
    # Chains
    burnin=20000, mcmc=5000, thin=5,
    # Starting values
    alpha_start=0, beta_start=0,
    lambda_start=0, W_start=0,
    V_alpha_start=1, 
    # Priors
    shape=0.5, rate=0.0005,
    mu_beta=0, V_beta=1.0E6,
    mu_lambda=0, V_lambda=10,
    # Various 
    seed=seed_number,
    verbose=1)
  return(my_jsdm)
}


parallel::mcMap(f = function(x,y,z) run_jsdm(x,y,z),
                y = list(y_valid, y_valid, y_valid, y_valid),
                x = list(dx_valid, dx_valid, dx_valid, dx_valid),
                z = as.list(seq(4)*16),
                mc.cores = 4
) -> pmod_jSDM_block_valid

names(pmod_jSDM_block_valid) <- paste0("chain_",seq_along(pmod_jSDM_block_valid))

str(pmod_jSDM_block_valid[1])

pmod_jSDM_block_valid %>%
  map2(.x = .,
       .y = as.list(paste0("./models/jsdm/pmod_jSDM_block_valid_",names(.),".rds.gz")),
       .f = ~ write_rds(x = .x, path = .y, compress = "gz")
  )

as.list(structure(.Data = c("chain_1","chain_2","chain_3","chain_4"), .Names = c("chain_1","chain_2","chain_3","chain_4"))) %>%
  map(.f = ~ read_rds(path = paste0("./models/jsdm/pmod_jSDM_block_valid_",.x,".rds.gz"))) -> pmod_jSDM_block_valid

str(pmod_jSDM_block_valid[1])



# pmod_jSDM_block_valid %>%
#   map2(.x = .,
#        .y = as.list(paste0("./models/jsdm/pmod_jSDM_block_valid_",names(.),".rds.gz")),
#        .f = ~ write_rds(x = .x, path = .y, compress = "gz")
#   )

#' POSTERIOR BETAS & DIAGNOSTICS

# map2(.f = ~ rbind(.x,.y),
#      .x = list(matrix(seq(5,8),nrow=2),matrix(seq(5,8),nrow=2)),
#      .y = list(matrix(seq(4),nrow=2),matrix(seq(4),nrow=2)))
# 
# pmap(.l = list(list(matrix(seq(5,8),nrow=2),matrix(seq(5,8),nrow=2)),
#                list(matrix(seq(4),nrow=2),matrix(seq(4),nrow=2)),
#                list(matrix(seq(9,12),nrow=2),matrix(seq(9,12),nrow=2))),
#      .f = rbind)

map(.x = pmod_jSDM_block_valid, .f = ~ pluck(.x, "mcmc.sp")) %>%
  pmap(.l = ., .f = rbind) %>%
  map(.x = ., .f = ~ as_tibble(.x)) %>%
  bind_rows(.id = "asv_number") %>%
  mutate(asv_id = structure(.Names = paste0("sp_",seq(ncol(y_valid))), .Data = colnames(y_valid))[asv_number]) %>%
  group_by(asv_number) %>%
  mutate(chain = rep(seq(4), each = 1000)) %>%
  ungroup() %>%
  left_join(asv_key, by = "asv_id") %>%
  select(contains("asv"), contains("seqvar"), contains("blast"), best, chain, everything()) %>%
  identity() -> valid_posterior_betas

valid_posterior_betas

map(.x = as.list(structure(.Names = paste0("chain_",seq(4)), .Data = seq(4))),
    .f = ~ filter(valid_posterior_betas, chain == .x)) %>%
  map2(.x = ., .y = names(.), .f = ~ write_rds(x = .x, path = paste0("./models/jsdm/valid_posterior_betas_",.y,".rds.gz"), compress = "gz"))




#' convergence diagnostics for posterior betas

valid_posterior_betas %>%
  gather(key = "beta", value = "posterior_sample", -asv_number, -asv_id, -seqvar_id, -blast_family, -blast_genus, -best, -chain) %>%
  select(chain, beta, asv_number, posterior_sample) %>%
  mutate(chain = paste0("chain_",chain)) %>%
  group_by(chain, beta, asv_number) %>%
  mutate(iter = seq_along(posterior_sample)) %>%
  ungroup() %>%
  spread(key = chain, value = posterior_sample) %>%
  select(-iter) %>%
  group_by(beta, asv_number) %>%
  nest() %>%
  summarise(rhat = map_dbl(.x = data, .f = ~ rstan::Rhat(as.matrix(.x)))) %>%
  identity() -> valid_beta_rhat

valid_beta_rhat

valid_beta_rhat %>%
  qplot(data = ., x = asv_number, y = rhat, color = beta, facets = ~ beta)



#' examine posterior betas

valid_posterior_betas %>%
  gather(key = "beta", value = "posterior_draws", -asv_number, -asv_id, -seqvar_id, -blast_family, -blast_genus, -best, -chain) %>%
  group_by(asv_number, asv_id, seqvar_id, blast_family, blast_genus, best, beta) %>%
  summarise_at(.vars = vars(posterior_draws), .funs = list("hpdi_low" = ~ rethinking::HPDI(.x, prob = 0.95)[1],
                                                           "hpdi_high" = ~ rethinking::HPDI(.x, prob = 0.95)[2],
                                                           "pi_low" = ~ rethinking::PI(.x, prob = 0.95)[1],
                                                           "pi_high" = ~ rethinking::PI(.x, prob = 0.95)[2],
                                                           "median" = ~ quantile(.x, 0.5))) %>%
  ungroup() %>%
  identity() -> valid_posterior_asv_summary

valid_posterior_asv_summary

valid_posterior_asv_summary %>%
  filter(grepl("Corynebacterium|Staphylococcus",blast_genus)) %>%
  ggplot(data = .) +
  geom_segment(aes(x = hpdi_low, xend = hpdi_high, y = seqvar_id, yend = seqvar_id)) +
  geom_point(aes(x = median, y = seqvar_id)) +
  facet_wrap(facets = ~ beta, scales = "free_x")



#' POSTERIOR RESIDUAL CORRELATION

pmod_jSDM_block_valid %>%
  map(get_residual_cor) -> residcor_jSDM_block_valid

residcor_jSDM_block_valid %>%
  map(.x = ., .f = ~ pluck(.x, "cor.median") %>%
        as_tibble() %>%
        mutate(asv_x = colnames(y_valid)) %>%
        gather(key = "asv_y", value = "cor", -asv_x) %>%
        mutate(asv_y = structure(.Data = colnames(y_valid), .Names = paste0("V",seq(ncol(y_valid))))[asv_y])
  ) %>%
  bind_rows(.id = "chain") -> resid_cor_valid
resid_cor_valid


resid_cor_valid %>%
  write_csv("./models/jsdm/valid_residual_cor.csv.gz")


resid_cor_valid %>%
  group_by(asv_x,asv_y) %>%
  summarise(interchain_cor_variance = var(cor)) %>%
  arrange(desc(interchain_cor_variance))


resid_cor_valid %>%
  qplot(data = ., x = asv_x, y = asv_y, geom = "tile", fill = cor) +
  facet_wrap(~ chain) +
  scale_fill_viridis_c(option = "magma", limits = c(-1,1)) +
  coord_equal() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 270, vjust = 0.5, hjust = 0))


resid_cor_valid %>%
  filter(asv_x %in% pull(filter(asv_key, grepl("aureus", best)),asv_id)) %>%
  mutate(genus_y = structure(.Data = asv_key$blast_genus, .Names = asv_key$asv_id)[asv_y],
         asv_genus = paste0(genus_y,"\n(",gsub("asv_","",asv_y),")")) %>%
  arrange(cor) %>%
  print(n=100)














#' ######################################
#' 
#' jSDM pooled data - parallel chains
#' 
#' ######################################

x_pooled %>%
  select(#age,
    #copd,
    #asthma,
    #ild,
    #lymphleuk,
    #dm,
    #cirrhosis,
    #chf,
    contains("before0_7")) %>%
  #mutate(age = age = scale(age)[,1]) %>%
  mutate_if(.predicate = ~ is.numeric(.x), .funs = ~ scale(.x)[,1]) %>%
  mutate_if(.predicate = ~ is.logical(.x), .funs = ~ as.numeric(.x)) %>%
  identity() -> dx_pooled
dx_pooled


run_jsdm <- function(x,y,seed_number) {
  my_jsdm <- jSDM_probit_block (
    # Response variable 
    presence_site_sp = as.matrix(y), 
    # Explanatory variables 
    site_suitability = ~.,   
    site_data = as.data.frame(x), n_latent=2,
    # Chains
    burnin=20000, mcmc=5000, thin=5,
    # Starting values
    alpha_start=0, beta_start=0,
    lambda_start=0, W_start=0,
    V_alpha_start=1, 
    # Priors
    shape=0.5, rate=0.0005,
    mu_beta=0, V_beta=1.0E6,
    mu_lambda=0, V_lambda=10,
    # Various 
    seed=seed_number,
    verbose=1)
  return(my_jsdm)
}


parallel::mcMap(f = function(x,y,z) run_jsdm(x,y,z),
                y = list(y_pooled, y_pooled, y_pooled, y_pooled),
                x = list(dx_pooled, dx_pooled, dx_pooled, dx_pooled),
                z = as.list(seq(4)*16),
                mc.cores = 4
) -> pmod_jSDM_block_pooled

names(pmod_jSDM_block_pooled) <- paste0("chain_",seq_along(pmod_jSDM_block_pooled))

str(pmod_jSDM_block_pooled[1])

pmod_jSDM_block_pooled %>%
  map2(.x = .,
       .y = as.list(paste0("./models/jsdm/pmod_jSDM_block_pooled_",names(.),".rds.gz")),
       .f = ~ write_rds(x = .x, path = .y, compress = "gz")
  )

as.list(structure(.Data = c("chain_1","chain_2","chain_3","chain_4"), .Names = c("chain_1","chain_2","chain_3","chain_4"))) %>%
  map(.f = ~ read_rds(path = paste0("./models/jsdm/pmod_jSDM_block_pooled_",.x,".rds.gz"))) -> pmod_jSDM_block_pooled

str(pmod_jSDM_block_pooled[1])



# pmod_jSDM_block_pooled %>%
#   map2(.x = .,
#        .y = as.list(paste0("./models/jsdm/pmod_jSDM_block_pooled_",names(.),".rds.gz")),
#        .f = ~ write_rds(x = .x, path = .y, compress = "gz")
#   )

#' POSTERIOR BETAS & DIAGNOSTICS

# map2(.f = ~ rbind(.x,.y),
#      .x = list(matrix(seq(5,8),nrow=2),matrix(seq(5,8),nrow=2)),
#      .y = list(matrix(seq(4),nrow=2),matrix(seq(4),nrow=2)))
# 
# pmap(.l = list(list(matrix(seq(5,8),nrow=2),matrix(seq(5,8),nrow=2)),
#                list(matrix(seq(4),nrow=2),matrix(seq(4),nrow=2)),
#                list(matrix(seq(9,12),nrow=2),matrix(seq(9,12),nrow=2))),
#      .f = rbind)

map(.x = pmod_jSDM_block_pooled, .f = ~ pluck(.x, "mcmc.sp")) %>%
  pmap(.l = ., .f = rbind) %>%
  map(.x = ., .f = ~ as_tibble(.x)) %>%
  bind_rows(.id = "asv_number") %>%
  mutate(asv_id = structure(.Names = paste0("sp_",seq(ncol(y_pooled))), .Data = colnames(y_pooled))[asv_number]) %>%
  group_by(asv_number) %>%
  mutate(chain = rep(seq(4), each = 1000)) %>%
  ungroup() %>%
  left_join(asv_key, by = "asv_id") %>%
  select(contains("asv"), contains("seqvar"), contains("blast"), best, chain, everything()) %>%
  identity() -> pooled_posterior_betas

pooled_posterior_betas

map(.x = as.list(structure(.Names = paste0("chain_",seq(4)), .Data = seq(4))),
    .f = ~ filter(pooled_posterior_betas, chain == .x)) %>%
  map2(.x = ., .y = names(.), .f = ~ write_rds(x = .x, path = paste0("./models/jsdm/pooled_posterior_betas_",.y,".rds.gz"), compress = "gz"))




#' convergence diagnostics for posterior betas

pooled_posterior_betas %>%
  gather(key = "beta", value = "posterior_sample", -asv_number, -asv_id, -seqvar_id, -blast_family, -blast_genus, -best, -chain) %>%
  select(chain, beta, asv_number, posterior_sample) %>%
  mutate(chain = paste0("chain_",chain)) %>%
  group_by(chain, beta, asv_number) %>%
  mutate(iter = seq_along(posterior_sample)) %>%
  ungroup() %>%
  spread(key = chain, value = posterior_sample) %>%
  select(-iter) %>%
  group_by(beta, asv_number) %>%
  nest() %>%
  summarise(rhat = map_dbl(.x = data, .f = ~ rstan::Rhat(as.matrix(.x)))) %>%
  identity() -> pooled_beta_rhat

pooled_beta_rhat

pooled_beta_rhat %>%
  qplot(data = ., x = asv_number, y = rhat, color = beta, facets = ~ beta)



#' examine posterior betas

pooled_posterior_betas %>%
  gather(key = "beta", value = "posterior_draws", -asv_number, -asv_id, -seqvar_id, -blast_family, -blast_genus, -best, -chain) %>%
  group_by(asv_number, asv_id, seqvar_id, blast_family, blast_genus, best, beta) %>%
  summarise_at(.vars = vars(posterior_draws), .funs = list("hpdi_low" = ~ rethinking::HPDI(.x, prob = 0.95)[1],
                                                           "hpdi_high" = ~ rethinking::HPDI(.x, prob = 0.95)[2],
                                                           "pi_low" = ~ rethinking::PI(.x, prob = 0.95)[1],
                                                           "pi_high" = ~ rethinking::PI(.x, prob = 0.95)[2],
                                                           "median" = ~ quantile(.x, 0.5))) %>%
  ungroup() %>%
  identity() -> pooled_posterior_asv_summary

pooled_posterior_asv_summary

pooled_posterior_asv_summary %>%
  filter(grepl("Corynebacterium|Staphylococcus",blast_genus)) %>%
  ggplot(data = .) +
  geom_segment(aes(x = hpdi_low, xend = hpdi_high, y = seqvar_id, yend = seqvar_id)) +
  geom_point(aes(x = median, y = seqvar_id)) +
  facet_wrap(facets = ~ beta, scales = "free_x")



#' POSTERIOR RESIDUAL CORRELATION

pmod_jSDM_block_pooled %>%
  map(get_residual_cor) -> residcor_jSDM_block_pooled

residcor_jSDM_block_pooled %>%
  map(.x = ., .f = ~ pluck(.x, "cor.median") %>%
        as_tibble() %>%
        mutate(asv_x = colnames(y_pooled)) %>%
        gather(key = "asv_y", value = "cor", -asv_x) %>%
        mutate(asv_y = structure(.Data = colnames(y_pooled), .Names = paste0("V",seq(ncol(y_pooled))))[asv_y])
  ) %>%
  bind_rows(.id = "chain") -> resid_cor_pooled
resid_cor_pooled


resid_cor_pooled %>%
  write_csv("./models/jsdm/pooled_residual_cor.csv.gz")


resid_cor_pooled %>%
  group_by(asv_x,asv_y) %>%
  summarise(interchain_cor_variance = var(cor)) %>%
  arrange(desc(interchain_cor_variance))


resid_cor_pooled %>%
  qplot(data = ., x = asv_x, y = asv_y, geom = "tile", fill = cor) +
  facet_wrap(~ chain) +
  scale_fill_viridis_c(option = "magma", limits = c(-1,1)) +
  coord_equal() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 270, vjust = 0.5, hjust = 0))


resid_cor_pooled %>%
  filter(asv_x %in% pull(filter(asv_key, grepl("aureus", best)),asv_id)) %>%
  mutate(genus_y = structure(.Data = asv_key$blast_genus, .Names = asv_key$asv_id)[asv_y],
         asv_genus = paste0(genus_y,"\n(",gsub("asv_","",asv_y),")")) %>%
  arrange(cor) %>%
  print(n=100)






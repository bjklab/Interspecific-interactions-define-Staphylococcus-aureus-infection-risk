#' ########################################
#' CLR binomial model plots
#' ########################################
#' 
#' depends: output from m_clr_binomial_models_replacement.R script
#' 



#' ########################################
#' COMBINE CLR COHORTS: LRTI ~ Sa ASVs
#' ########################################

CLR_lrti_binomial_Sa_cohorts <- bind_rows(mutate(read_csv("./models/binomial/CLRR_staph_lrti_model-AN_Sa_ASVs_only-reghs_scaled_reads_POOLED_seqvar.csv"), cohort = 'Pooled'),
                                      mutate(read_csv("./models/binomial/CLRR_staph_lrti_model-AN_Sa_ASVs_only-reghs_scaled_reads_PILOT_seqvar.csv"), cohort = 'Pilot'),
                                      mutate(read_csv("./models/binomial/CLRR_staph_lrti_model-AN_Sa_ASVs_only-reghs_scaled_reads_VALIDATION_seqvar.csv"), cohort = 'Validation')) %>%
  filter(grepl('beta',parameter))
CLR_lrti_binomial_Sa_cohorts


CLR_lrti_binomial_Sa_cohorts %>%
  mutate(cohort = factor(cohort, levels = c('Pilot','Validation','Pooled')),
         typeS10 = typeS_10_neg | typeS_10_pos,
         typeS5 = typeS_5_neg | typeS_5_pos,
         typeS25 = typeS_25_neg | typeS_25_pos,
         typeScol = factor(ifelse(typeS5, "< 0.05", ifelse(typeS10, "< 0.1", ifelse(typeS25, "< 0.25", ">= 0.25"))), levels = c("< 0.05", "< 0.1", "< 0.25", ">= 0.25"))) %>% #count(typeScol)
  group_by(seqvar_id) %>%
  mutate(checktypeS = sum(typeS25, na.rm = TRUE)) %>%
  ungroup() %>%
  #filter(checktypeS > 0) %>%
  mutate(blast_genus = gsub(" 1","",blast_genus),
         seqvar_lab = paste0(blast_genus,"\n(",seqvar_id,")")) %>%
  mutate_at(.vars = vars(mean, se_mean, sd, contains("%")), .funs = ~ exp(.x)) %>% # convert posterior samples to OR scale
  ggplot(data = .) +
  geom_vline(xintercept = 1, linetype = 2) +
  geom_segment(aes(x = `25%`, xend = `75%`, y = seqvar_id, yend = seqvar_id, colour = typeScol)) +
  geom_point(aes(x = `50%`, y = seqvar_id, colour = typeScol)) +
  facet_wrap(~cohort, nrow = 1) +
  #scale_colour_viridis_d(option = "C") +
  scale_colour_manual(values = c('< 0.05' = '#cb181d', '< 0.1' = '#fb6a4a', '< 0.25' = '#fcae91', '>= 0.25' = 'darkgrey'), drop = FALSE) +
  #scale_x_continuous(breaks = seq(from = -1, to = 5, by = 1), trans = "log2") +
  scale_x_continuous(breaks = c(0.5, 1, 2, 3, 4, 5), labels = c("0.5", "1", "2", "3", "4", "5"), trans = "log10") +
  #scale_x_log10() +
  theme_bw() +
  theme(strip.background = element_blank(),
        axis.text.x = element_text(colour = 'black'),
        axis.text.y = element_text(colour = 'black', size = 6),
        plot.title = element_text(hjust = 0),
        plot.subtitle = element_text(hjust = 0, size = 10),
        legend.position = 'bottom') +
  labs(x = expression(paste("Odds Ratio ",italic("S. aureus")," LRTI ","(posterior median and 50% CI)")),
       y = expression(paste(italic("S. aureus"), " Anterior Nares ASVs (CLR transformed)")),
       colour = "Type S Error") -> p_CLR_lrti_binomial_Sa_cohorts

p_CLR_lrti_binomial_Sa_cohorts

# ggsave(plot = p_CLR_lrti_binomial_Sa_cohorts, filename = "./figs/supp/p_CLR_lrti_binomial_Sa_cohorts.png", height = 4, width = 8, units = "in", dpi = 600)
# ggsave(plot = p_CLR_lrti_binomial_Sa_cohorts, filename = "./figs/supp/p_CLR_lrti_binomial_Sa_cohorts.pdf", height = 4, width = 8, units = "in")


# test version with viridis plasma color scale:
p_CLR_lrti_binomial_Sa_cohorts +
  scale_colour_manual(values = c('< 0.05' = '#7301A8FF', '< 0.1' = '#BD3786FF', '< 0.25' = '#ED7953FF', '>= 0.25' = 'darkgrey'), drop = FALSE)






#' ########################################
#' COMBINE CLR COHORTS: LRTI ~ Sa + other microbiome ASVs + antibiotics (7-days prior)
#' ########################################

CLR_lrti_binomial_Sa_ASV_abx_cohorts <- bind_rows(mutate(read_csv("./models/binomial/CLRR_staph_lrti_model-AN_Sa_ASVs_microbiome_abx-reghs_scaled_reads_POOLED_seqvar.csv"), cohort = 'Pooled'),
                                              mutate(read_csv("./models/binomial/CLRR_staph_lrti_model-AN_Sa_ASVs_microbiome_abx-reghs_scaled_reads_PILOT_seqvar.csv"), cohort = 'Pilot'),
                                              mutate(read_csv("./models/binomial/CLRR_staph_lrti_model-AN_Sa_ASVs_microbiome_abx-reghs_scaled_reads_VALIDATION_seqvar.csv"), cohort = 'Validation')) %>%
  filter(grepl('beta',parameter))
CLR_lrti_binomial_Sa_ASV_abx_cohorts


CLR_lrti_binomial_Sa_ASV_abx_cohorts %>%
  mutate(cohort = factor(cohort, levels = c('Pilot','Validation','Pooled')),
         typeS10 = typeS_10_neg | typeS_10_pos,
         typeS5 = typeS_5_neg | typeS_5_pos,
         typeS25 = typeS_25_neg | typeS_25_pos,
         typeScol = factor(ifelse(typeS5, "< 0.05", ifelse(typeS10, "< 0.1", ifelse(typeS25, "< 0.25", ">= 0.25"))), levels = c("< 0.05", "< 0.1", "< 0.25", ">= 0.25"))) %>% #count(typeScol)
  group_by(seqvar_id) %>%
  mutate(checktypeS = sum(typeS10) > 0) %>%
  ungroup() %>%
  filter(checktypeS) %>% # filter to type S error < 10%
  mutate(blast_genus = gsub(" 1","",blast_genus),
         seqvar_lab = paste0(blast_genus,"\n(",seqvar_id,")")) %>%
  mutate_at(.vars = vars(mean, se_mean, sd, contains("%")), .funs = ~ exp(.x)) %>% # convert posterior samples to OR scale
  ggplot(data = .) +
  geom_vline(xintercept = 1, linetype = 2) +
  geom_segment(aes(x = `25%`, xend = `75%`, y = seqvar_lab, yend = seqvar_lab, colour = typeScol)) +
  geom_point(aes(x = `50%`, y = seqvar_lab, colour = typeScol)) +
  facet_wrap(~cohort, nrow = 1) +
  scale_colour_manual(values = c('< 0.05' = '#cb181d', '< 0.1' = '#fb6a4a', '< 0.25' = '#fcae91', '>= 0.25' = 'darkgrey'), drop = FALSE) +
  #scale_x_continuous(breaks = seq(from = -1, to = 5, by = 1), trans = "log2") +
  scale_x_continuous(breaks = c(0.5, 1, 2, 3, 4, 5), labels = c("0.5", "1", "2", "3", "4", "5"), trans = "log10") +
  #scale_x_log10() +
  theme_bw() +
  theme(strip.background = element_blank(),
        axis.text.x = element_text(colour = 'black'),
        axis.text.y = element_text(colour = 'black', size = 6),
        plot.title = element_text(hjust = -0.5),
        plot.subtitle = element_text(hjust = 0, size = 10),
        legend.position = "bottom") +
  labs(x = expression(paste("Odds Ratio ",italic("S. aureus")," LRTI ","(posterior median and 50% CI)")),
       y = expression(paste("Non-",italic("S. aureus")," Anterior Nares ASVs (CLR transformed)")),
       colour = "Type S Error") -> p_CLR_lrti_binomial_Sa_ASV_abx_cohorts

p_CLR_lrti_binomial_Sa_ASV_abx_cohorts

# ggsave(plot = p_CLR_lrti_binomial_Sa_ASV_abx_cohorts, filename = "./figs/supp/p_CLR_lrti_binomial_Sa_ASV_abx_cohorts.png", height = 5, width = 8, units = "in", dpi = 600)
# ggsave(plot = p_CLR_lrti_binomial_Sa_ASV_abx_cohorts, filename = "./figs/supp/p_CLR_lrti_binomial_Sa_ASV_abx_cohorts.pdf", height = 5, width = 8, units = "in")


# test version with viridis plasma color scale:
p_CLR_lrti_binomial_Sa_ASV_abx_cohorts +
  scale_colour_manual(values = c('< 0.05' = '#7301A8FF', '< 0.1' = '#BD3786FF', '< 0.25' = '#ED7953FF', '>= 0.25' = 'darkgrey'), drop = FALSE)





#' #######################################
#' FIGURE: compose ORIGINAL two plots for Figure 2 (using cowplot package)
#' #######################################

ggdraw(plot_grid(p_CLR_lrti_binomial_Sa_cohorts + theme(legend.position = "none"), p_CLR_lrti_binomial_Sa_ASV_abx_cohorts, nrow = 2, rel_heights = c(4,16), labels = c("A","B"))) -> p_CLR_lrti_model_combined
p_CLR_lrti_model_combined

# ggsave(plot = p_CLR_lrti_model_combined, filename = "./figs/supp/p_CLR_lrti_model_combined.png", height = 20, width = 8, units = "in", dpi = 600)
# ggsave(plot = p_CLR_lrti_model_combined, filename = "./figs/supp/p_CLR_lrti_model_combined.pdf", height = 20, width = 8, units = "in")






#' ################################################
#' ################################################
#' integrate CLR-transformed and untransformed data
#' ################################################
#' ################################################

#' ########################################
#' COMBINE COHORTS: LRTI ~ Sa ASVs
#' ########################################

lrti_binomial_Sa_cohorts <- bind_rows(mutate(read_csv("./models/binomial/staph_lrti_model-AN_Sa_ASVs_only-reghs_scaled_reads_POOLED_seqvar.csv"), cohort = 'Pooled'),
                                      mutate(read_csv("./models/binomial/staph_lrti_model-AN_Sa_ASVs_only-reghs_scaled_reads_PILOT_seqvar.csv"), cohort = 'Pilot'),
                                      mutate(read_csv("./models/binomial/staph_lrti_model-AN_Sa_ASVs_only-reghs_scaled_reads_VALIDATION_seqvar.csv"), cohort = 'Validation')) %>%
  filter(grepl('beta',parameter))
lrti_binomial_Sa_cohorts


lrti_binomial_Sa_cohorts %>%
  mutate(cohort = factor(cohort, levels = c('Pilot','Validation','Pooled')),
         typeS10 = typeS_10_neg | typeS_10_pos,
         typeS5 = typeS_5_neg | typeS_5_pos,
         typeS25 = typeS_25_neg | typeS_25_pos,
         typeScol = factor(ifelse(typeS5, "< 0.05", ifelse(typeS10, "< 0.1", ifelse(typeS25, "< 0.25", ">= 0.25"))), levels = c("< 0.05", "< 0.1", "< 0.25", ">= 0.25"))) %>% #count(typeScol)
  group_by(seqvar_id) %>%
  mutate(checktypeS = sum(typeS25, na.rm = TRUE)) %>%
  ungroup() %>%
  #filter(checktypeS > 0) %>%
  mutate(blast_genus = gsub(" 1","",blast_genus),
         seqvar_lab = paste0(blast_genus,"\n(",seqvar_id,")")) %>%
  mutate_at(.vars = vars(mean, se_mean, sd, contains("%")), .funs = ~ exp(.x)) %>% # convert posterior samples to OR scale
  ggplot(data = .) +
  geom_vline(xintercept = 1, linetype = 2) +
  geom_segment(aes(x = `25%`, xend = `75%`, y = seqvar_id, yend = seqvar_id, colour = typeScol)) +
  geom_point(aes(x = `50%`, y = seqvar_id, colour = typeScol)) +
  facet_wrap(~cohort, nrow = 1) +
  #scale_colour_viridis_d(option = "C") +
  scale_colour_manual(values = c('< 0.05' = '#cb181d', '< 0.1' = '#fb6a4a', '< 0.25' = '#fcae91', '>= 0.25' = 'darkgrey'), drop = FALSE) +
  #scale_x_continuous(breaks = seq(from = -1, to = 5, by = 1), trans = "log2") +
  scale_x_continuous(breaks = c(0.5, 1, 2, 3, 4, 5), labels = c("0.5", "1", "2", "3", "4", "5"), trans = "log10") +
  #scale_x_log10() +
  theme_bw() +
  theme(strip.background = element_blank(),
        axis.text.x = element_text(colour = 'black'),
        axis.text.y = element_text(colour = 'black', size = 6),
        plot.title = element_text(hjust = 0),
        plot.subtitle = element_text(hjust = 0, size = 10),
        legend.position = 'bottom') +
  labs(x = expression(paste("Odds Ratio ",italic("S. aureus")," LRTI ","(posterior median and 50% CI)")),
       y = expression(paste(italic("S. aureus"), " Anterior Nares ASVs")),
       colour = "Type S Error") -> p_lrti_binomial_Sa_cohorts

p_lrti_binomial_Sa_cohorts

# test version with viridis plasma color scale:
p_lrti_binomial_Sa_cohorts +
  scale_colour_manual(values = c('< 0.05' = '#7301A8FF', '< 0.1' = '#BD3786FF', '< 0.25' = '#ED7953FF', '>= 0.25' = 'darkgrey'), drop = FALSE)





lrti_binomial_Sa_ASV_abx_cohorts <- bind_rows(mutate(read_csv("./models/binomial/staph_lrti_model-AN_Sa_ASVs_microbiome_abx-reghs_scaled_reads_POOLED_seqvar.csv"), cohort = 'Pooled'),
                                              mutate(read_csv("./models/binomial/staph_lrti_model-AN_Sa_ASVs_microbiome_abx-reghs_scaled_reads_PILOT_seqvar.csv"), cohort = 'Pilot'),
                                              mutate(read_csv("./models/binomial/staph_lrti_model-AN_Sa_ASVs_microbiome_abx-reghs_scaled_reads_VALIDATION_seqvar.csv"), cohort = 'Validation')) %>%
  filter(grepl('beta',parameter))
lrti_binomial_Sa_ASV_abx_cohorts


lrti_binomial_Sa_ASV_abx_cohorts %>%
  mutate(cohort = factor(cohort, levels = c('Pilot','Validation','Pooled')),
         typeS10 = typeS_10_neg | typeS_10_pos,
         typeS5 = typeS_5_neg | typeS_5_pos,
         typeS25 = typeS_25_neg | typeS_25_pos,
         typeScol = factor(ifelse(typeS5, "< 0.05", ifelse(typeS10, "< 0.1", ifelse(typeS25, "< 0.25", ">= 0.25"))), levels = c("< 0.05", "< 0.1", "< 0.25", ">= 0.25"))) %>% #count(typeScol)
  group_by(seqvar_id) %>%
  mutate(checktypeS = sum(typeS10) > 0) %>%
  ungroup() %>%
  filter(checktypeS) %>% # filter to type S error < 10%
  mutate(blast_genus = gsub(" 1","",blast_genus),
         seqvar_lab = paste0(blast_genus,"\n(",seqvar_id,")")) %>%
  mutate_at(.vars = vars(mean, se_mean, sd, contains("%")), .funs = ~ exp(.x)) %>% # convert posterior samples to OR scale
  ggplot(data = .) +
  geom_vline(xintercept = 1, linetype = 2) +
  geom_segment(aes(x = `25%`, xend = `75%`, y = seqvar_lab, yend = seqvar_lab, colour = typeScol)) +
  geom_point(aes(x = `50%`, y = seqvar_lab, colour = typeScol)) +
  facet_wrap(~cohort, nrow = 1) +
  scale_colour_manual(values = c('< 0.05' = '#cb181d', '< 0.1' = '#fb6a4a', '< 0.25' = '#fcae91', '>= 0.25' = 'darkgrey'), drop = FALSE) +
  #scale_x_continuous(breaks = seq(from = -1, to = 5, by = 1), trans = "log2") +
  scale_x_continuous(breaks = c(0.5, 1, 2, 3, 4, 5), labels = c("0.5", "1", "2", "3", "4", "5"), trans = "log10") +
  #scale_x_log10() +
  theme_bw() +
  theme(strip.background = element_blank(),
        axis.text.x = element_text(colour = 'black'),
        axis.text.y = element_text(colour = 'black', size = 6),
        plot.title = element_text(hjust = -0.5),
        plot.subtitle = element_text(hjust = 0, size = 10),
        legend.position = "bottom") +
  labs(x = expression(paste("Odds Ratio ",italic("S. aureus")," LRTI ","(posterior median and 50% CI)")),
       y = expression(paste("Non-",italic("S. aureus")," Anterior Nares ASVs")),
       colour = "Type S Error") -> p_lrti_binomial_Sa_ASV_abx_cohorts

p_lrti_binomial_Sa_ASV_abx_cohorts

# test version with viridis plasma color scale:
p_lrti_binomial_Sa_ASV_abx_cohorts +
  scale_colour_manual(values = c('< 0.05' = '#7301A8FF', '< 0.1' = '#BD3786FF', '< 0.25' = '#ED7953FF', '>= 0.25' = 'darkgrey'), drop = FALSE)

# which ASVs included in binomial model plots

lrti_binomial_Sa_cohorts %>%
  pull(seqvar_id) %>%
  unique() -> included_asv_sa_binomial
included_asv_sa_binomial


lrti_binomial_Sa_ASV_abx_cohorts %>%
  mutate(cohort = factor(cohort, levels = c('Pilot','Validation','Pooled')),
         typeS10 = typeS_10_neg | typeS_10_pos,
         typeS5 = typeS_5_neg | typeS_5_pos,
         typeS25 = typeS_25_neg | typeS_25_pos,
         typeScol = factor(ifelse(typeS5, "< 0.05", ifelse(typeS10, "< 0.1", ifelse(typeS25, "< 0.25", ">= 0.25"))), levels = c("< 0.05", "< 0.1", "< 0.25", ">= 0.25"))) %>% #count(typeScol)
  group_by(seqvar_id) %>%
  mutate(checktypeS = sum(typeS10) > 0) %>%
  ungroup() %>%
  filter(checktypeS) %>% # filter to type S error < 10%
  pull(seqvar_id) %>%
  unique() -> included_asv_binomial
included_asv_binomial




#' ########################################
#' ########################################
#' SIX-PANEL PLOTS
#' Sa only and Sa-ASV-abx
#' #######################################
#' #######################################

#'#####################################################
#' Sa models: CRL transformed and untransformed
#' ####################################################

mutate(lrti_binomial_Sa_cohorts, CLR = FALSE) %>%
  bind_rows(mutate(CLR_lrti_binomial_Sa_cohorts, CLR = TRUE)) %>%
  mutate(CLR = ifelse(CLR, "CLR Transformed", "Scaled"),
         CLR = factor(CLR, levels = c("Scaled","CLR Transformed"))) %>%
  mutate(cohort = factor(cohort, levels = c('Pilot','Validation','Pooled')),
         typeS10 = typeS_10_neg | typeS_10_pos,
         typeS5 = typeS_5_neg | typeS_5_pos,
         typeS25 = typeS_25_neg | typeS_25_pos,
         typeScol = factor(ifelse(typeS5, "< 0.05", ifelse(typeS10, "< 0.1", ifelse(typeS25, "< 0.25", ">= 0.25"))), levels = c("< 0.05", "< 0.1", "< 0.25", ">= 0.25"))) %>% #count(typeScol)
  group_by(seqvar_id) %>%
  mutate(checktypeS = sum(typeS25) > 0) %>%
  ungroup() %>%
  #filter(checktypeS) %>% # filter to type S error < 10%
  mutate(blast_genus = gsub(" 1","",blast_genus),
         seqvar_lab = paste0(blast_genus,"\n(",seqvar_id,")")) %>%
  mutate_at(.vars = vars(mean, se_mean, sd, contains("%")), .funs = ~ exp(.x)) %>% # convert posterior samples to OR scale
  ggplot(data = .) +
  geom_vline(xintercept = 1, linetype = 2) +
  geom_segment(aes(x = `25%`, xend = `75%`, y = seqvar_lab, yend = seqvar_lab, colour = typeScol)) +
  geom_point(aes(x = `50%`, y = seqvar_lab, colour = typeScol)) +
  facet_wrap(~cohort + CLR, nrow = 1) +
  scale_colour_manual(values = c('< 0.05' = '#cb181d', '< 0.1' = '#fb6a4a', '< 0.25' = '#fcae91', '>= 0.25' = 'darkgrey'), drop = FALSE) +
  #scale_x_continuous(breaks = seq(from = -1, to = 5, by = 1), trans = "log2") +
  scale_x_continuous(breaks = c(0.5, 1, 2, 3, 4, 5), labels = c("0.5", "1", "2", "3", "4", "5"), trans = "log10") +
  #scale_x_log10() +
  theme_bw() +
  theme(strip.background = element_blank(),
        axis.text.x = element_text(colour = 'black'),
        axis.text.y = element_text(colour = 'black', size = 6),
        plot.title = element_text(hjust = -0.5),
        plot.subtitle = element_text(hjust = 0, size = 10),
        legend.position = "bottom") +
  labs(x = expression(paste("Odds Ratio ",italic("S. aureus")," LRTI ","(posterior median and 50% CI)")),
       y = expression(paste(italic("S. aureus")," Anterior Nares ASVs")),
       colour = "Type S Error") -> p_lrti_binomial_Sa_CLR_off_on

p_lrti_binomial_Sa_CLR_off_on

# ggsave(plot = p_lrti_binomial_Sa_CLR_off_on, filename = "./figs/supp/p_lrti_binomial_Sa_CLR_off_on.png", height = 5, width = 9, units = "in", dpi = 600)
# ggsave(plot = p_lrti_binomial_Sa_CLR_off_on, filename = "./figs/supp/p_lrti_binomial_Sa_CLR_off_on.pdf", height = 5, width = 9, units = "in")


# test version with viridis plasma color scale:
p_lrti_binomial_Sa_CLR_off_on +
  scale_colour_manual(values = c('< 0.05' = '#7301A8FF', '< 0.1' = '#BD3786FF', '< 0.25' = '#ED7953FF', '>= 0.25' = 'darkgrey'), drop = FALSE)




#'#####################################################
#' Sa-ASV-abx models: CRL transformed and untransformed
#' ####################################################

mutate(lrti_binomial_Sa_ASV_abx_cohorts, CLR = FALSE) %>%
  bind_rows(mutate(CLR_lrti_binomial_Sa_ASV_abx_cohorts, CLR = TRUE)) %>%
  mutate(CLR = ifelse(CLR, "CLR Transformed", "Scaled"),
         CLR = factor(CLR, levels = c("Scaled","CLR Transformed"))) %>%
  mutate(cohort = factor(cohort, levels = c('Pilot','Validation','Pooled')),
         typeS10 = typeS_10_neg | typeS_10_pos,
         typeS5 = typeS_5_neg | typeS_5_pos,
         typeS25 = typeS_25_neg | typeS_25_pos,
         typeScol = factor(ifelse(typeS5, "< 0.05", ifelse(typeS10, "< 0.1", ifelse(typeS25, "< 0.25", ">= 0.25"))), levels = c("< 0.05", "< 0.1", "< 0.25", ">= 0.25"))) %>% #count(typeScol)
  group_by(seqvar_id) %>%
  mutate(checktypeS = sum(typeS10) > 0) %>%
  ungroup() %>%
  filter(checktypeS) %>% # filter to type S error < 10%
  filter(seqvar_id %in% included_asv_binomial) %>% # filter to seqvars for comparison with binomial model output
  mutate(blast_genus = gsub(" 1","",blast_genus),
         seqvar_lab = paste0(blast_genus,"\n(",seqvar_id,")")) %>%
  mutate_at(.vars = vars(mean, se_mean, sd, contains("%")), .funs = ~ exp(.x)) %>% # convert posterior samples to OR scale
  ggplot(data = .) +
  geom_vline(xintercept = 1, linetype = 2) +
  geom_segment(aes(x = `25%`, xend = `75%`, y = seqvar_lab, yend = seqvar_lab, colour = typeScol)) +
  geom_point(aes(x = `50%`, y = seqvar_lab, colour = typeScol)) +
  facet_wrap(~cohort + CLR, nrow = 1) +
  scale_colour_manual(values = c('< 0.05' = '#cb181d', '< 0.1' = '#fb6a4a', '< 0.25' = '#fcae91', '>= 0.25' = 'darkgrey'), drop = FALSE) +
  #scale_x_continuous(breaks = seq(from = -1, to = 5, by = 1), trans = "log2") +
  scale_x_continuous(breaks = c(0.5, 1, 2, 3, 4, 5), labels = c("0.5", "1", "2", "3", "4", "5"), trans = "log10") +
  #scale_x_log10() +
  theme_bw() +
  theme(strip.background = element_blank(),
        axis.text.x = element_text(colour = 'black'),
        axis.text.y = element_text(colour = 'black', size = 6),
        plot.title = element_text(hjust = -0.5),
        plot.subtitle = element_text(hjust = 0, size = 10),
        legend.position = "bottom") +
  labs(x = expression(paste("Odds Ratio ",italic("S. aureus")," LRTI ","(posterior median and 50% CI)")),
       y = expression(paste("Non-",italic("S. aureus")," Anterior Nares ASVs")),
       colour = "Type S Error") -> p_lrti_binomial_Sa_ASV_abx_CLR_off_on

p_lrti_binomial_Sa_ASV_abx_CLR_off_on

# ggsave(plot = p_lrti_binomial_Sa_ASV_abx_CLR_off_on, filename = "./figs/supp/p_lrti_binomial_Sa_ASV_abx_CLR_off_on.png", height = 8.5, width = 9, units = "in", dpi = 600)
# ggsave(plot = p_lrti_binomial_Sa_ASV_abx_CLR_off_on, filename = "./figs/supp/p_lrti_binomial_Sa_ASV_abx_CLR_off_on.pdf", height = 8.5, width = 9, units = "in")


# test version with viridis plasma color scale:
p_lrti_binomial_Sa_ASV_abx_CLR_off_on +
  scale_colour_manual(values = c('< 0.05' = '#7301A8FF', '< 0.1' = '#BD3786FF', '< 0.25' = '#ED7953FF', '>= 0.25' = 'darkgrey'), drop = FALSE)



#' #######################################
#' 12-panel FIGURE: compose CLR+/- plots (using cowplot package)
#' #######################################

ggdraw(plot_grid(p_lrti_binomial_Sa_CLR_off_on + theme(legend.position = "none"), p_lrti_binomial_Sa_ASV_abx_CLR_off_on, nrow = 2, rel_heights = c(3.5,5.5), labels = c("A","B"))) -> p_lrti_model_CLR_off_on_combined
p_lrti_model_CLR_off_on_combined

# ggsave(plot = p_lrti_model_CLR_off_on_combined, filename = "./figs/supp/p_lrti_model_CLR_off_on_combined.png", height = 9, width = 9, units = "in", dpi = 600)
# ggsave(plot = p_lrti_model_CLR_off_on_combined, filename = "./figs/supp/p_lrti_model_CLR_off_on_combined.svg", height = 9, width = 9, units = "in", system_fonts = list(sans = "Roboto"))
# ggsave(plot = p_lrti_model_CLR_off_on_combined, filename = "./figs/supp/p_lrti_model_CLR_off_on_combined.pdf", height = 9, width = 9, units = "in")





#
###
#####
###
#





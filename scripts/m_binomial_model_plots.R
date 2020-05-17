#' ########################################
#' binomial model plots
#' ########################################
#' 
#' depends: output from m_binomial_models.R script
#' 

library(tidyverse)
library(cowplot)

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
         typeScol = ifelse(typeS5, "< 0.05", ifelse(typeS10, "< 0.1", ifelse(typeS25, "< 0.25", ">= 0.25")))) %>% #count(typeScol)
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
  scale_colour_manual(values = c('< 0.05' = '#cb181d', '< 0.1' = '#fb6a4a', '< 0.25' = '#fcae91', '>= 0.25' = 'darkgrey')) +
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

# ggsave(plot = p_lrti_binomial_Sa_cohorts, filename = "./figs/main/p_lrti_binomial_Sa_cohorts.png", height = 4, width = 8, units = "in", dpi = 600)
# ggsave(plot = p_lrti_binomial_Sa_cohorts, filename = "./figs/main/p_lrti_binomial_Sa_cohorts.pdf", height = 4, width = 8, units = "in")


# test version with viridis plasma color scale:
p_lrti_binomial_Sa_cohorts +
  scale_colour_manual(values = c('< 0.05' = '#7301A8FF', '< 0.1' = '#BD3786FF', '< 0.25' = '#ED7953FF', '>= 0.25' = 'darkgrey'))





#' ########################################
#' COMBINE COHORTS: LRTI ~ Sa + other microbiome ASVs
#' ########################################
# 
# lrti_binomial_Sa_ASV_cohorts <- bind_rows(mutate(read_csv("./models/binomial/staph_lrti_model-AN_Sa_ASVs_microbiome-reghs_scaled_reads_POOLED_seqvar.csv"), cohort = 'Pooled'),
#                                           mutate(read_csv("./models/binomial/staph_lrti_model-AN_Sa_ASVs_microbiome-reghs_scaled_reads_PILOT_seqvar.csv"), cohort = 'Pilot'),
#                                           mutate(read_csv("./models/binomial/staph_lrti_model-AN_Sa_ASVs_microbiome-reghs_scaled_reads_VALIDATION_seqvar.csv"), cohort = 'Validation')) %>%
#   filter(grepl('beta',parameter))
# lrti_binomial_Sa_ASV_cohorts
# 
# 
# lrti_binomial_Sa_ASV_cohorts %>%
#   mutate(cohort = factor(cohort, levels = c('Pilot','Validation','Pooled')),
#          typeS10 = typeS_10_neg | typeS_10_pos,
#          typeS5 = typeS_5_neg | typeS_5_pos,
#          typeS25 = typeS_25_neg | typeS_25_pos,
#          typeScol = ifelse(typeS5, "< 0.05", ifelse(typeS10, "< 0.1", ifelse(typeS25, "< 0.25", ">= 0.25")))) %>% #count(typeScol)
#   group_by(seqvar_id) %>%
#   mutate(checktypeS = sum(typeS10) > 0) %>%
#   ungroup() %>%
#   filter(checktypeS) %>% # filter to at least one cohort with type S error < 10%
#   mutate(blast_genus = gsub(" 1","",blast_genus),
#          seqvar_lab = paste0(blast_genus,"\n(",seqvar_id,")")) %>%
#   mutate_at(.vars = vars(mean, se_mean, sd, contains("%")), .funs = ~ exp(.x)) %>% # convert posterior samples to OR scale
#   #
#   ggplot(data = .) +
#   geom_vline(xintercept = 1, linetype = 2) +
#   geom_segment(aes(x = `25%`, xend = `75%`, y = seqvar_lab, yend = seqvar_lab, colour = typeScol)) +
#   geom_point(aes(x = `50%`, y = seqvar_lab, colour = typeScol)) +
#   facet_wrap(~cohort, nrow = 1) +
#   scale_colour_manual(values = c('< 0.05' = '#cb181d', '< 0.1' = '#fb6a4a', '< 0.25' = '#fcae91', '>= 0.25' = 'darkgrey')) +
#   #scale_x_continuous(breaks = seq(from = -1, to = 5, by = 1), trans = "log2") +
#   scale_x_continuous(breaks = c(0.5, 1, 2, 3, 4, 5), labels = c("0.5", "1", "2", "3", "4", "5"), trans = "log10") +
#   #scale_x_log10() +
#   theme_bw() +
#   theme(strip.background = element_blank(),
#         axis.text.x = element_text(colour = 'black'),
#         axis.text.y = element_text(colour = 'black', size = 6),
#         plot.title = element_text(hjust = -0.5),
#         plot.subtitle = element_text(hjust = 0, size = 10),
#         legend.position = "bottom") +
#   labs(x = expression(paste("Odds Ratio ",italic("S. aureus")," LRTI ","(posterior median and 50% CI)")),
#        y = expression(paste("Non-",italic("S. aureus")," Anterior Nares ASVs")),
#        colour = "Type S Error") -> p_lrti_binomial_Sa_ASV_cohorts
# 
# p_lrti_binomial_Sa_ASV_cohorts
# 
# ggsave(plot = p_lrti_binomial_Sa_ASV_cohorts, filename = "./figs/main/p_lrti_binomial_Sa_ASV_cohorts.png", height = 5, width = 8, units = "in", dpi = 600)
# ggsave(plot = p_lrti_binomial_Sa_ASV_cohorts, filename = "./figs/main/p_lrti_binomial_Sa_ASV_cohorts.pdf", height = 5, width = 8, units = "in")





#' ########################################
#' COMBINE COHORTS: LRTI ~ Sa + other microbiome ASVs + antibiotics (7-days prior)
#' ########################################

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
         typeScol = ifelse(typeS5, "< 0.05", ifelse(typeS10, "< 0.1", ifelse(typeS25, "< 0.25", ">= 0.25")))) %>% #count(typeScol)
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
  scale_colour_manual(values = c('< 0.05' = '#cb181d', '< 0.1' = '#fb6a4a', '< 0.25' = '#fcae91', '>= 0.25' = 'darkgrey')) +
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

# ggsave(plot = p_lrti_binomial_Sa_ASV_abx_cohorts, filename = "./figs/main/p_lrti_binomial_Sa_ASV_abx_cohorts.png", height = 5, width = 8, units = "in", dpi = 600)
# ggsave(plot = p_lrti_binomial_Sa_ASV_abx_cohorts, filename = "./figs/main/p_lrti_binomial_Sa_ASV_abx_cohorts.pdf", height = 5, width = 8, units = "in")


# test version with viridis plasma color scale:
p_lrti_binomial_Sa_ASV_abx_cohorts +
  scale_colour_manual(values = c('< 0.05' = '#7301A8FF', '< 0.1' = '#BD3786FF', '< 0.25' = '#ED7953FF', '>= 0.25' = 'darkgrey'))






#' ########################################
#' ADD SPECIES LABELS -- COMBINE COHORTS: LRTI ~ Sa + other microbiome ASVs + antibiotics (7-days prior)
#' ########################################

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
         typeScol = ifelse(typeS5, "< 0.05", ifelse(typeS10, "< 0.1", ifelse(typeS25, "< 0.25", ">= 0.25")))) %>% #count(typeScol)
  group_by(seqvar_id) %>%
  mutate(checktypeS = sum(typeS10) > 0) %>%
  ungroup() %>%
  filter(checktypeS) %>% # filter to type S error < 10%
  #count(best)
  mutate(blast_genus = gsub(" 1","",blast_genus),
         seqvar_lab = paste0(blast_genus,"\n(",seqvar_id,")"),
         best_lab = paste0(gsub(" 1| ATCC 49725| W23144| PA99","",best),"\n(",seqvar_id,")")
         ) %>%
  mutate_at(.vars = vars(mean, se_mean, sd, contains("%")), .funs = ~ exp(.x)) %>% # convert posterior samples to OR scale
  ggplot(data = .) +
  geom_vline(xintercept = 1, linetype = 2) +
  geom_segment(aes(x = `25%`, xend = `75%`, y = best_lab, yend = best_lab, colour = typeScol)) +
  geom_point(aes(x = `50%`, y = best_lab, colour = typeScol)) +
  facet_wrap(~cohort, nrow = 1) +
  scale_colour_manual(values = c('< 0.05' = '#cb181d', '< 0.1' = '#fb6a4a', '< 0.25' = '#fcae91', '>= 0.25' = 'darkgrey')) +
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

# ggsave(plot = p_lrti_binomial_Sa_ASV_abx_cohorts, filename = "./figs/main/p_lrti_binomial_Sa_ASV_abx_cohorts.png", height = 5, width = 8, units = "in", dpi = 600)
# ggsave(plot = p_lrti_binomial_Sa_ASV_abx_cohorts, filename = "./figs/main/p_lrti_binomial_Sa_ASV_abx_cohorts.pdf", height = 5, width = 8, units = "in")


# test version with viridis plasma color scale:
p_lrti_binomial_Sa_ASV_abx_cohorts +
  scale_colour_manual(values = c('< 0.05' = '#7301A8FF', '< 0.1' = '#BD3786FF', '< 0.25' = '#ED7953FF', '>= 0.25' = 'darkgrey'))




#' #######################################
#' FIGURE: compose ORIGINAL two plots for Figure 2 (using cowplot package)
#' #######################################

ggdraw(plot_grid(p_lrti_binomial_Sa_cohorts + theme(legend.position = "none"), p_lrti_binomial_Sa_ASV_abx_cohorts, nrow = 2, rel_heights = c(3.5,5.5), labels = c("A","B"))) -> p_lrti_model_combined
p_lrti_model_combined

# ggsave(plot = p_lrti_model_combined, filename = "./figs/main/p_lrti_model_combined.png", height = 9, width = 8, units = "in", dpi = 600)
# ggsave(plot = p_lrti_model_combined, filename = "./figs/main/p_lrti_model_combined.svg", height = 9, width = 8, units = "in")
# ggsave(plot = p_lrti_model_combined, filename = "./figs/main/p_lrti_model_combined.pdf", height = 9, width = 8, units = "in")





#
###
#####
###
#



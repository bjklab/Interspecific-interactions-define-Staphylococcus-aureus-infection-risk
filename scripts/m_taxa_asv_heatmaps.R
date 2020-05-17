#' ########################################
#' heatmaps of family-level (all ASVs) & ASV-level (Sa & Corynebacterium)
#' ########################################
#' 
#' depends: asv_micu_complete
#' 

library(tidyverse)
library(cowplot)
set.seed(16)

list(read_csv("./tabs/AN_asv_summary.csv.gz",
              col_types = cols(`respiratory_cx_pathogen` = col_character(),
                               `blood_cx_pathogen` = col_character())),
     read_csv("./tabs/ET_asv_summary.csv.gz",
              col_types = cols(`respiratory_cx_pathogen` = col_character(),
                               `blood_cx_pathogen` = col_character()))) %>%
  bind_rows() -> asv_micu_complete
asv_micu_complete



#' ########################################
#' POOLED COHORTS
#' ########################################
#' subject breaks
asv_micu_complete %>%
  select(subject_id, subject_day, specimen_type, seqvar_id, read_count, blast_family, blast_genus, best) %>%
  filter(!is.na(seqvar_id) & specimen_type == "AN") %>%
  select(subject_id,subject_day) %>%
  distinct() %>%
  arrange(subject_id, subject_day) %>%
  filter(grepl("pilot",subject_id)) %>%
  count(subject_id) %>%
  pull(n) %>%
  cumsum() + 0.5 -> pilot_breaks
pilot_breaks <- pilot_breaks[1:(length(pilot_breaks)-1)]
pilot_breaks


asv_micu_complete %>%
  select(subject_id, subject_day, specimen_type, seqvar_id, read_count, blast_family, blast_genus, best) %>%
  filter(!is.na(seqvar_id) & specimen_type == "AN") %>%
  select(subject_id,subject_day) %>%
  distinct() %>%
  arrange(subject_id, subject_day) %>%
  filter(grepl("validation",subject_id)) %>%
  count(subject_id) %>%
  pull(n) %>%
  cumsum() + 0.5 -> valid_breaks
valid_breaks <- valid_breaks[1:(length(valid_breaks)-1)]
valid_breaks


subject_breaks <- tibble(cohort = c(rep("Pilot",length(pilot_breaks)),rep("Validation",length(valid_breaks))),
                         linebreaks = c(pilot_breaks,valid_breaks))


#' family-level AN heatmap
asv_micu_complete %>%
  select(subject_id, subject_num, subject_day, specimen_type, seqvar_id, read_count, blast_family, blast_genus, best) %>%
  filter(!is.na(seqvar_id) & specimen_type == "AN") %>%
  distinct() %>%
  group_by(subject_id, subject_day) %>%
  mutate(total_reads = sum(read_count, na.rm = TRUE)) %>%
  ungroup() %>%
  group_by(subject_id, subject_day, seqvar_id) %>%
  mutate(read_count = sum(read_count, na.rm = TRUE),
         seqvar_prop = read_count / total_reads) %>%
  ungroup() %>%
  group_by(subject_id, subject_day, blast_family) %>%
  mutate(family_prop = sum(read_count, na.rm = TRUE) / total_reads) %>%
  ungroup() %>%
  mutate(cohort = ifelse(grepl("pilot",subject_id), "Pilot", "Validation"),
         spec_label = paste0(cohort," ", stringr::str_pad(string = as.numeric(subject_num), width = 2, side = "left", pad = "0"), " Day ",stringr::str_pad(string = as.numeric(subject_day), width = 2, side = "left", pad = "0"))) %>%
  select(spec_label, blast_family, family_prop, cohort) %>%
  distinct() %>%
  filter(!is.na(blast_family)) %>%
  mutate(blast_family = factor(blast_family, levels = sort(unique(blast_family), decreasing = TRUE))) %>%
  ggplot(data = .) +
  geom_tile(aes(x = spec_label, y = blast_family, fill = family_prop)) +
  facet_wrap(~ cohort, ncol = 1, scales = "free_x") +
  #geom_segment(data = subject_breaks, aes(x = linebreaks, xend = linebreaks, y = -Inf, yend = Inf), colour = "white", size = 0.5) +
  geom_vline(data = subject_breaks, aes(xintercept = linebreaks), colour = "white", size = 0.1) +
  scale_fill_viridis_c(option = "magma", trans = "log10", na.value = "black", limits = c(0.00001,1), breaks = c(1,0.1,0.01,0.001,0.0001)) +
  scale_y_discrete(expand = c(0,0)) +
  scale_x_discrete(expand = c(0,0)) +
  theme_bw() +
  theme(strip.background = element_blank(),
        axis.text.x = element_text(angle = 270, size = 3, vjust = 0.5, hjust = 0),
        axis.text.y = element_text(vjust = 0.5, size = 6),
        legend.title = element_text(size = 5),
        legend.text = element_text(size = 5)) +
  labs(x = "", y = "", fill = "Proportional\nAbundance\nof Bacterial\nFamily in\nAnterior\nNares") -> p_family_AN_heat
p_family_AN_heat

# ggsave(plot = p_family_AN_heat, filename = "./figs/supp/p_anterior_nares_family_heatmap.pdf", height = 8, width = 12, units = "in")
# ggsave(plot = p_family_AN_heat, filename = "./figs/supp/p_anterior_nares_family_heatmap.png", height = 8, width = 12, units = "in", dpi = 600)



#' family-level ET heatmap
asv_micu_complete %>%
  select(subject_id, subject_day, specimen_type, seqvar_id, read_count, blast_family, blast_genus, best) %>%
  filter(!is.na(seqvar_id) & specimen_type == "ET") %>%
  select(subject_id,subject_day) %>%
  distinct() %>%
  arrange(subject_id, subject_day) %>%
  filter(grepl("pilot",subject_id)) %>%
  count(subject_id) %>%
  pull(n) %>%
  cumsum() + 0.5 -> et_pilot_breaks
et_pilot_breaks <- et_pilot_breaks[1:(length(et_pilot_breaks)-1)]
et_pilot_breaks


asv_micu_complete %>%
  select(subject_id, subject_day, specimen_type, seqvar_id, read_count, blast_family, blast_genus, best) %>%
  filter(!is.na(seqvar_id) & specimen_type == "ET") %>%
  select(subject_id,subject_day) %>%
  distinct() %>%
  arrange(subject_id, subject_day) %>%
  filter(grepl("validation",subject_id)) %>%
  count(subject_id) %>%
  pull(n) %>%
  cumsum() + 0.5 -> et_valid_breaks
et_valid_breaks <- et_valid_breaks[1:(length(et_valid_breaks)-1)]
et_valid_breaks


et_subject_breaks <- tibble(cohort = c(rep("Pilot",length(et_pilot_breaks)),rep("Validation",length(et_valid_breaks))),
                         linebreaks = c(et_pilot_breaks,et_valid_breaks))


asv_micu_complete %>%
  select(subject_id, subject_num, subject_day, specimen_type, seqvar_id, read_count, blast_family, blast_genus, best) %>%
  filter(!is.na(seqvar_id) & specimen_type == "ET") %>%
  distinct() %>%
  group_by(subject_id, subject_day) %>%
  mutate(total_reads = sum(read_count, na.rm = TRUE)) %>%
  ungroup() %>%
  group_by(subject_id, subject_day, seqvar_id) %>%
  mutate(read_count = sum(read_count, na.rm = TRUE),
         seqvar_prop = read_count / total_reads) %>%
  ungroup() %>%
  group_by(subject_id, subject_day, blast_family) %>%
  mutate(family_prop = sum(read_count, na.rm = TRUE) / total_reads) %>%
  ungroup() %>%
  mutate(cohort = ifelse(grepl("pilot",subject_id), "Pilot", "Validation"),
         spec_label = paste0(cohort," ", stringr::str_pad(string = as.numeric(subject_num), width = 2, side = "left", pad = "0"), " Day ",stringr::str_pad(string = as.numeric(subject_day), width = 2, side = "left", pad = "0"))) %>%
  select(spec_label, blast_family, family_prop, cohort) %>%
  distinct() %>%
  filter(!is.na(blast_family)) %>%
  mutate(blast_family = factor(blast_family, levels = sort(unique(blast_family), decreasing = TRUE))) %>%
  ggplot(data = .) +
  geom_tile(aes(x = spec_label, y = blast_family, fill = family_prop)) +
  facet_wrap(~ cohort, ncol = 1, scales = "free_x") +
  #geom_segment(data = et_subject_breaks, aes(x = linebreaks, xend = linebreaks, y = -Inf, yend = Inf), colour = "white", size = 0.5) +
  geom_vline(data = et_subject_breaks, aes(xintercept = linebreaks), colour = "white", size = 0.1) +
  scale_fill_viridis_c(option = "magma", trans = "log10", na.value = "black", limits = c(0.00001,1), breaks = c(1,0.1,0.01,0.001,0.0001)) +
  scale_y_discrete(expand = c(0,0)) +
  scale_x_discrete(expand = c(0,0)) +
  theme_bw() +
  theme(strip.background = element_blank(),
        axis.text.x = element_text(angle = 270, size = 5, vjust = 0.5, hjust = 0),
        axis.text.y = element_text(vjust = 0.5, size = 6),
        legend.title = element_text(size = 5),
        legend.text = element_text(size = 5)) +
  labs(x = "", y = "", fill = "Proportional\nAbundance\nof Bacterial\nFamily in\nEndotracheal\nAspirate") -> p_family_ET_heat
p_family_ET_heat

# ggsave(plot = p_family_ET_heat, filename = "./figs/supp/p_endotracheal_family_heatmap.pdf", height = 8, width = 12, units = "in")
# ggsave(plot = p_family_ET_heat, filename = "./figs/supp/p_endotracheal_family_heatmap.png", height = 8, width = 12, units = "in", dpi = 600)






#' ASV-level heatmap: Staphylococcus
asv_micu_complete %>%
  select(subject_id, subject_num, subject_day, specimen_type, seqvar_id, read_count, blast_family, blast_genus, best) %>%
  filter(!is.na(seqvar_id) & specimen_type == "AN") %>%
  distinct() %>%
  group_by(subject_id, subject_day) %>%
  mutate(total_reads = sum(read_count, na.rm = TRUE)) %>%
  ungroup() %>%
  group_by(subject_id, subject_day, seqvar_id) %>%
  mutate(read_count = sum(read_count, na.rm = TRUE),
         seqvar_prop = read_count / total_reads) %>%
  ungroup() %>%
  mutate(cohort = ifelse(grepl("pilot",subject_id), "Pilot", "Validation"),
         spec_label = paste0(cohort," ", stringr::str_pad(string = as.numeric(subject_num), width = 2, side = "left", pad = "0"), " Day ",stringr::str_pad(string = as.numeric(subject_day), width = 2, side = "left", pad = "0"))) %>%
  select(subject_id, spec_label, seqvar_id, seqvar_prop, cohort, blast_genus) %>%
  distinct() %>%
  filter(grepl("staphylococcus",tolower(blast_genus))) %>%
  ggplot(data = .) +
  geom_tile(aes(x = spec_label, y = seqvar_id, fill = seqvar_prop)) +
  facet_wrap(~ cohort, ncol = 1, scales = "free_x") +
  #geom_segment(data = subject_breaks, aes(x = linebreaks, xend = linebreaks, y = -Inf, yend = Inf), colour = "white", size = 0.5) +
  geom_vline(data = subject_breaks, aes(xintercept = linebreaks), colour = "white", size = 0.1) +
  scale_fill_viridis_c(option = "magma", trans = "log10", na.value = "black", limits = c(0.00001,1), breaks = c(1,0.1,0.01,0.001,0.0001)) +
  scale_y_discrete(expand = c(0,0)) +
  scale_x_discrete(expand = c(0,0)) +
  theme_bw() +
  theme(strip.background = element_blank(),
        axis.text.x = element_text(angle = 270, size = 3, vjust = 0.5, hjust = 1),
        axis.text.y = element_text(vjust = 0.5, size = 3),
        legend.title = element_text(size = 5),
        legend.text = element_text(size = 5)) +
  labs(x = "", y = "", fill = "Proportional\nAbundance of\nStaphylococcus\nASVs in\nAnterior Nares") -> p_staph_heat
p_staph_heat

# ggsave(plot = p_staph_heat, filename = "./figs/supp/p_anterior_nares_staph_asv_heatmap.pdf", height = 8, width = 12, units = "in")
# ggsave(plot = p_staph_heat, filename = "./figs/supp/p_anterior_nares_staph_asv_heatmap.png", height = 8, width = 12, units = "in", dpi = 600)





#' ASV-level heatmap: Corynebacterium
asv_micu_complete %>%
  select(subject_id, subject_num, subject_day, specimen_type, seqvar_id, read_count, blast_family, blast_genus, best) %>%
  filter(!is.na(seqvar_id) & specimen_type == "AN") %>%
  distinct() %>%
  group_by(subject_id, subject_day) %>%
  mutate(total_reads = sum(read_count, na.rm = TRUE)) %>%
  ungroup() %>%
  group_by(subject_id, subject_day, seqvar_id) %>%
  mutate(read_count = sum(read_count, na.rm = TRUE),
         seqvar_prop = read_count / total_reads) %>%
  ungroup() %>%
  mutate(cohort = ifelse(grepl("pilot",subject_id), "Pilot", "Validation"),
         spec_label = paste0(cohort," ", stringr::str_pad(string = as.numeric(subject_num), width = 2, side = "left", pad = "0"), " Day ",stringr::str_pad(string = as.numeric(subject_day), width = 2, side = "left", pad = "0"))) %>%
  select(subject_id, spec_label, seqvar_id, seqvar_prop, cohort, blast_genus) %>%
  distinct() %>%
  filter(grepl("corynebacterium",tolower(blast_genus))) %>%
  ggplot(data = .) +
  geom_tile(aes(x = spec_label, y = seqvar_id, fill = seqvar_prop)) +
  facet_wrap(~ cohort, ncol = 1, scales = "free_x") +
  #geom_segment(data = subject_breaks, aes(x = linebreaks, xend = linebreaks, y = -Inf, yend = Inf), colour = "white", size = 0.5) +
  geom_vline(data = subject_breaks, aes(xintercept = linebreaks), colour = "white", size = 0.1) +
  scale_fill_viridis_c(option = "magma", trans = "log10", na.value = "black", limits = c(0.00001,1), breaks = c(1,0.1,0.01,0.001,0.0001)) +
  scale_y_discrete(expand = c(0,0)) +
  scale_x_discrete(expand = c(0,0)) +
  theme_bw() +
  theme(strip.background = element_blank(),
        axis.text.x = element_text(angle = 270, size = 3, vjust = 0.5, hjust = 1),
        axis.text.y = element_text(vjust = 0.5, size = 3),
        legend.title = element_text(size = 5),
        legend.text = element_text(size = 5)) +
  labs(x = "", y = "", fill = "Proportional\nAbundance of\nCorynebacterium\nASVs in\nAnterior Nares") -> p_coryne_heat
p_coryne_heat

# ggsave(plot = p_coryne_heat, filename = "./figs/supp/p_anterior_nares_coryne_asv_heatmap.pdf", height = 8, width = 12, units = "in")
# ggsave(plot = p_coryne_heat, filename = "./figs/supp/p_anterior_nares_coryne_asv_heatmap.png", height = 8, width = 12, units = "in", dpi = 600)



#' combined heatmap
ggdraw(plot_grid(p_family_AN_heat, p_family_ET_heat, p_staph_heat, p_coryne_heat, nrow = 2, ncol = 2, rel_heights = c(8,8), rel_widths = c(12,12),
                 labels = c("A. Anterior Nares Bacterial Community",
                            "B. Endotracheeal Bacterial Community",
                            "C. Anterior Nares Staphylococcus ASVs",
                            "D. Anterior Nares Corynebacterium ASVs"),
                 hjust = -0.1)) -> p_heatmap_combined
p_heatmap_combined

# ggsave(plot = p_heatmap_combined, filename = "./figs/supp/p_heatmap_combined.png", height = 16, width = 24, units = "in", dpi = 600)
# ggsave(plot = p_heatmap_combined, filename = "./figs/supp/p_heatmap_combined.svg", height = 16, width = 24, units = "in")
# ggsave(plot = p_heatmap_combined, filename = "./figs/supp/p_heatmap_combined.pdf", height = 16, width = 24, units = "in")





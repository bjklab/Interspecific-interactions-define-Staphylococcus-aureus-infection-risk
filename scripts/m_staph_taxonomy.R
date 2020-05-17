#' ########################################
#' Staph taxonomy analysis
#' ########################################
#' 
#' depends: asv_micu_complete & micu_cohort
#' 

library(tidyverse)
library(patchwork)
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
  filter(specimen_type %in% c("AN","ET")) %>% # AN and ET specimens only
  group_by(seqvar_id) %>%
  filter(sum(read_count, na.rm = TRUE) > 0) %>% # filter out any ASVs not present in this specimen set
  ungroup() %>%
  select(embl_taxonomy, blast_genus, best, seqvar_id,subject_id) %>%
  #mutate(study_phase = str_to_title(str_extract(string = subject_id, pattern = "pilot|validation"))) %>%
  select(-subject_id) %>%
  distinct() %>%
  filter(grepl("Staphylococcus",embl_taxonomy) | grepl("Staphylococcus",blast_genus) | grepl("Staphylococcus", best)) %>%
  mutate(staph_best_species = str_extract(string = best, pattern = "aureus|capitis|epidermidis|haemolyticus|lugdunensis|pettenkoferi|sacchorolyticus|schleiferi"),
         staph_best_species = paste0("Staphylococcus ",staph_best_species),
         staph_best_species = gsub(pattern = "NA", replacement = "- species unresolved", x = staph_best_species)) %>%
  #count(staph_best_species,study_phase) %>%
  count(staph_best_species) %>%
  mutate(staph_best_species = factor(staph_best_species, levels = rev(staph_best_species))) %>%
  qplot(data = ., y = staph_best_species, yend = staph_best_species, x = 0, xend = n, geom = "segment") +
  theme_bw() +
  theme(strip.background = element_blank(),
        text = element_text(size=10)) +
  labs(x = expression(paste("Count of All ",italic("Staphylococcus")," amplicon sequence variants",sep = "")),
       #title = expression(paste("ASVs Assigned to ",italic("Staphylococcus"),sep = "")),
       y = "High Confidence Taxonomic Assignment") -> p_count


asv_micu_complete %>%
  filter(subject_id %in% micu_cohort$subject_id) %>% # included subjects only
  filter(specimen_type %in% c("AN","ET")) %>% # AN and ET specimens only
  group_by(seqvar_id) %>%
  filter(sum(read_count, na.rm = TRUE) > 0) %>% # filter out any ASVs not present in this specimen set
  ungroup() %>%
  select(embl_taxonomy, blast_genus, best, seqvar_id,subject_id) %>%
  #mutate(study_phase = str_to_title(str_extract(string = subject_id, pattern = "pilot|validation"))) %>%
  select(-subject_id) %>%
  distinct() %>%
  filter(grepl("Staphylococcus",embl_taxonomy) | grepl("Staphylococcus",blast_genus) | grepl("Staphylococcus", best)) %>%
  mutate(staph_best_species = str_extract(string = best, pattern = "aureus|capitis|epidermidis|haemolyticus|lugdunensis|pettenkoferi|sacchorolyticus|schleiferi"),
         staph_best_species = paste0("Staphylococcus ",staph_best_species),
         staph_best_species = gsub(pattern = "NA", replacement = "- species unresolved", x = staph_best_species)) %>%
  #count(staph_best_species,study_phase) %>%
  count(staph_best_species) %>%
  mutate(staph_best_species = factor(staph_best_species, levels = rev(staph_best_species)),
         n = n / sum(n)) %>%
  qplot(data = ., y = staph_best_species, yend = staph_best_species, x = 0, xend = n, geom = "segment") +
  theme_bw() +
  theme(strip.background = element_blank(),
        text = element_text(size=10)) +
  labs(x = expression(paste("Proportion of All ",italic("Staphylococcus")," amplicon sequence variants",sep = "")),
       #title = expression(paste("ASVs Assigned to ",italic("Staphylococcus"),sep = "")),
       y = "High Confidence Taxonomic Assignment") -> p_prop

p_Sa_taxonomy <- p_count + p_prop + plot_layout(ncol = 1)

# ggsave(plot = p_Sa_taxonomy, filename = "./figs/supp/p_Sa_taxonomy.pdf", height = 6, width = 6, units = "in")
# ggsave(plot = p_Sa_taxonomy, filename = "./figs/supp/p_Sa_taxonomy.svg", height = 6, width = 6, units = "in", system_fonts = list(sans = "Roboto"))
# ggsave(plot = p_Sa_taxonomy, filename = "./figs/supp/p_Sa_taxonomy.png", height = 6, width = 6, units = "in", dpi = 600)


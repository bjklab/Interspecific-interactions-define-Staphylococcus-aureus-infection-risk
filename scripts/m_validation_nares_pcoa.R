#' ########################################
#' evaluate whether AN communities differ by missing/not-missing ET specimens
#' ########################################
#' 
#' depends: asv_micu_complete
#' 


library(tidyverse)
library(tidybayes)
library(bayesplot)
library(ggnewscale)
library(ggsci)
library(cowplot)
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
#' VALIDATION COHORT (no missingness in PILOT)
#' ########################################



asv_micu_complete %>%
  filter(subject_id %in% micu_cohort$subject_id) %>%
  filter(study == 'micuvalidation') %>%
  select(seqvar_id, specimen_id, study, specimen_type, read_count, subject_id, subject_day) %>%
  distinct() %>%
  group_by(subject_id) %>%
  mutate(subject_ET_yn = sum(specimen_type == "ET") > 0) %>%
  ungroup() %>%
  select(specimen_id, subject_ET_yn) %>%
  distinct() -> subject_ET_check

subject_ET_check




#' ############################################
#' DISTANCE: Weighted Jaccard
#' ############################################

asv_micu_complete %>%
  filter(subject_id %in% micu_cohort$subject_id) %>%
  filter(study == 'micuvalidation') %>%
  filter(specimen_type == "AN") %>%
  select(seqvar_id, specimen_id, read_count) %>%
  distinct() %>%
  complete(seqvar_id, specimen_id, fill = list('read_count' = 0)) %>%
  spread(key = seqvar_id, value = read_count) %>%
  column_to_rownames(var = "specimen_id") %>%
  as.matrix() %>%
  .[rowSums(.) > 0,] %>%
  .[,colSums(.) > 0] %>% 
  vegan::vegdist(x = ., method = 'jaccard', binary = FALSE) -> weighted_jaccard_dist

str(weighted_jaccard_dist)



#' ############################################
#' DISTANCE: unweighted Jaccard
#' ############################################

asv_micu_complete %>%
  filter(subject_id %in% micu_cohort$subject_id) %>%
  filter(study == 'micuvalidation') %>%
  filter(specimen_type == "AN") %>%
  select(seqvar_id, specimen_id, read_count) %>%
  distinct() %>%
  complete(seqvar_id, specimen_id, fill = list('read_count' = 0)) %>%
  spread(key = seqvar_id, value = read_count) %>%
  column_to_rownames(var = "specimen_id") %>%
  as.matrix() %>%
  .[rowSums(.) > 0,] %>%
  .[,colSums(.) > 0] %>% 
  vegan::vegdist(x = ., method = 'jaccard', binary = TRUE) -> unweighted_jaccard_dist

str(unweighted_jaccard_dist)






#' ############################################
#' DISTANCE: weighted normalized UniFrac
#' ############################################

asv_micu_complete %>%
  filter(subject_id %in% micu_cohort$subject_id) %>%
  filter(study == 'micuvalidation') %>%
  filter(specimen_type == "AN") %>%
  select(seqvar_id, specimen_id, read_count) %>%
  distinct() %>%
  complete(seqvar_id, specimen_id, fill = list('read_count' = 0)) %>%
  spread(key = seqvar_id, value = read_count) %>%
  column_to_rownames(var = "specimen_id") %>%
  as.matrix() %>%
  .[rowSums(.) > 0,] %>%
  .[,colSums(.) > 0] %>% #str()
  phyloseq::otu_table(object = ., taxa_are_rows = FALSE) %>% #str()
  phyloseq::phyloseq(., phyloseq::phy_tree(ape::read.tree("./data/rooted_tree.nwk"))) %>% #str()
  phyloseq::UniFrac(physeq = ., weighted = TRUE, normalized = TRUE) -> weighted_unifrac_dist

str(weighted_unifrac_dist)




#' ############################################
#' DISTANCE: unweighted UniFrac
#' ############################################


asv_micu_complete %>%
  filter(subject_id %in% micu_cohort$subject_id) %>%
  filter(study == 'micuvalidation') %>%
  filter(specimen_type == "AN") %>%
  select(seqvar_id, specimen_id, read_count) %>%
  distinct() %>%
  complete(seqvar_id, specimen_id, fill = list('read_count' = 0)) %>%
  spread(key = seqvar_id, value = read_count) %>%
  column_to_rownames(var = "specimen_id") %>%
  as.matrix() %>%
  .[rowSums(.) > 0,] %>%
  .[,colSums(.) > 0] %>% #str()
  phyloseq::otu_table(object = ., taxa_are_rows = FALSE) %>% #str()
  phyloseq::phyloseq(., phyloseq::phy_tree(ape::keep.tip(ape::read.tree("./data/rooted_tree.nwk"),colnames(.)))) %>% #str()
  phyloseq::UniFrac(physeq = ., weighted = FALSE, normalized = TRUE) -> unweighted_unifrac_dist

str(unweighted_unifrac_dist)




#' ############################################
#' PERMANOVA/adonis: all distance metrics
#' ############################################

# vegan::adonis(weighted_unifrac_dist ~ subject_ET_yn, data = as.data.frame(filter(subject_ET_check, specimen_id %in% labels(weighted_unifrac_dist))))$aov.tab %>%
#   broom::tidy()
# vegan::adonis2(weighted_unifrac_dist ~ subject_ET_yn, data = as.data.frame(filter(subject_ET_check, specimen_id %in% labels(weighted_unifrac_dist)))) %>%
#   broom::tidy()


et_bias_check <- list("Weighted Jaccard" = weighted_jaccard_dist, "Unweighted Jaccard" = unweighted_jaccard_dist, "Weighted UniFrac" = weighted_unifrac_dist, "Unweighted UniFrac" = unweighted_unifrac_dist)

et_bias_check %>%
  map2_dfr(.x = ., .y = names(.), .f = ~ vegan::adonis2(.x ~ subject_ET_yn, data = as.data.frame(filter(subject_ET_check, specimen_id %in% labels(.x)))) %>%
        broom::tidy() %>%
          mutate(distance_metric = .y)) -> et_bias_adonis

et_bias_adonis %>%
  filter(term == "subject_ET_yn") %>%
  select(distance_metric, R2) -> et_bias_adonis_summary

et_bias_adonis_summary










#' ############################################
#' PCoA Plots
#' ############################################

pcoa_melt <- function(dist_mat, metadata = subject_ET_check) {
  pcoa = ape::pcoa(D = dist_mat)
  pcoa %>%
    pluck('vectors') %>%
    as.data.frame() %>%
    as_tibble(rownames = NA) %>%
    rownames_to_column(var = 'specimen_id') %>%
    distinct() -> p_coords
  pcoa %>%
    pluck('values') %>%
    as.data.frame() %>%
    as_tibble() %>%
    mutate(axis = paste0('Axis.',seq_along(Relative_eig))) -> p_eig
  p_coords %>%
    gather(key = axis, value = coord, -specimen_id) %>%
    left_join(metadata, by = 'specimen_id') %>%
    left_join(p_eig, by = 'axis') %>%
    distinct() %>%
    mutate(coord_label = paste0("PCoA ",axis," (",round(Relative_eig,3)*100,"% variance)")) %>%
    filter(axis == "Axis.1" | axis == "Axis.2") %>%
    select(specimen_id, coord, coord_label, subject_ET_yn) -> p_check
  
  return(p_check)
}




pcoa_melt(weighted_jaccard_dist)


pcoa_melt(weighted_jaccard_dist) %>%
  spread(key = "coord_label", value = "coord") %>%
  ggplot(data = ., aes(x = `PCoA Axis.1 (11.4% variance)`, y = `PCoA Axis.2 (5.7% variance)`, colour = subject_ET_yn)) +
  geom_point() +
  #stat_ellipse(type = 't', level = 0.95) +
  ggsci::scale_color_nejm() +
  theme_bw() +
  labs(colour = 'ET specimen\nalso collected\nfrom subject',
       title = "Principal Coordinate Analysis (weighted Jaccard distances)",
       subtitle = '187 Validation Cohort AN Specimens from 76 Subjects'
       ) 



pcoa_melt(weighted_unifrac_dist)


pcoa_melt(weighted_unifrac_dist) %>%
  spread(key = "coord_label", value = "coord") %>%
  ggplot(data = ., aes(x = `PCoA Axis.1 (86% variance)`, y = `PCoA Axis.2 (10.7% variance)`, colour = subject_ET_yn)) +
  geom_point() +
  #stat_ellipse(type = 't', level = 0.95) +
  ggsci::scale_color_nejm() +
  theme_bw() +
  labs(colour = 'ET specimen\nalso collected\nfrom subject',
       title = "Principal Coordinate Analysis (weighted UniFrac distances)",
       subtitle = '187 Validation Cohort AN Specimens from 76 Subjects') 





pcoa_melt(unweighted_jaccard_dist)


pcoa_melt(unweighted_jaccard_dist) %>%
  spread(key = "coord_label", value = "coord") %>%
  ggplot(data = ., aes(x = `PCoA Axis.1 (9.1% variance)`, y = `PCoA Axis.2 (5.9% variance)`, colour = subject_ET_yn)) +
  geom_point() +
  #stat_ellipse(type = 't', level = 0.95) +
  ggsci::scale_color_nejm() +
  theme_bw() +
  labs(colour = 'ET specimen\nalso collected\nfrom subject',
       title = "Principal Coordinate Analysis (unweighted Jaccard distances)",
       subtitle = '187 Validation Cohort AN Specimens from 76 Subjects') 



pcoa_melt(unweighted_unifrac_dist)


pcoa_melt(unweighted_unifrac_dist) %>%
  spread(key = "coord_label", value = "coord") %>%
  ggplot(data = ., aes(x = `PCoA Axis.1 (15.3% variance)`, y = `PCoA Axis.2 (9.2% variance)`, colour = subject_ET_yn)) +
  geom_point() +
  #stat_ellipse(type = 't', level = 0.95) +
  ggsci::scale_color_nejm() +
  theme_bw() +
  labs(colour = 'ET specimen\nalso collected\nfrom subject',
       title = "Principal Coordinate Analysis (unweighted UniFrac distances)",
       subtitle = '187 Validation Cohort AN Specimens from 76 Subjects') 








#' save PCoA plots

pcoa_melt(weighted_jaccard_dist) %>%
  spread(key = "coord_label", value = "coord") %>%
  ggplot(data = ., aes(x = `PCoA Axis.1 (11.4% variance)`, y = `PCoA Axis.2 (5.7% variance)`, colour = subject_ET_yn)) +
  geom_point() +
  #stat_ellipse(type = 't', level = 0.95) +
  ggsci::scale_color_nejm() +
  theme_bw() +
  labs(colour = 'ET specimen\nalso collected\nfrom subject',
       title = "Principal Coordinate Analysis (weighted Jaccard distances)",
       subtitle = '187 Validation Cohort AN Specimens from 76 Subjects') -> p_validation_AN_pcoa_WJ

p_validation_AN_pcoa_WJ

#ggsave(plot = p_validation_AN_pcoa_WJ, filename = "./figs/supp/p_validation_AN_pcoa.pdf", height = 4, width = 6, units = "in")
#ggsave(plot = p_validation_AN_pcoa_WJ, filename = "./figs/supp/p_validation_AN_pcoa.png", height = 4, width = 6, units = "in", dpi = 600)


#' save all distances

list("Weighted Jaccard" = pcoa_melt(weighted_jaccard_dist),
      "Unweighted Jaccard" = pcoa_melt(unweighted_jaccard_dist),
      "Weighted UniFrac" = pcoa_melt(weighted_unifrac_dist),
      "Unweighted UniFrac" = pcoa_melt(unweighted_unifrac_dist)) %>%
  bind_rows(.id = "distance_metric") %>%
  mutate(coord_simple = stringr::str_extract(string = coord_label, pattern = "Axis\\.1|Axis\\.2")) %>%
  select(-coord_label) %>%
  spread(key = "coord_simple", value = "coord") %>%
  left_join(et_bias_adonis_summary, by = "distance_metric") %>%
  mutate(r2_label = paste0("R2 = ",signif(R2, digits = 2))) %>%
  group_by(distance_metric) %>%
  mutate(xpos = quantile(x = Axis.1, 0.2),
         ypos = quantile(x = Axis.2, 0.8)) %>%
  ungroup() %>%
  ggplot(data = ., aes(x = Axis.1, y = Axis.2, colour = subject_ET_yn)) +
  geom_point() +
  geom_text(data = et_bias_adonis_summary, aes(label = paste0("R^2 == ",signif(R2, digits = 3))), x = 0, y = 0.4, colour = "black", parse = TRUE) +
  facet_wrap(~ distance_metric) +
  #stat_ellipse(type = 't', level = 0.95) +
  ggsci::scale_color_nejm() +
  theme_bw() +
  theme(legend.position = "bottom") +
  theme(strip.background = element_blank()) +
  labs(colour = 'ET specimen also collected from subject',
       title = "Principal Coordinate Analysis",
       subtitle = "187 Validation Cohort AN Specimens from 76 Subjects") -> p_all_dist
  
p_all_dist

# ggsave(plot = p_all_dist, filename = "./figs/supp/p_validation_AN_pcoa_R2.pdf", height = 7, width = 6, units = "in")
# ggsave(plot = p_all_dist, filename = "./figs/supp/p_validation_AN_pcoa_R2.svg", height = 7, width = 6, units = "in", system_fonts = list(sans = "Roboto"))
# ggsave(plot = p_all_dist, filename = "./figs/supp/p_validation_AN_pcoa_R2.png", height = 7, width = 6, units = "in", dpi = 600)


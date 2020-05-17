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






#' #######################################
#' FIGURE: ET and AN S.aureus proportional abundance
#' #######################################



asv_micu_complete %>%
  filter(subject_id %in% micu_cohort$subject_id) %>% # included subjects only
  mutate(Sa_prop = sumstaph / specimen_read_total) %>%
  select(subject_id,subject_day,specimen_type,Sa_prop) %>%
  distinct() %>%
  filter(specimen_type %in% c("AN","ET")) -> et_an_comparison

et_an_comparison %>%
  mutate(specimen_type = replace(specimen_type,specimen_type == "AN", "Anterior Nares"),
         specimen_type = replace(specimen_type,specimen_type == "ET", "Endotracheal")) %>%
  ggplot(data = .) +
  geom_density(aes(x = Sa_prop, fill = specimen_type), colour = "black", alpha = 0.6) +
  theme_bw() +
  ggsci::scale_fill_nejm() +
  theme(legend.position = 'bottom') +
  scale_x_log10(labels = scales::scientific) +
  labs(x = expression(paste("Proportional Abundance of ",italic("S. aureus"))),
       y = "Density",
       fill = "Sampling Site") -> p_et_v_an_Sa_prop_density
p_et_v_an_Sa_prop_density

# ggsave(plot = p_et_v_an_Sa_prop_density, filename = "./figs/main/p_et_v_an_Sa_prop_density.png", height = 4, width = 4, units = "in", dpi = 600)
# ggsave(plot = p_et_v_an_Sa_prop_density, filename = "./figs/main/p_et_v_an_Sa_prop_density.pdf", height = 4, width = 4, units = "in")



et_an_comparison %>%
  mutate(specimen_type = replace(specimen_type,specimen_type == "AN", "Anterior Nares"),
         specimen_type = replace(specimen_type,specimen_type == "ET", "Endotracheal")) %>%
  ggplot(data = .) +
  geom_histogram(aes(x = Sa_prop, fill = specimen_type), colour = "black", alpha = 0.6) +
  facet_wrap(facets = ~ specimen_type, ncol = 1, scales = "free_y") +
  theme_bw() +
  ggsci::scale_fill_nejm() +
  theme(legend.position = 'none',
        strip.background = element_blank()) +
  scale_x_log10(labels = scales::scientific) +
  labs(x = expression(paste("Proportional Abundance of ",italic("S. aureus"))),
       y = "Count",
       fill = "Sampling Site") -> p_et_v_an_Sa_prop_histogram
p_et_v_an_Sa_prop_histogram

# ggsave(plot = p_et_v_an_Sa_prop_histogram, filename = "./figs/main/p_et_v_an_Sa_prop_histogram.png", height = 4, width = 4, units = "in", dpi = 600)
# ggsave(plot = p_et_v_an_Sa_prop_histogram, filename = "./figs/main/p_et_v_an_Sa_prop_histogram.pdf", height = 4, width = 4, units = "in")






#' #######################################
#' FIGURE: ET vs AN S.aureus proportional abundance
#' #######################################
et_an_comparison %>%
  group_by(subject_id,subject_day) %>%
  filter(paste(sort(specimen_type), collapse = " ") == "AN ET") %>% # only include same-day same-subject pairs
  spread(key = specimen_type, value = Sa_prop) -> et_an_same_day
et_an_same_day


et_an_same_day %>%
  ggplot(data = .) +
  geom_point(aes(x = AN, y = ET)) +
  geom_smooth(aes(x = AN, y = ET), method = "lm", colour = "black") +
  scale_x_log10(labels = scales::scientific) +
  scale_y_log10(labels = scales::scientific) +
  #coord_equal() +
  theme_bw() +
  labs(x = expression(paste("Anterior Nares ",italic("S. aureus")," Abundance")),
       y = expression(paste("Endotracheal ",italic("S. aureus")," Abundance")),
       title = "B") -> p_et_v_an_Sa_prop_biplot
p_et_v_an_Sa_prop_biplot

# ggsave(plot = p_et_v_an_Sa_prop_biplot, filename = "./figs/main/p_et_v_an_Sa_prop_biplot.png", height = 4, width = 4, units = "in", dpi = 600)
# ggsave(plot = p_et_v_an_Sa_prop_biplot, filename = "./figs/main/p_et_v_an_Sa_prop_biplot.pdf", height = 4, width = 4, units = "in")

et_an_same_day %>%
  ungroup() %>%
  mutate(AN = log(AN),
         ET = log(ET)) %>%
  filter(abs(AN) != Inf & abs(ET) != Inf) -> log_et_an_same_day

rstanarm::stan_glm(formula = ET ~ AN,
                   family = gaussian(),
                   data = log_et_an_same_day,
                   algorithm = "sampling",
                   iter = 1000,
                   chains = 8,
                   verbose = TRUE,
                   seed = 16) -> model_et_an_same_day

lm_colors <- colorspace::sequential_hcl(5, palette = "Light Grays")

log_et_an_same_day %>%
  tidyr::expand(AN = seq(min(AN, na.rm = TRUE), max(AN, na.rm = TRUE), length.out = 100)) %>%
  tidybayes::add_fitted_draws(model_et_an_same_day) %>% 
  tidybayes::median_qi(.width = c(.95, .8, .5)) %>% 
  ungroup() %>%
  mutate_at(.vars = vars(AN,.value,.lower,.upper), .funs = ~ exp(.x)) %>%
  ggplot() + 
  aes(x = AN, y = .value) + 
  geom_lineribbon(aes(fill = factor(.width), ymin = .lower, ymax = .upper, color = .point)) + 
  geom_point(data = et_an_same_day, aes(x = AN, y = ET)) +
  scale_color_manual(aesthetics = c("fill", "color"),
                     # I have no idea how this ordering works. Takes some trial and error
                     breaks = c("0.5", "0.8", "0.95", "median"),
                     values = c(lm_colors[2:4], lm_colors[1])) + 
  guides(fill = guide_legend(title = "Posterior CI",
                             override.aes = list(color = c(NA, NA, NA, lm_colors[1]),
                                                 fill = c(lm_colors[2:4], NA))),
         color = FALSE) +
  scale_x_log10(labels = scales::scientific) +
  scale_y_log10(labels = scales::scientific) +
  #coord_equal() +
  theme_bw() +
  theme(legend.position = 'bottom') +
  labs(x = expression(paste("Anterior Nares ",italic("S. aureus")," Abundance")),
       y = expression(paste("Endotracheal ",italic("S. aureus")," Abundance"))) -> p_et_v_an_Sa_prop_model_biplot
p_et_v_an_Sa_prop_model_biplot

# ggsave(plot = p_et_v_an_Sa_prop_model_biplot, filename = "./figs/main/p_et_v_an_Sa_prop_model_biplot.png", height = 4, width = 4, units = "in", dpi = 600)
# ggsave(plot = p_et_v_an_Sa_prop_model_biplot, filename = "./figs/main/p_et_v_an_Sa_prop_model_biplot.pdf", height = 4, width = 4, units = "in")



#' #######################################
#' FIGURE: ASV-level ET vs AN S.aureus proportional abundance
#' #######################################
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
  count(seqvar_id)

et_an_asv_same_day %>%
  group_by(seqvar_id) %>%
  filter(sum(AN, na.rm = TRUE) + sum(ET, na.rm = TRUE) > 0) %>%  # remove ASVs not found in this specimen set
  ungroup() -> et_an_asv_same_day_filter
et_an_asv_same_day_filter


et_an_asv_same_day_filter %>%
  ggplot(data = .) +
  geom_point(aes(x = AN, y = ET, colour = seqvar_id)) +
  geom_smooth(aes(x = AN, y = ET), method = "lm", colour = "black") +
  scale_x_log10(labels = scales::scientific) +
  scale_y_log10(labels = scales::scientific) +
  #coord_equal() +
  scale_color_manual(values = rainbow(10)) +
  theme_bw() +
  labs(x = expression(paste("Anterior Nares ",italic("S. aureus"), " ASV Abundance")),
       y = expression(paste("Endotracheal ",italic("S. aureus"), " ASV Abundance")),
       colour = expression(paste(italic("S. aureus")," Amplicon Sequence Variant"))) -> p_et_v_an_asv_Sa_prop_biplot
p_et_v_an_asv_Sa_prop_biplot

# ggsave(plot = p_et_v_an_asv_Sa_prop_biplot, filename = "./figs/main/p_et_v_an_asv_Sa_prop_biplot.png", height = 4, width = 6, units = "in", dpi = 600)
# ggsave(plot = p_et_v_an_asv_Sa_prop_biplot, filename = "./figs/main/p_et_v_an_asv_Sa_prop_biplot.pdf", height = 4, width = 6, units = "in")




et_an_asv_same_day %>%
  group_by(seqvar_id) %>%
  filter(sum(AN, na.rm = TRUE) + sum(ET, na.rm = TRUE) > 0) %>%  # remove ASVs not found in this specimen set
  ungroup() %>%
  mutate(AN = log(AN),
         ET = log(ET)) %>%
  filter(abs(AN) != Inf & abs(ET) != Inf) -> log_et_an_asv_same_day

rstanarm::stan_glm(formula = ET ~ AN,
                   family = gaussian(),
                   data = log_et_an_asv_same_day,
                   algorithm = "sampling",
                   iter = 1000,
                   chains = 8,
                   verbose = TRUE,
                   seed = 16) -> model_et_an_asv_same_day

lm_colors <- colorspace::sequential_hcl(5, palette = "Light Grays")

log_et_an_asv_same_day %>%
  tidyr::expand(AN = seq(min(AN, na.rm = TRUE), max(AN, na.rm = TRUE), length.out = 100)) %>%
  tidybayes::add_fitted_draws(model_et_an_same_day) %>% 
  tidybayes::median_qi(.width = c(.95, .8, .5)) %>% 
  ungroup() %>%
  mutate_at(.vars = vars(AN,.value,.lower,.upper), .funs = ~ exp(.x)) %>%
  ggplot() + 
  aes(x = AN, y = .value) + 
  geom_lineribbon(aes(fill = factor(.width), ymin = .lower, ymax = .upper)) + 
  scale_color_manual(aesthetics = c("fill", "color"),
                     # I have no idea how this ordering works. Takes some trial and error
                     breaks = c("0.5", "0.8", "0.95", "median"),
                     values = c(lm_colors[2:4], lm_colors[1]), guide = FALSE) + 
  #guides(fill = guide_legend(title = "Posterior Intervals",
  #                           override.aes = list(color = c(NA, NA, NA, lm_colors[1]),
  #                                               fill = c(lm_colors[2:4], NA))),
  #       color = FALSE) +
  new_scale_color() +
  geom_point(data = et_an_asv_same_day_filter, aes(x = AN, y = ET, colour = seqvar_id)) +
  scale_colour_manual(values = rainbow(10)) +
  scale_x_log10(labels = scales::scientific) +
  scale_y_log10(labels = scales::scientific) +
  #coord_equal() +
  theme_bw() +
  #theme(legend.position = 'bottom') +
  labs(x = expression(paste("Anterior Nares ",italic("S. aureus"), " ASV Abundance")),
       y = expression(paste("Endotracheal ",italic("S. aureus"), " ASV Abundance")),
       colour = expression(paste(italic("S. aureus")," Amplicon Sequence Variant"))) -> p_et_v_an_asv_Sa_prop_model_biplot
p_et_v_an_asv_Sa_prop_model_biplot

# ggsave(plot = p_et_v_an_asv_Sa_prop_model_biplot, filename = "./figs/main/p_et_v_an_asv_Sa_prop_model_biplot.png", height = 4, width = 6, units = "in", dpi = 600)
# ggsave(plot = p_et_v_an_asv_Sa_prop_model_biplot, filename = "./figs/main/p_et_v_an_asv_Sa_prop_model_biplot.pdf", height = 4, width = 6, units = "in")




#' #######################################
#' FIGURE: compose all three for Figure 1 (using cowplot package)
#' #######################################

cowplot::get_legend(p_et_v_an_asv_Sa_prop_model_biplot) -> p_asv_legend

ggdraw(plot_grid(plot_grid(p_et_v_an_Sa_prop_histogram, p_et_v_an_Sa_prop_model_biplot, nrow = 1, rel_widths = c(1,1), labels = c("A","B")),
                 plot_grid(p_et_v_an_asv_Sa_prop_model_biplot + theme(legend.position = 'none'), p_asv_legend, nrow = 1, rel_widths = c(1.2,0.8), labels = c("C","")),
                 ncol = 1, rel_heights = c(1,1))) -> p_et_v_an_combined
p_et_v_an_combined
       
# ggsave(plot = p_et_v_an_combined, filename = "./figs/main/p_et_v_an_combined.png", height = 8, width = 8, units = "in", dpi = 600)
# ggsave(plot = p_et_v_an_combined, filename = "./figs/main/p_et_v_an_combined.svg", height = 8, width = 8, units = "in")
# ggsave(plot = p_et_v_an_combined, filename = "./figs/main/p_et_v_an_combined.pdf", height = 8, width = 8, units = "in")




#
###
#####
###
#








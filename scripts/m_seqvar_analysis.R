#' ########################################
#' binomial models relating Sa LRTI ~ AN microbiome (day 0)
#' ########################################
#' 
#' depends: binomial model results & trees
#' 
#' visualize tree with Staph and Corynebacterium effects
#' see: https://aschuerch.github.io/posts/2017-04-24-blog-post-1
#' for plot design


library(tidyverse)
library(cowplot)
library(ggnewscale)
library(ape)
library(vegan)
library(ggtree)

lrti_binomial_model_effects <- bind_rows(mutate(read_csv("./models/binomial/staph_lrti_model-AN_Sa_ASVs_microbiome_abx-reghs_scaled_reads_POOLED_seqvar.csv"), cohort = 'Pooled'),
                                              mutate(read_csv("./models/binomial/staph_lrti_model-AN_Sa_ASVs_microbiome_abx-reghs_scaled_reads_PILOT_seqvar.csv"), cohort = 'Pilot'),
                                              mutate(read_csv("./models/binomial/staph_lrti_model-AN_Sa_ASVs_microbiome_abx-reghs_scaled_reads_VALIDATION_seqvar.csv"), cohort = 'Validation')) %>%
  filter(grepl('beta',parameter))
lrti_binomial_model_effects



lrti_binomial_model_effects %>%
  #filter(cohort == "Pooled") %>%
  filter(grepl("staphylococcus",tolower(blast_genus))) %>%
  select(seqvar_id,`50%`,best,contains("typeS_10"),contains("typeS_5"),cohort) %>%
  distinct() %>%
  mutate(median_OR = exp(`50%`)) -> staph_effects
staph_effects



staph_effects %>%
  filter(typeS_10_neg|typeS_10_pos)



staph_effects %>%
  filter(typeS_5_neg|typeS_5_pos)



staph_effects %>%
  filter(cohort == "Pooled") -> staph_effects_pooled



ape::read.tree("./data/unrooted_tree.nwk") %>%
  ape::keep.tip(staph_effects_pooled$seqvar_id) %>%
  ggtree(., layout = 'rectangular') %<+% staph_effects_pooled + 
  geom_tiplab(aes(fill = median_OR),#aes(fill = factor(hit)),
              color = "black", # color for label font
              geom = "label",  # labels not text
              label.padding = unit(0.15, "lines"), # amount of padding around the labels
              label.size = 0) + # size of label border
  scale_fill_viridis_c(option = "C", direction = -1) +
  theme(legend.position = c(0.9,0.8),
        plot.title = element_text(hjust = 0.5) #, #legend.title = element_blank(), # no title #legend.key = element_blank()) # no keys
        ) +
  #theme(plot.margin = unit(c(1,1,1,1),"in")) +
  xlim_tree(0.2) +
  labs(fill = "Median\nPosterior\nOdds Ratio\nSa LRTI",
       title = expression(paste("Anterior Nares ",italic("Staphylococcus")," ASVs: Effect on Risk of ",italic("S. aureus")," LRTI (pooled cohorts)"))) -> p_staph_effects
p_staph_effects





lrti_binomial_model_effects %>%
  #filter(cohort == "Pooled") %>%
  filter(grepl("corynebacterium",tolower(blast_genus))) %>%
  select(seqvar_id,`50%`,best,contains("typeS_10"),contains("typeS_5"),cohort) %>%
  distinct() %>%
  mutate(median_OR = exp(`50%`)) -> coryn_effects
coryn_effects



coryn_effects %>%
  filter(typeS_10_neg|typeS_10_pos)



coryn_effects %>%
  filter(typeS_5_neg|typeS_5_pos)



coryn_effects %>%
  filter(cohort == "Pooled") -> coryn_effects_pooled




ape::read.tree("./data/unrooted_tree.nwk") %>%
  ape::keep.tip(coryn_effects_pooled$seqvar_id) %>%
  ggtree(., layout = 'rectangular') %<+% coryn_effects_pooled + 
  geom_tiplab(aes(fill = median_OR),#aes(fill = factor(hit)),
              color = "black", # color for label font
              geom = "label",  # labels not text
              label.padding = unit(0.15, "lines"), # amount of padding around the labels
              label.size = 0) + # size of label border
  scale_fill_viridis_c(option = "D", direction = -1) +
  theme(legend.position = c(0.1,0.8),
        plot.title = element_text(hjust = 0.5)#, #legend.title = element_blank(), # no title #legend.key = element_blank()) # no keys
  ) +
  #theme(plot.margin = unit(c(1,1,1,1),"in")) +
  xlim_tree(0.2) +
  labs(fill = "Median\nPosterior\nOdds Ratio\nSa LRTI",
       title = expression(paste("Anterior Nares ",italic("Corynebacterium")," ASVs: Effect on Risk of ",italic("S. aureus")," LRTI (pooled cohorts)"))) -> p_coryn_effects
p_coryn_effects



#' #######################################
#' FIGURE: compose all two plots (using cowplot package)
#' #######################################

ggdraw(plot_grid(p_staph_effects, p_coryn_effects, nrow = 2, rel_heights = c(11,11), labels = c("A","B"))) -> p_tree_combined
p_tree_combined

# ggsave(plot = p_tree_combined, filename = "./figs/supp/p_tree_combined.png", height = 22, width = 10, units = "in", dpi = 600)
# ggsave(plot = p_tree_combined, filename = "./figs/supp/p_tree_combined.svg", height = 22, width = 10, units = "in", system_fonts = list(sans = "Roboto"))
# ggsave(plot = p_tree_combined, filename = "./figs/supp/p_tree_combined.pdf", height = 22, width = 10, units = "in")




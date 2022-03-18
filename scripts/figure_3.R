## Rscript of plots in the paper
library(ggplot2)
library(tidyr)
library(readr)
library(dplyr)
library(ggsci)
library(cowplot)

# tandem duplication
count_clust_tdg = read_tsv("../data_scripts/TDG_cluster_count.tab", col_names = c("species", "cl_count"))
count_clust_tdg_corals = count_clust_tdg %>% filter(species == "Acropora millepora" |
                                                      species == "Acropora digitifera" |
                                                      species == "Montipora capitata" |
                                                      species == "Galaxea fascicularis" |
                                                      species == "Porites lobata" |
                                                      species == "Porites lutea" |
                                                      species == "Porites evermanni" |
                                                      species == "Pocillopora verrucosa" |
                                                      species == "Pocillopora damicornis" |
                                                      species == "Pocillopora meandrina" |
                                                      species == "Stylophora pistillata" |
                                                      species == "Goniastrea aspera" |
                                                      species == "Orbicella faveolata" |
                                                      species == "Fungia sp."
)

count_clust_tdg_corals = count_clust_tdg_corals %>% 
  mutate(clade = case_when(
    species %in% c("Acropora millepora", "Acropora digitifera", "Montipora capitata",
                   "Galaxea fascicularis", "Porites lobata", "Porites lutea", "Porites rus",
                   "Porites evermanni") ~ "Complex",
    species %in% c("Pocillopora verrucosa", "Pocillopora damicornis", "Pocillopora meandrina",
                   "Stylophora pistillata", "Goniastrea aspera", "Orbicella faveolata", "Fungia spp") ~ "Robust"
  ))

sum_clust_tdg_corals = count_clust_tdg_corals %>%
  group_by(species) %>%
  summarise(sum = sum(cl_count)) %>%
  mutate(clade = case_when(
    species %in% c("Acropora millepora", "Acropora digitifera", "Montipora capitata",
                   "Galaxea fascicularis", "Porites lobata", "Porites lutea", "Porites rus",
                   "Porites evermanni") ~ "Complex",
    species %in% c("Pocillopora verrucosa", "Pocillopora damicornis", "Pocillopora meandrina",
                   "Stylophora pistillata", "Goniastrea aspera", "Orbicella faveolata", "Fungia spp") ~ "Robust"
  ))


## Fig. 3 plots
# TDG
# TDG count
sum_clust_tdg_corals = sum_clust_tdg_corals %>%
  filter(species %in% c("Porites lobata", "Porites lutea", "Porites evermanni",
                        "Porites rus", "Pocillopora verrucosa",
                        "Pocillopora damicornis", "Pocillopora meandrina"))

p7 = ggplot(sum_clust_tdg_corals, aes(x = factor(species, levels = c("Porites lutea", "Porites rus", "Porites evermanni", "Porites lobata",
                                                                     "Pocillopora verrucosa", "Pocillopora damicornis", "Pocillopora meandrina")),
                                      y = sum
)) +
  geom_bar(aes(fill=clade), stat = "identity") +
  geom_label(aes(label = sum), position = position_stack(), size = 4) +
  scale_fill_futurama() +
  coord_flip() +
  theme_classic() +
  theme(axis.title.x=element_text(size=15), axis.text.x=element_text(size = 15),
        axis.title.y=element_blank(), axis.text.y=element_text(size = 15, face = "italic"),
        legend.title = element_blank(), legend.text = element_text(size=15)) + #, legend.position=c(0.8,0.9)
  scale_y_continuous("Number of TDG", breaks = seq(0,15000,2000))

# TDG cluster sizes
count_clust_tdg_corals = count_clust_tdg_corals %>%
  filter(species %in% c("Porites lobata", "Porites lutea", "Porites evermanni",
                        "Porites rus", "Pocillopora verrucosa",
                        "Pocillopora damicornis", "Pocillopora meandrina"))
p8 = ggplot(count_clust_tdg_corals, aes(x = factor(species, levels = c("Porites lobata", "Porites evermanni", "Porites lutea", "Porites rus",
                                                                       "Pocillopora meandrina", "Pocillopora verrucosa", "Pocillopora damicornis")),
                                        y = cl_count, color = clade)) +
  geom_boxplot() +
  scale_colour_futurama() +
  coord_flip() +
  theme_classic() +
  theme(axis.title.x=element_text(size=15), axis.text.x=element_text(size = 15),
        axis.title.y=element_blank(), axis.text.y=element_blank(),
        legend.title = element_blank(), legend.text = element_text(size=15), legend.position=c(0.8,0.9)) +
  scale_y_continuous("log10(number of TDG per cluster)", trans="log10", breaks = c(2,3,4,5,10,20,30,40,50,70), limits = c(2,80))


fig3=ggdraw() +
  draw_plot(p7, 0, 0, 0.57, 1) +
  draw_plot(p8, 0.57, 0, 0.43, 1) +
  draw_plot_label(c("A","B"), c(0,0.55), c(1,1), size = 15)
save_plot("fig3.svg", fig3, base_width = 15, base_height = 5)

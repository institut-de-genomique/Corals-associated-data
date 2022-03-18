library(ggplot2)
library(tidyr)
library(readr)
library(dplyr)
library(ggsci)
library(cowplot)

####################
## data
genome = data.frame(species = c("Porites lobata", "Porites lutea", "Porites evermanni", "Porites rus", "Acropora digitifera", "Acropora millepora", "Galaxea fascicularis", "Montipora capitata",
                                "Pocillopora meandrina", "Pocillopora damicornis", "Pocillopora verrucosa", "Pocillopora acuta", "Orbicella faveolata", "Stylophora pistillata", "Fungia sp.", "Goniastrea aspera",
                                "Aiptasia", "Nematostella vectensis",
                                "Amplexidiscus fenestrafer", "Discosoma sp.",
                                "Dendronephthya gigantea",
                                "Hydra vulgaris", "Clytia hemisphaerica",
                                "Aurelia aurita (ATL)"),
                    clade = c("Complex", "Complex", "Complex", "Complex", "Complex", "Complex", "Complex", "Complex",
                              "Robust", "Robust", "Robust", "Robust", "Robust", "Robust", "Robust", "Robust",
                              "Actinaria", "Actinaria",
                              "Corallimorpharia", "Corallimorpharia",
                              "Octocorallia",
                              "Hydrozoa", "Hydrozoa",
                              "Scyphozoa"),
                    size = c(646, 552, 604, 470, 416, 475, 334, 886,
                             347, 234, 381, 352, 485, 400, 606, 765,
                             255, 356,
                             370, 444,
                             286,
                             854, 445,
                             377),
                    kmer_estim = c(551, 610, 522, NA, 433, 443, 513, 652,
                                   313, NA, 332, NA, 448, 569, 642, 785,
                                   NA,NA,
                                   NA,NA,
                                   NA,
                                   NA,NA,
                                   NA),
                    n50 = c(2.154, 0.661, 0.171, 0.137, 1.097, 1.091, 0.088, 0.541,
                            4.754, 0.326, 0.334, 0.147, 1.162, 0.457, 0.323, 0.519,
                            0.441, 0.473,
                            0.510, 0.772,
                            1.446,
                            1.002, 0.366,
                            1.043
                    ),
                    tech = c("LR","SR", "SR", "SR","LR","LR","SR","SR",
                             "LR","SR","SR", "SR", "SR","SR","SR","SR",
                             "SR","SR",
                             "SR","SR",
                             "LR",
                             "SR","SR",
                             "SR"
                    ),
                    gene = c(42872, 31126, 40389, 39511, 26275, 28188, 22418, 63227,
                             32095, 26077, 27439, 79505, 32587, 27384, 38209, 35901,
                             26042, 24780,
                             21372, 23199,
                             28741,
                             36059, 26727,
                             38007)
)
genome$clade = factor(genome$clade, levels = c("Complex", "Robust", "Actinaria", "Corallimorpharia", "Octocorallia", "Hydrozoa", "Scyphozoa"))

# corals only
coral_genomes = genome %>%
  filter(clade=="Complex" | clade=="Robust")

# busco metrics (genome assembly) using metazoan dataset
busco_corail = read_tsv("../data_scripts/busco_all_metrics.tsv", col_names = c("species", "complete", "complete_single", "complete_dupli", "frag", "miss"))
busco_corail = busco_corail %>%
  mutate(tag = paste("S:", complete_single, ", D:", complete_dupli, ", F", frag, ", M", miss, sep = " ")) %>%
  select(-complete) %>%
  rename(Single = complete_single, Duplicated = complete_dupli, Fragmented = frag, Missed = miss) %>%
  gather("metric", "value", 2:5)

# interpro domains
ipr_count = read_tsv("../data_scripts/species_domains.count", col_names = c("species", "count", "total", "ratio"))
ipr_coral_count = ipr_count %>%
  filter(species %in% c("Acropora millepora", "Acropora digitifera", "Montipora capitata",
                        "Galaxea fascicularis", "Porites lobata", "Porites lutea", "Porites evermanni",
                        "Porites rus", "Pocillopora verrucosa", "Pocillopora acuta",
                        "Pocillopora damicornis", "Pocillopora meandrina", "Stylophora pistillata",
                        "Goniastrea aspera", "Orbicella faveolata", "Fungia spp")) %>%
  mutate(clade = case_when(
    species %in% c("Acropora millepora", "Acropora digitifera", "Montipora capitata",
                   "Galaxea fascicularis", "Porites lobata", "Porites lutea", "Porites rus",
                   "Porites evermanni") ~ "Complex",
    species %in% c("Pocillopora verrucosa", "Pocillopora damicornis", "Pocillopora meandrina", "Pocillopora acuta",
                   "Stylophora pistillata", "Goniastrea aspera", "Orbicella faveolata", "Fungia spp") ~ "Robust"))


## Ratio of genes included in OGs grouped by species
# Load gene count by OG and species
og_gene_count = read.table("../data_scripts/Orthogroups.GeneCount.tsv", header = TRUE)

# Load number of genes by species
gene_species = read.table("../data_scripts/gene_count_species", col.names=c("species", "total"))
gene_species$species = as.character(gene_species$species)

# Reshape og table
og_count = og_gene_count %>%
  gather("species", "count", 2:26) %>% # lines into colums
  filter(count != Total) %>% # discard OG containing only 1 species
  group_by(species) %>% # group by species names
  summarise(genes_in_og = sum(count)) # sum by species the number of genes in OGs

# merge the two tables, only with species described in og count
merge_count = inner_join(og_count, gene_species, by="species")

# ratio of genes contained in OG
merge_count$ratio = (merge_count$genes_in_og / merge_count$total) * 100

# clades
merge_count$species = gsub("_"," ",merge_count$species,fixed=TRUE)
merge_count = merge_count %>%
  mutate(clade = case_when(
    species %in% c("Acropora millepora", "Acropora digitifera", "Montipora capitata",
                   "Galaxea fascicularis", "Porites lobata", "Porites lutea", "Porites rus",
                   "Porites evermanni") ~ "Complex",
    species %in% c("Pocillopora verrucosa", "Pocillopora damicornis", "Pocillopora meandrina",
                   "Stylophora pistillata", "Goniastrea aspera", "Orbicella faveolata", "Fungia spp") ~ "Robust",
    species %in% c("Discosoma sp", "Amplexidiscus fenestrafer",
                   "Nematostella vectensis", "Aiptasia sp") ~ "Sea Anemones",
    species %in% c("Dendronephthya gigantea") ~ "Octocorallia",
    species %in% c("Aurelia aurita PAC", "Aurelia aurita ATL", "Morbakka virulenta",
                   "Hydra vulgaris", "Clytia hemisphaerica") ~ "Jellyfish"
  ))

# keep corals only and sort by clade and values
merge_count = merge_count %>%
  filter(clade %in% c("Complex", "Robust")) %>%
  #arrange(ratio)
  arrange(clade, ratio)

# add P. acuta (not included in orthofinder analysis)
acuta_count = data.frame(species = "Pocillopora acuta", genes_in_og = 0, total = 79505, ratio = 0, clade = "Robust")
merge_count = bind_rows(merge_count, acuta_count)



####### plots

## Fig 1
# genome size
p1 = ggplot(coral_genomes, aes(x = factor(species, levels = c("Acropora millepora", "Acropora digitifera", "Montipora capitata",
                                                              "Galaxea fascicularis", "Porites lobata", "Porites lutea", "Porites evermanni",
                                                              "Porites rus", "Pocillopora verrucosa", "Pocillopora acuta",
                                                              "Pocillopora damicornis", "Pocillopora meandrina", "Stylophora pistillata",
                                                              "Goniastrea aspera", "Orbicella faveolata", "Fungia sp."
)))) +
  geom_bar(aes(y = kmer_estim), stat = "identity", width = 0.35, color="#0BA498", fill="#0BA498", position = position_nudge(x=-0.3)) +
  geom_bar(aes(y = size, fill= clade), stat = "identity", position = position_nudge(x=0.1), width = 0.45) +
  geom_label(aes(y = size, label = size), size = 4, nudge_x = 0.3) +
  scale_fill_futurama() +
  coord_flip() +
  theme_classic() +
  theme(legend.title=element_blank(),
        axis.title.x=element_text(size=15), axis.text.x=element_text(size = 15),
        axis.title.y=element_blank(), axis.text.y=element_text(size = 15, face = "italic"),
        legend.position = "none") +
  scale_y_continuous("Genome assembly size (Mb)", breaks = seq(0,1000,200), limits = c(0,1000))

# N50
p2 = ggplot(coral_genomes, aes(x = factor(species, levels = c("Acropora millepora", "Acropora digitifera", "Montipora capitata",
                                                              "Galaxea fascicularis", "Porites lobata", "Porites lutea", "Porites evermanni",
                                                              "Porites rus", "Pocillopora verrucosa", "Pocillopora acuta",
                                                              "Pocillopora damicornis", "Pocillopora meandrina", "Stylophora pistillata",
                                                              "Goniastrea aspera", "Orbicella faveolata", "Fungia sp."
)))) +
  geom_bar(aes(y = n50, fill= clade), stat = "identity") +
  geom_label(aes(y = n50, label = n50), size = 4) +
  scale_fill_futurama() +
  coord_flip() +
  theme_classic() +
  theme(legend.title=element_blank(),
        axis.title.x=element_text(size=15), axis.text.x=element_text(size = 15),
        axis.title.y=element_blank(), axis.text.y=element_text(size = 15, face = "italic"),
        legend.text = element_text(size=15), legend.position = c(0.8, 0.2), legend.box = "horizontal") +
  scale_y_continuous("N50 (Mb)", limits = c(0,5))

# Gene count
p3 = ggplot(coral_genomes, aes(x = factor(species, levels = c("Acropora millepora", "Acropora digitifera", "Montipora capitata",
                                                              "Galaxea fascicularis", "Porites lobata", "Porites lutea", "Porites evermanni",
                                                              "Porites rus", "Pocillopora verrucosa", "Pocillopora acuta",
                                                              "Pocillopora damicornis", "Pocillopora meandrina", "Stylophora pistillata",
                                                              "Goniastrea aspera", "Orbicella faveolata", "Fungia sp."
)),
y = gene
)) +
  
  geom_bar(aes(fill=clade), stat = "identity") +
  geom_label(aes(y = gene, label = gene), size = 4) +
  scale_fill_futurama() +
  coord_flip() +
  theme_classic() +
  theme(legend.title=element_blank(),
        axis.title.x=element_text(size=15), axis.text.x=element_text(size = 15),
        axis.title.y=element_blank(), axis.text.y=element_text(size = 15, face = "italic"),
        legend.text = element_text(size=15)) +
  scale_y_continuous("Number of genes", breaks = seq(0,80000,20000), limits = c(0,80000))

# busco
p4 = ggplot(busco_corail, aes(x = factor(species, levels = c("Acropora millepora", "Acropora digitifera", "Montipora capitata",
                                                             "Galaxea fascicularis", "Porites lobata", "Porites lutea", "Porites evermanni",
                                                             "Porites rus", "Pocillopora verrucosa", "Pocillopora acuta",
                                                             "Pocillopora damicornis", "Pocillopora meandrina", "Stylophora pistillata",
                                                             "Goniastrea aspera", "Orbicella faveolata", "Fungia spp"
)),
y = value,
fill= factor(metric, levels = rev(c("Single", "Duplicated", "Fragmented", "Missed"))))) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = value), position = position_stack(vjust = 0.5), size = 4.5, color="#101010") +
  scale_fill_manual(values = c("#F21A00", "#E1AF00", "#78B7C5", "#3B9AB2")) +
  coord_flip() +
  theme_classic() +
  theme(legend.title=element_blank(),
        axis.title.x=element_text(size=15), axis.text.x=element_text(size = 15),
        axis.title.y=element_blank(), axis.text.y=element_text(size = 15, face = "italic"),
        legend.text = element_text(size=15), legend.position = "top"
  ) +
  scale_y_continuous("% of Metazoan BUSCO genes", breaks = seq(0,100,10))

# domains
p5 = ggplot(ipr_coral_count, aes(x = factor(species, levels = c("Acropora millepora", "Acropora digitifera", "Montipora capitata",
                                                                "Galaxea fascicularis", "Porites lobata", "Porites lutea", "Porites evermanni",
                                                                "Porites rus", "Pocillopora verrucosa", "Pocillopora acuta",
                                                                "Pocillopora damicornis", "Pocillopora meandrina", "Stylophora pistillata",
                                                                "Goniastrea aspera", "Orbicella faveolata", "Fungia spp"
)),
y = ratio)) +
  geom_bar(aes(fill=clade), stat = "identity") +
  geom_label(aes(y = ratio, label = ratio), size = 4) +
  scale_fill_futurama() +
  coord_flip() +
  theme_classic() +
  theme(legend.title=element_blank(),
        axis.title.x=element_text(size=15), axis.text.x=element_text(size = 15),
        axis.title.y=element_blank(), axis.text.y=element_blank(),
        legend.position = "none") +
  scale_y_continuous("% of genes with domains", breaks = seq(0,100,20), limits = c(0,100))

# nb of genes per OG
p6 = ggplot(merge_count, aes(x= factor(species, levels = c("Acropora millepora", "Acropora digitifera", "Montipora capitata",
                                                           "Galaxea fascicularis", "Porites lobata", "Porites lutea", "Porites evermanni",
                                                           "Porites rus", "Pocillopora verrucosa", "Pocillopora acuta",
                                                           "Pocillopora damicornis", "Pocillopora meandrina", "Stylophora pistillata",
                                                           "Goniastrea aspera", "Orbicella faveolata", "Fungia spp"
)),
y=ratio)) +
  geom_bar(aes(fill=clade), stat = "identity") +
  geom_label(aes(y = ifelse(species!="Pocillopora acuta",ratio,NA), label = round(ratio,1)), size = 4) +
  scale_fill_futurama() +
  coord_flip() +
  theme_classic() +
  theme(axis.title.x=element_text(size=15), axis.text.x=element_text(size = 15),
        axis.title.y=element_blank(), axis.text.y=element_blank(),
        legend.position = "none" ) +
  scale_y_continuous("% of genes in OGs", limits = c(0,100), breaks = seq(0,100,10))


fig1_1=ggdraw() +
  #             x, y, width, heigh
  draw_plot(p1, 0, 0.4, 1, 0.6) +
  draw_plot(p2, 0, 0, 1, 0.4) +
  #                                            X                       Y
  draw_plot_label(c("B","C"), c(0,0), c(1,0.5), size = 15)
save_plot("fig1_1.svg", fig1_1, base_width = 8, base_height = 10)

fig1_2=ggdraw() +
  #             x, y, width, heigh
  draw_plot(p3, 0, 0.6, 0.48, 0.4) +
  draw_plot(p5, 0.48, 0.6, 0.26, 0.4) +
  draw_plot(p6, 0.72, 0.6, 0.26, 0.4) +
  draw_plot(p4, 0, 0, 1, 0.6) +
  #                                            X                       Y
  draw_plot_label(c("D","E","F","G"), c(0,0.42,0.7,0), c(1,1,1,0.6), size = 15)
save_plot("fig1_2.svg", fig1_2, base_width = 15.9, base_height = 10)


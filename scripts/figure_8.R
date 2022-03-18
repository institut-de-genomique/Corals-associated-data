## Rscript of plots in the paper
library(ggplot2)
library(tidyr)
library(readr)
library(dplyr)
library(ggsci)
library(cowplot)
library(gplots)


## Fig 8
# Correlation table of TDG abundance profils across metatranscriptomic samples
exp_corr_pid_tdg = read_tsv("../data_scripts/pmea_all_sample.genes_sameTDG.corr_pid.tsv", col_names = c("id1", "id2", "exp_cor", "exp_pval", "pid", "ds"))

exp_corr_pid_tdg = exp_corr_pid_tdg %>%
  mutate(gr_pid = case_when(
    pid < 40 ~ "<40",
    pid >= 40 & pid < 50 ~ "40-<50",
    pid >= 50 & pid < 60 ~ "50-<60",
    pid >= 60 & pid < 70 ~ "60-<70",
    pid >= 70 & pid < 80 ~ "70-<80",
    pid >= 80 & pid < 90 ~ "80-<90",
    pid >= 90 ~ "90-100"
  ))

exp_corr_pid_tdg = exp_corr_pid_tdg %>%
  mutate(gr_ds = case_when(
    ds < 0.1 ~ "<0.1",
    ds >= 0.1 & ds < 0.2 ~ "0.1-<0.2",
    ds >= 0.2 & ds < 0.3 ~ "0.2-<0.3",
    ds >= 0.3 & ds < 0.4 ~ "0.3-<0.4",
    ds >= 0.4 & ds < 0.5 ~ "0.4-<0.5",
    ds >= 0.5 & ds < 0.6 ~ "0.5-<0.6",
    ds >= 0.6 & ds < 0.7 ~ "0.6-<0.7",
    ds >= 0.7 & ds < 0.8 ~ "0.7-<0.8",
    ds >= 0.8 & ds < 0.9 ~ "0.8-<0.9",
    ds >= 0.9 & ds < 1 ~ "0.9-<1",
    ds >= 1 ~ ">1"
  ))

p9=ggplot(exp_corr_pid_tdg, aes(x=factor(gr_ds,
                                         levels = c("<0.1","0.1-<0.2","0.2-<0.3","0.3-<0.4","0.4-<0.5","0.5-<0.6","0.6-<0.7","0.7-<0.8","0.8-<0.9","0.9-<1",">1"),
                                         ordered = TRUE),
                                y=exp_cor)) +
  geom_boxplot(varwidth = TRUE) +
  theme_classic() +
  theme(axis.title.x=element_text(size=30), axis.text.x=element_text(size = 30, angle = 30, vjust=0.5),
        axis.title.y=element_text(size=30), axis.text.y=element_text(size = 30)) +
  scale_x_discrete("dS") +
  scale_y_continuous("correlation of abundances", breaks = seq(-1,1,0.2), labels = c(-1,-0.8,-0.6,-0.4,-0.2,0,0.2,0.4,0.6,0.8,1)) #, breaks = seq(-1,1,0.2), labels = c(-1,-0.8,-0.6,-0.4,-0.2,0,0.2,0.4,0.6,0.8,1)
ggsave("fig_8B.svg", p9, width=18, height=8)



# Fig 8C
og_fonc = read_tsv("../data_scripts/OG.ordred.Categ.tsv", col_names = c("OG", "fonc", "rank"))
og_tpm = read_tsv("../data_scripts/og_amplif_from_supp_table_gene_tpm.tsv", col_names = c("sample", "gene", "OG", "TPM"))

og_gene_mean_tpm = og_tpm %>%
  group_by(sample,OG) %>%
  summarise(TPM = mean(TPM))

mat = inner_join(og_gene_mean_tpm, og_fonc, by = "OG") %>%
  filter(! is.na(fonc)) %>%
  spread(sample, TPM)
mat = mat[order(mat$rank),]
mtx = as.matrix(mat[, -1:-3])
rownames(mtx) = mat$OG

png("fig_8C.png", width=800, height = 800)
heatmap.2(mtx, trace = "none", scale = "row", Colv = F, Rowv = F, col = bluered(100))
dev.off()

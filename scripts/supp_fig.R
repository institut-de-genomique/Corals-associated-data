## Rscript of plots in the paper
library(ggplot2)
library(tidyr)
library(readr)
library(dplyr)
library(ggsci)
library(cowplot)
library(gplots)

# Fig S18
og_non_amplif_tpm = read_tsv("../data_scripts/og_non_amplif_tpm_all.tsv", col_names = c("sample", "gene", "OG", "TPM"))
og_gene_mean_tpm = og_non_amplif_tpm %>%
  group_by(sample,OG) %>%
  summarise(TPM = mean(TPM))

mat = og_gene_mean_tpm %>%
  filter(TPM>0) %>%
  spread(sample, TPM)
mat = na.omit(mat)
mat_samp = mat %>% sample_n(size = 200) # random sample of 200 OG
mtx = as.matrix(mat_samp[, -1])
rownames(mtx) = mat_samp$OG
png("fig_S18.png", width=800, height = 800)
heatmap.2(mtx, trace = "none", scale = "row", Colv = F, Rowv = F, col = bluered(100))
dev.off()

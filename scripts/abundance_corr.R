library(readr)
library(tidyr)
library(dplyr)
library(tibble)
library(corrplot)
library(Hmisc)


### data
# fichier de TPM/FPKM des mapping metaT sur tous les genes
abundances = read_tsv("../data_scripts/pmea_all_sample.genes.results", col_names = TRUE)

### calcul des correlations des niveaux d'abondances des genes
mtx_abund = abundances %>% select(sample, gene_id, TPM) %>%
  spread(gene_id, TPM) %>%
  select(-sample)
mtx_abund = as.data.frame(mtx_abund)

# correlation
cor_mtx_abund = cor(mtx_abund)

# pvalue
cormat = rcorr(cor_mtx_abund)

# cleaning
rm(abundances)
rm(mtx_abund)
rm(cor_mtx_abund)

# formatage en tab et print
cor_r = rownames_to_column(as.data.frame(cormat$r), var = "row")
cor_r = gather(cor_r, column, cor, -1)
cor_p = rownames_to_column(as.data.frame(cormat$P), var = "row")
cor_p = gather(cor_p, column, p, -1)
cor_p_matrix = left_join(cor_r, cor_p, by = c("row", "column"))
write_tsv(cor_p_matrix, "pmea_all_sample.genes.corr.tsv")


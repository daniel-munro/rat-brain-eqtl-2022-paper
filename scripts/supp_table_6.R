library(tidyverse)

genes <- read_tsv("data/genes.txt", col_types = "cc----------")

sig <- read_tsv("data/coloc/SMR_sig.tsv", col_types = "cccdcddd") |>
    left_join(genes, by = "gene_id", relationship = "many-to-one") |>
    relocate(trait, .after = tissue) |>
    relocate(gene_name, .after = gene_id) |>
    arrange(trait, gene_name, p_SMR)

write_tsv(sig, "tables/Supp_Table_6.txt")

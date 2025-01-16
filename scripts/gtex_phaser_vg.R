# phASER-derived Vg estimates from Marcela.

library(tidyverse)

vg <- list.files("data/anevah/gtex", pattern = "*.txt", recursive = TRUE,
                 full.names = TRUE) |>
    read_tsv(col_types = "cd-", id = "path") |>
    mutate(tissue = str_match(path, "gtex/(\\w+)/Common")[, 2], .after = path) |>
    select(-path)

kept <- read_csv("data/anevah/Vg_estimates_phaser_tuned.csv",
                 col_types = cols(.default = "c")) |>
    pivot_longer(everything(), names_to = "tissue", values_to = "GENE_ID") |>
    filter(!is.na(GENE_ID))

vg_filt <- vg |>
    inner_join(kept, by = c("tissue", "GENE_ID")) |>
    pivot_wider(GENE_ID, names_from = tissue, values_from = Vg_GeneWise)

write_tsv(vg_filt, "data/anevah/gtex_phaser_vg.tsv.gz")

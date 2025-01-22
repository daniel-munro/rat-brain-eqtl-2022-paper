library(tidyverse)

expr_rat <- read_tsv("data/expression/medianGeneExpression.txt.gz",
                     col_types = "cddddddddd") |>
    pivot_longer(-geneId, names_to = "tissue", values_to = "tpm") |>
    filter(tissue %in% c("IL", "LHb", "NAcc", "OFC", "PL")) |>
    rename(gene_id = geneId) |>
    summarise(n_expr = sum(tpm > 1), .by = gene_id) |>
    filter(n_expr >= 1)

expr_gtex <- read_tsv(
    "data/gtex/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct.gz",
    skip = 2,
    col_types = cols(Name = "c", Description = "-", .default = "d")
) |>
    pivot_longer(-Name, names_to = "tissue", values_to = "tpm") |>
    filter(str_sub(tissue, 1, 5) == "Brain") |>
    mutate(gene_id = str_replace(Name, "\\..+$", "")) |>
    summarise(n_expr = sum(tpm > 1), .by = gene_id) |>
    filter(n_expr >= 1)

eqtl_rat <- read_tsv("data/eqtls/top_assoc.txt", col_types = "cc-------------d--") |>
    filter(gene_id %in% expr_rat$gene_id) |>
    summarise(eqtl_rat = any(qval < 0.05), .by = gene_id)

eqtl_gtex <- list.files("data/gtex/GTEx_Analysis_v8_eQTL", pattern = "*Brain_*",
                        full.names = TRUE) |>
    read_tsv(id = "file",
             col_types = cols(gene_id = "c", gene_name = "c", qval = "d", .default = "-")) |>
    mutate(gene_id = str_replace(gene_id, "\\..+$", "")) |>
    filter(gene_id %in% expr_gtex$gene_id) |>
    summarise(eqtl_gtex = any(qval <= 0.05),
              .by = c(gene_id, gene_name))

pli <- read_tsv("data/human/forweb_cleaned_exac_r03_march16_z_data_pLI.txt.gz",
                col_types = "-c-----------------d")

d <- read_tsv("data/gtex/orthologs.txt", col_types = "cc") |>
    mutate(gene_id_human = str_replace(gene_id_human, "\\..+$", "")) |>
    inner_join(eqtl_rat, by = c("gene_id_rat" = "gene_id")) |>
    inner_join(eqtl_gtex, by = c("gene_id_human" = "gene_id")) |>
    inner_join(pli, by = c("gene_name" = "gene")) |>
    mutate(
        eqtl = case_when(
            eqtl_rat & !eqtl_gtex ~ "HS rat\nonly",
            !eqtl_rat & eqtl_gtex ~ "GTEx\nonly",
            eqtl_rat & eqtl_gtex ~ "Both",
            !eqtl_rat & !eqtl_gtex ~ "Neither",
        ) |>
            fct_relevel("Both", "GTEx\nonly", "HS rat\nonly")
    )

d_stats <- d |>
    summarise(mean_pLI = mean(pLI),
              SE_pLI = sd(pLI) / sqrt(length(pLI)),
              n = n(),
              .by = eqtl)

##############################
## Point plot: figure panel ##
##############################

d_stats |>
    mutate(group = str_glue("{eqtl}\n({n})")) |>
    ggplot(aes(x = group, y = mean_pLI, ymin = mean_pLI - SE_pLI,
               ymax = mean_pLI + SE_pLI)) +
    geom_pointrange(size = 0.75, fatten = 1.5, show.legend = FALSE) +
    theme_bw() +
    theme(
        panel.grid = element_blank(),
        axis.text.x = element_text(hjust = 0.5),
    ) +
    xlab("Ortholog pairs grouped by\ncis-eQTL status") +
    ylab("pLI (mean ± SE)")

ggsave("figures/figure5/figure5c.png", width = 2.5, height = 2.5)

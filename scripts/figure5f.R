library(tidyverse)

rat <- tibble(tissue = c("IL", "LHb", "NAcc", "OFC", "PL")) |>
    reframe(
        read_tsv(str_glue("data/anevah/output/Vg.{tissue}.tsv.gz"), col_types = "cd"),
        .by = tissue
    ) |>
    group_by(GENE_ID) |>
    summarise(SDg_rat = sqrt(mean(Vg)))

gtex <- read_tsv("data/anevah/gtex_phaser_vg.tsv.gz", col_types = cols(GENE_ID = "c", .default = "d")) |>
    pivot_longer(-GENE_ID, names_to = "tissue", values_to = "Vg") |>
    filter(str_sub(tissue, 1, 3) == "BRN",
           !is.na(Vg)) |>
    summarise(SDg_GTEx = sqrt(mean(Vg)), .by = GENE_ID)

gene_map <- read_tsv("data/gtex/orthologs.txt", col_types = "cc") |>
    mutate(gene_id_human = str_replace(gene_id_human, "\\..+$", ""))

d <- gene_map |>
    inner_join(rat, by = c("gene_id_rat" = "GENE_ID")) |>
    inner_join(gtex, by = c("gene_id_human" = "GENE_ID")) |>
    mutate(SDg_rat = pmin(SDg_rat, 0.5),
           SDg_GTEx = pmin(SDg_GTEx, 0.5))
write_tsv(d, "data/gtex/orthologs/ortholog_SDg.tsv")

stats <- d |>
    summarise(
        pairs = n(),
        pearson = cor(SDg_rat, SDg_GTEx),
        pearson_p = cor.test(SDg_rat, SDg_GTEx)$p.value,
        spearman = cor(SDg_rat, SDg_GTEx, method = "spearman"),
        spearman_p = cor.test(SDg_rat, SDg_GTEx, method = "spearman")$p.value,
    )
lab <- str_c(str_glue(" r = {format(stats$pearson, digits = 2)}"),
             str_glue("P = {format(stats$pearson_p, digits = 2)}"),
             sep = "\n")

d |>
    ggplot(aes(x = SDg_rat, y = SDg_GTEx)) +
    geom_point(size = 0.5, alpha = 0.3) +
    annotate(geom = "text", x = 0.2, y = 0.45, hjust = 0, label = lab, color = "#5555cc") +
    coord_fixed() +
    expand_limits(x = 0, y = 0) +
    scale_x_continuous(expand = c(0, 0.01)) +
    scale_y_continuous(expand = c(0, 0.01)) +
    theme_bw() +
    theme(panel.grid = element_blank()) +
    xlab(expression(SD^G*" (Rat)")) +
    ylab(expression(SD^G*" (Human)"))

ggsave("figures/figure5/figure5f.png", width = 2.5, height = 2.5)

###########
## Stats ##
###########

with(d, cor.test(SDg_rat, SDg_GTEx, method = "spearman"))

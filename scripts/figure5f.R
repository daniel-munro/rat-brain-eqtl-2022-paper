library(tidyverse)

# # Copied median TPM table from RatGTEx for convenience:
# expressed <- read_tsv("data/expression/medianGeneExpression.txt.gz", col_types = "c--ddd-dd-") |>
#     rename(GENE_ID = geneId) |>
#     pivot_longer(-GENE_ID, names_to = "tissue", values_to = "TPM") |>
#     filter(!is.na(TPM),
#            TPM > 10) |>
#     select(-TPM)
## Decided not to use, didn't have much effect.

rat <- tibble(tissue = c("IL", "LHb", "NAcc", "OFC", "PL")) |>
    group_by(tissue) |>
    summarise(
        read_tsv(str_glue("data/anevah/output/Vg.{tissue}.tsv.gz"), col_types = "cd"),
        .groups = "drop"
    ) |>
    # inner_join(expressed, by = c("tissue", "GENE_ID")) |>
    group_by(GENE_ID) |>
    summarise(SDg_rat = sqrt(mean(Vg)))
    # summarise(SDg_rat = sqrt(1 / mean(1 / Vg))) # Harmonic mean
    ## Randomly choose one value:
    # slice_sample(n = 1) |>
    # mutate(SDg_rat = sqrt(Vg)) |>
    # ungroup()

gtex <- read_tsv("data/anevah/gtex_phaser_vg.tsv.gz", col_types = cols(GENE_ID = "c", .default = "d")) |>
    pivot_longer(-GENE_ID, names_to = "tissue", values_to = "Vg") |>
    filter(str_sub(tissue, 1, 3) == "BRN",
           !is.na(Vg)) |>
    group_by(GENE_ID) |>
    summarise(SDg_GTEx = sqrt(mean(Vg)))
    # summarise(SDg_GTEx = sqrt(1 / mean(1 / Vg))) # Harmonic mean
    ## Randomly choose one value:
    # slice_sample(n = 1) |>
    # mutate(SDg_GTEx = sqrt(Vg)) |>
    # ungroup()

gene_map <- read_tsv("data/gtex/orthologs.txt", col_types = "cc") |>
    mutate(gene_id_human = str_replace(gene_id_human, "\\..+$", ""))

d <- gene_map |>
    inner_join(rat, by = c("gene_id_rat" = "GENE_ID")) |>
    inner_join(gtex, by = c("gene_id_human" = "GENE_ID")) |>
    mutate(SDg_rat = pmin(SDg_rat, 0.5),
           SDg_GTEx = pmin(SDg_GTEx, 0.5))
# write_tsv(d, "data/gtex/orthologs/ortholog_SDg.tsv")

# deming <- d |>
#     summarise({
#         coef <- deming::deming(SDg_GTEx ~ SDg_rat)$coefficients
#         tibble(intercept = coef[1],
#                slope = coef[2])
#     })
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
    # geom_abline(aes(intercept = intercept, slope = slope), data = deming, color = "#8888cc") +
    annotate(geom = "text", x = 0.2, y = 0.45, hjust = 0, label = lab, color = "#5555cc") +
    coord_fixed() +
    expand_limits(x = 0, y = 0) +
    scale_x_continuous(expand = c(0, 0.01)) +
    scale_y_continuous(expand = c(0, 0.01)) +
    theme_bw() +
    theme(panel.grid = element_blank()) +
    xlab(expression(SD^G*" (Rat)")) +
    ylab(expression(SD^G*" (Human)"))
    # ggtitle(str_glue("Pearson's r = {with(mean_ortho, cor(SDg_GTEx, SDg_rat) |> round(4))}"))

ggsave("figures/figure5/figure5f.png", width = 2.5, height = 2.5)

###########
## Stats ##
###########

with(d, cor.test(SDg_rat, SDg_GTEx, method = "spearman"))

library(tidyverse)

load_gtex <- function(dir) {
    list.files(dir, full.names = TRUE) |>
        read_tsv(col_types = "c-d-d")
}

rat <- tibble(tissue = c("IL", "LHb", "NAcc", "OFC", "PL")) |>
    group_by(tissue) |>
    summarise(
        read_tsv(str_glue("data/gemma/{tissue}.cis_pve.txt"), col_types = "cd-"),
        .groups = "drop"
    ) |>
    filter(!is.na(pve)) |>
    mutate(pve = pmax(0, pmin(pve, 1))) |>
    group_by(gene_id) |>
    summarise(mean_h2_rat = mean(pve))

gene_map <- read_tsv("data/gtex/orthologs.txt", col_types = "cc") |>
    mutate(gene_id_human = str_replace(gene_id_human, "\\..+$", ""))

# gtex_genes <- read_tsv("data/gtex/gtex_genes.txt", col_names = c("gene_id", "Gene"),
#                        col_types = "cc")
human_genes <- read_tsv("data/gtex/human_genes.txt", col_types = "c-c") |>
    rename(gene_id = `Gene stable ID`,
           Gene = `Gene name`)

gtex <- tibble(dir = list.dirs("data/gtex/heritability/GTEx_v8_h2")) |>
    filter(str_detect(dir, "Brain")) |>
    group_by(dir) |>
    summarise(load_gtex(dir), .groups = "drop") |>
    filter(!is.na(h2cis)) |>
    mutate(h2cis = pmax(0, pmin(h2cis, 1))) |>
    left_join(human_genes, by = "Gene") |>
    group_by(dir, gene_id) |>
    summarise(h2cis = mean(h2cis), .groups = "drop") |>
    group_by(gene_id) |>
    summarise(mean_h2_human = mean(h2cis), .groups = "drop")

d <- gene_map |>
    inner_join(rat, by = c("gene_id_rat" = "gene_id")) |>
    inner_join(gtex, by = c("gene_id_human" = "gene_id"))
# write_tsv(d, "data/gtex/orthologs/ortholog_h2.tsv")

# deming <- d |>
#     slice_sample(n = 2000) |>
#     summarise({
#         coef <- deming::deming(mean_h2_human ~ mean_h2_rat)$coefficients
#         tibble(intercept = coef[1],
#                slope = coef[2])
#     })
stats <- d |>
    summarise(
        pairs = n(),
        pearson = cor(mean_h2_rat, mean_h2_human),
        pearson_p = cor.test(mean_h2_rat, mean_h2_human)$p.value,
        spearman = cor(mean_h2_rat, mean_h2_human, method = "spearman"),
        spearman_p = cor.test(mean_h2_rat, mean_h2_human, method = "spearman")$p.value,
    )
lab <- str_c(
    #str_glue(" r = {format(stats$pearson, digits = 2)}"),
    str_glue(" r = {round(stats$pearson, 2) |> format(nsmall = 2)}"),
    str_glue("P = {format(stats$pearson_p, digits = 2)}"),
    sep = "\n"
)
d |>
    ggplot(aes(x = mean_h2_rat, y = mean_h2_human)) +
    geom_point(size = 0.5, alpha = 0.3) +
    # geom_abline(aes(intercept = intercept, slope = slope), data = deming, color = "#8888cc") +
    annotate(geom = "text", x = 0.55, y = 0.9, hjust = 0, label = lab, color = "#5555cc") +
    coord_fixed() +
    scale_x_continuous(expand = c(0, 0.01)) +
    scale_y_continuous(expand = c(0, 0.01)) +
    expand_limits(x = 0:1, y = 0:1) +
    theme_bw() +
    theme(panel.grid = element_blank()) +
    xlab(expression("Mean "*h^2*" (Rat)")) +
    ylab(expression("Mean "*h^2*" (Human)"))

ggsave("figures/figure5/figure5e.png", width = 2.7, height = 2.5)

with(d, cor.test(mean_h2_rat, mean_h2_human, method = "spearman")) |> unlist()

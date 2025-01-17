library(tidyverse)

CI <- function(values, fun, R = 1e4) {
    b <- boot::boot(values, function(d, i) fun(d[i]), R)
    ci <- boot::boot.ci(b, conf = 0.95, type = "norm")
    ci$normal[2:3] |> set_names(c("CI_low", "CI_high"))
}

d <- read_tsv("data/gtex/orthologs/ortholog_SDg.tsv", col_types = "ccdd")

gsets <- read_tsv("data/gtex/orthologs/gene_sets.tsv", col_types = "cc")

d_sets <- d |>
    inner_join(gsets, by = "gene_id_human") |>
    select(-gene_id_rat,-gene_id_human) |>
    pivot_longer(
        c(SDg_rat, SDg_GTEx),
        names_prefix = "SDg_",
        names_to = "organism",
        values_to = "SDg"
    ) |>
    mutate(organism = c(rat = "Rat", GTEx = "Human")[organism]) |>
    group_by(set, organism) |>
    summarise(
        SDg_median = median(SDg),
        CI = CI(SDg, median) |> list(),
        n = n(),
        .groups = "drop"
    ) |>
    unnest_wider(CI) |>
    mutate(set = set |>
               fct_reorder2(SDg_median, organism,
                            function (x, y) median(x[y == "Human"]),
                            .desc = FALSE))
median_all <- d |>
    select(SDg_rat, SDg_GTEx) |>
    pivot_longer(c(SDg_rat, SDg_GTEx), names_prefix = "SDg_",
                 names_to = "organism", values_to = "SDg") |>
    mutate(organism = c(rat = "Rat", GTEx = "Human")[organism]) |>
    group_by(organism) |>
    summarise(median = median(SDg))

size_labels <- d_sets |>
    filter(organism == "Human") |>
    mutate(label = str_glue("({n})"))

d_sets |>
    ggplot(aes(x = set, y = SDg_median, ymin = CI_low, ymax = CI_high, color = organism)) +
    geom_hline(aes(yintercept = median, color = organism), data = median_all, alpha = 0.5) +
    geom_pointrange(fatten = 1.5, position = position_dodge(width = 0.5)) +
    geom_text(aes(color = NULL, label = label), data = size_labels,
              y = 0.002, angle = 90, hjust = 0, size = 3) +
    scale_color_manual(values = c("#619CFFFF", "#F8766DFF")) +
    expand_limits(y = 0.004) +
    theme_bw() +
    theme(
        axis.text.x = element_text(hjust = 1, vjust = 1, angle = 30),
        panel.grid = element_blank(),
        legend.title = element_blank(),
        legend.position = "left",
        legend.margin = margin(0, -5, 0, 0, unit = "pt"),
    ) +
    xlab("Human gene sets") +
    ylab(expression(SD^G))

ggsave("figures/figure5/figure5g.png", width = 7, height = 2.9)

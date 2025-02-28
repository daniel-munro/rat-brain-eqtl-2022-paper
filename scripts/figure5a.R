library(tidyverse)

top_assoc <- read_tsv("data/eqtls/top_assoc.txt", col_types = "ccciiciiccdddddddd")
egenes <- filter(top_assoc, qval < 0.05)

sample_size <- read_tsv("data/samples.txt", col_types = "cc") |>
    count(brain_region) |>
    deframe()

genes_tested <- read_tsv("data/genes.txt", col_types = "c------lllll") |>
    reframe(tribble(
        ~tissue, ~expressedGeneCount,
        "IL", sum(in_expr_IL),
        "LHb", sum(in_expr_LHb),
        "NAcc", sum(in_expr_NAcc),
        "OFC", sum(in_expr_OFC),
        "PL", sum(in_expr_PL),
    ))

rat_tissues <- egenes |>
    count(tissue, name = "eGeneCount") |>
    mutate(group = "HS rat brain",
           sample_size = sample_size[tissue]) |>
    left_join(genes_tested, by = "tissue")

gtex_tissues <- jsonlite::fromJSON("data/gtex/tissueInfo.json")[[1]] |>
    filter(hasEGenes) |>
    select(tissue = tissueSiteDetailId,
           sample_size = rnaSeqAndGenotypeSampleCount,
           eGeneCount,
           expressedGeneCount) |>
    mutate(group = if_else(str_detect(tissue, "Brain"), "Human brain", "Human non-brain"))

gtex_subsam <- gtex_tissues |>
    filter(group == "Human brain") |>
    summarise(
        read_tsv(str_glue("data/gtex/subsample/{tissue}.cis_qtl.txt.gz"),
                 col_types = cols(qval = "d", .default = "-")) |>
            summarise(eGeneCount = sum(qval <= 0.05)),
        .by = tissue
    ) |>
    mutate(sample_size = 81L,
           group = "Human brain\n(subsampled)")

#################################################
## Number of eGenes in relation to sample size ##
#################################################

colors <- RColorBrewer::brewer.pal(3, "Set1") # Darker, swap blue and green so brain is blue
bind_rows(rat_tissues, gtex_tissues) |>
    ggplot(aes(x = sample_size, y = eGeneCount, color = group, shape = group)) +
    geom_point(size = 1) +
    scale_color_manual(values = colors) +
    geom_boxplot(aes(shape = NULL), data = gtex_subsam, color = colors[2],
                 width = 0.4, show.legend = FALSE) +
    scale_x_sqrt(breaks = 100 * 1:7,
                 labels = c(100, 200, 300, 400, "", 600, "")) +
    scale_y_sqrt(breaks = 2000 * 1:10,
                 labels = c(2000, 4000, 6000, 8000, 10000, "", 14000, "", 18000, "")) +
    xlab("Sample size") +
    ylab("eGene count") +
    theme_bw() +
    theme(
        legend.title = element_blank(),
        legend.position = c(0.7, 0.15),
        legend.background = element_rect(color = "black", linewidth = 0.1),
        legend.spacing.y = unit(-2, "pt"),
        legend.key.width = unit(2, "pt"),
        legend.key.height = unit(8, "pt"),
        legend.margin = margin(t = 5, b = 3, l = 3, r = 3),
        legend.text = element_text(size = 8),
        panel.grid = element_blank(),
    )

ggsave("figures/figure5/figure5a.png", width = 2.7, height = 2.5)

bind_rows(rat_tissues, gtex_tissues) |>
    mutate(eGenes_per_sample = eGeneCount / sample_size,
           eGenes_per_sample_per_1k_genes = eGenes_per_sample / (expressedGeneCount / 1000)) |>
    group_by(group) |>
    skimr::skim()

with(gtex_tissues, cor(sample_size, eGeneCount))

gtex_subsam |>
    summarise(mean = mean(eGeneCount),
              sd = sd(eGeneCount))

###########
## Stats ##
###########

gtex_tissues |>
    filter(group == "Human brain") |>
    summarise(mean_n = mean(sample_size),
              sd_n = sd(sample_size),
              range_n = list(range(sample_size)))

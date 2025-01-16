library(tidyverse)
library(patchwork)

probs <- function(prob) {
    names(dimnames(prob)) <- c("individual", "Strain", "SNP")
    cubelyr::as.tbl_cube(prob) |>
        as_tibble() |>
        mutate(Strain = strains[Strain])
}

strains <- c(
    A = "ACI/N", B = "BN/SsN", C = "BUF/N", D = "F344/N",
    E = "M520/N", F = "MR/N", G = "WN/N", H = "WKY/N"
)

d <- tibble(chrom = 1:20) |>
    group_by(chrom) |>
    summarise(
        readRDS(str_glue("data/haplotypes/haplotype_probs_chr{chrom}.rds")) |>
            probs(),
        .groups = "drop"
    )

mean_probs <- d |>
    group_by(SNP, Strain) |>
    summarise(prob = mean(prob), .groups = "drop") |>
    separate(SNP, c("chrom", "pos"), sep = ":", convert = TRUE) |>
    mutate(pos = pos / 1e6,
           chrom = factor(chrom, levels = str_c("chr", 1:20)))

mean_probs |>
    # filter(chrom %in% str_c("chr", 1:6)) |>
    ggplot(aes(x = pos, y = prob, fill = Strain)) +
    facet_wrap(~ chrom, ncol = 2, dir = "v", strip.position = "left") +
    geom_area() +
    scale_fill_brewer(type = "qual", palette = 6) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    theme_bw() +
    xlab("Position (Mb)") +
    ylab("Estimated Proportion") +
    theme(
        axis.text.y = element_blank(),
        # panel.grid.major.y = element_blank(),
        # panel.grid.minor.y = element_blank(),
        panel.grid = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.x = element_line(size = 0.25),
        panel.border = element_blank(),
        strip.background = element_blank(),
        legend.position = "top",
    )

ggsave("figures/figureS1/figureS1_raw.png", width = 11, height = 8, bg = "white")

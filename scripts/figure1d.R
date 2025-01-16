library(VariantAnnotation)
library(tidyverse)

d <- info(readVcf("data/genotype/P50.rnaseq.88.unpruned.vcf.gz")) |>
    as_tibble() |>
    mutate(AC = as.integer(AC),
           AF = AC / AN,
           MAF = pmin(AF, 1 - AF))

# d |>
#     ggplot(aes(x = MAF)) +
#     geom_histogram(bins = 30, boundary = 0) +
#     scale_y_continuous(breaks = c(50000, 100000, 150000),
#                        labels = c(50, 100, 150),
#                        expand = c(0, 0)) +
#     scale_x_continuous(expand = c(0, 0)) +
#     theme_minimal() +
#     theme(panel.grid = element_blank()) +
#     ylab("SNPs (thousands)")

gtex <- tibble(all = read_lines("data/gtex/GTEx_AF.subsample.txt.gz")) |>
    mutate(AF = str_match(all, "AF=([0-9.]+);")[, 2] |> as.double(),
           MAF = pmin(AF, 1 - AF))

# # Overlapping:
# bind_rows(
#     d |> mutate(data = "HS rats"),
#     gtex |> mutate(data = "GTEx")
# ) |>
#     ggplot(aes(x = MAF, fill = data)) +
#     geom_histogram(bins = 30, boundary = 0, position = "identity", alpha = 0.5) +
#     scale_y_continuous(breaks = c(1e6, 2e6, 3e6),
#                        labels = c(1, 2, 3),
#                        expand = c(0, 0)) +
#     scale_x_continuous(expand = c(0, 0)) +
#     theme_minimal() +
#     theme(panel.grid = element_blank(),
#           legend.title = element_blank(),
#           legend.position = c(0.6, 0.5)) +
#     ylab("SNPs (millions)")

# Density instead of histograms:
max_r <- max(density(d$MAF, adjust = 5)$y)
max_g <- max(density(gtex$MAF, adjust = 5)$y)
bind_rows(
    d |> mutate(data = "HS rats"),
    gtex |> mutate(data = "GTEx"),
) |>
    mutate(data = fct_rev(data)) |>
    ggplot(aes(x = MAF, y = ..scaled.., color = data)) +
    geom_density(adjust = 5, size = 0.8, outline.type = "full", show.legend = FALSE) +
    scale_x_continuous(expand = c(0.02, 0)) +
    scale_y_continuous(expand = c(0.02, 0), breaks = c(0, 50, 100, 150) / max_g,
                       sec.axis = sec_axis(~ ., breaks = c(0, 1, 2) / max_r)) +
    scale_color_manual(values = c("#F8766DFF", "#619CFFFF")) + # from: (scales::hue_pal())(3)[c(1, 3)]
    theme_bw() +
    theme(
        panel.grid = element_blank(),
        # panel.border = element_blank(),
        # axis.line.x = element_line(),
        # axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
    ) +
    ylab("") # Line up with LD decay panel

ggsave("figures/figure1/figure1d.png", width = 2.6, height = 2.5, dpi = 300)

# # Separate panels:
# bind_rows(
#     d |> mutate(data = "88 HS rats"),
#     gtex |> mutate(data = "GTEx (subsampled SNPs)")
# ) |>
#     ggplot(aes(x = MAF)) +
#     facet_wrap(~data, ncol = 1, scales = "free_y") +
#     geom_histogram(bins = 30, boundary = 0) +
#     # scale_y_continuous(breaks = c(2e4, 4e4, 6e4, 1e6, 2e6, 3e6),
#     #                    labels = c(0.02, 0.04, 0.06, 1, 2, 3),
#     scale_y_continuous(breaks = c(5e4, 1e5, 1.5e5, 2e5, 1e6, 2e6, 3e6),
#                        labels = c(0.05, 0.1, 0.15, 0.2, 1, 2, 3),
#     expand = c(0, 0)) +
#     scale_x_continuous(expand = c(0, 0)) +
#     theme_minimal() +
#     theme(panel.grid = element_blank()) +
#     ylab("SNPs (millions)")
# 
# ggsave("genotypes/MAF2.png", width = 3, height = 2.5, dpi = 300)

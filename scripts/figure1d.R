library(VariantAnnotation)
library(tidyverse)

d <- info(readVcf("data/genotype/P50.rnaseq.88.unpruned.vcf.gz")) |>
    as_tibble() |>
    mutate(AC = as.integer(AC),
           AF = AC / AN,
           MAF = pmin(AF, 1 - AF))

gtex <- tibble(all = read_lines("data/gtex/GTEx_AF.subsample.txt.gz")) |>
    mutate(AF = str_match(all, "AF=([0-9.]+);")[, 2] |> as.double(),
           MAF = pmin(AF, 1 - AF))

# Density plot
max_r <- max(density(d$MAF, adjust = 5)$y)
max_g <- max(density(gtex$MAF, adjust = 5)$y)
bind_rows(
    d |> mutate(data = "HS rats"),
    gtex |> mutate(data = "GTEx"),
) |>
    mutate(data = fct_rev(data)) |>
    ggplot(aes(x = MAF, y = after_stat(scaled), color = data)) +
    geom_density(adjust = 5, size = 0.8, outline.type = "full", show.legend = FALSE) +
    scale_x_continuous(expand = c(0.02, 0)) +
    scale_y_continuous(expand = c(0.02, 0), breaks = c(0, 50, 100, 150) / max_g,
                       sec.axis = sec_axis(~ ., breaks = c(0, 1, 2) / max_r)) +
    scale_color_manual(values = c("#F8766DFF", "#619CFFFF")) + # from: (scales::hue_pal())(3)[c(1, 3)]
    theme_bw() +
    theme(
        panel.grid = element_blank(),
        axis.text.y = element_blank(),
    ) +
    ylab("") # Line up with LD decay panel

ggsave("figures/figure1/figure1d.png", width = 2.6, height = 2.5, dpi = 300)

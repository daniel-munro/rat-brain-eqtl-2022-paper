library(VariantAnnotation)
library(tidyverse)

N_BINS <- 100
MAX_DIST <- 1e7

ld_rat <- read_tsv("data/genotype/LD.txt.gz", col_types = "ciiid")
ld_human <- read_tsv("data/gtex/GTEx_LD.txt.gz", col_types = "ciiid")

ld <- bind_rows(
    ld_rat |> mutate(data = "HS rats", .before = 1),
    ld_human |> mutate(data = "GTEx", .before = 1)
) |>
    mutate(data = fct_rev(data))

ld |>
    slice_sample(n = 1e4) |>
    ggplot(aes(x = distance / 1e6, y = r2, color = data)) +
    geom_point(alpha = 0.5) +
    stat_smooth() +
    theme_bw()

ld_mean <- ld |>
    mutate(dist_bin = cut_interval(distance, N_BINS, labels = FALSE)) |>
    group_by(data, dist_bin) |>
    summarise(r2_mean = mean(r2), .groups = "drop") |>
    mutate(bin_mid = ((dist_bin - 0.5) / N_BINS) * MAX_DIST / 1e6)

ld_mean |>
    ggplot(aes(x = bin_mid, y = r2_mean, color = data)) +
    geom_line(size = 0.8, show.legend = FALSE) +
    scale_x_continuous(breaks = c(0, 2, 4, 6, 8, 10),
                       expand = c(0.01, 0.01)) +
    scale_y_continuous(expand = c(0.002, 0.002)) +
    scale_color_manual(values = (scales::hue_pal())(3)[c(1, 3)]) +
    expand_limits(x = c(0, 1), y = c(0, 1)) +
    theme_bw() +
    theme(
        panel.grid = element_blank(),
    ) +
    xlab("Distance (Mb)") +
    ylab(expression("Average LD "*(r^2)))

ggsave("figures/figure1/figure1c.png", width = 2.8, height = 2.5)

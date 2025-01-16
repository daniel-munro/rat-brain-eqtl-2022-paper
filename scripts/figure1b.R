library(tidyverse)

N_CHROMS <- 20
SEGS_PER_CHROM <- 50 # On average

set.seed(20210928)
segs <- N_CHROMS * SEGS_PER_CHROM
d <- tibble(start = runif(segs),
            chrom = sample(1:N_CHROMS, segs, replace = TRUE)) |>
    bind_rows(tibble(start = 0, chrom = 1:N_CHROMS)) |>
    arrange(chrom, start) |>
    group_by(chrom) |>
    mutate(end = lead(start, default = 1)) |>
    ungroup() |>
    mutate(strain = sample(1:8, n(), replace = TRUE) |> as.factor())

ggplot(d, aes(xmin = chrom, xmax = chrom + 0.9, ymin = start, ymax = end,
              fill = strain)) +
    geom_rect(show.legend = FALSE) +
    scale_fill_brewer(type = "qual", palette = 6) +
    theme_minimal() +
    theme(axis.text = element_blank(),
          panel.grid = element_blank())

ggsave("figures/figure1/figure1b_mosaics.png", width = 4, height = 3)

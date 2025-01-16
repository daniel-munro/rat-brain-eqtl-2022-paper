# Kinship matrix for 75 NAcc samples and for all chromosomes (i.e. not LOCO) generated with:
# gemma -gk -g geno.txt.gz -p pheno.txt -outdir kinship -o all
# This is a centered relatedness matrix

library(tidyverse)

d <- read_tsv("data/gemma/kinship/all.cXX.txt",
              col_types = cols(.default = "d"),
              col_names = FALSE) |>
    mutate(rat1 = 1:n()) |>
    pivot_longer(-rat1,
                 names_to = "rat2",
                 values_to = "kinship") |>
    mutate(rat2 = rat2 |> str_sub(2) |> as.integer())

# Centered relatedness values:

d |>
    filter(rat1 <= rat2) |>
    ggplot(aes(x = kinship)) +
    ggtitle("Centered kinship values, including self-kinships") +
    geom_histogram(bins = 30) +
    annotate("rect", xmin = 0.105, xmax = 0.18, ymin = -20, ymax = 100,
             color = "red", fill = NA) +
    annotate("text", label = "siblings?", x = 0.1425, y = 150, color = "red") +
    annotate("rect", xmin = 0.27, xmax = 0.38, ymin = -20, ymax = 100,
             color = "red", fill = NA) +
    annotate("text", label = "self", x = 0.325, y = 150, color = "red")

d |>
    filter(rat1 <= rat2) |>
    mutate(Relationship = if_else(rat1 == rat2, "Self", "Kinship")) |>
    ggplot(aes(x = kinship, fill = Relationship)) +
    geom_histogram(bins = 100) +
    theme_bw() +
    theme(
        panel.grid = element_blank(),
        legend.position = c(0.7, 0.7),
    ) +
    xlab("Centered relatedness") +
    ylab("Rat pairs within cohort")

ggsave("figures/figureS3/figureS3a.png", width = 3, height = 3)

# Scaled relative to self = 1:
d |>
    group_by(rat1) |>
    mutate(kinship = kinship / (kinship[rat1 == rat2])) |>
    ungroup() |>
    filter(rat1 <= rat2) |>
    ggplot(aes(x = kinship)) +
    ggtitle("Kinship relative to self") +
    geom_histogram(bins = 30)


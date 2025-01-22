suppressPackageStartupMessages(library(tidyverse))
library(UpSetR)

eqtls <- read_tsv("data/eqtls/eqtls_indep.txt", col_types = "ccciiciiccdddddid")
x <- eqtls |>
    distinct(tissue, gene_id) |>
    with(split(gene_id, tissue))

colors <- viridis::viridis(5)
t <- list("IL", "LHb", "NAcc", "OFC", "PL")
queries <- list(
    list(query = intersects, params = t[1], color = colors[1], active = TRUE),
    list(query = intersects, params = t[2], color = colors[1], active = TRUE),
    list(query = intersects, params = t[3], color = colors[1], active = TRUE),
    list(query = intersects, params = t[4], color = colors[1], active = TRUE),
    list(query = intersects, params = t[5], color = colors[1], active = TRUE),
    
    list(query = intersects, params = t[c(1, 2)], color = colors[2], active = TRUE),
    list(query = intersects, params = t[c(1, 3)], color = colors[2], active = TRUE),
    list(query = intersects, params = t[c(1, 4)], color = colors[2], active = TRUE),
    list(query = intersects, params = t[c(1, 5)], color = colors[2], active = TRUE),
    list(query = intersects, params = t[c(2, 3)], color = colors[2], active = TRUE),
    list(query = intersects, params = t[c(2, 4)], color = colors[2], active = TRUE),
    list(query = intersects, params = t[c(2, 5)], color = colors[2], active = TRUE),
    list(query = intersects, params = t[c(3, 4)], color = colors[2], active = TRUE),
    list(query = intersects, params = t[c(3, 5)], color = colors[2], active = TRUE),
    list(query = intersects, params = t[c(4, 5)], color = colors[2], active = TRUE),
    
    list(query = intersects, params = t[c(1, 2, 3)], color = colors[3], active = TRUE),
    list(query = intersects, params = t[c(1, 2, 4)], color = colors[3], active = TRUE),
    list(query = intersects, params = t[c(1, 2, 5)], color = colors[3], active = TRUE),
    list(query = intersects, params = t[c(1, 3, 4)], color = colors[3], active = TRUE),
    list(query = intersects, params = t[c(1, 3, 5)], color = colors[3], active = TRUE),
    list(query = intersects, params = t[c(1, 4, 5)], color = colors[3], active = TRUE),
    list(query = intersects, params = t[c(2, 3, 4)], color = colors[3], active = TRUE),
    list(query = intersects, params = t[c(2, 3, 5)], color = colors[3], active = TRUE),
    list(query = intersects, params = t[c(2, 4, 5)], color = colors[3], active = TRUE),
    list(query = intersects, params = t[c(3, 4, 5)], color = colors[3], active = TRUE),
    
    list(query = intersects, params = t[c(1, 2, 3, 4)], color = colors[4], active = TRUE),
    list(query = intersects, params = t[c(1, 2, 3, 5)], color = colors[4], active = TRUE),
    list(query = intersects, params = t[c(1, 2, 4, 5)], color = colors[4], active = TRUE),
    list(query = intersects, params = t[c(1, 3, 4, 5)], color = colors[4], active = TRUE),
    list(query = intersects, params = t[c(2, 3, 4, 5)], color = colors[4], active = TRUE),
    
    list(query = intersects, params = t[c(1, 2, 3, 4, 5)], color = colors[5], active = TRUE)
)

png("figures/figure4/figure4c.png", width = 8, height = 5.5, units = "in", res = 300)
upset(
    fromList(x),
    queries = queries,
    show.numbers = FALSE,
    sets.x.label = "eGenes",
    mainbar.y.label = "Intersection Size (eGenes)",
    mb.ratio = c(0.6, 0.4),
    text.scale = 2,
    point.size = 4,
    line.size = 1,
)
dev.off()

###################
## Overlap stats ##
###################

d <- eqtls |>
    distinct(tissue, gene_id) |>
    mutate(eqtl = TRUE) |>
    complete(gene_id, tissue, fill = list(eqtl = FALSE)) |>
    pivot_wider(id_cols = gene_id, names_from = tissue, values_from = eqtl)

with(d, sum(IL & PL & OFC & !LHb & !NAcc)) # Highest 3-tissue eGene set
with(d, sum(NAcc & PL & OFC & !LHb & !IL)) # Next-highest 3-tissue eGene set

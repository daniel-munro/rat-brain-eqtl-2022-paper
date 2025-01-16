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

# png("stats/eGene_upset.png", width = 5, height = 3, units = "in", res = 300)
# pdf("stats/eGene_upset.pdf", width = 5, height = 3, onefile = FALSE)
## To fix obscured intersection bars, save as SVG, raise upper plot in Affinity, and export.
# svglite::svglite("stats/eGene_upset.svg", width = 5, height = 2.5)
png("figures/figure4/figure4c.png", width = 8, height = 5.5, units = "in", res = 300)
upset(
    fromList(x),
    # main.bar.color = c(rep(colors[1], 5),
    #                    rep(colors[2], 10),
    #                    rep(colors[3], 10),
    #                    rep(colors[4], 5),
    #                    rep(colors[5], 1)),
    # matrix.color = c("blue", "red"))
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

# length(setdiff(intersect(intx[["IL"]], x[["PL"]]), x[["OFC"]])), unique(c(x[["LHb"]], x[["NAcc"]]))))
d <- eqtls |>
    distinct(tissue, gene_id) |>
    mutate(eqtl = TRUE) |>
    complete(gene_id, tissue, fill = list(eqtl = FALSE)) |>
    pivot_wider(gene_id, names_from = tissue, values_from = eqtl)

with(d, sum(IL & PL & OFC & !LHb & !NAcc)) # Highest 3-tissue eGene set
with(d, sum(NAcc & PL & OFC & !LHb & !IL)) # Next-highest 3-tissue eGene set

library(tidyverse)
# library(umap)
library(patchwork)

# samples <- read_tsv("data/samples.txt", col_types = "cc")
# 
# expr <- read.csv(
#     "~/Dropbox (Scripps Research)/HS-RNASeq/quantitation/EnsemblGene_v2/ensembl-gene_raw-tpm.txt",
#     sep = "\t",
#     check.names = FALSE
# )
# expr <- expr[, colnames(expr) %in% samples$library]

# Actually, I think we should use upper-quartile-scaled values. From my understanding all processing up to these files was done as a combined dataset.

expr <- tibble(tissue = c("IL", "LHb", "NAcc", "OFC", "PL")) |>
    group_by(tissue) |>
    summarise(
        read_tsv(str_glue("data/expression/ensembl-gene_log2_{tissue}.bed.gz"),
                 col_types = cols(`#chr` = "-", start = "-", end = "-",
                                  gene_id = "c", .default = "d")) |>
            pivot_longer(-gene_id, names_to = "rat_id", values_to = "expr"),
        .groups = "drop"
    ) |>
    group_by(gene_id) |>
    filter(n_distinct(tissue) == 5) |>
    ungroup() |>
    pivot_wider(c(rat_id, tissue), names_from = gene_id, values_from = expr)
tmp <- str_c(expr$rat_id, expr$tissue, sep="_")
expr <- expr |>
    select(-rat_id, -tissue) |>
    as.data.frame()
rownames(expr) <- tmp
# expr <- expr[rowSums(expr) > 0, ]

# um <- umap(expr)
# d <- tibble(sample = rownames(um$layout),
#             UMAP1 = um$layout[, 1],
#             UMAP2 = um$layout[, 2]) |>
#     separate(sample, c("rat_id", "tissue"), sep = "_")

pca <- prcomp(expr, scale. = TRUE)
d <- pca$x[, 1:4] |>
    as_tibble(rownames = "sample") |>
    separate(sample, c("rat_id", "tissue"), sep = "_")

eig <- pca$sdev ^ 2
pve <- (100 * eig / sum(eig)) |> round()

##########################
## Single panel (PC1+2) ##
##########################

tissue_labels <- tribble(
    ~tissue, ~PC1, ~PC2,
    "IL",    5,    20,
    "LHb",   75,   5,
    "NAcc",  15,   -45,
    "OFC",   20,   50,
    "PL",    10,   35,
)

d |>
    ggplot(aes(x = PC1, y = PC2, color = tissue)) +
    geom_point(size = 0.5, show.legend = FALSE) +
    scale_color_manual(values = c("#377eb8", "#4daf4a", "#e41a1c", "#ff7f00", "#984ea3")) +
    # scale_color_manual(values = c("#567beb", "#57bf45", "#ffa600", "#e2656b", "#ffef3d")) +
    coord_fixed() +
    geom_text(aes(label = tissue), data = tissue_labels, fontface = "bold", show.legend = FALSE) +
    # labs(color = "Tissue") +
    xlab(str_glue("PC1 ({pve[1]}%)")) +
    ylab(str_glue("PC2 ({pve[2]}%)")) +
    theme_minimal()

# ggsave("expression/expr_cluster_fig.png", width = 4, height = 3.5, dpi = 300, bg = "white")

#################################
## Double panel (PC1+2, PC3+4) ##
#################################

# tissue_labels <- tribble(
#     ~tissue, ~PC1, ~PC2, ~PC3, ~PC4,
#     "IL",    5,    20,   90,   45,
#     "LHb",   75,   5,    130,  0,
#     "NAcc",  15,   -45,  -120, -20,
#     "OFC",   20,   50,   100,  -50,
#     "PL",    10,   35,   120,  30,
# )
# tissue_labels <- tribble(
#     ~tissue, ~PC1, ~PC2, ~PC3, ~PC4,
#     "IL",    5,    20,   90,   45,
#     "LHb",   65,   5,    130,  0,
#     "NAcc",  15,   -45,  -120, -15,
#     "OFC",   20,   50,   100,  -50,
#     "PL",    10,   35,   120,  30,
# )
## For smaller squares:
tissue_labels <- tribble(
    ~tissue, ~PC1, ~PC2, ~PC3, ~PC4,
    "IL",    5,    20,   90,   45,
    "LHb",   65,   5,    130,  0,
    "NAcc",  15,   -45,  -120, -15,
    "OFC",   20,   50,   100,  -50,
    "PL",    10,   35,   120,  30,
)

p1 <- d |>
    ggplot(aes(x = PC1, y = PC2, color = tissue)) +
    geom_point(size = 0.5, show.legend = FALSE) +
    scale_color_manual(values = c("#377eb8", "#4daf4a", "#e41a1c", "#ff7f00", "#984ea3")) +
    # coord_fixed() +
    geom_text(aes(label = tissue), data = tissue_labels, fontface = "bold", show.legend = FALSE) +
    xlab(str_glue("PC1 ({pve[1]}%)")) +
    ylab(str_glue("PC2 ({pve[2]}%)")) +
    theme_bw() +
    theme(panel.grid = element_blank())

p2 <- d |>
    ggplot(aes(x = PC3, y = PC4, color = tissue)) +
    geom_point(size = 0.5, show.legend = FALSE) +
    scale_color_manual(values = c("#377eb8", "#4daf4a", "#e41a1c", "#ff7f00", "#984ea3")) +
    # coord_fixed() +
    expand_limits(x = c(min(tissue_labels$PC3) - 20, max(tissue_labels$PC3) + 15)) +
    geom_text(aes(label = tissue), data = tissue_labels, fontface = "bold", show.legend = FALSE) +
    xlab(str_glue("PC3 ({pve[3]}%)")) +
    ylab(str_glue("PC4 ({pve[4]}%)")) +
    theme_bw() +
    theme(panel.grid = element_blank())

p1 / p2

## For stacked coord_fixed:
# ggsave("expression/expr_cluster_fig.png", width = 4, height = 5.5, dpi = 300, bg = "white")
# ## For side-by-side squares:
# ggsave("expression/expr_cluster_fig.png", width = 6, height = 3, dpi = 300, bg = "white")
## For stacked squares:
ggsave("figures/figure1/figure1ef.png", width = 2.8, height = 5, dpi = 300, bg = "white")

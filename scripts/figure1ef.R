library(tidyverse)
library(patchwork)

# Using upper-quartile-scaled values. From my understanding all processing up to these files was done as a combined dataset.

expr <- tibble(tissue = c("IL", "LHb", "NAcc", "OFC", "PL")) |>
    reframe(
        read_tsv(str_glue("data/expression/ensembl-gene_log2_{tissue}.bed.gz"),
                 col_types = cols(`#chr` = "-", start = "-", end = "-",
                                  gene_id = "c", .default = "d")) |>
            pivot_longer(-gene_id, names_to = "rat_id", values_to = "expr"),
        .by = tissue
    ) |>
    filter(n_distinct(tissue) == 5, .by = gene_id) |>
    pivot_wider(id_cols = c(rat_id, tissue), names_from = gene_id, values_from = expr)
tmp <- str_glue("{expr$rat_id}_{expr$tissue}")
expr <- expr |>
    select(-rat_id, -tissue) |>
    as.data.frame()
rownames(expr) <- tmp

pca <- prcomp(expr, scale. = TRUE)
d <- pca$x[, 1:4] |>
    as_tibble(rownames = "sample") |>
    separate_wider_delim(sample, "_", names = c("rat_id", "tissue"))

eig <- pca$sdev ^ 2
pve <- (100 * eig / sum(eig)) |> round()

#################################
## Double panel (PC1+2, PC3+4) ##
#################################

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
    geom_text(aes(label = tissue), data = tissue_labels, fontface = "bold", show.legend = FALSE) +
    xlab(str_glue("PC1 ({pve[1]}%)")) +
    ylab(str_glue("PC2 ({pve[2]}%)")) +
    theme_bw() +
    theme(panel.grid = element_blank())

p2 <- d |>
    ggplot(aes(x = PC3, y = PC4, color = tissue)) +
    geom_point(size = 0.5, show.legend = FALSE) +
    scale_color_manual(values = c("#377eb8", "#4daf4a", "#e41a1c", "#ff7f00", "#984ea3")) +
    expand_limits(x = c(min(tissue_labels$PC3) - 20, max(tissue_labels$PC3) + 15)) +
    geom_text(aes(label = tissue), data = tissue_labels, fontface = "bold", show.legend = FALSE) +
    xlab(str_glue("PC3 ({pve[3]}%)")) +
    ylab(str_glue("PC4 ({pve[4]}%)")) +
    theme_bw() +
    theme(panel.grid = element_blank())

p1 / p2
ggsave("figures/figure1/figure1ef.png", width = 2.8, height = 5, dpi = 300, bg = "white")

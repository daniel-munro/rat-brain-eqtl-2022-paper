library(tidyverse)

# Look at gene strand to orient TSS distance.
genes <- read_tsv("data/genes.txt", col_types = "c----c------")

eqtls <- read_tsv("data/eqtls/eqtls_indep.txt", col_types = "ccciiciiccdddddid") |>
    left_join(genes, by = "gene_id") |>
    mutate(tss_distance = if_else(strand == "+", pos - tss, tss - pos) / 1e6)

#####################
## GTEx comparison ##
#####################
# Compare top eQTL per gene, not all independent eQTLs, since non-primary eQTLs
# are probably further away and there are many more of them in GTEx.

gtex_top <- tibble(
    file = list.files("data/gtex/GTEx_Analysis_v8_eQTL",
    full.names = TRUE
)) |>
    group_by(file) |>
    summarise(
        read_tsv(
            file,
            col_types = cols(
                gene_id = "c", qval = "d", pval_beta = "-",
                log2_aFC = "d", tss_distance = "i", .default = "-"
            )
        ),
        .groups = "drop"
    ) |>
    filter(qval <= 0.05) |> # GTEx site says to do <=
    mutate(
        tissue = str_match(file, "eQTL/(.+)\\.v8")[, 2],
        group = if_else(str_detect(tissue, "Brain"), "GTEx brain", "GTEx other"),
        abs_tss_distance = abs(tss_distance)
    )

egenes <- eqtls |>
    filter(rank == 1)

tss <- bind_rows(
    gtex_top |>
        filter(group == "GTEx brain") |>
        select(tissue, group, abs_tss_distance),
    egenes |>
        mutate(group = "Rat brain",
               abs_tss_distance = abs(pos - tss)) |>
        select(tissue, group, abs_tss_distance)
)

tss |>
    mutate(group = fct_rev(group),
           abs_tss_distance = abs_tss_distance / 1e6) |>
    ggplot(aes(x = abs_tss_distance, color = group)) +
    geom_density(adjust = 0.5, size = 0.8, outline.type = "full", show.legend = FALSE) +
    scale_x_continuous(expand = c(0.02, 0)) +
    scale_y_continuous(expand = c(0.02, 0)) +
    scale_color_manual(values = c("#F8766DFF", "#619CFFFF")) +
    theme_bw() +
    theme(
        panel.grid = element_blank(),
    ) +
    xlab("TSS distance (Mb)") +
    ylab("Density") # Line up with LD decay panel

ggsave("figures/figure5/figure5b.png", width = 2.5, height = 2.5, dpi = 300)

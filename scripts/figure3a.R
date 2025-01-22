library(tidyverse)
library(patchwork)

genome_bin <- function(chrom, pos, width = 2e7) {
    chr_len <- c(282763074, 266435125, 177699992, 184226339, 173707219,
                 147991367, 145729302, 133307652, 122095297, 112626471,
                 90463843, 52716770, 114033958, 115493446, 111246239,
                 90668790, 90843779, 88201929, 62275575, 56205956)
    n_bins <- ceiling(chr_len / width)
    bins <- tibble(chrom = 1:20) |>
        reframe(tibble(chr_bin = 1:n_bins[chrom]), .by = chrom) |>
        mutate(bin = str_c(chrom, "_", chr_bin)) |>
        pull(bin)
    tibble(chrom, pos) |>
        mutate(bin = str_c(chrom, "_", ceiling(pos / width)) |>
                   factor(levels = bins)) |>
        pull(bin)
}

chrom_boundaries <- function(bins) {
    chrs <- str_split(levels(bins), "_", simplify = TRUE)[, 1]
    which(!duplicated(chrs))[-1] - 0.5
}

chrom_label_locs <- function(bins) {
    chrs <- str_split(levels(bins), "_", simplify = TRUE)[, 1]
    first <- which(!duplicated(chrs))
    last <- c(first[-1] - 1, length(chrs))
    levels(bins)[round((first + last) / 2)]
}

genes <- read_tsv("data/genes.txt", col_types = "c-c---i-----") |>
    rename(gene_chrom = chrom,
           gene_pos = tss)

eqtls <- read_tsv("data/tensorqtl/NQCT.trans_qtl_pairs.txt.gz", col_types = "ccd---") |>
    rename(gene_id = phenotype_id) |>
    separate_wider_delim(variant_id, ":", names = c("chrom", "pos"),
             cols_remove = FALSE) |>
    mutate(chrom = str_replace(chrom, "chr", "") |> as.integer(),
           pos = as.integer(pos)) |>
    left_join(genes, by = "gene_id")

binned <- eqtls |>
    mutate(gene_bin = genome_bin(gene_chrom, gene_pos),
           var_bin = genome_bin(chrom, pos)) |>
    summarise(min_pval = min(pval),
              .by = c(gene_chrom, gene_bin, chrom, var_bin))

binned |>
    mutate(log10_min_pval = pmin(-log10(min_pval), 10)) |>
    ggplot(aes(x = gene_bin, y = var_bin, fill = log10_min_pval)) +
    geom_tile() +
    geom_vline(xintercept = chrom_boundaries(binned$gene_bin),
               color = "black", linewidth = 0.1) +
    geom_hline(yintercept = chrom_boundaries(binned$var_bin),
               color = "black", linewidth = 0.1) +
    scale_x_discrete(drop = FALSE,
                     breaks = chrom_label_locs(binned$gene_bin),
                     labels = c("chr1", 2:10, "\n11", 12, "\n13", 14, "\n15", 16, "\n17", 18, "\n19", 20)) +
    scale_y_discrete(drop = FALSE,
                     breaks = chrom_label_locs(binned$gene_bin),
                     labels = c("chr1", 2:20)) +
    scale_fill_gradientn(colors = c("#cccccc", "#5555ff", "red"),
                         breaks = 5:10,
                         labels = c(5:9, "10+")) +
    coord_fixed() +
    theme_minimal() +
    theme(
        panel.grid = element_blank(),
        legend.position = "top",
    ) +
    xlab("Gene location") +
    ylab("Variant location") +
    labs(fill = expression(-log[10](P)))

ggsave("figures/figure3/figure3a.png", width = 5.2, height = 6, bg = "white")

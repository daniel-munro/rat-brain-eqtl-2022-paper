library(tidyverse)
library(patchwork)

chr_len <- c(282763074, 266435125, 177699992, 184226339, 173707219,
             147991367, 145729302, 133307652, 122095297, 112626471,
             90463843, 52716770, 114033958, 115493446, 111246239,
             90668790, 90843779, 88201929, 62275575, 56205956)
label_locs <- cumsum(c(0, chr_len[1:19])) + chr_len / 2
grid_locs <- cumsum(c(0, chr_len))

eqtls <- read_tsv("data/eqtls/eqtls_indep.txt", col_types = "cc---c-----------")
esnps <- read_tsv("data/tensorqtl/PQCT.cis_qtl_signif.txt.gz", col_types = "cc----ddd-") |>
    rename(gene_id = phenotype_id) |>
    filter(gene_id %in% eqtls$gene_id) |>
    mutate(variant_id = str_replace(variant_id, "chr", ""),
           logp = -log10(pval_nominal)) |>
    separate(variant_id, c("chrom", "pos"), sep = ":", convert = TRUE) |>
    mutate(gpos = pos + cumsum(c(0, chr_len))[chrom],
           gcolor = as.factor((chrom - 1) %% 2)) |>
    arrange(gpos)

gwas <- read_tsv("data/coloc/adiposity_GWAS/allChr_physiological_retrofat.assoc.txt",
                 col_types = "-c------------d",
                 col_names = c("variant_id", "p_score")) |>
    mutate(variant_id = str_replace(variant_id, "chr", ""),
           logp = -log10(p_score)) |>
    separate(variant_id, c("chrom", "pos"), sep = ":", convert = TRUE) |>
    mutate(gpos = pos + cumsum(c(0, chr_len))[chrom],
           gcolor = as.factor((chrom - 1) %% 2)) |>
    arrange(gpos)

smr <- read_tsv("data/coloc/SMR_sig.tsv", col_types = "cccdcddd") |>
    filter(tissue == "PL",
           trait == "retrofat") |>
    mutate(variant_id = str_replace(variant_id, "chr", ""),
           logp = -log10(p_SMR)) |>
    separate(variant_id, c("chrom", "pos"), sep = ":", convert = TRUE) |>
    mutate(gpos = pos + cumsum(c(0, chr_len))[chrom]) |>
    arrange(gpos)

ggplot(gwas, aes(x = gpos, y = logp, color = gcolor)) +
    geom_point(size = 0.5, show.legend = FALSE) +
    geom_point(data = smr, color = "#984ea3", shape = 5, stroke = 1) +
    expand_limits(x = c(0, sum(chr_len))) +
    scale_x_continuous(breaks = label_locs, labels = 1:20, expand = c(0.02, 0),
                       minor_breaks = grid_locs) +
    scale_color_manual(values = c("#111111", "#444444")) +
    theme_bw() +
    theme(
        panel.grid.major = element_blank(),
        panel.grid.minor.y = element_blank(),
    ) +
    xlab("Chromosome") +
    ylab(expression(-log[10](P[GWAS]*" or "*P[SMR])))

ggsave("figures/figure6/figure6a1.png", width = 7, height = 2)

p <- ggplot(esnps, aes(x = gpos, y = logp, color = gcolor)) +
    geom_point(size = 0.5, show.legend = FALSE) +
    expand_limits(x = c(0, sum(chr_len))) +
    expand_limits(y = c(0, max(esnps$logp) * 1.05)) +
    scale_x_continuous(breaks = label_locs, labels = 1:20, expand = c(0.02, 0),
                       minor_breaks = grid_locs) +
    scale_y_continuous(expand = c(0.01, 0)) +
    scale_color_manual(values = c("#984ea3", "#c698cd")) +
    theme_bw() +
    theme(
        panel.grid.major = element_blank(),
        panel.grid.minor.y = element_blank(),
        legend.box.margin = margin(0, 0, 0, -5, unit = "pt"),
    ) +
    xlab("Chromosome") +
    ylab(expression(-log[10](P[eQTL])))

ggsave("figures/figure6/figure6a2.png", p, width = 7, height = 2)

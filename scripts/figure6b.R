library(tidyverse)

chr_len <- c(282763074, 266435125, 177699992, 184226339, 173707219,
             147991367, 145729302, 133307652, 122095297, 112626471,
             90463843, 52716770, 114033958, 115493446, 111246239,
             90668790, 90843779, 88201929, 62275575, 56205956)
label_locs <- cumsum(c(0, chr_len[1:19])) + chr_len / 2
grid_locs <- cumsum(c(0, chr_len))

smr <- read_tsv("data/coloc/SMR.tsv", col_types = "cccdcddd") |>
    mutate(sig = p_SMR < 0.05 / n(),
           .by = c(tissue, trait)) |>
    mutate(variant_id = str_replace(variant_id, "chr", ""),
           logp = -log10(p_SMR)) |>
    separate_wider_delim(variant_id, ":", names = c("chrom", "pos")) |>
    mutate(chrom = as.integer(chrom),
           pos = as.integer(pos),
           gpos = pos + cumsum(c(0, chr_len))[chrom],
           gcolor = as.factor((chrom - 1) %% 2)) |>
    arrange(gpos)

thresh <- smr |>
    summarise(threshold = -log10(0.05 / n()),
              .by = c(tissue, trait)) |>
    summarise(threshold = median(threshold),
              .by = tissue) # Within tissue they're almost identical

smr |>
    filter(sig) |>
    ggplot(aes(x = gpos, y = logp, color = tissue, shape = tissue)) +
    # Workaround since only one color scale can be used:
    geom_point(aes(color = NULL), data = filter(smr, !sig, gcolor == 0), size = 0.5,
               color = "#999999", show.legend = FALSE) +
    geom_point(aes(color = NULL), data = filter(smr, !sig, gcolor == 1), size = 0.5,
               color = "#bbbbbb", show.legend = FALSE) +
    geom_point(size = 1) +
    geom_hline(aes(yintercept = threshold, color = tissue), data = thresh, size = 0.2, alpha = 0.5) +
    expand_limits(x = c(0, sum(chr_len))) +
    expand_limits(y = c(0, max(smr$logp) * 1.05)) +
    scale_x_continuous(breaks = label_locs, labels = 1:20, expand = c(0.02, 0),
                       minor_breaks = grid_locs) +
    scale_y_continuous(expand = c(0.01, 0)) +
    scale_color_manual(values = c("#377eb8", "#4daf4a", "#e41a1c", "#ff7f00", "#984ea3")) +
    scale_shape_manual(values = c(1:5)) +
    theme_bw() +
    theme(
        legend.position = "inside",
        legend.position.inside = c(0.98, 0.95),
        legend.justification = c(1, 1),
        legend.background = element_rect(color = "black", size = 0.2),
        legend.spacing.y = unit(5, "pt"),
        legend.key.size = unit(10, "pt"),
        panel.grid.major = element_blank(),
        panel.grid.minor.y = element_blank(),
    ) +
    xlab("Chromosome") +
    ylab(expression(-log[10](P[SMR]))) +
    labs(color = "Tissue", shape = "Tissue")

ggsave("figures/figure6/figure6b.png", width = 7, height = 2.7)

## Trait labels:

smr |>
    filter(sig) |>
    distinct(trait, chrom, pos) |>
    arrange(chrom, pos, trait) |>
    View()

## Stats:

genes <- read_tsv("data/genes.txt", col_types = "cc----------")

sig <- read_tsv("data/coloc/SMR_sig.tsv", col_types = "cccdcddd") |>
    left_join(genes, by = "gene_id", relationship = "many-to-one")

sig |>
    count(tissue, trait, sort = TRUE)

sig |>
    count(gene_name) |>
    filter(n >= 4) |>
    arrange(gene_name)

smr |>
    summarise(threshold = 0.05 / n(), .by = c(tissue, trait)) |>
    summarise(min = min(threshold),
              max = max(threshold))

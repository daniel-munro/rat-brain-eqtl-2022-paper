library(VariantAnnotation)
library(tidyverse)
library(patchwork)

# Look at gene strand to orient TSS distance.
genes <- read_tsv("data/genes.txt", col_types = "c----c------")

eqtls <- read_tsv("data/eqtls/eqtls_indep.txt", col_types = "ccciiciiccdddddid") |>
    filter(tissue == "NAcc") |>
    left_join(genes, by = "gene_id") |>
    mutate(tss_distance = if_else(strand == "+", pos - tss, tss - pos) / 1e6,
           maf = pmin(af, 1 - af))

#########
## aFC ##
#########

p1 <- eqtls |>
    mutate(log2_aFC = pmax(-4, pmin(log2_aFC, 4))) |>
    ggplot(aes(x = log2_aFC)) +
    geom_histogram(aes(y = ..density..), bins = 100) +
    scale_x_continuous(expand = c(0, 0), limits = c(-4.001, 4.001)) +
    scale_y_continuous(expand = c(0, 0), limits = c(0, 1.4)) +
    xlab(expression("aFC "*(log[2]))) +
    ylab("Density") +
    theme_bw() +
    theme(
        panel.grid = element_blank(),
    )

read_tsv("data/eqtls/eqtls_indep.txt", col_types = "ccciiciiccdddddid") |>
    summarise(fc_2 = mean(abs(log2_aFC) <= 1),
              fc_4 = mean(abs(log2_aFC) <= 2))

######################
## allele frequency ##
######################

all_snps <- info(readVcf("data/genotype/P50.rnaseq.88.unpruned.vcf.gz")) |>
    as_tibble() |>
    mutate(AC = as.integer(AC),
           AF = AC / AN,
           maf = pmin(AF, 1 - AF))
nbins <- 18
p2 <- eqtls |>
    ggplot(aes(x = maf)) +
    geom_histogram(aes(y = after_stat(density)), binwidth = 0.5 / nbins, boundary = 0.5) +
    stat_bin(aes(y = after_stat(density)),
             data = all_snps,
             geom = "step",
             binwidth = 0.5 / nbins,
             boundary = 0.5,
             color = "red",
             position = position_nudge(x = -0.5 * (0.5 / nbins))) +
    scale_x_continuous(expand = c(0, 0), limits = c(-0.001, 0.501)) +
    scale_y_continuous(expand = c(0, 0), limits = c(0, 3.7)) +
    xlab("Minor allele frequency") +
    ylab("Density") +
    theme_bw() +
    theme(
        panel.grid = element_blank(),
    )

##################
## TSS distance ##
##################

p3 <- eqtls |>
    ggplot(aes(x = tss_distance)) +
    geom_histogram(aes(y = after_stat(density)), bins = 50) +
    scale_x_continuous(expand = c(0, 0), limits = c(-1.001, 1.001)) +
    scale_y_continuous(expand = c(0, 0), limits = c(0, 1.7)) +
    xlab("Distance from TSS (Mb)") +
    ylab("Density") +
    theme_bw() +
    theme(
        panel.grid = element_blank(),
    )

p1 / plot_spacer() / p2 / plot_spacer() / p3 + plot_layout(heights = c(1, 0.1, 1, 0.1, 1))
ggsave("figures/figure3/figure3cde.png", width = 2.75, height = 7)

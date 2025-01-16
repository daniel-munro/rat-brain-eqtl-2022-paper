# Based on ../haplotype_skew.R

library(tidyverse)

probs <- function(prob) {
    names(dimnames(prob)) <- c("individual", "Strain", "SNP")
    cubelyr::as.tbl_cube(prob) |>
        as_tibble() |>
        mutate(Strain = strains[Strain])
}

strains <- c(
    A = "ACI/N", B = "BN/SsN", C = "BUF/N", D = "F344/N",
    E = "M520/N", F = "MR/N", G = "WN/N", H = "WKY/N"
)
# chr_len <- c(282763074, 266435125, 177699992)
# n_SNPs <- round(2000 * chr_len / chr_len[1])

d <- tibble(chrom = 1:3) |>
    group_by(chrom) |>
    summarise(
        readRDS(str_glue("data/haplotypes/haplotype_probs_chr{chrom}.rds")) |>
            probs(),
        .groups = "drop"
    ) |>
    group_by(SNP, Strain) |>
    summarise(prob = mean(prob), .groups = "drop") |>
    separate(SNP, c("chrom", "pos"), sep = ":", convert = TRUE) |>
    mutate(pos = pos / 1e6)

d |>
    mutate(Strain = str_glue("{Strain} ({round(deframe(tot)[Strain] * 100, 1)}%)")) |>
    ggplot(aes(x = pos, y = prob, fill = Strain)) +
    facet_wrap(~ chrom, ncol = 1, strip.position = "left") +
    geom_area() +
    scale_fill_brewer(type = "qual", palette = 6) +
    theme_bw() +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    theme(
        axis.text.y = element_blank(),
        # panel.grid.major.y = element_blank(),
        # panel.grid.minor.y = element_blank(),
        panel.grid = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.x = element_line(size = 0.25),
        panel.border = element_blank(),
        strip.background = element_blank(),
    ) +
    xlab("Position (Mb)") +
    ylab("Estimated Proportion")

ggsave("figures/figure2/figure2a.png", width = 6, height = 2.5)

# ###################################
# ## Overall haplotype proportions ##
# ###################################
# 
# library(VariantAnnotation)
# 
# # I saved 2000 loci per chromosome, so get original counts to weight the probabilities.
# counts <- VariantAnnotation::readVcf("data/genotype/P50.rnaseq.88.unpruned.vcf.gz") |>
#     VariantAnnotation::info() |>
#     rownames() |>
#     str_replace("chr", "") |>
#     str_replace(":.+$", "") |>
#     table() |>
#     enframe(name = "chrom", value = "count") |>
#     mutate(chrom = as.integer(chrom)) |>
#     # filter(chrom %in% 1:3) |>
#     mutate(weight = count / sum(count))
# 
# tot <- tibble(chrom = 1:20) |>
#     group_by(chrom) |>
#     summarise(
#         readRDS(str_glue("data/haplotypes/haplotype_probs_chr{chrom}.rds")) |>
#             probs(),
#         .groups = "drop"
#     ) |>
#     left_join(counts, by = "chrom") |>
#     group_by(Strain) |>
#     summarise(prob = weighted.mean(prob, weight), .groups = "drop")
# 
# tot |>
#     ggplot(aes(x = 1, y = prob, fill = Strain)) +
#     geom_col(show.legend = FALSE) +
#     coord_polar("y") +
#     scale_fill_brewer(type = "qual", palette = 6) +
#     theme_void()
# 
# ggsave("figures/figure2/figure2a2.png")

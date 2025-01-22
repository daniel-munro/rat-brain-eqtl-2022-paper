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

d <- tibble(chrom = 1:3) |>
    reframe(
        readRDS(str_glue("data/haplotypes/haplotype_probs_chr{chrom}.rds")) |>
            probs(),
        .by = chrom
    ) |>
    summarise(prob = mean(prob), .by = c(SNP, Strain)) |>
    separate_wider_delim(SNP, ":", names = c("chrom", "pos")) |>
    mutate(pos = as.integer(pos) / 1e6)

# I saved 2000 loci per chromosome, so get original counts to weight the probabilities.
counts <- VariantAnnotation::readVcf("data/genotype/P50.rnaseq.88.unpruned.vcf.gz") |>
    VariantAnnotation::info() |>
    rownames() |>
    str_replace("chr", "") |>
    str_replace(":.+$", "") |>
    table() |>
    enframe(name = "chrom", value = "count") |>
    mutate(chrom = as.integer(chrom)) |>
    mutate(weight = count / sum(count))

tot <- tibble(chrom = 1:20) |>
    reframe(
        readRDS(str_glue("data/haplotypes/haplotype_probs_chr{chrom}.rds")) |>
            probs(),
        .by = chrom
    ) |>
    left_join(counts, by = "chrom") |>
    summarise(prob = weighted.mean(prob, weight), .by = Strain)

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
        panel.grid = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.x = element_line(size = 0.25),
        panel.border = element_blank(),
        strip.background = element_blank(),
    ) +
    xlab("Position (Mb)") +
    ylab("Estimated Proportion")

ggsave("figures/figure2/figure2a.png", width = 6, height = 2.5)

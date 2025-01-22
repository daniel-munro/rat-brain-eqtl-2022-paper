library(VariantAnnotation)
library(tidyverse)

r2 <- function(m1, m2) {
    # for each row i, computes correlation between m1[i, ] and m2[i, ]
    a1 <- m1 - rowMeans(m1)
    a2 <- m2 - rowMeans(m2)
    (rowSums(a1 * a2) / sqrt(rowSums(a1 ^ 2) * rowSums(a2 ^ 2))) ^ 2
}

MAX_DIST <- 1e7

filename <- "data/genotype/P50.rnaseq.88.unpruned.vcf.gz"
gt <- VariantAnnotation::readGT(filename)
geno <- apply(gt, 2, function(x) c("0|0" = 0, "0|1" = 1, "1|0" = 1, "1|1" = 2)[x])
rownames(geno) <- rownames(gt)
rm(gt)

mac <- tibble(SNP = rownames(geno),
              AC = rowSums(geno)) |>
    mutate(AC = pmin(AC, (ncol(geno) * 2) - AC)) |>
    filter(AC / (ncol(geno) * 2) > 0.2) |>
    separate_wider_delim(SNP, ":", names = c("chrom", "pos"), cols_remove = FALSE) |>
    mutate(pos = as.integer(pos))

ld_pairs <- mac |>
    reframe({
        n_pairs <- round(0.01 * length(SNP) ^ 2)
        # print(str_glue("{unique(chrom)}, AC {unique(AC)}, {n_pairs} pairs"))
        sample1 <- sample.int(length(SNP), n_pairs, replace = TRUE)
        sample2 <- sample.int(length(SNP), n_pairs, replace = TRUE)
        tibble(SNP.x = SNP[sample1],
               pos.x = pos[sample1],
               SNP.y = SNP[sample2],
               pos.y = pos[sample2]) |>
            distinct() |>
            filter(pos.x < pos.y,
                   pos.y - pos.x <= MAX_DIST)
    }, .by = c(chrom, AC)) |>
    slice_sample(n = 1e6)

ld <- ld_pairs |>
    mutate(distance = abs(pos.y - pos.x),
           r2 = r2(geno[SNP.x, , drop = FALSE],
                   geno[SNP.y, , drop = FALSE])) |>
    select(-AC, -SNP.x, -SNP.y)

write_tsv(ld, "data/genotype/LD.txt.gz")

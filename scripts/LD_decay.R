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
# rng <- GRanges(1, IRanges(1, 1e7))
# gt <- VariantAnnotation::readGT(filename, param = rng)
gt <- VariantAnnotation::readGT(filename)
geno <- apply(gt, 2, function(x) c("0|0" = 0, "0|1" = 1, "1|0" = 1, "1|1" = 2)[x])
rownames(geno) <- rownames(gt)
rm(gt)

mac <- tibble(SNP = rownames(geno),
              AC = rowSums(geno)) |>
    mutate(AC = pmin(AC, (ncol(geno) * 2) - AC)) |>
    filter(AC / (ncol(geno) * 2) > 0.2) |>
    separate(SNP, c("chrom", "pos"), sep = ":", convert = TRUE, remove = FALSE)

ld_pairs <- mac |>
    # slice(1:10) |>
    group_by(chrom, AC) |>
    summarise({
        n_pairs <- round(0.01 * length(SNP) ^ 2)
        print(str_glue("{unique(chrom)}, AC {unique(AC)}, {n_pairs} pairs"))
        # crossing(pos.x = pos, pos.y = pos) |>
        #     slice_sample(prop = 0.01) |> # for efficiency
        # tibble(pos.x = sample(pos, 1e5, replace = TRUE), # for more efficiency
        #        pos.y = sample(pos, 1e5, replace = TRUE)) |>
        #     distinct() |>
        # crossing(SNP.x = SNP, SNP.y = SNP) |>
        #     slice_sample(prop = 0.001) |>
        #     separate(SNP.x, c("chrom.x", "pos.x"), sep = ":", convert = TRUE, remove = FALSE) |>
        #     separate(SNP.y, c("chrom.y", "pos.y"), sep = ":", convert = TRUE, remove = FALSE) |>
        sample1 <- sample.int(length(SNP), n_pairs, replace = TRUE)
        sample2 <- sample.int(length(SNP), n_pairs, replace = TRUE)
        # print(length(SNP))
        tibble(SNP.x = SNP[sample1],
               pos.x = pos[sample1],
               SNP.y = SNP[sample2],
               pos.y = pos[sample2]) |>
            distinct() |>
            filter(pos.x < pos.y,
                   pos.y - pos.x <= MAX_DIST)
    }, .groups = "drop") |>
    slice_sample(n = 1e6)

ld <- ld_pairs |>
    # group_by(chrom, AC) |>
    # mutate(SNP.x = str_c(chrom, ":", pos.x),
    #        SNP.y = str_c(chrom, ":", pos.y),
    mutate(distance = abs(pos.y - pos.x),
           r2 = r2(geno[SNP.x, , drop = FALSE],
                   geno[SNP.y, , drop = FALSE])) |>
    # str_glue("{unique(chrom)}, AC {unique(AC)}"))) |>
    # ungroup() |>
    select(-AC, -SNP.x, -SNP.y)

write_tsv(ld, "data/genotype/LD.txt.gz")

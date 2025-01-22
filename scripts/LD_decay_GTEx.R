library(VariantAnnotation)
library(tidyverse)

load_geno_chr <- function(chr, filename) {
    rng <- GRanges(chr, IRanges(1, 5e8))
    gt <- VariantAnnotation::readGT(filename, param = rng)
    geno <- apply(gt, 2, function(x) c("0|0" = 0, "0|1" = 1, "1|0" = 1, "1|1" = 2)[x])
    rownames(geno) <- rownames(gt)
    rm(gt)
    geno
}

r2 <- function(m1, m2) {
    # for each row i, computes correlation between m1[i, ] and m2[i, ]
    a1 <- m1 - rowMeans(m1)
    a2 <- m2 - rowMeans(m2)
    (rowSums(a1 * a2) / sqrt(rowSums(a1 ^ 2) * rowSums(a2 ^ 2))) ^ 2
}

MAX_DIST <- 1e7
filename <- "data/gtex/GTEx_MAF20.vcf.gz"

geno <- str_c("chr", c(1:22, "X")) |>
    lapply(load_geno_chr, filename)
geno <- do.call(rbind, geno)

mac <- tibble(SNP = rownames(geno),
              AC = rowSums(geno)) |>
    mutate(AC = pmin(AC, (ncol(geno) * 2) - AC)) |>
    filter(AC / (ncol(geno) * 2) > 0.2) |>
    separate_wider_delim(SNP, "_", names = c("chrom", "pos", "ref", "alt", "b"),
                         cols_remove = FALSE) |>
    mutate(pos = as.integer(pos))

ld_pairs <- mac |>
    reframe({
        print(str_glue("{unique(chrom)}, AC {unique(AC)}"))
        crossing(SNP.x = SNP, SNP.y = SNP) |>
            slice_sample(prop = 0.02) |>
            separate_wider_delim(
                SNP.x, "_", names = c("chrom.x", "pos.x", "ref.x", "alt.x", "b.x"),
                cols_remove = FALSE
            ) |>
            separate_wider_delim(
                SNP.y, "_", names = c("chrom.y", "pos.y", "ref.y", "alt.y", "b.y"),
                cols_remove = FALSE
            ) |>
            mutate(pos.x = as.integer(pos.x),
                   pos.y = as.integer(pos.y)) |>
            filter(pos.x < pos.y,
                   pos.y - pos.x <= MAX_DIST)
    }, .by = c(chrom, AC)) |>
    slice_sample(n = 1e6)

ld <- ld_pairs |>
    mutate(distance = abs(pos.y - pos.x),
           r2 = r2(geno[SNP.x, , drop = FALSE],
                   geno[SNP.y, , drop = FALSE])) |>
    select(chrom, pos.x, pos.y, distance, r2)

write_tsv(ld, "data/gtex/GTEx_LD.txt.gz")

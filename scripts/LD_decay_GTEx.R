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
filename <- "GTEx_MAF20.vcf.gz"

# geno <- str_c("chr", 1:3) |>
geno <- str_c("chr", c(1:22, "X")) |>
    lapply(load_geno_chr, filename)
geno <- do.call(rbind, geno)

mac <- tibble(SNP = rownames(geno),
              AC = rowSums(geno)) |>
    mutate(AC = pmin(AC, (ncol(geno) * 2) - AC)) |>
    filter(AC / (ncol(geno) * 2) > 0.2) |>
    separate(SNP, c("chrom", "pos", "ref", "alt", "b"), sep = "_", convert = TRUE,
             remove = FALSE)

ld_pairs <- mac |>
    # slice(1:10) |>
    # filter(AC %% 5 == 0) |>
    group_by(chrom, AC) |>
    summarise({
        print(str_glue("{unique(chrom)}, AC {unique(AC)}"))
        # crossing(pos.x = pos, pos.y = pos) |>
        #     slice_sample(prop = 0.01) |> # for efficiency
        ## For more efficiency:
        # tibble(SNP.x = sample(SNP, 1e5, replace = TRUE),
        #        SNP.y = sample(SNP, 1e5, replace = TRUE)) |>
        #     left_join(select(mac, SNP.x = SNP, pos.x = pos), by = "SNP.x") |>
        #     left_join(select(mac, SNP.y = SNP, pos.y = pos), by = "SNP.y") |>
        ## For even more efficiency:
        # n_pairs <- 
        # sample1 <- sample.int(length(SNP), 1e4, replace = TRUE)
        # sample2 <- sample.int(length(SNP), 1e4, replace = TRUE)
        # print(length(SNP))
        # tibble(SNP.x = SNP[sample1],
        #        pos.x = pos[sample1],
        #        SNP.y = SNP[sample2],
        #        pos.y = pos[sample2]) |>
        #     distinct() |>
        crossing(SNP.x = SNP, SNP.y = SNP) |>
            slice_sample(prop = 0.02) |>
            # left_join(select(mac, SNP.x = SNP, pos.x = pos), by = "SNP.x") |>
            # left_join(select(mac, SNP.y = SNP, pos.y = pos), by = "SNP.y") |>
            separate(SNP.x, c("chrom.x", "pos.x", "ref.x", "alt.x", "b.x"), sep = "_", convert = TRUE,
                     remove = FALSE) |>
            separate(SNP.y, c("chrom.y", "pos.y", "ref.y", "alt.y", "b.y"), sep = "_", convert = TRUE,
                     remove = FALSE) |>
            filter(pos.x < pos.y,
                   pos.y - pos.x <= MAX_DIST)
    }, .groups = "drop") |>
    slice_sample(n = 1e6)
# print(nrow(ld_pairs))

ld <- ld_pairs |>
    # group_by(chrom, AC) |>
    mutate(distance = abs(pos.y - pos.x),
           r2 = r2(geno[SNP.x, , drop = FALSE],
                   geno[SNP.y, , drop = FALSE])) |>
    # str_glue("{unique(chrom)}, AC {unique(AC)}"))) |>
    # ungroup() |>
    # select(-AC, -SNP.x, -SNP.y)
    select(chrom, pos.x, pos.y, distance, r2)

write_tsv(ld, "data/gtex/GTEx_LD.txt.gz")

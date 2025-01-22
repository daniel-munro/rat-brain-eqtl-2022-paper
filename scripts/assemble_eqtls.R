library(tidyverse)

load_tensorqtl <- function(tensorqtl_out) {
    read_tsv(tensorqtl_out, col_types = "ci----ci--dddd-ddd") |>
        rename(gene_id = phenotype_id) |>
        separate_wider_delim(variant_id, ":", names = c("chrom", "pos"),
                             cols_remove = FALSE) |>
        mutate(chrom = str_replace(chrom, "chr", ""),
               pos = as.integer(pos),
               tss = pos - tss_distance, .after = tss_distance) |>
        select(-tss_distance)
}

load_tensorqtl_ind <- function(tensorqtl_out) {
    read_tsv(tensorqtl_out, col_types = "ci----ci--dddd-cd") |>
        rename(gene_id = phenotype_id) |>
        separate_wider_delim(variant_id, ":", names = c("chrom", "pos"),
                 cols_remove = FALSE) |>
        mutate(chrom = str_replace(chrom, "chr", ""),
               pos = as.integer(pos),
               tss = pos - tss_distance, .after = tss_distance) |>
        select(-tss_distance)
}

load_tensorqtl_splice <- function(tensorqtl_out) {
    read_tsv(tensorqtl_out, col_types = "ci----ci--dddd-dcidd") |>
        separate_wider_delim(variant_id, ":", names = c("chrom", "pos"),
                 cols_remove = FALSE) |>
        mutate(chrom = str_replace(chrom, "chr", ""),
               pos = as.integer(pos),
               tss = pos - tss_distance, .after = tss_distance) |>
        select(-tss_distance)
}

load_tensorqtl_ind_splice <- function(tensorqtl_out) {
    read_tsv(tensorqtl_out, col_types = "ci----ci--dddd-dcii") |>
        separate_wider_delim(variant_id, ":", names = c("chrom", "pos"),
                 cols_remove = FALSE) |>
        mutate(chrom = str_replace(chrom, "chr", ""),
               pos = as.integer(pos),
               tss = pos - tss_distance, .after = tss_distance) |>
        select(-tss_distance)
}

load_afc <- function(afc_out) {
    read_tsv(afc_out, col_types = "cc--d--") |>
        rename(gene_id = pid,
               variant_id = sid) |>
        mutate(log2_aFC = signif(log2_aFC, 6))
}

load_ase <- function(ase_out, min_het_inds = 10) {
    read_tsv(ase_out, col_types = "cc--i---d-------------------") |>
        filter(var_het_n >= min_het_inds) |>
        select(gene_id = gene,
               variant_id = var_id,
               log2_aFC_ASE = var_het_afc) |>
        mutate(log2_aFC_ASE = signif(log2_aFC_ASE, 6))
}

genes <- read_tsv("data/genes.txt", col_types = "cc----------") |>
    deframe()

# Codes describe the eQTL pipeline version used:
# - 1st letter of tissue
# - Q for quantile-normalized expression
# - C for use of covariates
# - T for tensorQTL for eQTL mapping
codes <- c(IL = "IQCT", LHb = "LQCT", NAcc = "NQCT", OFC = "OQCT", PL = "PQCT")

alleles <- read_tsv("data/genotype/alleles.txt.gz", col_types = "ccc",
                    col_names = c("variant_id", "ref", "alt"))

afc <- tibble(tissue = names(codes)) |>
    reframe(
        load_afc(str_glue("data/afc/{codes[tissue]}.aFC.txt")),
        .by = tissue
    )

top_assoc <- tibble(tissue = names(codes)) |>
    reframe(
        load_tensorqtl(str_glue("data/tensorqtl/{codes[tissue]}.cis_qtl.txt.gz")),
        .by = tissue
    ) |>
    left_join(afc, by = c("tissue", "gene_id", "variant_id")) |>
    mutate(gene_name = genes[gene_id], .after = gene_id) |>
    relocate(tss, .after = gene_name) |>
    left_join(alleles, by = "variant_id") |>
    relocate(ref, alt, .after = pos)

write_tsv(top_assoc, "data/eqtls/top_assoc.txt")

eqtls_ind <- tibble(tissue = names(codes)) |>
    reframe(
        load_tensorqtl_ind(
            str_glue("data/tensorqtl/{codes[tissue]}.cis_independent_qtl.txt.gz")
        ),
        .by = tissue
    ) |>
    left_join(afc, by = c("tissue", "gene_id", "variant_id")) |>
    mutate(gene_name = genes[gene_id], .after = gene_id) |>
    relocate(tss, .after = gene_name) |>
    left_join(alleles, by = "variant_id") |>
    relocate(ref, alt, .after = pos)

write_tsv(eqtls_ind, "data/eqtls/eqtls_indep.txt")
write_tsv(eqtls_ind, "tables/Supp_Table_2.txt")

ase <- tibble(tissue = names(codes)) |>
    reframe(
        load_ase(str_glue("data/afc/{codes[tissue]}.ASE_aFC.txt")),
        .by = tissue
    )

write_tsv(ase, "data/eqtls/ASE_aFC.txt")

###########
## sQTLs ##
###########

tissues <- c("IL", "LHb", "NAcc", "OFC", "PL")

top_assoc_s <- tibble(tissue = tissues) |>
    reframe(
        load_tensorqtl_splice(str_glue("data/splice/{tissue}_splice.cis_qtl.txt.gz")),
        .by = tissue
    )

write_tsv(top_assoc_s, "data/splice/top_assoc_splice.txt")

sqtls <- tibble(tissue = tissues) |>
    reframe(
        load_tensorqtl_ind_splice(str_glue("data/splice/{tissue}_splice.cis_independent_qtl.txt.gz")),
        .by = tissue
    )

write_tsv(sqtls, "data/splice/sqtls_indep.txt")
write_tsv(sqtls, "tables/Supp_Table_5.txt")

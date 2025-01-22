library(tidyverse)

#############################################
## Supp Table 3: Top tissue-specific eQTLs ##
#############################################
# Get eGenes that are specific to one tissue and have the largest ratio of
# p-values of the significant and the next-strongest association.

eqtls <- read_tsv("data/eqtls/eqtls_indep.txt", col_types = "ccciiciiccdddddid")

top_assoc <- read_tsv("data/eqtls/top_assoc.txt", col_types = "ccciiciiccdddddddd")

spec_eqtl <- top_assoc |>
    mutate(n_sig = sum(qval < 0.05),
           n_tested = n(),
           .by = gene_id) |>
    filter(n_sig == 1,
           # Assume untested gene/tissue combos would be even less significant:
           n_tested >= 2) |>
    group_by(gene_id, gene_name) |>
    arrange(pval_beta) |>
    summarise(variant_id = variant_id[1],
              eQTL_tissue = tissue[1],
              eQTL_log10_pval = -log10(pval_beta[1]),
              next_log10_pval = -log10(pval_beta[2]),
              log10_pval_diff = eQTL_log10_pval - next_log10_pval,
              abs_log2_aFC_diff = abs(log2_aFC[1]) - abs(log2_aFC[2]),
              .groups = "drop") |>
    arrange(desc(log10_pval_diff))

spec_eqtl |>
    mutate(
        eQTL_log10_pval = round(eQTL_log10_pval, 6),
        next_log10_pval = round(next_log10_pval, 6),
        log10_pval_diff = round(log10_pval_diff, 6),
        abs_log2_aFC_diff = round(abs_log2_aFC_diff, 6)
    ) |>
    write_tsv("tables/Supp_Table_3.txt")

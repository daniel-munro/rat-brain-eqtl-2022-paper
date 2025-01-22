library(tidyverse)

##################################################
## Supp Table 4: Top tissue-specific expression ##
##################################################
# The only meaningful comparison I can think of is fraction of samples per tissue
# that are expressed, e.g. read counts >0, >10, etc.
# Sort by difference between 1st and 2nd-highest fraction.

min_reads <- 10

samples <- read_tsv("data/samples.txt", col_types = "cc")
genes <- read_tsv("data/genes.txt", col_types = "cc-----lllll") |>
    filter(in_expr_IL | in_expr_LHb | in_expr_NAcc | in_expr_OFC | in_expr_PL)
gene_names <- genes |>
    select(gene_id, gene_name) |>
    deframe()

# Use raw count matrix to include all genes in all tissues.
expr <- read.table("data/expression/ensembl-gene_raw-counts.txt",
                   check.names = FALSE) |>
    as_tibble(rownames = "gene_id") |>
    filter(gene_id %in% genes$gene_id) |>
    pivot_longer(-gene_id, names_to = "library", values_to = "expr") |>
    filter(library %in% samples$library) |>
    left_join(samples, by = "library", relationship = "many-to-one")
expr_frac <- expr |>
    summarise(frac_expr = mean(expr >= min_reads),
              frac_nonzero = mean(expr > 0),
              .by = c(brain_region, gene_id))

spec_expr <- expr_frac |>
    group_by(gene_id) |>
    arrange(desc(frac_expr)) |>
    summarise(top_expr_tissue = brain_region[1],
              top_expr_frac = frac_expr[1],
              next_expr_frac = frac_expr[2],
              expr_frac_diff = frac_expr[1] - frac_expr[2],
              .groups = "drop") |>
    filter(expr_frac_diff > 0.9) |>
    arrange(desc(expr_frac_diff)) |>
    mutate(gene_name = gene_names[gene_id], .after = gene_id)

spec_expr |>
    mutate(
        top_expr_frac = round(top_expr_frac, 6),
        next_expr_frac = round(next_expr_frac, 6),
        expr_frac_diff = round(expr_frac_diff, 6),
    ) |>
    write_tsv("tables/Supp_Table_4.txt")

########################################################################
## Determine read count cutoff for a gene being expressed in a tissue ##
########################################################################

egenes <- filter(top_assoc, qval < 0.05)

expr_cutoffs <- tibble(cutoff = c(1, 2, 5, 10, 20, 50, 100)) |>
    reframe(
        expr |>
            summarise(frac_expr = mean(expr >= cutoff),
                      .by = c(brain_region, gene_id)),
        .by = cutoff
    ) |>
    left_join(
        egenes |>
            select(brain_region = tissue, gene_id) |>
            mutate(is_eGene = TRUE),
        by = c("brain_region", "gene_id"),
        relationship = "many-to-one"
    ) |>
    replace_na(list(is_eGene = FALSE))

expr_cutoffs |>
    mutate(cutoff = as.factor(cutoff)) |>
    ggplot(aes(x = cutoff, y = frac_expr, fill = is_eGene)) +
    geom_violin(scale = "width") +
    geom_boxplot(alpha = 0.5, outlier.size = 0.25)

expr |>
    summarise(median_count = median(expr),
              .by = c(brain_region, gene_id)) |>
    left_join(
        egenes |>
            select(brain_region = tissue, gene_id) |>
            mutate(is_eGene = TRUE),
        by = c("brain_region", "gene_id"),
        relationship = "one-to-one"
    ) |>
    replace_na(list(is_eGene = FALSE)) |>
    ggplot(aes(x = median_count, fill = is_eGene)) +
    geom_density(alpha = 0.5) +
    scale_x_log10()

# I'm not seeing a good indicator for which cutoff to choose, so I think 10 is fine.

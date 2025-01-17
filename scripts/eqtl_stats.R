suppressPackageStartupMessages(library(tidyverse))

top_assoc <- read_tsv("data/eqtls/top_assoc.txt", col_types = "ccciiciiccdddddddd")

egenes <- filter(top_assoc, qval < 0.05)

eqtls <- read_tsv("data/eqtls/eqtls_indep.txt", col_types = "ccciiciiccdddddid")

# Number of significant eGenes per tissue:
top_assoc |>
    group_by(tissue) |>
    summarise(n = sum(qval < 0.05),
              total = n()) |>
    mutate(frac = n / total)

# Total unique eGenes:
n_distinct(egenes$gene_id)

# eGenes found in N tissues:
egenes |>
    group_by(gene_id) |>
    summarise(n_tissues = n(), .groups = "drop") |>
    count(n_tissues) |>
    mutate(frac = n / sum(n))

# Average number of eGenes with 2 and 3 independent eQTLs
eqtls |>
    count(tissue, gene_id, name = "eQTLs") |>
    count(tissue, eQTLs) |>
    group_by(eQTLs) |>
    summarise(mean_n = mean(n))

# Additional eQTLs per tissue
eqtls |>
    group_by(tissue) |>
    summarise(n = sum(rank > 1),
              perc_incr = n / sum(rank == 1),
              .groups = "drop") |>
    summarise(mean_n = mean(n),
              mean_perc_incr = mean(perc_incr))

# Total genes with 2+ eQTLs:
eqtls |>
    count(tissue, gene_id, name = "eQTLs") |>
    filter(eQTLs > 1) |>
    distinct(gene_id) |>
    count()

# Count samples:
samples <- read_tsv("data/samples.txt", col_types = "cc") |>
    separate(library, into = c("rat_id", "tissue"))
count(samples, brain_region)
n_distinct(samples$rat_id)

# Total cis-QTLs:
sqtls <- read_tsv("data/splice/sqtls_indep.txt", col_types = "cciccicdddddcii")
nrow(eqtls) + nrow(sqtls)

###########################
## Count expressed genes ##
###########################

expr <- eqtls |>
    distinct(tissue) |>
    group_by(tissue) |>
    summarise(
        read_tsv(
            str_glue("data/expression/ensembl-gene_log2_{tissue}.bed.gz"),
            col_types = cols(`#chr` = '-', start = '-', end = '-', gene_id = 'c', .default = 'd')
        ) |>
            pivot_longer(-gene_id, names_to = 'sample', values_to = 'expr'),
        .groups = "drop"
    )

# These tables are already filtered, so I'll just count the genes in them:
expr |>
    group_by(tissue, gene_id) |>
    filter(sum(expr > 0) == 0)

# No. expressed and fraction of expressed genes that had eQTLs:
expr |>
    distinct(tissue, gene_id) |>
    left_join(
        select(egenes, tissue, gene_id, qval),
        by = c("tissue", "gene_id")
    ) |>
    group_by(tissue) |>
    summarise(n_expr = n(),
              n_egene = sum(!is.na(qval)),
              frac_egene = mean(!is.na(qval)),
              .groups = "drop")

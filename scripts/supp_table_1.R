library(tidyverse)

############################################
## Supp Table 1: Seq/sample summary stats ## 
############################################

tissues <- c(Acbc = "NAcc", IL = "IL", LHB = "LHb", PL = "PL", VoLo = "OFC")
meta <- read_csv("data/qc/metadata/metadata_p50_hao_chen_2014.csv",
                 col_types = "ccccDDcDddDccciiddciiddcccccc") |>
    mutate(tissue = tissues[brain_region])

sample_counts <- meta |>
    group_by(tissue) |>
    summarise(
        n_initial_samples = n(),
        n_sequencing_QC_pass = sum(QC_reads == "pass"),
        n_geno_match_QC_pass = sum(QC_reads == "pass" & QC_geno_match == "pass"),
        n_expression_QC_pass = sum(QC_pass == "pass"),
    )

# Read counts etc. are just for final sample sets
other_stats <- meta |>
    filter(QC_pass == "pass") |>
    group_by(tissue) |>
    summarise(
        raw_reads_total = sum(raw_reads),
        raw_reads_mean = mean(raw_reads),
        raw_reads_SD = sd(raw_reads),
        unique_mapped_reads_mean = mean(uniq_mapped_reads),
        unique_mapped_reads_SD = sd(uniq_mapped_reads),
    )

# Counts of reads used in gene quantification for final sample sets
samples <- read_tsv("data/samples.txt", col_types = "cc")
expr_reads <- read.csv(
    "data/expression/ensembl-gene_raw-counts.txt",
    sep = "\t",
    check.names = FALSE
)
expr_reads <- expr_reads[, colnames(expr_reads) %in% samples$library]
expr_reads_stats <- expr_reads |>
    as_tibble(rownames = "gene_id") |>
    pivot_longer(-gene_id, names_to = "library", values_to = "count") |>
    mutate(rat_id = str_sub(library, 1, 10),
           tissue = str_sub(library, 12),
           tissue = tissues[tissue]) |> # For speed
    group_by(tissue, rat_id) |>
    summarise(count = sum(count),
              .groups = "drop") |>
    group_by(tissue) |>
    summarise(expr_quant_reads_mean = mean(count),
              expr_quant_reads_SD = sd(count),
              .groups = "drop")

# Expressed gene count
expr <- tibble(tissue = c("IL", "LHb", "NAcc", "OFC", "PL")) |>
    group_by(tissue) |>
    summarise(
        read_tsv(str_glue("data/expression/ensembl-gene_log2_{tissue}.bed.gz"),
                 col_types = cols(gene_id = "c", .default = "-")) |>
            count(name = "n_expressed_genes"),
        .groups = "drop"
    )

# Combine
d <- sample_counts |>
    left_join(other_stats, by = "tissue") |>
    left_join(expr_reads_stats, by = "tissue") |>
    left_join(expr, by = "tissue") |>
    pivot_longer(-tissue, names_to = "Stat") |>
    mutate(value = value |> round()) |>
    pivot_wider(id_cols = Stat, names_from = tissue, values_from = value)

write_tsv(d, "tables/Supp_Table_1.txt")

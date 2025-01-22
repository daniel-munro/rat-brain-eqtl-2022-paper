library(tidyverse)

## Save CSV files that can be opened in Excel and copied into GEO metadata files.

d <- read_csv("data/qc/metadata/metadata_p50_hao_chen_2014.csv",
              col_types = cols(.default = "c")) |>
    select(library, brain_region, sex, rat_batch, QC_reads, QC_geno_match,
           QC_expr_outlier, QC_expr_zeros, QC_pass, file_locations) |>
    arrange(brain_region) |>
    # There are identical file names in different batches, so make unique:
    mutate(file_locations = str_replace_all(file_locations, "/", "_")) |>
    mutate(id = str_c("Sample ", 1:n()), .before = 1, .by = brain_region) |>
    mutate(description = case_when(
        QC_reads == "reject" ~ "Removed: low read count",
        QC_geno_match == "reject" ~ "Removed: genotype mismatch",
        QC_expr_outlier == "reject" ~ "Removed: expression outlier",
        QC_expr_zeros == "reject" ~ "Removed: too many zero counts",
        TRUE ~ ""
    ), .before = file_locations)

write_csv(d, "data/qc/GEO/sample_table.csv")

###############
## Raw files ##
###############

fastq <- tibble(batchnum = 1:5) |>
    reframe(
        read_table(str_glue("data/qc/fastq_md5/md5_batch{batchnum}.txt"), col_types = "cc",
                   col_names = c("md5sum", "file")),
        .by = batchnum
    ) |>
    separate_wider_delim(file, "/", names = c("batch", "file"), cols_remove = TRUE) |>
    filter(str_detect(file, "Undetermined", negate = TRUE)) |>
    mutate(tissue = str_split(file, "_", simplify = TRUE)[, 2],
           type = "fastq") |>
    arrange(tissue, file) |>
    mutate(number = 1:n(), .before = 1, .by = tissue) |>
    # There are identical file names in different batches, so make unique:
    mutate(file = str_glue("{batch}_{file}")) |>
    select(number, tissue, file, type, md5sum)

write_csv(fastq, "data/qc/GEO/fastq_table.csv")

library(tidyverse)

date_clean <- function(dates) {
    map_chr(dates, ~as.character(as.Date(.x, tryFormats = c("%Y-%m-%d", "%m/%d/%Y"))))
}

init_samples <- read_tsv("data/library_info.txt", col_types = "cccicccccccddc") %>%
    # Dates in D/M/YYYY and YYYY-MM-DD formats are intermixed:
    mutate(BrainRegion_DissDate = date_clean(BrainRegion_DissDate),
           sequencing_batch = str_replace_all(sequencing_batch, "_", ";"),
           RNADate = date_clean(RNADate),
           NanoDropDate = date_clean(NanoDropDate),
           DateShipped = date_clean(DateShipped))
# Rename library/rat_id for swapped pair:
i <- which(init_samples$library == "000789FFF0_LHB")
j <- which(init_samples$library == "000789FFF9_LHB")
init_samples$library[i] <- "000789FFF9_LHB"
init_samples$library[j] <- "000789FFF0_LHB"
init_samples$rat_id[i] <- "000789FFF9"
init_samples$rat_id[j] <- "000789FFF0"

rats <- read_tsv("data/rat_info.txt", col_types = "ciccccccccTiiddcccc") %>%
    rename(rfid = rat_id,
           sex = sex_mf)

fastq <- read_tsv("data/qc/metadata/fastq_files_listing.txt", col_types = "cccci") %>%
    mutate(path = str_glue("batch{batch}/{filename}")) %>%
    group_by(library) %>%
    summarise(file_locations = str_c(path, collapse = ";"),
              .groups = "drop")
# Rename library/file_locations for swapped pair:
i <- which(fastq$library == "000789FFF0_LHB")
j <- which(fastq$library == "000789FFF9_LHB")
fastq$library[i] <- "000789FFF9_LHB"
fastq$library[j] <- "000789FFF0_LHB"
fastq$file_locations[i] <- str_replace_all(fastq$file_locations[i], "FFF0", "FFF9")
fastq$file_locations[j] <- str_replace_all(fastq$file_locations[j], "FFF9", "FFF0")

reads <- read_tsv("data/qc/metadata/read_stats.txt", col_types = cols(
    file = "c",
    `Number of input reads` = "i",
    `Uniquely mapped reads number` = "i",
    `Uniquely mapped reads %` = "c",
    `% of reads unmapped: too many mismatches` = "c",
    `% of reads unmapped: too short` = "c",
    `% of reads unmapped: other` = "c",
    .default = "-"
)) %>%
    rename(
        raw_reads = `Number of input reads`,
        uniq_mapped_reads = `Uniquely mapped reads number`,
        uniq_mapped_ratio = `Uniquely mapped reads %`,
        unmapped_mismatch_ratio = `% of reads unmapped: too many mismatches`,
        unmapped_short_ratio = `% of reads unmapped: too short`,
        unmapped_other_ratio = `% of reads unmapped: other`
    ) %>%
    mutate(across(ends_with("ratio"),
                  ~as.double(str_replace(.x, "%", "")) / 100),
           unmapped_reads_ratio = unmapped_mismatch_ratio + unmapped_short_ratio + 
               unmapped_other_ratio,
           library = str_extract(file, "\\w+_\\w+")) %>%
    select(-file,
           -unmapped_mismatch_ratio,
           -unmapped_short_ratio,
           -unmapped_other_ratio)

seq_problems <- read_lines("~/Dropbox/HS-RNASeq/info/problem_samples_sequencing.txt")
expr_outliers <- read_lines("data/qc/metadata/expr_outliers.txt")
expr_zeros <- c("000789FFF8_Acbc", "00078A0215_IL")

final_samples <- read_tsv("data/samples.txt", col_types = "cc") %>%
    arrange(brain_region, library)

# library, rfid, brain_region, etc.
d <- init_samples %>%
    select(library,
           rfid = rat_id,
           brain_region,
           sequencing_batch,
           brain_region_diss_date = BrainRegion_DissDate,
           RNA_date = RNADate,
           RNA_method = RNAMethod,
           NanoDrop_date = NanoDropDate,
           NanoDrop_conc = NanoDropConc,
           NanoDrop_ratio_260_280 = NanoDropRatio260.280,
           date_shipped = DateShipped) %>%
    filter(rfid != "Undetermined") # Excluded from expression quantitation

# sex, rat_batch, etc.
d <- d %>%
    left_join(
        rats %>%
            select(rfid,
                   sex,
                   rat_batch,
                   rat_dissection_datetime = datetime_dissected,
                   rat_dissection_age = dissection_age,
                   rat_body_weight_g = body_weight_g,
                   rat_length_wo_tail_cm = length_wo_tail_cm,
                   rat_length_w_tail_cm = length_w_tail_cm),
        by = "rfid"
    )

# file_locations
d <- d %>%
    left_join(fastq, by = "library")

# raw_reads, uniq_mapped_reads, uniq_mapped_ratio, unmapped_reads_ratio
d <- d %>%
    left_join(reads, by = "library")

# QC_reads
d <- d %>%
    mutate(QC_reads = if_else(library %in% seq_problems, "reject", "pass"))

# QC_geno_match
geno_dups <- c("00078A0224_Acbc" = "00078A0244",
               "00078A2667_Acbc" = "000789FFF8",
               "000789FFF8_IL" = "00078A07A2")
geno_nomatch <- c("00078A18A7_Acbc", "00078A18A7_IL", "00078A2667_IL",
                  "00078A18A7_LHB", "00078A18A7_PL", "00078A2667_PL",
                  "00078A18A7_VoLo", "00078A2667_VoLo")
d <- d %>%
    mutate(QC_geno_match = if_else(library %in% c(names(geno_dups), geno_nomatch),
                                   "reject", "pass"))

# geno_match_comment
d <- d %>%
    mutate(geno_match_comment = case_when(
        library == "000789FFF0_LHB" ~ "Original label before genotype mismatch fix: 000789FFF9_LHB",
        library == "000789FFF9_LHB" ~ "Original label before genotype mismatch fix: 000789FFF0_LHB",
        library %in% geno_nomatch ~ "No matching genotype found",
        TRUE ~ as.character(NA)
    )) %>%
    mutate(geno_match_comment = if_else(library %in% names(geno_dups),
                                        str_glue("Matches {geno_dups[library]} genotype and is therefore apparent duplicate"),
                                        geno_match_comment))

# QC_expr_outlier
d <- d %>%
    mutate(QC_expr_outlier = if_else(library %in% expr_outliers, "reject", "pass"))

# QC_expr_zeros
d <- d %>%
    mutate(QC_expr_zeros = if_else(library %in% expr_zeros, "reject", "pass"))

# QC_pass
d <- d %>%
    mutate(QC_pass = if_else(QC_reads == "pass" &
                                 QC_geno_match == "pass" &
                                 QC_expr_outlier == "pass" &
                                 QC_expr_zeros == "pass",
                             "pass", "reject"))

d %>%
    arrange(library) %>% # Because of swapped samples
    mutate(across(ends_with("ratio"), as.character)) %>%  # Avoid precision issue
    write_csv("data/qc/metadata/metadata_p50_hao_chen_2014.csv")

# Confirm QC_pass matches final_samples
final <- d %>%
    filter(QC_pass == "pass") %>%
    select(brain_region, library) %>%
    arrange(brain_region, library)
identical(final$library, final_samples$library)

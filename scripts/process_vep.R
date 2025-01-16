library(GenomicRanges)
library(tidyverse)

genes <- read_tsv("data/genes.txt", col_types = "c-c---illlll") %>%
    filter(in_expr_IL | in_expr_LHb | in_expr_NAcc | in_expr_OFC | in_expr_PL)
cis_rng <- with(genes, GRanges(chrom, IRanges(tss - 1e6, tss + 1e6)))

vep_all <- read_tsv("data/vep/vep.txt.gz", comment = "##", na = "-",
                    col_types = "cc-c-cc------c") %>%
    rename(variant_id = `#Uploaded_variation`)
vep_all_rng <- vep_all %>%
    separate(variant_id, c("chrom", "pos"), sep = ":", convert = TRUE) %>%
    mutate(chrom = str_replace(chrom, "chr", "")) %>%
    with(GRanges(chrom, IRanges(pos, pos)))

labels <- c(
    "5' UTR" = "5_prime_UTR_variant",
    "3' UTR" = "3_prime_UTR_variant",
    "Missense" = "missense_variant",
    "Synonymous" = "synonymous_variant",
    "Splice region"= "splice_region_variant",
    "Intron" = "intron_variant",
    "Noncoding transcript" = "non_coding_transcript_variant",
    "Noncoding transcript" = "non_coding_transcript_exon_variant",
    "Intergenic - upstream" = "upstream_gene_variant",
    "Intergenic - downstream" = "downstream_gene_variant",
    "Intergenic - other" = "intergenic_variant"
)

vep <- vep_all %>%
    filter(countOverlaps(vep_all_rng, cis_rng) > 0) %>%
    select(variant_id, Consequence, Gene) %>%
    separate_rows(Consequence, sep = ",") %>%
    group_by(Consequence) %>%
    filter(
        n() > 100, # Only drops 0.01% of annotations.
        Consequence != "NMD_transcript_variant" # Doesn't describe the variant itself
    ) %>%
    ungroup() %>%
    mutate(Consequence = Consequence %>%
               # str_replace("non_coding_transcript_exon_variant",
               #             "non_coding_transcript_variant") %>%
               fct_recode(!!!labels)) %>%
    distinct(variant_id, Consequence, Gene)

write_tsv(vep, "data/vep/processed_vep.txt.gz")

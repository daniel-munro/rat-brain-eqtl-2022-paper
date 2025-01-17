library(tidyverse)

human_genes <- read_tsv("data/gtex/gtex_genes.txt", col_types = "cc",
                        col_names = c("id", "name")) |>
    select(name, id) |>
    deframe()

# Subset to ortholog pairs with SDg used in the figure:
orthologs <- read_tsv("data/gtex/orthologs/ortholog_SDg.tsv", col_types = "-c--") |>
    distinct() |>
    pull()

gsets <- read_tsv("data/human/GeneSets.txt", col_types = cols(.default = "c"), comment = "#") |>
    # All but first 2 columns are gene symbols:
    mutate(across(-starts_with("Autosomal"), ~ human_genes[.x])) |>
    pivot_longer(everything(), names_to = "set", values_to = "gene_id_human") |>
    mutate(set = set |>
               str_replace("_zz", " [") |>
               str_replace("zz", "]") |>
               str_replace("PMID", "PMID:") |>
               str_replace_all("_", " ")) |>
    filter(!is.na(gene_id_human)) |>
    distinct()

##################
## GWAS Catalog ##
##################

gwas <- read_tsv("data/human/gwas_catalog_v1.0.2-associations_e105_r2022-02-21.tsv.gz",
                 col_types = cols(STUDY = "c", `DISEASE/TRAIT` = "c",
                                  `REPORTED GENE(S)` = "c", .default = "-")) |>
    rename(study = STUDY,
           trait = `DISEASE/TRAIT`,
           gene = `REPORTED GENE(S)`) |>
    filter(!is.na(gene),
           !(gene %in% c('NR', 'Intergenic', 'intergenic', 'NR x NR'))) |>
    separate_rows(gene, sep = ", ") |>
    mutate(trait = trait |>
               str_replace(" \\(amplitude at temporal datapoints\\)", "") |>
               str_replace("[Bb]ody mass index", "BMI") |>
               str_replace("Autism spectrum disorder", "ASD") |>
               str_replace("Electrocardiogram", "ECG"))

gwas_top_10 <- gwas |>
    count(trait, sort = TRUE) |>
    filter(!(trait %in% c(
        "ASD or schizophrenia",
        "Hip circumference adjusted for BMI",
        "Waist-to-hip ratio adjusted for BMI",
        "Waist circumference adjusted for BMI",
        "Waist-hip index",
        "Serum metabolite levels",
        "A body shape index"
    ))) |>
    slice(1:10) |>
    pull(trait)

gwas <- gwas |>
    mutate(gene_id_human = human_genes[gene]) |>
    filter(!is.na(gene_id_human))

gwas |>
    count(gene, sort = TRUE)

gwas |>
    count(study, sort = TRUE)
gwas |>
    count(trait, sort = TRUE)

#############
## Combine ##
#############

df <- bind_rows(
    gsets |>
        filter(gene_id_human %in% orthologs,
               str_detect(set, "GWAS", negate = TRUE)) |>
        group_by(set) |>
        filter(n() >= 20) |>
        ungroup(),
    gwas |>
        filter(gene_id_human %in% orthologs,
               trait %in% gwas_top_10) |>
        select(set = trait, gene_id_human) |>
        mutate(set = str_c(set, " [GWAS]"))
)

write_tsv(df, "data/gtex/orthologs/gene_sets.tsv")

suppressPackageStartupMessages(library(tidyverse))

ase <- read_tsv("data/eqtls/ASE_aFC.txt", col_types = "cccd")

eqtls <- read_tsv("data/eqtls/eqtls_indep.txt", col_types = "ccciiciiccdddddid") |>
    left_join(ase, by = c("tissue", "gene_id", "variant_id"),
              relationship = "one-to-one")

counts <- eqtls |>
    summarise(n_eQTLs = n(),
              n_eQTLs_with_ASE = sum(!is.na(log2_aFC_ASE)),
              .by = tissue)
corrs <- eqtls |>
    filter(!is.na(log2_aFC_ASE) & !is.na(log2_aFC)) |>
    summarise(
        R = cor(log2_aFC_ASE, log2_aFC),
        rho = cor(log2_aFC_ASE, log2_aFC, method = "spearman"),
        beta = deming::deming(log2_aFC ~ log2_aFC_ASE)$coefficients[2],
        .by = tissue
    ) |>
    left_join(counts, by = "tissue", relationship = "one-to-one") |>
    mutate(
        stats = str_c("r=", format(R, digits = 2),
                      " rho=", format(rho, digits = 2),
                      " β=", format(beta, digits = 2)),
        count = str_glue("n={n_eQTLs_with_ASE} (of {n_eQTLs})")
    )

eqtls |>
    filter(!is.na(log2_aFC_ASE),
           !is.na(log2_aFC)) |>
    ggplot(aes(x = log2_aFC_ASE, y = log2_aFC)) +
    facet_wrap(~ tissue, nrow = 1) +
    geom_point(size = 0.25, alpha = 0.5) +
    geom_text(aes(x = -7, y = 7, label = stats), data = corrs, hjust = "left") +
    geom_text(aes(x = -7, y = 5.2, label = count), data = corrs, hjust = "left") +
    expand_limits(x = c(-7, 7), y = c(-7, 7)) +
    theme_minimal() +
    xlab(expression(log[2]*"aFC (ASE)")) +
    ylab(expression(log[2]*"aFC (eQTL model)"))

ggsave("figures/FigureS2.png", width = 12, height = 3, bg = "white")

###########
## Stats ##
###########

skimr::skim(corrs)

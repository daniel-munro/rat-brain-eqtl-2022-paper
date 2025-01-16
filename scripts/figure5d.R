# Based on rat_human_aFC.Rmd

library(tidyverse)

eqtls <- read_tsv("data/eqtls/eqtls_indep.txt", col_types = "ccciiciiccdddddid") |>
    filter(rank == 1) |>
    group_by(gene_id) |>
    summarise(mean_abs_rat = mean(abs(log2_aFC)))

# # Get homologs
# # (https://www.r-bloggers.com/2016/10/converting-mouse-to-human-gene-names-with-biomart-package/)
# human <- biomaRt::useMart("ensembl", dataset = "hsapiens_gene_ensembl")
# rat <- biomaRt::useMart("ensembl", dataset = "rnorvegicus_gene_ensembl")
# gene_map <- biomaRt::getLDS(
#     attributes = c("ensembl_gene_id"),
#     # filters = "ensembl_gene_id",
#     # values = unique(eqtls$gene_id),
#     mart = rat,
#     attributesL = c("ensembl_gene_id_version"),
#     martL = human,
#     uniqueRows = TRUE
# ) |>
#     rename(gene_id_rat = Gene.stable.ID,
#            gene_id_human = Gene.stable.ID.version)
# write_tsv(gene_map, "data/gtex/orthologs.txt")
gene_map <- read_tsv("data/gtex/orthologs.txt", col_types = "cc") |>
    mutate(gene_id_human = str_replace(gene_id_human, "\\..+$", ""))

gtex <- list.files("data/gtex/GTEx_Analysis_v8_eQTL", pattern = "*Brain_*",
                   full.names = TRUE) |>
    read_tsv(id = "file",
             col_types = cols(gene_id = "c", qval = "d", pval_beta = "-",
                              log2_aFC = "d", .default = "-")) |>
    filter(qval <= 0.05) |>  # GTEx site says to do <=
    mutate(gene_id = str_replace(gene_id, "\\..+$", "")) |>
    group_by(gene_id) |>
    summarise(mean_abs_human = mean(abs(log2_aFC)),
              .groups = "drop")

# All homolog pairs with eQTLs in rat and human brain:
d <- gene_map |>
    inner_join(eqtls, by = c("gene_id_rat" = "gene_id")) |>
    inner_join(gtex, by = c("gene_id_human" = "gene_id"))
# write_tsv(d, "data/gtex/orthologs/ortholog_aFC.tsv")

stats <- d |>
    summarise(
        pairs = n(),
        pearson = cor(mean_abs_rat, mean_abs_human),
        pearson_p = cor.test(mean_abs_rat, mean_abs_human)$p.value,
        spearman = cor(mean_abs_rat, mean_abs_human, method = "spearman"),
        spearman_p = cor.test(mean_abs_rat, mean_abs_human, method = "spearman")$p.value,
    )

# deming <- d |>
#     summarise({
#         coef <- deming::deming(mean_abs_human ~ mean_abs_rat)$coefficients
#         tibble(intercept = coef[1],
#                slope = coef[2])
#     })
# lab <- str_c(str_glue("rho={format(stats$spearman, digits = 2)}"),
#              str_glue("   P={format(stats$spearman_p, digits = 2)}"),
#              sep = "\n")
lab <- str_c(str_glue(" r = {format(stats$pearson, digits = 2)}"),
             str_glue("P = {format(stats$pearson_p, digits = 2)}"),
             sep = "\n")
d |>
    mutate(mean_abs_rat = pmin(mean_abs_rat, 4),
           mean_abs_human = pmin(mean_abs_human, 4)) |>
    ggplot(aes(x = mean_abs_rat, y = mean_abs_human)) +
    geom_point(size = 0.5, alpha = 0.3) +
    # geom_abline(aes(intercept = intercept, slope = slope), data = deming,
    #             color = "#8888cc") +
    annotate(geom = "text", x = 2.2, y = 3.5, hjust = 0, label = lab, color = "#5555cc") +
    coord_fixed() +
    # scale_x_continuous(expand = c(0, 0), limits = c(0, max(d$mean_abs_rat) + 0.2)) +
    # scale_y_continuous(expand = c(0, 0), limits = c(0, max(d$mean_abs_human) + 0.2)) +
    scale_x_continuous(expand = c(0, 0.1)) +
    scale_y_continuous(expand = c(0, 0.1)) +
    xlab(expression("Mean |"*log[2]*"aFC| (Rat)")) +
    ylab(expression("Mean |"*log[2]*"aFC| (Human)")) +
    theme_bw() +
    theme(
        panel.grid = element_blank(),
    )

ggsave("figures/figure5/figure5d.png", width = 2.5, height = 2.5, bg = "white")

#####################################
## Correlation with non-brain GTEx ##
#####################################

gtex_nonbrain <- list.files("data/gtex/GTEx_Analysis_v8_eQTL", full.names = TRUE)
gtex_nonbrain <- gtex_nonbrain[str_detect(gtex_nonbrain, "Brain", negate = TRUE)] |>
    read_tsv(id = "file",
             col_types = cols(gene_id = "c", qval = "d", pval_beta = "-",
                              log2_aFC = "d", .default = "-")) |>
    filter(qval <= 0.05) |>  # GTEx site says to do <=
    mutate(gene_id = str_replace(gene_id, "\\..+$", "")) |>
    group_by(gene_id) |>
    summarise(mean_abs_human = mean(abs(log2_aFC)),
              .groups = "drop")

# All homolog pairs with eQTLs in rat and human brain:
d_nonbrain <- gene_map |>
    inner_join(eqtls, by = c("gene_id_rat" = "gene_id")) |>
    inner_join(gtex_nonbrain, by = c("gene_id_human" = "gene_id"))
# write_tsv(d, "analysis/orthologs/ortholog_aFC.tsv")

stats_nonbrain <- d_nonbrain |>
    summarise(
        pairs = n(),
        pearson = cor(mean_abs_rat, mean_abs_human),
        pearson_p = cor.test(mean_abs_rat, mean_abs_human)$p.value,
        spearman = cor(mean_abs_rat, mean_abs_human, method = "spearman"),
        spearman_p = cor.test(mean_abs_rat, mean_abs_human, method = "spearman")$p.value,
    )
stats_nonbrain

#########################################
## Do high-aFC eQTLs also have low AF? ##
#########################################

gtex2 <- list.files("data/gtex/GTEx_Analysis_v8_eQTL", pattern = "*Brain_*",
                   full.names = TRUE) |>
    read_tsv(id = "file",
             col_types = cols(gene_id = "c", qval = "d", maf = "d",
                              log2_aFC = "d", .default = "-")) |>
    filter(qval <= 0.05) |>  # GTEx site says to do <=
    group_by(gene_id) |>
    summarise(mean_maf = mean(maf),
              mean_abs_human = mean(abs(log2_aFC)),
              .groups = "drop")
gtex2 |>
    ggplot(aes(x = mean_maf, y = mean_abs_human)) +
    geom_point(size = 0.2, alpha = 0.5) +
    geom_smooth

library(tidyverse)
library(patchwork)

labels <- c("5' UTR", "3' UTR", "Missense", "Synonymous", "Splice region", "Intron",
            "Noncoding transcript", "Intergenic - upstream", "Intergenic - downstream",
            "Intergenic - other")

esnps <- read_tsv("data/eqtls/eqtls_indep.txt", col_types = "ccciiciiccdddddid") |>
    distinct(tissue, variant_id)

ssnps <- read_tsv("data/splice/sqtls_indep.txt", col_types = "cciccicdddddcii") |>
    distinct(tissue, variant_id)

vep <- read_tsv("data/vep/processed_vep.txt.gz", col_types = "cc-") |>
    mutate(
        esnp_IL = variant_id %in% esnps$variant_id[esnps$tissue == "IL"],
        esnp_LHb = variant_id %in% esnps$variant_id[esnps$tissue == "LHb"],
        esnp_NAcc = variant_id %in% esnps$variant_id[esnps$tissue == "NAcc"],
        esnp_OFC = variant_id %in% esnps$variant_id[esnps$tissue == "OFC"],
        esnp_PL = variant_id %in% esnps$variant_id[esnps$tissue == "PL"],
        ssnp_IL = variant_id %in% ssnps$variant_id[ssnps$tissue == "IL"],
        ssnp_LHb = variant_id %in% ssnps$variant_id[ssnps$tissue == "LHb"],
        ssnp_NAcc = variant_id %in% ssnps$variant_id[ssnps$tissue == "NAcc"],
        ssnp_OFC = variant_id %in% ssnps$variant_id[ssnps$tissue == "OFC"],
        ssnp_PL = variant_id %in% ssnps$variant_id[ssnps$tissue == "PL"],
    )

snp_counts <- vep |>
    distinct(variant_id, .keep_all = TRUE) |>
    with(c(esnp_IL = sum(esnp_IL),
           esnp_LHb = sum(esnp_LHb),
           esnp_NAcc = sum(esnp_NAcc),
           esnp_OFC = sum(esnp_OFC),
           esnp_PL = sum(esnp_PL),
           ssnp_IL = sum(ssnp_IL),
           ssnp_LHb = sum(ssnp_LHb),
           ssnp_NAcc = sum(ssnp_NAcc),
           ssnp_OFC = sum(ssnp_OFC),
           ssnp_PL = sum(ssnp_PL)))

enrich <- vep |>
    group_by(Consequence) |>
    summarise(esnp_IL = sum(esnp_IL),
              esnp_LHb = sum(esnp_LHb),
              esnp_NAcc = sum(esnp_NAcc),
              esnp_OFC = sum(esnp_OFC),
              esnp_PL = sum(esnp_PL),
              ssnp_IL = sum(ssnp_IL),
              ssnp_LHb = sum(ssnp_LHb),
              ssnp_NAcc = sum(ssnp_NAcc),
              ssnp_OFC = sum(ssnp_OFC),
              ssnp_PL = sum(ssnp_PL),
              n_total = n(),
              .groups = "drop") |>
    pivot_longer(esnp_IL:ssnp_PL, names_to = "type_tissue", values_to = "n_xsnp") |>
    mutate(
        frac_xsnp = n_xsnp / snp_counts[type_tissue],
        frac_total = n_total / n_distinct(vep$variant_id),
        log2_enrich = log2(frac_xsnp / frac_total)
    ) |>
    separate(type_tissue, c("type", "tissue")) |>
    mutate(Consequence = fct_relevel(Consequence, labels)) |>
    arrange(Consequence)

p1 <- enrich |>
    filter(type == "esnp") |>
    group_by(Consequence) |>
    summarise(mean_enr = mean(log2_enrich),
              sd_enr = sd(log2_enrich),
              .groups = "drop") |>
    mutate(Consequence = fct_rev(Consequence)) |>
    ggplot(aes(x = Consequence, y = mean_enr, ymin = mean_enr - sd_enr,
               ymax = mean_enr + sd_enr, color = type)) +
    geom_pointrange(fatten = 1) +
    geom_hline(yintercept = 0, lty = 2) +
    coord_flip() +
    theme_bw() +
    theme(panel.grid = element_blank()) +
    xlab(NULL) +
    ylab(expression(log[2]*" fold enrichment in eSNPs"))

p2 <- enrich |>
    group_by(Consequence, type) |>
    summarise(mean_frac = mean(frac_xsnp),
              sd_frac = sd(frac_xsnp),
              .groups = "drop") |>
    mutate(Consequence = fct_rev(Consequence),
           type = type |> fct_recode(eSNPs = "esnp", sSNPs = "ssnp") |> fct_rev()) |>
    ggplot(aes(x = Consequence, y = mean_frac, ymin = mean_frac - sd_frac,
               ymax = mean_frac + sd_frac, fill = type)) +
    geom_col(position = "dodge") +
    geom_linerange(position = position_dodge(width = 1)) +
    coord_flip() +
    scale_y_continuous(breaks = c(0, 0.2, 0.4)) +
    guides(fill = guide_legend(reverse = TRUE)) +
    theme_bw() +
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          panel.grid = element_blank(),
          legend.position = c(0.6, 0.75),
          legend.key.size = unit(10, "pt"),
          legend.spacing.x = unit(5, "pt"),
          legend.title = element_blank()) +
    xlab(NULL) +
    ylab("Proportion\nof e/sSNPs")

p1 + p2 + plot_layout(widths = c(3, 2))
ggsave("figures/figure4/figure4g.png", width = 6, height = 2.2)

library(tidyverse)

genes <- read_tsv("data/genes.txt", col_types = "c-c---i-----") |>
    rename(gene_chrom = chrom)

d_lm <- read_tsv("data/gemma/qq/Acbc.lm.assoc.txt.gz", col_types = "c-c--------d")

d_lmm <- read_tsv("data/gemma/qq/Acbc.lmm.assoc.txt.gz", col_types = "c-c----------d--")

d <- bind_rows(
    d_lm |> mutate(Method = "LM"),
    d_lmm |> mutate(Method = "LMM")
) |>
    separate_wider_delim(rs, ":", names = c("chrom", "pos")) |>
    mutate(chrom = str_replace(chrom, "chr", ""),
           pos = as.integer(pos)) |>
    left_join(genes, by = "gene_id", relationship = "many-to-one") |>
    mutate(type = case_when(
        gene_chrom == chrom & abs(pos - tss) < 1e6 ~ "cis",
        gene_chrom == chrom & abs(pos - tss) < 5e6 ~ "mid",
        TRUE ~ "trans"
    )) |>
    filter(type != "mid")

ggplot(d, aes(x = p_wald)) +
    facet_grid(rows = vars(type), cols = vars(Method), scales = "free") +
    geom_histogram(bins = 100, boundary = 0)

tmp <- d |>
    mutate(type = fct_recode(type,
                             "TSS distance < 1 Mb" = "cis",
                             "TSS distance > 5 Mb" = "trans")) |>
    reframe(
        p = p_wald,
        p_ref = ecdf(p)(p),
        .by = c(Method, type)
    )

ggplot(tmp, aes(x = -log10(p_ref), y = -log10(p), color = Method)) +
    facet_wrap(~ type) +
    geom_abline(slope = 1, intercept = 0, linewidth = 0.5, lty = 2) +
    geom_point(size = 0.5) +
    geom_text(data = tibble(p_ref = 10^-3.5, p = 10^-2.8, Method = NA,
                            type = "TSS distance < 1 Mb"),
              label = "y = x", color = "black", show.legend = FALSE) +
    scale_color_manual(values = c("#ff5555", "#5555ff")) +
    guides(color = guide_legend(override.aes = list(size = 2))) +
    xlab(expression("Expected "*-log[10](P))) +
    ylab(expression("Observed "*-log[10](P))) +
    theme_bw() +
    theme(
        panel.grid = element_blank(),
        strip.background = element_blank(),
        legend.position = c(0.65, 0.7),
        legend.background = element_rect(fill = "white", color = "#CCCCCC"),
    )

ggsave("figures/figureS3/figureS3d.png", width = 5, height = 3, bg = "white")

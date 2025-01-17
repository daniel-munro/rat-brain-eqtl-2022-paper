library(tidyverse)
library(patchwork)

top_lm <- read_tsv("data/gemma/NAcc.lm.assoc.txt.gz",
                   col_types = "c-c------ddd") |>
    group_by(gene_id) |>
    sample_n(1) |>
    ungroup() |>
    mutate(z = beta / se)

top_lmm <- read_tsv("data/gemma/NAcc.lmm.assoc.txt.gz",
                    col_types = "c-c-----dd---d--") |>
    group_by(gene_id) |>
    sample_n(1) |>
    ungroup() |>
    mutate(z = beta / se)

pve <- read_tsv("data/gemma/NAcc.grm_pve.pve.txt", col_types = "cd-") |>
    mutate(pve = pmax(0, pmin(pve, 1)))

top <- full_join(
    top_lm |> select(gene_id, p_lm = p_wald, z_lm = z),
    top_lmm |> select(gene_id, p_lmm = p_wald, z_lmm = z),
    by = "gene_id"
) |>
    left_join(pve, by = "gene_id")

########################################
## Scatter plot of LMM vs LM z-scores ##
########################################

p1_r <- with(top, cor(abs(z_lm), abs(z_lmm), method = "spearman")) |> signif(4)
top |>
    slice_sample(prop = 1L, replace = FALSE) |>
    ggplot(aes(x = abs(z_lm), y = abs(z_lmm), color = pve)) +
    geom_abline(slope = 1, intercept = 0, color = "#cccccc") +
    annotate("text", x = 30, y = 29, hjust = 0, label = "x = y") +
    geom_point(size = 0.25) +
    coord_fixed() +
    scale_color_viridis_c(begin = 0, end = 0.8) +
    annotate("text", x = 2, y = 33, hjust = 0,
             label = str_glue("rho = {p1_r}\nN = {nrow(top)}")) +
    theme_bw() +
    theme(
        panel.grid = element_blank(),
        legend.position = c(0.85, 0.3),
        legend.box.background = element_rect(color = "black"),
    ) +
    labs(color = expression(PVE[GRM])) +
    xlab("|Z-score|, fixed effects model") +
    ylab("|Z-score|, linear mixed model")

ggsave("figures/figureS3/figureS3b.png", width = 3.5, height = 3.5, device = png)

###########################################
## eGene overlap at different thresholds ##
###########################################

overlap <- tibble(threshold = 10 ^ seq(from = -10, to = -2, length.out = 500)) |>
    group_by(threshold) |>
    summarise(`LM only` = sum(top$p_lm < threshold & top$p_lmm >= threshold),
              `LMM only` = sum(top$p_lm >= threshold & top$p_lmm < threshold),
              both = sum(top$p_lm < threshold & top$p_lmm < threshold),
              .groups = "drop") |>
    pivot_longer(-threshold, names_to = "group", values_to = "eGenes") |>
    mutate(group = fct_relevel(group, "LM only", "both", "LMM only"))

overlap |>
    ggplot(aes(x = threshold, y = eGenes, fill = group)) +
    geom_col(position = "fill", width = 0.02) +
    scale_x_log10(expand = c(0, 0)) +
    scale_y_continuous(breaks = c(0, 0.5, 1), expand = c(0, 0)) +
    scale_fill_manual(values = c("#ff5555", "#aa55aa", "#5555ff")) +
    theme_bw() +
    theme(panel.grid = element_blank(),
          panel.border = element_blank()) +
    xlab("P-value threshold") +
    ylab("Proportion of eGenes") +
    labs(fill = "Below\nthreshold in")

ggsave("figures/figureS3/figureS3c.png", width = 6, height = 2)

#############################
## Stats related to figure ##
#############################

top |>
    summarise(Spearman_p = cor(abs(z_lm), abs(z_lmm), method = "spearman"))

tibble(threshold = c(1e-9, 1e-6, 1e-3)) |>
    group_by(threshold) |>
    summarise(both = sum(top$p_lm < threshold & top$p_lmm < threshold),
              total = sum(top$p_lm < threshold | top$p_lmm < threshold),
              .groups = "drop") |>
    mutate(fraction = both / total)

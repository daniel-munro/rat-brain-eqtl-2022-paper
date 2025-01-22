library(tidyverse)
library(patchwork)

top_assoc <- read_tsv("data/eqtls/top_assoc.txt", col_types = "ccciiciiccdddddddd")

egenes <- filter(top_assoc, qval < 0.05)

eqtls <- read_tsv("data/eqtls/eqtls_indep.txt", col_types = "ccciiciiccdddddid")

# Num eGenes per tissue
p1 <- egenes |>
    mutate(tissue = fct_infreq(tissue) |> fct_rev()) |> # to match upset plot order
    ggplot(aes(x = tissue)) +
    geom_bar() +
    scale_y_continuous(expand = c(0, 0), limits = c(0, max(table(egenes$tissue)) * 1.02)) +
    xlab(NULL) +
    ylab("eGenes") +
    theme_bw() +
    theme(axis.text.x = element_text(hjust = 1, vjust = 0.5, angle = 90),
          panel.grid = element_blank())

# Num eGenes found in N tissues
tmp <- egenes |>
    summarise(n_tissues = n(), .by = gene_id) |>
    count(n_tissues)
p2 <- ggplot(tmp, aes(x = n_tissues, y = n, fill = n_tissues)) +
    geom_col(show.legend = FALSE) +
    scale_y_continuous(expand = c(0, 0), limits = c(0, max(tmp$n) * 1.02)) +
    scale_fill_viridis_c() +
    xlab("No. tissues") +
    ylab("eGenes") +
    theme_bw() +
    theme(panel.grid = element_blank())

######################################
## Combine with sQTL plots to align ##
######################################

sqtls <- read_tsv("data/splice/sqtls_indep.txt", col_types = "cciccicdddddcii")

sgenes <- sqtls |>
    distinct(tissue, group_id) |>
    rename(gene_id = group_id)

# Num sGenes per tissue
p3 <- sgenes |>
    mutate(tissue = fct_infreq(tissue) |> fct_rev()) |> # to match upset plot order
    ggplot(aes(x = tissue)) +
    geom_bar() +
    scale_y_continuous(expand = c(0, 0), limits = c(0, max(table(sgenes$tissue)) * 1.02)) +
    xlab(NULL) +
    ylab("sGenes") +
    theme_bw() +
    theme(axis.text.x = element_text(hjust = 1, vjust = 0.5, angle = 90),
          panel.grid = element_blank())

# Num sGenes found in N tissues
tmp2 <- sgenes |>
    summarise(n_tissues = n(), .by = gene_id) |>
    count(n_tissues)
p4 <- ggplot(tmp2, aes(x = n_tissues, y = n, fill = n_tissues)) +
    geom_col(show.legend = FALSE) +
    scale_y_continuous(expand = c(0, 0), limits = c(0, max(tmp2$n) * 1.02)) +
    scale_fill_viridis_c() +
    xlab("No. tissues") +
    ylab("sGenes") +
    theme_bw() +
    theme(panel.grid = element_blank())

p1 + p2 + p3 + p4 +
    theme(plot.margin = margin(0.7, 0.7, 0.7, 0.7, unit = "cm")) +
    plot_layout(ncol = 2)
ggsave("figures/figure4/figure4abde.png", width = 3.5, height = 6.2)

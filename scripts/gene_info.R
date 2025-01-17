library(tidyverse)

genes <- read_tsv("data/Rnor_6.0_anno/Rattus_norvegicus.Rnor_6.0.99.genes.bed",
                  col_types = "ciic-c---c",
                  col_names = c("chrom", "start", "end", "gene_id", "strand", "etc")) %>%
    mutate(gene_name = str_match(etc, 'gene_name "([^"]+)"')[, 2],
           tss = if_else(strand == "+", start, end)) %>%
    select(gene_id, gene_name, chrom, start, end, strand, tss)

# Record whether gene is included in each expression table:
for (tissue in c("IL", "LHb", "NAcc", "OFC", "PL")) {
    expr <- read_tsv(str_glue("data/expression/ensembl-gene_log2_{tissue}.bed.gz"),
                     col_types = cols(gene_id = "c", .default = "-"))
    genes[[str_c("in_expr_", tissue)]] <- genes$gene_id %in% expr$gene_id
}

write_tsv(genes, "data/genes.txt")

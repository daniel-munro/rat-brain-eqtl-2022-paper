library(VariantAnnotation)
library(tidyverse)

################
## SNP counts ##
################

# ac <- readInfo("data/genotype/P50.rnaseq.88.unpruned.vcf.gz", "AC")
d <- info(readVcf("data/genotype/P50.rnaseq.88.unpruned.vcf.gz")) |>
    as_tibble() |>
    mutate(AC = as.integer(AC),
           MAF = pmin(AC / AN, 1 - AC / AN))

summary(d$MAF)
nrow(d)

read_lines("data/genotype/imputing/observed.snplist.txt") |>
    length()

#################
## Read counts ##
#################

meta <- read_csv("data/qc/metadata/metadata_p50_hao_chen_2014.csv",
                 col_types = "ccccDDcDddDccciiddciiddcccccc")

sum(meta$QC_pass == "pass")

meta |>
    summarise(raw_reads_mean = mean(raw_reads),
              raw_reads_sd = sd(raw_reads))

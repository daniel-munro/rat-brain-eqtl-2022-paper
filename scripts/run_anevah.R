# Based on https://github.com/PejLab/anevah/blob/master/vignettes/anevah.Rmd
# I can no longer install the dependency PejLab/bln due to
# `RcppExports.cpp:4:10: fatal error: 'Rcpp.h' file not found`, so now I'm loading
# the functions.

library(tidyverse)
library(pracma) # For eps()

source("tools/bln/R/bln.R")
source("tools/anevah/R/anevah.R")
source("tools/anevah/R/bln.R")
source("tools/anevah/R/Vg.R")

run_anevah <- function(tissue) {
    ref_counts <- read.delim(
        str_glue("data/anevah/input/anevah_{tissue}_ref_counts.tsv"),
        header = TRUE,
        row.names = 1,
        as.is = TRUE, #nrow = 100
    )
    alt_counts <- read.delim(
        str_glue("data/anevah/input/anevah_{tissue}_alt_counts.tsv"),
        header = TRUE,
        row.names = 1,
        as.is = TRUE, #nrow = 100
    )
    anevah(ref_counts, alt_counts)
}

for (tissue in c("IL", "LHb", "NAcc", "OFC", "PL")) {
    out <- run_anevah(tissue)
    write_tsv(out, str_glue("data/anevah/output/Vg.{tissue}.tsv.gz"))
}

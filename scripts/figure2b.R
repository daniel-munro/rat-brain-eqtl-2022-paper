library(tidyverse)

########################
## Run the simulation ##
########################

next_generation <- function(genos, mating = "circular") {
    if (mating == "circular") {
        x <- genos |>
            mutate(father = c(male[2:length(male)], male[1]))
    } else {
        x <- genos |>
            mutate(father = male[sample(length(male), replace = FALSE)])
    }
    x |>
        select(cage, mother = female, father) |>
        rowwise() |>
        mutate(female = list(c(sample(mother, 1), sample(father, 1))),
               male = list(c(sample(mother, 1), sample(father, 1)))) |>
        ungroup() |>
        select(cage, female, male)
}

geno_sim <- function(pairs = 50, generations = 90, mating = "circular") {
    genos1 <- tibble(cage = 1:pairs) |>
        group_by(cage) |>
        mutate(female = list(sample(LETTERS[1:8], 2, replace = TRUE)),
               male = list(sample(LETTERS[1:8], 2, replace = TRUE))) |>
        ungroup()
    genos <- list(genos1)
    for (gen in 2:generations) {
        genos <- c(genos, list(next_generation(genos[[length(genos)]], mating)))
    }
    bind_rows(genos, .id = "generation") |>
        mutate(generation = as.integer(generation))
}

genos <- tibble(permutation = 1:200) |>
    reframe(geno_sim(generations = 90), .by = permutation)

counts <- genos |>
    reframe(
        tibble(strain = c(unlist(female), unlist(male))) |>
            count(strain),
        .by = c(permutation, generation)
    )

write_tsv(counts, "data/haplotypes/breeding_sim_circ.tsv")

genos2 <- tibble(permutation = 1:200) |>
    reframe(geno_sim(generations = 90, mating = "random"), .by = permutation)

counts2 <- genos2 |>
    reframe(
        tibble(strain = c(unlist(female), unlist(male))) |>
            count(strain),
        .by = c(permutation, generation)
    )

write_tsv(counts2, "data/haplotypes/breeding_sim_rand.tsv")

#######################################
## Load real haplotype probabilities ##
#######################################

probs <- function(prob, n_SNPs = NULL) {
    names(dimnames(prob)) <- c("individual", "strain", "SNP")
    if (!is.null(n_SNPs)) {
        SNP_reps <- round(seq(from = 1, to = dim(prob)[3], length.out = n_SNPs))
        prob <- prob[, , SNP_reps]
    }
    cubelyr::as.tbl_cube(prob) |>
        as_tibble()
}

chr_len <- c(282763074, 266435125, 177699992, 184226339, 173707219,
             147991367, 145729302, 133307652, 122095297, 112626471,
             90463843, 52716770, 114033958, 115493446, 111246239,
             90668790, 90843779, 88201929, 62275575, 56205956)
n_SNPs <- round(2000 * chr_len / chr_len[1])

d <- tibble(chrom = 1:20) |>
    reframe(
        readRDS(str_glue("data/haplotypes/haplotype_probs_chr{chrom}.rds")) |>
            probs(n_SNPs = n_SNPs[chrom]) |>
            summarise(prob = mean(prob), .by = c(SNP, strain)),
        .by = chrom
    ) |>
    mutate(SNP = as.integer(as.factor(SNP))) |>
    filter(SNP %in% sample(unique(SNP), 200, replace = FALSE))

#####################
## Compare entropy ##
#####################

sim_circ <- read_tsv("data/haplotypes/breeding_sim_circ.tsv", col_types = "iici")
sim_rand <- read_tsv("data/haplotypes/breeding_sim_rand.tsv", col_types = "iici")

counts <- bind_rows(
    d |>
        mutate(data = "Observed proportions",
               generation = 80) |>
        select(-chrom),
    sim_circ |>
        mutate(data = "Simulated - circular mating") |>
        rename(SNP = permutation) |>
        mutate(prob = n / sum(n), .by = c(generation, SNP)) |>
        select(-n),
    sim_rand |>
        mutate(data = "Simulated - random mating") |>
        rename(SNP = permutation) |>
        mutate(prob = n / sum(n), .by = c(generation, SNP)) |>
        select(-n)
)

entropies <- counts |>
    summarise(entropy = entropy::entropy.empirical(prob, unit = "log2"),
              .by = c(data, SNP, generation))

entropies |>
    filter(generation <= 80,
           generation %% 10 == 0) |>
    mutate(generation = as.factor(generation),
           data = fct_relevel(data,
                              "Simulated - circular mating",
                              "Simulated - random mating")) |>
    ggplot(aes(x = generation, y = entropy, fill = data)) +
    geom_boxplot(outlier.size = 0.5) +
    expand_limits(y = 0) +
    scale_fill_manual(values = c("#37a9bf", "#376dbf", "#b03884")) +
    theme_bw() +
    theme(
        legend.position = "inside",
        legend.position.inside = c(0.21, 0.25),
        legend.title = element_blank(),
        panel.grid = element_blank(),
    ) +
    xlab("Generation") +
    ylab("Shannon entropy (base 2)")

ggsave("figures/figure2/figure2b.png", width = 5.5, height = 3)

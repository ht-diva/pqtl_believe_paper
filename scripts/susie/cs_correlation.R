

# November 11, 2026
# Dariush Ghasemi

# This script will identify pairwise LD (r2) between
# credible sets for 2644 pQTLs with seconday signal.


# Loci with >1 CS to explore LD correlation
cs_secondary <- cs_believe_sum %>% filter(nset > 1)


find_ld <- function(seqid, locus) {
  
  my_rds <- glue(
    path_susie, "results/believe_bgen/susierss/cs_fitness/",
    seqid, "_", gsub("chr", "", locus), "_annot.rds"
  ) %>%
    readRDS()
  
  my_rds$sld %>%
    as.data.frame() %>%
    tibble::rownames_to_column("cs_id1") %>%
    tidyr::pivot_longer(
      cols = -cs_id1,
      names_to = "cs_id2",
      values_to = "LD"
    ) %>%
    mutate(
      seqid = seqid,
      locus = locus
    ) %>%
    select(seqid, locus, cs_id1, cs_id2, LD)
  
}


cs_ld_secondary <- map2_df(
  cs_secondary$seqid,
  cs_secondary$locus,
  find_ld
  )



cs_ld_secondary %>%
  dplyr::filter(
    cs_id1 != cs_id2,
    LD >= 0.5
    ) %>% 
  distinct(seqid, locus)

# Total loci with >1 CS: 2644 of 8,924
# Number of loci with correlated sets
#   - LD > 0.5: 81
#   - LD > 0.6: 40
#   - LD > 0.7: 23
#   - LD > 0.8: 13
#   - LD > 0.9: 13

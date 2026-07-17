
# ======================================== # 
# Take inputs from susie_believe.R script.
# ======================================== # 

# distribution across chroms
cs_believe %>%
  filter(cs_id == "no_credible") %>%
  count(chr)

# pQTLs without CS
pqtls_nocs <- cs_believe %>%
  filter(cs_id == "no_credible") %>%
  right_join(lb_believe_annot[, c(1, 4:13, 43:45)], join_by(seqid, locus))


#-----------------#
# Loci without CS
pqtls_nocs %>%
  mutate(
    MAF = ifelse(EAF < 0.5, EAF, 1 - EAF),
    MAF_cat = case_when(MAF > 0.01 ~ "MAF > 0.01", TRUE ~ "MAF <= 0.01")
  ) %>%
  ggplot(aes(x = chr, y = MLOG10P, size = loci_width)) + #shape = MAF_cat
  geom_jitter(fill = "orange", color = "steelblue", shape = 21) + #, size = 3
  scale_x_continuous(breaks  = c(1:12)) +
  theme_light() + labs(x = "Chromosomal position")

pqtls_nocs %>%
  ggplot(aes(x = BETA, y = MLOG10P,
             xmin = BETA - SE, xmax = BETA + SE)) + 
  geom_pointrange(fill = "orange", color = "steelblue", shape = 21) +
  geom_vline(xintercept = 0, linetype = 2) +
  theme_light() 


pqtls_nocs %>%
  ggplot(aes(x = factor(loci_cat, levels = level_order))) +
  geom_histogram(fill = "steelblue", color = "orange", stat = "count") +
  scale_y_continuous(breaks  = c(1:7)) +
  theme_light() + labs(x = "Locus size", y = "#Loci")

pqtls_nocs %>%
  mutate(cs_id = ifelse(is.na(cs_id), "with_credible", cs_id)) %>%
  ggplot(aes(x = cs_id, y = MLOG10P)) +
  geom_violin(fill = "steelblue", trim = FALSE) +
  theme_light() + labs(x = "Locus")

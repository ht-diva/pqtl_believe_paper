

path_susie <- "/scratch/dariush.ghasemi/projects/pqtl_susie/"
path_cs_believe <- glue(path_susie, "results/believe_nomir/susierss/cs_list")

path_pips <- glue(path_susie, "results/believe_test/susierss/cs_summary")
path_pips_combined <- glue(path_susie, "results/believe_test/susierss/combined_cssums.tsv")

path_pip_cslist <- glue(path_susie, "results/believe_estVarF/susierss/combined_credibles.tsv")


#----------------------------#
# -----      CS List    -----
#----------------------------#

level_order <- c('1-SNP', '1bp-100Kbp', '100-500Kbp', '500K-1Mbp', '1-2Mbp', '2-5Mbp', '>5Mbp')


# Credible sets
cs_believe <- fread(path_pip_cslist)


# Files with credible sets
cs_lists_believe <- list.files(
  path = path_cs_believe,
  pattern = ".cslist",
  full.names = T)

cslist_class <- c(rep("character",3), rep("numeric",4), "character")

# combine all cs files
cs_believe <- map_dfr(cs_lists_believe, 
                      ~ fread(.x, colClasses = cslist_class))


#----------------------------#
#----     CS Sumstats    ----
#----------------------------#

# Sumstats of credible sets variants along with PIP
list_cssums <- list.files(
  path = path_pips,
  pattern = "seq.(\\d)+.(\\d)+_(\\d)+_(\\d)+_(\\d)+.cssum$",
  full.names = T
)

combined_cssums <- map_dfr(list_cssums, function(path) {
  filename <- basename(path) %>% gsub('.cssum', '', .)
  fread(path) %>% 
    mutate(seqid_locus = filename) %>% 
    relocate(seqid_locus)
}
)

# Save combined sumstats with PIP for Giulia Pontali
fwrite(combined_cssums, file = path_pips_combined, sep = "\t")


#----------------------------#
# ----   Missing jobs    ----
#----------------------------#

cs_believe %>% distinct(seqid, locus) %>% dim()

# Believe study loci whose job got this error:
# Error: Estimating residual variance failed: the estimated value is negative
lb_believe_annot %>%
  dplyr::mutate(locus = str_c("chr", chr, "_", start, "_", end)) %>%
  #select(seqid, locus) %>%
  left_join(
    cs_believe %>% summarize(nset = n(), .by = c(seqid, locus)),
    join_by(seqid, locus)
    ) %>%
  dplyr::filter(is.na(nset)) %>%
  #select(- locus, - nset) %>%
  #write.csv("/scratch/dariush.ghasemi/projects/pqtl_susie/config/loci/believe_negative_res_var.csv", quote = F, row.names = F)
  count(cis_or_trans, loci_cat) %>% 
  tidyr::spread(value = n, "loci_cat") %>%
  DT::datatable()

#----------------------------#
# -----      Counts     -----
#----------------------------#

# CS summary for modelling
cs_believe_sum <- cs_believe %>%
  filter(cs_id != "no_credible") %>%
  select(- cs_snps) %>%
  summarize(
    n_cs   = n_distinct(cs_id),
    n_snps = sum(ncs),
    cs_power = mean(cs_log10bf),
    cs_avgr2 = mean(cs_avg_r2),
    cs_minr2 = mean(cs_min_r2),
    .by = c(seqid, locus)
    )

### For exploration

# number of pQTLs without CS
cs_believe %>%
  filter(cs_id == "no_credible") %>%
  distinct(seqid, locus) %>%
  # cis/trans annotation
  left_join(
    lb_believe_annot %>% select(seqid, locus, cis_or_trans),
    join_by(seqid, locus)
  ) %>%
  count(cis_or_trans)

# number of pQTLs (loci) with at least one credible set (CS)
cs_believe %>%
  filter(cs_id != "no_credible") %>%
  distinct(seqid, locus)

# number of CS per locus
cs_believe %>%
  filter(cs_id != "no_credible") %>%
  summarize(
    ncs = n(),
    .by = c(seqid, locus)
  ) %>%
  count(ncs, name = "n_pQTLs")

# pQTLs with >1 CS
cs_believe %>%
  filter(cs_id != "no_credible") %>% #head(200) %>%
  select(- cs_snps) %>%
  mutate(n_cs = n_distinct(cs_id), .by = c(seqid, locus)) %>%
  mutate(secondary = n_cs > 1,
         onesnp_cs = ncs == 1) %>%
  #filter(onesnp_cs) %>% distinct(seqid, locus) # 2470 loci with 1-SNP CS
  distinct(seqid, locus, secondary, onesnp_cs) %>% # 10,319 rows
  left_join(
    lb_believe_annot %>% select(seqid, locus, cis_or_trans),
    join_by(seqid, locus)
  ) %>%
  #count(cis_or_trans, secondary) # CS=1: 362 cis, 5513 trans; CS>1: 1143 cis, 1875 trans
  count(cis_or_trans, secondary, onesnp_cs) %>%
  spread(secondary, n)


# append number of CS per locus
cs_count <- cs_believe %>%
  filter(cs_id != "no_credible") %>%
  mutate(n_credsets = n_distinct(cs_id), .by = c(seqid, locus)) %>%
  distinct(seqid, locus, n_credsets)

cs_count %>%
  left_join(
    lb_believe_annot %>% select(seqid, locus, cis_or_trans),
    join_by(seqid, locus)
    ) %>%
  #count(cis_or_trans)
  summarize(sum(n_credsets), .by = cis_or_trans)


#----------------------------#
# -----      Plots      -----
#----------------------------#

# histogram of no. loci per CS size
cs_count %>%
  distinct(seqid, locus, n_credsets) %>% # keep number of CS per locus
  ggplot(aes(x = n_credsets)) +
  geom_histogram(
    stat = "count", fill = "#b9936c",
    position = position_dodge(0.9, preserve = 'single')
  ) +
  stat_count(
    aes(label=..count.., y=..count.. + 170),
    geom = 'text', color = 'gold', #'#bd5734'
    position = position_dodge(.9),
  ) +
  scale_x_continuous(
    breaks = c(1:10), expand = c(0,0) # remove space between plot and y-axis
    ) +
  scale_y_continuous(
    breaks = c(seq(0, 6000, 500)), 
    limits = c(0, 6400), expand = c(0,0) # remove space between plot and x-axis
  ) +
  labs(
    x = "\nNumber of credible sets per pQTL (excluding 11 without CS)",
    y = "Number of loci\n"
  ) +
  theme_classic() +
  theme(
    axis.title = element_text(size = 12, face = 2, color = "steelblue2"),
    axis.text  = element_text(size = 11, color = "skyblue2"),
    axis.ticks.length = unit(2, "mm"),
    #panel.grid.minor = element_line(color = "white"),
    #panel.grid.major = element_line(color = "white")
    plot.background  = element_rect(fill = "transparent", color = NA),
    panel.background = element_rect(fill = "transparent", color = NA)
    )


ggsave("03-Jun-26_believe_susie_cs_count.png", bg = "transparent",
       plot = last_plot(), width = 8, height = 6, dpi = 300)



#----------------------------#
# -----   Diagnostics   -----
#----------------------------#

cs_believe %>%
  ggplot(aes(y = cs_log10bf)) +
  geom_boxplot() +
  theme_minimal()

ggsave("22-Jan-26_believe_susie_cs_log10bf.jpg",
       plot = last_plot(), width = 6, height = 8, dpi = 300)


cs_believe %>%
  mutate(
    ncs = ifelse(cs_id == "no_credible", 0, ncs)
    ) %>%
  ggplot(aes(y = ncs)) +
  geom_boxplot() +
  # geom_jitter(
  #   stat = "count", #bins = 18,
  #   position = position_dodge(0.9, preserve = 'single')
  # ) +
  # stat_count(
  #   aes(label=..count.., y=..count.. + 130), 
  #   geom = 'text', color = '#bd5734',
  #   position = position_dodge(.9),
  # ) +
  labs(#x = "No. of Loci",
       y = "No. of SNPs in Credible Sets") +
  theme_minimal()



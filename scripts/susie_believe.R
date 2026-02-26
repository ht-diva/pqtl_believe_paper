

path_susie <- "/scratch/dariush.ghasemi/projects/pqtl_susie/"
path_cs_believe <- glue(path_susie, "results/believe/susierss/cs_list")
path_pips <- glue(path_susie, "results/believe_test/susierss/cs_summary")
path_pips_combined <- glue(path_susie, "results/believe_test/susierss/combined_cssums.tsv")

#----------------------------#
# -----   SuSiE Output  -----
#----------------------------#

# files with credible sets
cs_lists_believe <- list.files(path = path_cs_believe, pattern = ".cslist", full.names = T)

# combine all cs files
cs_believe <- map(cs_lists_believe, function(path) fread(path, colClasses = c(
  rep("character",3), rep("numeric",4), "character"))) %>%
  bind_rows()

#----------------------------#
# -----      Counts     -----
#----------------------------#

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


# append number of CS per locus
cs_believe_count <- cs_believe %>%
  filter(cs_id != "no_credible") %>%
  group_by(seqid, locus) %>%
  mutate(n_credsets = n_distinct(cs_id)) %>%
  ungroup()



#----------------------------#
# -----      Plots      -----
#----------------------------#

# histogram of no. loci per CS size
cs_believe_count %>%
  distinct(seqid, locus, n_credsets) %>% # keep number of CS per locus
  ggplot(aes(x = n_credsets)) +
  geom_histogram(
    stat = "count", fill = "#b9936c",
    position = position_dodge(0.9, preserve = 'single')
  ) +
  stat_count(
    aes(label=..count.., y=..count.. + 170),
    geom = 'text', color = '#bd5734',
    position = position_dodge(.9),
  ) +
  scale_x_continuous(
    breaks = c(1:10), expand = c(0,0) # remove space between plot and y-axis
    ) +
  scale_y_continuous(
    breaks = c(seq(0, 6500, 500)), 
    limits = c(0, 6900), expand = c(0,0) # remove space between plot and x-axis
  ) +
  labs(
    x = "\nNumber of credible sets (excluding 21 without CS)",
    y = "Number of loci\n"
  ) +
  theme_light() +
  theme(
    axis.title = element_text(size = 12, face = 2),
    axis.text  = element_text(size = 11),
    axis.ticks.length = unit(2, "mm"),
    panel.grid.minor = element_line(color = "white"),
    panel.grid.major = element_line(color = "white")
  )


ggsave("22-Jan-26_believe_susie_cs_count.jpg",
       plot = last_plot(), width = 8, height = 6, dpi = 150)


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


#----------#
# loci with missing susie fine-mapping
list_locuseq <- list.files(
  path = path_cslist,
  pattern = "seq.(\\d)+.(\\d)+_(\\d)+_(\\d)+_(\\d)+.cslist$",
  recursive = F, full.names = F
)

cslist_done <- tibble(list_locuseq) %>%
  transmute(
    locuseq = str_remove(list_locuseq, ".cslist"),
    seqid = str_split_fixed(locuseq, "_", 2)[,1],
    locus = str_split_fixed(locuseq, "_", 2)[,2]
  )

#----------#
ld_matrices <- list.files(
  path = glue(path_susie, "results/believe/tmp/"),
  pattern = "seq.(\\d)+.(\\d)+_(\\d)+_(\\d)+_(\\d)+_ld.matrix$",
  full.names = T
) %>%
  tibble() %>%
  transmute(
    ld_size_gb = file.size(.)/(1e9),
    memory  = (ld_size_gb * 4.5 + 4)*1024,
    locuseq = basename(.) %>% str_remove("_ld.matrix")
    )

write.csv(ld_matrices, "ld_matrix_size_required_mem.csv")


#----------#
list_report <- list.files(
  path = glue(path_susie, "results/believe/susierss/cs_report"),
  pattern = "seq.(\\d)+.(\\d)+_(\\d)+_(\\d)+_(\\d)+.report$",
  full.names = T
)

# find other regions with literally similar SNPs size
cs_report <- map_dfr(list_report[1:1000], fread)
cs_report %>% filter(
  #between(n_snp_pgen, 9200, 9300),
  #between(n_snp_pgen, 14000, 15000),
  between(n_snp_pgen, 50000, 56000),
  )


plot(cs_report$n_snp_pgen, cs_report$n_snp_gwas)
hist(cs_report$n_snp_gwas, nclass = 30)


list_report <- list.files(
  path = glue(path_susie, "results/believe/susierss/cs_report"),
  pattern = "seq.(\\d)+.(\\d)+_(\\d)+_(\\d)+_(\\d)+.report$",
  full.names = T
)

# 9320 SNPs
ld_matrices %>% filter(locuseq %in% locuseq2test)
# A tibble: 2 × 3
# ld_size_gb memory locuseq                         
# <dbl>  <dbl> <chr>                           
# 1      2.23  15527. seq.22834.21_1_196246962_198992…
# 2      0.942  8917. seq.4440.15_1_156706526_1580416…  -- failed
# 3                   seq.10754.113_6_27895071_32982423 (55,112snps) -- done
# 4      29G  147Gb  seq.13689.2_6_28042655_32982423 (54,359snps) -- failed


# [pqtl_susie]dariush.ghasemi@cnode33$ seff 26846416
# Job ID: 26846416
# Cluster: slurm
# User/Group: dariush.ghasemi/diangelantonio
# State: COMPLETED (exit code 0)
# Cores: 1
# CPU Utilized: 08:42:24
# CPU Efficiency: 98.34% of 08:51:14 core-walltime
# Job Wall-clock time: 08:51:14
# Memory Utilized: 168.63 GB
# Memory Efficiency: 74.07% of 227.68 GB





ld_nsnps <- ld_matrices %>%
  filter(locuseq != "seq.4440.15_1_156706526_158041681") %>%
  mutate(
    snpfile = paste0(path_susie, "results/believe/tmp/", locuseq, "_ld.headers")
    ) %>%
  rowwise() %>%
  mutate(nsnps = fread(snpfile) %>% nrow()) %>%
  select(-snpfile)

plot(ld_nsnps$ld_size_gb, ld_nsnps$nsnps)
  

fread(glue(path_susie, "cs_list.believe")) %>%
  transmute(locuseq = basename(cs_list) %>% str_remove(".cslist")) %>%
  inner_join(ld_matrices, join_by(locuseq)) %>% View()

locuseq2test <- c(
  "seq.4440.15_1_156706526_158041681",
  "seq.22834.21_1_196246962_198992432"
)

# loci to rerun
lb_believe %>%
  mutate(locuseq = paste0(seqid, "_", chr, "_", start, "_", end)) %>%
  filter(!locuseq %in% cslist_done$locuseq) %>%
  filter(loci_cat != ">5Mbp") %>%
  select(-locuseq) %>% #count(loci_cat)
  write.csv(glue(path_susie, "config/believe_loci_failed.csv"), row.names = F, quote = F)



#----------------------------#
#----  Sumstat with PIP  ----
#----------------------------#

list_cssums <- list.files(
  path = path_pips,
  pattern = "seq.(\\d)+.(\\d)+_(\\d)+_(\\d)+_(\\d)+.cssum$",
  full.names = T
)

combined_cssums <- map_dfr(list_cssums, fread)

# save combined sumstats with PIP for Giulia Pontali
fwrite(combined_cssums, file = path_pips_combined, sep = "\t")



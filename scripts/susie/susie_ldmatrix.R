


ld_matrices <- list.files(
  path = glue(path_susie, "results/believe/tmp/"),
  pattern = "seq.(\\d)+.(\\d)+_(\\d)+_(\\d)+_(\\d)+_ld.matrix$",
  full.names = T
) 



#----------#
ld_matrices %>%
  tibble() %>%
  transmute(
    ld_size_gb = file.size(.)/(1e9),
    memory  = (ld_size_gb * 4.5 + 4)*1024,
    locuseq = basename(.) %>% str_remove("_ld.matrix")
  )

write.csv(ld_matrices, "ld_matrix_size_required_mem.csv")

#----------#

# find other regions with literally similar SNPs size
cs_report %>% filter(
  #between(n_snp_pgen, 9200, 9300),
  #between(n_snp_pgen, 14000, 15000),
  between(n_snp_pgen, 50000, 56000),
)


plot(cs_report$n_snp_pgen, cs_report$n_snp_gwas)
hist(cs_report$n_snp_gwas, nclass = 30)


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



#----------#

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





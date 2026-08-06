
library(tidyverse)

# Inputs
path_gnh_quant <- "/scratch/dariush.ghasemi/projects/pqtl_coloc/results/believe_cis_quant/tmp/gwas"
path_gnh_binar <- "/scratch/dariush.ghasemi/projects/pqtl_coloc/results/believe_cis_binary/tmp/gwas"

# test GWAS
path_gnh_ldlc <- "/scratch/dariush.ghasemi/projects/pqtl_coloc/tiledb/results/APOE_19_44819866_45114528/_genesandhealth_v010_quantitative_traits_median_values_07d4a34298.csv.gz"

test_gwas <- fread(path_gnh_ldlc)

which.max(test_gwas$MLOG10P)
test_gwas[which.max(test_gwas$MLOG10P), ]


cis_hits_quant <- fread("G&H_hits_cis_pqtls_quant.csv")

cis_hits_quant %>%
  filter(MLOG10P > 10) %>%
  mutate(locuseq = dirname(path_sumstat) %>% str_remove(".*gwas/")) %>%
  distinct(locuseq, meta_trait_desc, pass_thr, path_sumstat) %>%
  distinct(locuseq)

# List GWASs
gwas_list <- list.files(
  #path = path_gnh_quant, # 161,142 quant GWASs for 1506 cis pQTLs
  path = path_gnh_binar,  # 796,674 binary GWASs
  pattern = ".csv.gz",
  recursive = TRUE,
  full.names = TRUE
)

# APOE
# 19_44819866_45114528

find_hit <- function(gwas_file, threshold = 4){
  
  gwas  <- fread(gwas_file)
  
  # Flag if any variant passes threshold
  gwas$pass_thr <- any(gwas$MLOG10P > threshold)
  
  # Add file path
  gwas$path_sumstat <- gwas_file
  
  # Maximum value
  max_log10p <- max(gwas$MLOG10P, na.rm = TRUE)
  
  # Accounts for multiple rows
  gwas[gwas$MLOG10P == max_log10p, ]
  
}


map_dfr(gwas_list[c(1,2, 108,109, 215,216)], find_hit)

# List of hits for each quant trait in G&H GWAS at cis loci
cis_hits_quant <- map_dfr(gwas_list, find_hit)
cis_hits_binar <- map_dfr(gwas_list, find_hit)

# 10:45 quant started
# 19:00 quant ended
# 00:20 binary started
# failed OOM error

write.csv(cis_hits_quant, "G&H_hits_cis_pqtls_quant.csv", quote = T, row.names = F)

cis_hits_quant %>% 
  #count(pass_thr)
  mutate(locuseq = dirname(path_sumstat) %>% str_remove(".*gwas/")) %>%
  distinct(locuseq, meta_trait_desc, pass_thr, path_sumstat) %>%
  #summarise(has_hit = any(pass_thr), .by = "locuseq") %>%
  count(locuseq, pass_thr) %>%
  spread(pass_thr, n, sep = "_") %>%
  arrange(- pass_thr_TRUE) %>% #head(10)
  summarize(sum(pass_thr_FALSE))

# From total 161,142 region-trait (1506 cis pQTLs * 107 quant trait GWASs)
# 78,068 region-trait GWAS had at least one significant variant (p > 1e-4)
# 83,074 region-trait GWAS were non-significant. 

"/scratch/dariush.ghasemi/projects/pqtl_susie/results/believe_estVarF/tmp/seq.6912.6_2_89242311_89246140_sumstat.csv"
"/scratch/dariush.ghasemi/projects/pqtl_susie/results/believe_estVarF/tmp/seq.6364.7_12_6288176_7250542_sumstat.csv"
"/scratch/dariush.ghasemi/projects/pqtl_susie/results/believe_estVarF/tmp/seq.11372.2_5_177339176_177415473_sumstat.csv"

"/scratch/dariush.ghasemi/projects/pqtl_susie/results/believe_estVarF/tmp/seq.23690.228_5_177388765_177415473_sumstat.csv" %>%
  fread(header = FALSE, col.names = headers, sep = "\t", data.table = FALSE) %>% 
  View()

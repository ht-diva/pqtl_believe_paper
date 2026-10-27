
library(tidyverse)
library(glue)
library(data.table)


my_locuseq <- "seq.9253.52_9_132740402_133566789" # lambda 0.2645
my_locuseq <- "seq.8309.12_3_45982252_52219282"
my_locuseq <- "seq.3206.4_11_10500217_10698007" # from config/believe_test.csv
my_locuseq <- "seq.11563.51_3_165718150_166124518" # NO CS

my_chr <- str_split_fixed(my_locuseq, "_", 4)[[2]]
my_beg <- str_split_fixed(my_locuseq, "_", 4)[[3]]
my_end <- str_split_fixed(my_locuseq, "_", 4)[[4]]

path_susie <- "/scratch/dariush.ghasemi/projects/pqtl_susie/"
path_workspace <- "results/believe_nomir/tmp/"
path_workspace <- "results/believe_bgen_test/tmp/"
path_fitness <- "results/believe_bgen/susierss/cs_fitness/"

path_sumstat   <- glue(path_susie, path_workspace, my_locuseq, "_sumstat.csv")
path_ld_header <- glue(path_susie, path_workspace, my_locuseq, "_ld.headers")
path_ld_matrix <- glue(path_susie, path_workspace, my_locuseq, "_ld.matrix")

# SuSiE model fitness
path_susie_obj <- glue(path_susie, path_fitness, my_locuseq, "_fit.rds")

# LD computed based on hardcall genotype (bed/bim/fam)
path_ld_bed <- glue(path_susie, "/plink_ld/ld/", my_locuseq, "_ld.headers")
path_sumstat_bgen <- "/exchange/healthds/pQTL/BELIEVE/Working_shared/bgeno_assoc_phen_seq.9253.52_res/bgeno_assoc_phen_seq.9253.52_res.gwaslab.tsv.gz"

# Original genotype before harmonization
path_geno_pgen <- glue(path_susie, "plink_ld/no_harm/", my_locuseq, ".pgen")
path_geno_psam <- str_replace(path_geno_pgen, ".pgen", ".psam")
path_geno_pvar <- str_replace(path_geno_pgen, ".pgen", ".pvar")

# LD with PLINK without harmonizing genotype
path_ld_header_noharm <- glue(path_susie, "plink_ld/no_harm/", my_locuseq, ".unphased.vcor1.vars")
path_ld_matrix_noharm <- glue(path_susie, "plink_ld/no_harm/", my_locuseq, ".unphased.vcor1")



# BGEN GWAS created by Solene, fastGWA accounted dosages
path_gwas_bgen <- "/group/diangelantonio/users/Solene/pQTL/BELIEVE/gcta_dosage/bgeno_assoc_phen_seq.9253.52_res.fastGWA"
path_gwas_regenie <- "/exchange/healthds/pQTL/pQTL_workplace/regenie_unrelated_samples/test_targets_1_2_output/test_targets_1_2/results/gwas/seq.9253.52_res.gwas.regenie.gz"
path_gwas_noqc <- "/exchange/healthds/pQTL/pQTL_workplace/regenie_unqced/test_targets_1_2_output/test_targets_1_2/results/gwas/seq.9253.52_res.gwas.regenie.gz"
path_geno  <- "/center/healthds/pQTL/BELIEVE/Genetic_QC_files/HDS_BELIEVE_final_HDS.pgen"


# Believe
headers <- c("CHR", "POS", "SNPID", "EA", "NEA", "EAF", "BETA", "SE", "P", "MLOG10P", "N", "Z")


# Params to reproduce SuSiE pipe
compute_ld_from_X <- FALSE
susie_L <- 10
susie_iter <- 1000
susie_min_abs_cor <- 0.5
susie_est_resvar <- FALSE  # TRUE if using in-sample LD



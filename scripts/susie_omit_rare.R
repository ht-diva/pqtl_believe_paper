
library(tidyverse)


my_locuseq <- "seq.9253.52_9_132740402_133566789" # lambda 0.2645
path_susie <- "/scratch/dariush.ghasemi/projects/pqtl_susie/"

path_sumstat   <- glue(path_susie, "results/believe_nomir/tmp/", my_locuseq, "_sumstat.csv")
path_ld_header <- glue(path_susie, "results/believe_nomir/tmp/", my_locuseq, "_ld.headers")
path_ld_matrix <- glue(path_susie, "results/believe_nomir/tmp/", my_locuseq, "_ld.matrix")


# Believe
headers <- c("CHR", "POS", "SNPID", "EA", "NEA", "EAF", "BETA", "SE", "P", "MLOG10P", "Z")


# Read GWAS subset
sumstat <- fread(path_sumstat, header = F, col.names = headers, data.table = F)

# Eleminate variants with MAF <= 0.01
sumstat_norare <- sumstat %>%
  dplyr::mutate(MAF = ifelse(EAF <= 0.5, EAF, 1 - EAF)) %>%
  dplyr::filter(MAF > 0.01)


#-------------------------------#
# -----     LD symmetry    -----
#-------------------------------#

# load LD matrix created with Plink2 --r-unphased 'matrix' 'ref-based'
ld_header <- fread(path_ld_header, header = FALSE, col.names = "SNP")
R <- fread(path_ld_matrix, header = F, data.table = F) %>% as.matrix()
rownames(R) <- colnames(R) <- ld_header$SNP


# Keep common variants (MAF > 0.01) in LD matrix
R_norare <- R[sumstat_norare$SNPID, sumstat_norare$SNPID]


# take values
n_believe <- 9216
betas     <- sumstat_norare$BETA
se_betas  <- sumstat_norare$SE
z_scores  <- betas / se_betas


# The estimated λ is
lambda <- susieR::estimate_s_rss(z = z_scores, R = R_norare, n = n_believe)


# Lambda value slightly reduced: 0.2645 -> 0.2602


# Base diagnostic plot 
# Compute conditional Z-score
condz <- susieR::kriging_rss(z = z_scores, R = R_norare, n = n_believe)

# Plot + lambda
condz$plot + 
  labs(
    title = paste("SeqID-Locus= ", my_locuseq),
    subtitle = paste("λ =", signif(lambda, 4))
  )


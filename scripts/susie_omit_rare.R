
library(tidyverse)


my_locuseq <- "seq.9253.52_9_132740402_133566789" # lambda 0.2645
path_susie <- "/scratch/dariush.ghasemi/projects/pqtl_susie/"

path_sumstat   <- glue(path_susie, "results/believe_nomir/tmp/", my_locuseq, "_sumstat.csv")
path_ld_header <- glue(path_susie, "results/believe_nomir/tmp/", my_locuseq, "_ld.headers")
path_ld_matrix <- glue(path_susie, "results/believe_nomir/tmp/", my_locuseq, "_ld.matrix")

# LD computed based on hardcall genotype (bed/bim/fam)
path_ld_header <- "/scratch/dariush.ghasemi/projects/pqtl_susie/plink_ld/ld/seq.9253.52_9_132740402_133566789_ld.headers"
path_ld_matrix <- gsub(".headers", ".matrix", path_ld_header)

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
betas     <- sumstat$BETA
se_betas  <- sumstat$SE
z_scores  <- (betas / se_betas)


# The estimated λ is
lambda_hardcall <- susieR::estimate_s_rss(z = z_scores, R = R, n = n_believe)
lambda_dosage   <- susieR::estimate_s_rss(z = z_scores, R = R, n = n_believe)



# Lambda value slightly reduced: 0.2645 -> 0.2602

#------------#

#estimated effect from fitted model
beta_hat <- colSums(res_rss$alpha*res_rss$mu)
z2 <- z_scores
z2[unique(unlist(res_rss$sets$cs))] = 0
z3 <- as.vector(sqrt(n_believe - 1) * R_norare%*%beta_hat)
res_z <- z_scores - z3

susieR::estimate_s_rss(z = res_z, R = R_norare, n = n_believe)

# Remove 10 CS variants from sumstat and R matrix
unique(unlist(res_rss$sets$cs))
# Will tell us how influenced the lambda is from those single huge signals

cs_snps <- unique(unlist(res_rss$sets$cs))

lambda_noncs <- estimate_s_rss(
  z_scores[- cs_snps],
  R_norare[- cs_snps, - cs_snps],
  n_believe)


# Base diagnostic plot 
# Compute conditional Z-score
susieR::kriging_rss(z = res_z, R = R_norare, n = n_believe)

#------------#

condz_dosage <- susieR::kriging_rss(
  z = z_scores,
  R = R,
  n = n_believe
  )


# Plot + lambda
qq_dosage <- condz_dosage$plot + 
  labs(
    #title = paste("SeqID-Locus= ", my_locuseq),
    subtitle = paste("λ =", signif(lambda_dosage, 4), "| dosage")
  )


# Plots for GWAS-LD with MAF filter
#qq_hardcall_norare <- qq_hardcall
#qq_dosage_norare   <- qq_dosage

qq_joint_norare <- gridExtra::grid.arrange(
  qq_dosage, qq_hardcall,
  qq_dosage_norare, qq_hardcall_norare,
  nrow = 2, ncol = 2)


ggsave("10-Jun-26_lambda_dosage_vs_hardcall.png", 
       plot = qq_joint_norare, width = 7.5, height = 6.5, dpi = 200)


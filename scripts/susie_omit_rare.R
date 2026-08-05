
library(tidyverse)


my_locuseq <- "seq.9253.52_9_132740402_133566789" # lambda 0.2645
path_susie <- "/scratch/dariush.ghasemi/projects/pqtl_susie/"

path_sumstat   <- glue(path_susie, "results/believe_nomir/tmp/", my_locuseq, "_sumstat.csv")
path_ld_header <- glue(path_susie, "results/believe_nomir/tmp/", my_locuseq, "_ld.headers")
path_ld_matrix <- glue(path_susie, "results/believe_nomir/tmp/", my_locuseq, "_ld.matrix")

# LD computed based on hardcall genotype (bed/bim/fam)
path_ld_bed <- glue(path_susie, "/plink_ld/ld/", my_locuseq, "_ld.headers")
path_sumstat_bgen <- "/exchange/healthds/pQTL/BELIEVE/Working_shared/bgeno_assoc_phen_seq.9253.52_res/bgeno_assoc_phen_seq.9253.52_res.gwaslab.tsv.gz"


#-------------------------------#
# -----     Read GWAS      -----
#-------------------------------#

# Believe
headers <- c("CHR", "POS", "SNPID", "EA", "NEA", "EAF", "BETA", "SE", "P", "MLOG10P", "Z")

# Read GWAS subset
sumstat <- fread(path_sumstat, header = F, col.names = headers, data.table = F)
sumstat_bgen <- fread(path_sumstat_bgen)


# Eleminate variants with MAF <= 0.01
sumstat_norare <- sumstat %>%
  dplyr::mutate(MAF = ifelse(EAF <= 0.5, EAF, 1 - EAF)) %>%
  dplyr::filter(MAF > 0.01)


# Subset of GWAS for the locus +/-100 Kbp
sumstat_dosbase <- sumstat_bgen %>% 
  filter(CHR == 9, POS <= 133566789 + 100000, POS >= 132740402 - 100000)


#-------------------------------#
# -----    Compare GWAS    -----
#-------------------------------#

# Check SNP order
summary(match(sumstat$SNPID, sumstat_dosbase$SNPID))


par(mfrow = c(2,2))

# Visualize freq, effect size, and test statistic
plot(sumstat$EAF, sumstat_dosbase$EAF,
     main = "GWAS EAF", xlab = "with PGEN", ylab = "with BGEN")

plot(sumstat$BETA, sumstat_dosbase$BETA,
     main = "GWAS BETA", xlab = "with PGEN", ylab = "with BGEN")

plot(sumstat$Z, sumstat_dosbase$Z,
     main = "GWAS Z", xlab = "with PGEN", ylab = "with BGEN")

plot(sumstat$MLOG10P, sumstat_dosbase$MLOG10P,
     main = "GWAS MLOG10P", xlab = "with PGEN", ylab = "with BGEN")

# Add a 45-degree bisecting line (intercept = 0, slope = 1)
abline(a = 0, b = 1, col = "red", lwd = 2, lty = 2)


#-------------------------------#
# -----       Read LD      -----
#-------------------------------#

# load LD matrix created from PGEN
ld_header <- fread(path_ld_header, header = FALSE, col.names = "SNP")
R <- fread(path_ld_matrix, header = F, data.table = F) %>% as.matrix()
rownames(R) <- colnames(R) <- ld_header$SNP

# load LD matrix created from BED
ld_header_bed <- fread(path_ld_bed, header = FALSE, col.names = "SNP")
R_bed <- fread(gsub(".headers", ".matrix", path_ld_bed), header = F, data.table = F) %>% as.matrix()
rownames(R_bed) <- colnames(R_bed) <- ld_header_bed$SNP

# LD plot
#heatmap(R_bed)


# Keep common variants (MAF > 0.01) in LD matrix
R_norare <- R[sumstat_norare$SNPID, sumstat_norare$SNPID]


# take values
n_believe <- 9216
betas     <- sumstat$BETA
se_betas  <- sumstat$SE
z_scores  <- (betas / se_betas)

#-------------------------------#
# -----       Diag LD      -----
#-------------------------------#

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


#-------------------------------#
# ----  Hardcall vs Dosage  ----
#-------------------------------#

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


#-------------------------------#
# -----    PGEN vs. BGEN   -----
#-------------------------------#

# GWAS with PGEN vs. GWAS with BGEN
lambda_bgen_ld_pgen <- susieR::estimate_s_rss(
  z = sumstat_dosbase$BETA/sumstat_dosbase$SE, R=R, n=9216
)

lambda_bgen_ld_bed  <- susieR::estimate_s_rss(
  z = sumstat_dosbase$BETA/sumstat_dosbase$SE, R=R_bed, n=9216
  )

# Kriging plot
krig_bgen_ld_pgen <- susieR::kriging_rss(
  z = sumstat_dosbase$BETA/sumstat_dosbase$SE, R=R, n=9216
  )

krig_bgen_ld_bed  <- susieR::kriging_rss(
  z = sumstat_dosbase$BETA/sumstat_dosbase$SE, R=R_bed, n=9216
  )

# Label and store plot
qq_bgen_ld_pgen <- krig_bgen_ld_pgen$plot + 
  labs(subtitle = paste("λ =", signif(lambda_bgen_ld_pgen, 4), "| GWAS with BGEN - LD with dosage"))

qq_bgen_ld_bed  <- krig_bgen_ld_bed$plot  +
  labs(subtitle = paste("λ =", signif(lambda_bgen_ld_bed,  4), "| GWAS with BGEN - LD with hardcall"))

# Join plot
qq_joint_bgen <- gridExtra::grid.arrange(
  qq_bgen_ld_pgen, qq_bgen_ld_bed,
  nrow = 1, ncol = 2)


ggsave("11-Jun-26_lambda_pgen_vs_bgen.png", 
       plot = qq_joint_bgen, width = 11.5, height = 6.5, dpi = 200)


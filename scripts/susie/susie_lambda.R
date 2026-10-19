#=============================#
# June 4, 2026
# Dariush Ghasemi

# This script investigate high- 
# lambda issue for BELIEVE's pQTLs.
#=============================#



#----------------------------#
# ----   LD with PLINK   ----
#----------------------------#

# load LD matrix created from PGEN
ld_header <- fread(path_ld_header_noharm, header = FALSE, col.names = "SNP")
R_plink   <- fread(path_ld_matrix_noharm, header = F, data.table = F) %>% as.matrix()
rownames(R_plink) <- colnames(R_plink) <- ld_header$SNP

dim(R_plink)

#----------------------------#
# ----   LD with Dosage  ----
#----------------------------#

# Read psam & pvar
#geno_pgen <- pgenlibr::NewPgen(path_geno_pgen)
#psam_df <- read.delim(path_geno_psam, header = TRUE, comment.char = "")
#pvar_df <- read.delim(header = TRUE, comment.char = "")

# Extract dosages for all of variants
#dosage <- pgenlibr::ReadList(geno_pgen, 1:nrow(pvar_df), meanimpute = FALSE)

# Add variant/sample IDs as column/row names
#colnames(dosage) <- pvar_df$ID
#rownames(dosage) <- psam_df$IID

# Cross variants to avoid allele mismatch
#common_snps <- intersect(sumstat$SNP, ld_header$SNP)
length(common_snps)

#X <- dosage[, common_snps] %>% as.matrix()

# Compute LD (r) correlation matrix
#R <- cor(X, use = "pairwise")

dim(R)

#heatmap(R_plink[1:1000, 1:1000], Rowv = NA, Colv = NA)
#heatmap(R[1:1000, 1:1000], Rowv = NA, Colv = NA)


#----------------------------#
# ----        GWAS       ----
#----------------------------#

gwas <- fread(path_gwas)
#gwas_unrelated <- fread(path_gwas_regenie)
gwas_noqc <- fread(path_gwas_noqc)


sumstat <- gwas_noqc %>%
  dplyr::filter(
    #CHR == my_chr, between(POS, my_beg, my_end)
    CHROM == my_chr, between(GENPOS, my_beg, my_end)
  )

dim(sumstat)

sumstat <- ld_header %>%
  inner_join(
    sumstat %>% mutate(Z = BETA/SE),
    join_by(SNP == ID)
  )

dim(sumstat)


#----------------------------#
# ----       Lambda      ----
#----------------------------#

common_snps <- intersect(sumstat$SNP, ld_header$SNP)
length(common_snps)

n_believe <- min(sumstat$N)


# The estimated λ is
lambda <- susieR::estimate_s_rss(
  z = sumstat$Z,
  R = R,
  #R = R_plink[common_snps, common_snps],
  #n = n_believe,
  #method = "null-pseudomle", # default = "null-mle", or "null-partialmle"
  )


# Without specifying N, the estimated λ was:

#  - method = "null-mle":        0.05005966
#  - method = "null-pseudomle":  0.08724696
#  - method = "null-partialmle": 0.00076204


# plot for unrelated GWAS
plt_kriging_unadj <- susieR::kriging_rss(z = sumstat$Z, R = R) #, n = n_believe

plt_kriging$plot +
  labs(
    #title = paste0(my_locuseq, "\nGWAS on 6772 unrelated individuals in BELIEVE"),
    title = paste0(my_locuseq, "\nGWAS without QC"),
    subtitle = paste("λ =", signif(lambda, 4))
  )

#----------#
par(mfrow = c(1, 2))

plot(
  plt_kriging_unadj$conditional_dist$z, plt_kriging_unadj$conditional_dist$z_std_diff,
  xlab= "Cond_Z", ylab = "Cond_Z_std_diff",
  main = "Un-adjusted", #"Adjusted for N",
  cex = 1, pch = 1 #lab = "Adjusted"
  )

abline(h = 0)
abline(v = 0)

dev.off()


#----------------------------#
# ----   Lambda vs. Z    ----
#----------------------------#


find_lambda <- function(z_thre){
  
  sumstat2 <- ld_header %>%
    inner_join(
      sumstat %>% dplyr::mutate(Z = BETA/SE) %>% dplyr::filter(between(Z, - z_thre, z_thre)),
      join_by(SNP)
    )
  
  common_snps <- intersect(sumstat2$SNP, ld_header$SNP)
  nsnps_overlap <- length(common_snps)
  
  
  lambda <- susieR::estimate_s_rss(z = sumstat2$Z, R = R_plink[common_snps, common_snps], n = n_believe)
  
  tribble(
    ~threshold, ~nsnps_overlap, ~lambda,
    z_thre, nsnps_overlap, lambda
    )

}

list_lambda <- map_dfr(seq(80, 5, -5), find_lambda)

list_lambda %>% DT::datatable()#knitr::kable(options = 3)

list_lambda %>%
  ggplot(aes(x = threshold, y = lambda))+
  geom_point() +
  geom_line() +
  scale_x_continuous(breaks = seq(0,80, 5))+
  scale_y_continuous(breaks = seq(0,.25, .025))+
  ggtitle(my_locuseq) +
  labs(x = "GWAS with |Z| < threshold") +
  theme_light()+
  theme(axis.text = element_text(size = 12))

ggsave("13-Jul-26_zscore_vs_lambda_locus_9_132740402_133566789.png",
       height=5.5, width=7.5, dpi = 200)

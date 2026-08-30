
# Install susieR version 0.16.6
remotes::install_github("stephenslab/susieR")

#-------------------------------#
# -----     SuSiE GWAS     -----
#-------------------------------#

path_gnh_ld <- "/scratch/dariush.ghasemi/projects/pqtl_susie/ld/ld_ten_loci/"

fread("/scratch/dariush.ghasemi/projects/pqtl_coloc/results/believe_cis_quant/tmp/gwas/seq.3583.54_22_50495131_50746518/seq.3583.54_22_50495131_50746518_genesandhealth_v010_quantitative_traits_median_values_f5ff31e8c6.csv.gz")

# GWAS from Genes&Health study
my_seqid <- "seq.9175.48"
my_locus <- "21_40688572_40801365"
my_locuseq <- str_c(my_seqid, "_", my_locus)


### GWAS results
path_coloc <- "/scratch/dariush.ghasemi/projects/pqtl_coloc/results/believe_cis_quant/tmp/gwas/"


# 107 non-binary traits GWAS
path_gwas <- list.files(
  glue(path_coloc, my_locuseq),
  pattern = "*.csv.gz", recursive = TRUE, full.names = TRUE
  )

# Import
sumstat <- fread(path_gwas[1])

sumstat <- sumstat %>%
  dplyr::filter(SNPID %in% ld_headers$SNP)
  #filter(POS >= 40688572, POS <= 40801365)

dim(sumstat)
n_distinct(sumstat$SNPID)


# Regional plot
sumstat %>%
  ggplot(aes(x = POS, y = MLOG10P)) +
  geom_point(size = 3, fill = "#7e4a35", shape = 21) +
  labs(x = "Genomic Position (hg38)") + 
  #ggtitle(paste0(seqid_locus, "(size: ", my_locus$loci_cat,")")) + 
  theme_light() +
  theme(
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    axis.ticks.length = unit(2.5,"mm")
  )


### LD matrix
path_ld_matrix <- glue(path_gnh_ld, my_locus, "_ld.matrix")
path_ld_header <- gsub(".matrix", ".header", path_ld_matrix)


# load LD matrix created with Plink2 --r-unphased 'matrix' 'ref-based'
ld_headers <- fread(path_ld_header, header = FALSE, col.names = "SNP")
R <- fread(path_ld_matrix) %>% as.matrix()
rownames(R) <- colnames(R) <- ld_headers$SNP


dim(R)
R[1:5, 1:5]


# Check SNPs overlap
intersect(sumstat$SNPID, ld_headers$SNP) %>% length()

# Select only variants in region sumstat
R <- R[sumstat$SNPID, sumstat$SNPID]

heatmap(R, Rowv = NA, Colv = NA)

# N in GWAS
n_GnH <- min(sumstat$meta_total_samples, na.rm = TRUE)

susieR::estimate_s_rss(
  z = sumstat$BETA/sumstat$SE, R=R, n=9216
)

susieR::kriging_rss(
  z = sumstat$BETA/sumstat$SE,
  R = R, n = n_GnH
)

compute_ld_from_X <- FALSE
susie_L <- 10
susie_iter <- 1000
susie_min_abs_cor <- 0.5
susie_est_resvar <- FALSE  # TRUE if using in-sample LD

# Check point
identical(sumstat$SNPID, rownames(R))


err_handling <- function(e) { stop("❌ SuSiE failed: ", e$message) }

res_rss_score <- tryCatch(
  susie_rss(
    #bhat = sumstat$BETA,
    #shat = sumstat$SE,
    z = z_scores,
    z_method = "score",
    n = n_believe,
    R = R,
    L = susie_L,
    max_iter = susie_iter,
    min_abs_corr = susie_min_abs_cor,
    estimate_residual_variance = susie_est_resvar,
  ),
  error = err_handling
)


# SuSiE failed: The correlation matrix (7170 by 7170) is 
# not a positive semidefinite matrix. The smallest eigenvalue 
# is -1.26221817719385e-05. You can bypass this by "check_R = FALSE" 
# which instead sets negative eigenvalues to 0 to allow for continued computations.


# Extract results
full_res <- summary(res_rss)
cs   <- full_res$cs    # containing CS impurity indices
vars <- full_res$vars  # containing CS Posterior Inclusion Probabilities



# SuSiE Plot
susieR::susie_plot(
  res_rss,
  y = "PIP",
  b = betas,
  xlab = "Variants",
  add_bar = FALSE,
  add_legend = TRUE,
  main = paste("SeqID_locus:", my_locuseq) #, "\nmethod=score"
)



# list of the entire SNPs with PIP
snps_pip <- vars %>%
  transmute(
    cs_id = cs,
    SNPID = sumstat$SNPID[variable],
    PIP = variable_prob
  )

# subset of GWAS results for CS variants
cs_summary <- sumstat %>%
  left_join(snps_pip, by = "SNPID") %>%
  # append CS characteristics to sumstat
  left_join(cs[1:4], join_by(cs_id == cs)) %>%
  filter(cs_id > 0)

# list of CS variants
cs_list <- cs_summary %>%
  summarize(
    cs_snps = paste(SNPID, collapse = ","),
    .by = cs_id
  ) %>%
  full_join(cs, join_by(cs_id == cs))


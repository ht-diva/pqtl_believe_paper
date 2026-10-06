
library(dplyr)
library(stringr)
library(data.table)
library(glue)
library(susieR)



# locus with no credibles
my_locuseq <- "seq.9940.35_2_241163697_241740802" #"seq.8922.4_3_165436169_165755715"

# input file is too small
my_locuseq <- "seq.21660.4_9_117083803_117084672" #"seq.23771.17_5_176836532_176842474"

# VTN locus
my_locuseq <- "seq.13125.45_17_25283026_29975846" #"seq.15333.11_17_25854347_28884532"

# LD is not positive semi-definite
my_locuseq <- "seq.22003.4_3_165435919_165755916"

# locus with incorrectly indexed credible
my_locuseq <- "seq.2638.12_17_6958690_7205230"

# locus failed in interval error out of bound
my_locuseq <- "seq.23030.4_19_51923541_52759987"

# ❌ SuSiE failed: Estimating residual variance failed: the estimated value is negative
#my_locuseq <- "seq.22377.27_2_89119570_94182751"
my_locuseq <- "seq.22405.61_2_89119570_94877043"

# Believe Test locus
my_locuseq <- "seq.18373.13_22_43928847_43998522"

# inputs
path_proj <- "/scratch/dariush.ghasemi/projects/pqtl_susie/results/"

path_cs_rds <- glue(path_proj, "susierss/cs_fitness/", my_locuseq, "_fit.rds")
path_cs_rds <- "cs_fitness/seq.15333.11_17_25854347_28884532_fit.rds"
path_cs_rds <- "susierss/cs_fitness/seq.23030.4_19_51923541_52759987_fit.rds"
example_fit <- glue(path_proj, "believe_test/susierss/cs_fitness/seq.11836.144_21_41449719_41533425_fit.rds")
example_sumstat <- glue(path_proj, "believe_test/tmp/seq.11836.144_21_41449719_41533425_sumstat.csv")


# LD files
path_sumstat   <- glue(path_proj, "believe_test/tmp/", my_locuseq, "_sumstat.csv")
path_ld_header <- glue(path_proj, "believe_negvar_89loci/tmp/", my_locuseq, "_ld.headers")
path_ld_matrix <- glue(path_proj, "believe_negvar_89loci/tmp/", my_locuseq, "_ld.matrix")


#-------------------------------#
#-----     GWAS sumstat    -----
#-------------------------------#

# input files
#path_pgen <- glue(path_proj, "tmp/", my_locuseq, "_dosage.pgen")


# Interval GWAS column names
#headers <- c("CHR", "POS", "SNPID", "EA", "NEA", "EAF", "N", "BETA", "SE", "MLOG10P", "CHISQ")

# Believe
headers <- c("CHR", "POS", "SNPID", "EA", "NEA", "EAF", "BETA", "SE", "P", "MLOG10P", "Z")


# Read GWAS sumstat
#sumstat <- fread(example_sumstat, col.names = sumstat_header)
#sumstat <- fread(path_sumstat, sep = "\t", data.table = FALSE)
sumstat <- fread(path_sumstat, header = FALSE, col.names = headers, sep = "\t", data.table = FALSE)

#colnames(sumstat)[which(names(sumstat) == "##CHR")] <- "CHR"

n_distinct(sumstat$SNPID)
n_distinct(sumstat$POS)
nrow(sumstat)


#-------------------------------#
#-----    Count Variants   -----
#-------------------------------#

# Per-row variant type
sumstat2 <- sumstat %>%
  mutate(
    is_indel = nchar(EA) != 1 | nchar(NEA) != 1,
    is_snp   = nchar(EA) == 1 & nchar(NEA) == 1
  )

# Count indels and SNPs
n_indels <- sum(sumstat2$is_indel)
n_snps   <- sum(sumstat2$is_snp)

# Site-level allele counts
site_counts <- sumstat2 %>%
  group_by(CHR, POS) %>%
  summarise(
    n_alleles = n_distinct(c(EA, NEA)),
    .groups = "drop"
  )

# Count bi-allelic and multi-allelic sites
n_biallelic    <- sum(site_counts$n_alleles == 2)
n_multiallelic <- sum(site_counts$n_alleles > 2)

list(
  n_snps = n_snps,
  n_indels = n_indels,
  n_biallelic = n_biallelic,
  n_multiallelic = n_multiallelic
)


# take values
betas    <- sumstat$BETA
se_betas <- sumstat$SE
z_scores <- betas / se_betas
#n        <- min(sumstat$N, na.rm = TRUE)
n <- 9216


#-------------------------------#
# -----     LD symmetry    -----
#-------------------------------#

# load compressed LD matrix
#fread("/scratch/dariush.ghasemi/projects/pqtl_susie/plink_ld/ld/eq.8280.238_17_26592946_26739127_ld.matrix.zst")


# load LD matrix created with Plink2 --r-unphased 'matrix' 'ref-based'
ld_header <- fread(path_ld_header, header = FALSE, col.names = "SNP")
R <- fread(path_ld_matrix, header = F, data.table = F) %>% as.matrix()
rownames(R) <- colnames(R) <- ld_header$SNP

dim(R)
nrow(sumstat)

# heatmap of LD
png("22-Oct-25_heatmap_ld_with_cor.png", width=8, height=6, units = "in", res = 300)

heatmap(ld, main = 'Plink-based LD')
heatmap(R,  main = 'Correlation-based LD')

dev.off()

isSymmetric(R) # TRUE

#-----------#
# Check if SNPIDs are in the same order in LD and in GWAS
# if LD SNPID does not match with GWAS SNPID, it returns NA
match(
  rownames(R_norare),
  sumstat_norare$SNPID,
  #nomatch = 0
  ) %>% summary()

#-----------#

# The estimated λ is
susieR::estimate_s_rss(z = z_scores, R=R, n=n)

# Base diagnostic plot
condz = susieR::kriging_rss(z = z_scores, R=R, n=n)

condz$plot +
  labs(
    subtitle = paste("λ =", signif(lambda, 4))
  )


e <- eigen(R, symmetric = TRUE, only.values = TRUE)$values
cat("Minimum eigenvalue:", min(e), "\n")
cat("Number negative:", sum(e < 0), "\n")
cat("Condition number:", max(e) / min(abs(e)), "\n")


# Capture error message
warn_txt <- NA_character_

withCallingHandlers(
  {
    lambda <- susieR::estimate_s_rss(z = betas / se_betas, R = R, n = n)
  }, 
  message = function(m) {                 # Let it continue to the sink
    warn_txt <<- gsub("\033\\[[0-9;]*m|\\n", "", conditionMessage(m)) # omit HTML coloring warning
    message(conditionMessage(m))          # re-emit original (with color) → goes to log
    invokeRestart("muffleMessage")        # suppress duplicate print
  }
)

message("✅ The estimated λ is ", lambda)

lambda <- tryCatch({
  
  value <- susieR::estimate_s_rss(z = z_scores, R = R, n = n)
  
  list(
    lambda = value,
    warnme = NA_character_
  )
  
}, message = function(w) {
  
  message("⚠️ ", conditionMessage(w))
  
  return(
    list(
      lambda = value,
      warnme = gsub("\033\\[[0-9;]*m", "", conditionMessage(w))
    )
  )
  })



#-------------------------------#
# -----     Run SuSiE      -----
#-------------------------------#

err_handling <- function(e) { stop("❌ SuSiE failed: ", e$message) }

res_rss <- tryCatch(
  susie_rss(
    bhat = betas,
    shat = se_betas,
    n = n_believe,
    R = R_norare,
    L = 10,
    max_iter = 1000,
    min_abs_corr = 0.5,
    estimate_residual_variance = FALSE
  ),
  error = err_handling
)


#----------------------------#
#----  Examp. SuSiE obj  ----
#----------------------------#

# load summary of susie_rss 
cs_rds <- readRDS(glue(path_proj, path_cs_rds))
my_rss <- readRDS(example_fit)

# full model summary
full_res <- summary(my_rss)

# Credible sets
#cs <- susie_get_cs(res_rss, X = NULL) # issue #257: generates more credible sets as it applies NO impurity filter

cs   <- full_res$cs    # containing CS impurity indices
vars <- full_res$vars  # containing CS Posterior Inclusion Probabilities

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
  left_join(cs[1:4], join_by(cs_id == cs)) %>%
  filter(cs_id > 0)


#-------------------------------#
# -----   Check credible   -----
#-------------------------------#


# CS correlation
get_cs_correlation(res_rss, Xcorr = R_norare) %>% corrplot::corrplot()


# retrieve credibles
my_cs <- susie_get_cs(res_rss, X = NULL) # issue #257: generates more credible sets as it applies NO impurity filter
my_cs <- susie_get_cs(res_rss, Xcorr = ld, min_abs_corr = 0.5) # generates correct credible sets passing impurity filter.

my_cs$cs  # this way resulted in 3 CS, due to lack of impurity filter

full_res <- summary(res_rss) # this way gave rise to 2 CS
cs <- full_res$cs
cs

vars <- full_res$vars

# list of all snps with PIP
snps_pip <- vars %>%
  transmute(
    cs_id = cs,
    SNPID = sumstat[variable, "SNPID"],
    PIP = variable_prob
    )

cs_summary <- sumstat %>%
  left_join(snps_pip, by = "SNPID") %>%
  filter(cs_id > 0)

cs_summary %>%
  summarize(cs_snps = paste(SNPID, collapse = ","), .by = cs_id) %>%
  full_join(cs, join_by(cs_id == cs)) %>%
  mutate(
    seqid = tag_seqid,
    locus = tag_locus,
    ncs = str_count(variable, ",") + 1
    ) %>%
  select(seqid, locus, cs_id, cs_log10bf, cs_avg_r2, cs_min_r2, ncs, cs_snps)

  

# this does NOT return whole credibles
test_res <- summary(cs_rds)

lapply(
  seq_along(cs$cs),
  function(k) {cs$cs[[k]]}
  )

cs_details <- test_res$cs %>%
  data.frame() %>%
  mutate(
    seqid = tag_seqid,
    locus = tag_locus,
    cs_id = paste0("CS", cs)
  ) %>%
  select(- cs, - variable)


cs_summary <- lapply(unlist(cs$cs), function(k) {
  
  # nm is like "L1", "L2", "L10", etc.
  #cs_num <- as.integer(sub("^L", "", k))   # original CS index
  cs_num <- as.integer(k)   # original CS index
  idx <- cs$cs[[k]]                        # indices of variants (relative to input)
  
  # guard against NULL/empty/invalid indices
  if (is.null(idx) || length(idx) == 0) return(NULL)
  idx <- as.integer(idx)
  #idx <- idx[is.finite(idx) & idx >= 1 & idx <= nrow(sumstat)]
  if (length(idx) == 0) return(NULL)
  
  data.table(
    cs_id  = paste0("CS", cs_num),
    cs_num = cs_num
  )
  
}) %>% rbindlist()


cs_summary %>%
  summarize(
    cs_snps = paste0(unique(SNP), collapse = ","), 
    .by = c("seqid", "locus", "cs_id") # keep seqid and locus in CS list
  ) %>%
  full_join(cs_details, ., join_by(seqid, locus, cs_id)) %>%
  select(cs_id, cs_log10bf, cs_avg_r2, cs_min_r2)


# 
res_rss$sets$lbf
res_rss$lbf_variable


# library("ggpubr")
# library("cowplot")
# library("gridExtra")


# Create a side-by-side layout
par(mfrow = c(1, 2))

susie_plot(betas/se_betas, y = "z", b=betas, add_legend = TRUE)


# plot of posterior priorities
susieR::susie_plot(
  #res_rss, 
  res_susie_ild,
  y="PIP", 
  b=betas, 
  xlab="Variants",
  add_bar = FALSE,
  add_legend = TRUE,
  #main = 'Correlation-based LD'
  #main = 'Plink-based LD'
  main = paste("SeqID-Locus:", my_locuseq)
)

susie_plot_iteration(res_rss, L=10)

susie_get_posterior_mean(res_rss, prior_tol = 1e-09) %>% head()
susie_get_posterior_sd(res, prior_tol = 1e-09)
susieR::susie_get_niter(res_rss)

coloc::plot_dataset(d = rds1_annot, susie_obj = TRUE)



#-------------------------------#
# -----  Coloc with SuSiE  -----
#-------------------------------#

rds1_annot <- coloc::annotate_susie(res_rss, sumstat$SNPID, R)
rds2_annot <- coloc::annotate_susie(res_rss, sumstat$SNPID, R)


# run coloc with susie
library(coloc)
data(coloc_test_data)

res_coloc <- coloc.susie(
  dataset1 = rds1_annot,
  dataset2 = rds2_annot,
  back_calculate_lbf = FALSE,
  susie.args = list()    # only needed if you want coloc to re-run susie via runsusie
)

res_coloc$summary



# Annotate and save full model fitness with LD matrix for coloc
res_rss_annot <- coloc::annotate_susie(res_rss, sumstat$SNPID, R)

saveRDS(res_rss_annot, file = out_cs_annot)
message("✅ Saved LD-annotated SuSiE full summary for coloc: ", out_cs_annot)

#----------------------------------------#
# -----    Create and Save Plots   ------
#----------------------------------------#


# date
start_time <- Sys.time()

end_time <- Sys.time()

# Compute elapsed time
elapsed_time <- end_time - start_time

as.numeric(elapsed_time, units = "secs")
cat("Elapsed time:", round(as.numeric(elapsed_time, units="mins"), 3), "minutes\n")


library(tictoc)

tic("simulation")
x <- replicate(100, mean(rnorm(1e6)))
toc()


#-------------------------------#
# -----     LocusZoom      -----
#-------------------------------#


cs_sumstat <- fread(glue(path_proj, "susierss/cs_summary/", my_locuseq, ".cssum"))

annot_sumstat <- sumstat %>% left_join(cs_sumstat)

annot_sumstat %>% count(cs_id)

annot_sumstat %>%
  mutate(cs_idf = str_replace_na(cs_id, "no_cs") %>% fct_rev()) %>%
  ggplot(aes(x = POS, y = MLOG10P, color = factor(cs_id), label = cs_id)) +
  geom_point(size = 3, alpha = 0.8) + #fill = "white", shape = 21,
  scale_y_continuous(breaks = seq(0,1500, 100), limits = c(0,1500)) +
  guides(color = guide_legend(nrow = 1, position = "bottom", title = "Credible set")) +
  coord_cartesian(clip = "off") +
  theme_bw() +
  ggrepel::geom_label_repel(fill = NA, xlim = c(-Inf, Inf), ylim = c(-Inf, Inf))
  #ggrepel::geom_text_repel(
    #aes(label = cs_id,), 
    #color = 'white', 
  #  size = 3.5,
    #box.padding = unit(0.35, "lines"),
  #  box.padding = 0.5,
    #point.padding = unit(0.3, "lines")
  #  ) +


ggsave(
  paste0("04-Nov-25_regional_plot_", my_locuseq, "_with_credibles.png"),
  dpi = 300, height = 5.5, width = 7.5
  )





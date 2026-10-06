
library(glue)
library(data.table)
library(stringr)
library(dplyr)
library(Rmpfr)
library(coloc)
library(purrr)


# coloc results
path_coloc_res <- "/scratch/dariush.ghasemi/projects/pqtl_coloc/results/believe_cis/combined_coloc_results_quant.tsv"
path_coloc <- "/scratch/dariush.ghasemi/projects/pqtl_coloc/results/believe/tmp/coloc"
path_pwas  <- "/scratch/dariush.ghasemi/projects/pqtl_coloc/results/believe/tmp/pwas/"
path_gwas  <- "/scratch/dariush.ghasemi/projects/pqtl_coloc/results/believe/tmp/gwas/"

#locuseq <- "seq.7056.16_22_17604191_17674295"
#locuseq <- "seq.6255.74_20_2674174_3022373"
locuseq <- "seq.21547.6_7_75024_284666"
pwas_file <- glue(path_pwas, locuseq, "_sumstat.csv.gz")

# Read colocalization results
coloc_res <- fread(path_coloc_res)


# Handle missings when computing p-value from beta&sd
safe_pnorm <- function(b, se, p=FALSE) {
  
  # Ensure the vectors are of the same length
  if(length(b) != length(se)) {
    stop("Beta and SE must be of the same length")
  }
  
  k  <- length(b)
  b  <- as.numeric(b)
  se <- as.numeric(se)
  
  # Initialize result vector with NA values
  result <- rep(NA, k)
  
  # Identify non-missing and non-zero indices
  i <- which(!is.na(b) & !is.na(se) & se != 0)
  
  # compute z-score and take absolute, raise digits with mpfr, apply pnorm for non-missing values
  z_score <- b[i] / se[i]
  z_mpfr <- Rmpfr::mpfr(- abs(z_score), 120)
  p_mpfr <- 2 * pnorm(z_mpfr)
  mlog10p <- - log10(p_mpfr)
  
  # print p-value in character format and mlog10p in numeric
  if(p==TRUE){
    # reformat to mpfr character, then to numeric (don't set digits for MLOG10P)
    mlog10p_mpfr <- Rmpfr::formatMpfr(p_mpfr, scientific = TRUE, digits = 6)
    result[i] <- mlog10p_mpfr
  } else {
    mlog10p_mpfr <- Rmpfr::formatMpfr(mlog10p, scientific = TRUE)
    result[i] <- as.numeric(mlog10p_mpfr)
  }
  
  return(result)
}


prepare4coloc <- function(data){
  
  temp  <- data |>
    dplyr::rename(position = POS, beta = BETA) |>
    dplyr::distinct(position, .keep_all = TRUE) |> # remove duplicate SNPs
    dplyr::rename_with(~gsub("meta_total_samples", "N", .x)) |>
    dplyr::mutate(
      snp = paste0(CHR, ":", position),
      varbeta = SE^2,
      pvalues = safe_pnorm(beta, SE, p = TRUE),
      MAF = ifelse(EAF < 0.5, EAF, 1- EAF),
      sdY = coloc:::sdY.est(varbeta, MAF, N)
      ) |>
    dplyr::select(position, snp, beta, varbeta, MAF, pvalues, sdY, dplyr::any_of("SNPID"))
  
  temp$type <- "quant"
  odata <- as.list(na.omit(temp))
  odata$type <- unique(odata$type)
  odata$sdY <- unique(odata$sdY)
  
  return(odata)
}

#-------------------------------#
# -----     Read files     -----
#-------------------------------#

head_interval <- c("CHR", "POS", "SNPID", "EA", "NEA", "EAF", "N", "BETA", "SE", "MLOG10P", "CHISQ")

sum_locuseq <- fread(pwas_file)#, header = F, col.names = head_interval)


# files with credible sets
sums_lists <- list.files(
  path = glue(path_gwas, locuseq),
  pattern = ".csv.gz", full.names = T
)

sum_test <- fread(sums_lists[1])

#-------------------------------#
# -----     Check files    -----
#-------------------------------#

# check if there is any duplicates
sum_locuseq |> duplicated() |> table()
sum_test |> duplicated() |> table()

# check for duplicate SNPs
sum_test[duplicated(sum_test$POS), ]
sum_locuseq[duplicated(sum_locuseq$POS), ]
sum_locuseq |> dplyr::filter(POS == "3007694")

# removing indels
sum_locuseq |>
  group_by(CHR, POS) |>
  filter(n() <= 1) |>
  ungroup()

sum_locuseq |> prepare4coloc() |> plot_dataset()
sum_test |> prepare4coloc() |> plot_dataset()

sum_locuseq |>
  dplyr::rename(position = POS, beta = BETA) |>
  dplyr::distinct(position, .keep_all = TRUE) |> # remove duplicate SNPs
  dplyr::rename_with(~gsub("meta_total_samples", "N", .x))|>
  dplyr::mutate(
    snp = paste0(CHR, ":", POS),
    varbeta = SE*SE,
    pvalues = safe_pnorm(beta, SE, p = TRUE),
    MAF = ifelse(EAF < 0.5, EAF, 1- EAF),
    sdY = coloc:::sdY.est(varbeta, MAF, N)
  ) |>
  select(CHR, POS, snp, beta, varbeta, EAF, MAF, N, pvalues, MLOG10P)

#test_df <- fread(sums_lists[2]) |> prepare4coloc()
#seq_df  <- sum_locuseq |> prepare4coloc()


#-------------------------------#
# -----     Run Coloc      -----
#-------------------------------#

run_coloc <- function(pfile, gfile){
  
  pwas <- fread(pfile)#, header = F, col.names = head_interval
  gwas <- fread(gfile)
  
  annot_pwas <- prepare4coloc(pwas)
  annot_gwas <- prepare4coloc(gwas)

  # run coloc standard
  res <- coloc::coloc.abf(annot_pwas, annot_gwas)

  res_h4 <- res$summary %>% t() %>% as.data.frame() %>% select(nsnps, PP.H4.abf)

  seqid <- pfile %>% basename() %>% stringr::str_remove("_sumstat.csv.gz")
  pheno <- unique(gwas$meta_trait_desc)

  res_final <- data.frame(
    "seqid" = seqid,
    "phenotype" = pheno
  ) %>%
    cbind(res_h4)

  return(res_final)
}


traits2test <- expand.grid(pwas_file, sums_lists[19], stringsAsFactors = FALSE)

res_combin <- map2_df(traits2test$Var1, traits2test$Var2, run_coloc)



#-------------------------------#
# ----  Check SNPs overlap  ----
#-------------------------------#

check_resolution <- function(pfile, gfile){
  
  pwas <- fread(pfile, header = F, col.names = head_interval)
  gwas <- fread(gfile)
  
  pwas <- pwas |> 
    dplyr::distinct(POS, .keep_all = TRUE) |>
    dplyr::mutate(snp = paste0(CHR, ":", POS))
  
  gwas <- gwas |> 
    dplyr::distinct(POS, .keep_all = TRUE) |>
    dplyr::mutate(snp = paste0(CHR, ":", POS))
  
  n_overlap <- length(intersect(pwas$snp, gwas$snp))
  
  return(n_overlap)
  }

map2_dbl(traits2test$Var1[1:10], traits2test$Var2[1:10], check_resolution)


#-------------------------------#
# -----    Check tileDB    -----
#-------------------------------#

# GWAS files in G&H
gwas_lists <- list.files(
  path = path_gwas,
  pattern = ".csv.gz", full.names = T, recursive = T
)

gwas_downloaded <- gwas_lists %>%
  tibble(.name_repair = ~ "full_path") %>%
  transmute(
    locuseq = dirname(full_path) %>% str_replace("^*.+/gwas/", ""),
    filename = basename(full_path) %>% str_replace("^*.+values_", "")
    ) %>%
  summarise(
    n_traits = n(),
    .by = "locuseq"
  )

# number of incomplete GWASes
gwas_downloaded %>%
  count(n_traits < 107) %>%
  kableExtra::kable(format = "pipe")

hist(gwas_downloaded$n_traits, nclass = 30)

# characteristics of incomplete GWASes
lb_believe %>%
  mutate(locuseq = paste0(seqid, "_", chr, "_", start, "_", end)) %>%
  select(locuseq, loci_width, loci_cat) %>%
  right_join(gwas_downloaded, join_by(locuseq)) %>%
  mutate(job = ifelse(n_traits < 107, "incomplete", "completed")) %>%
  count(job, loci_cat) %>% #spread(job, n, fill = 0)
  ggplot(aes(
    x = factor(loci_cat, levels = level_order),
    y = n,
    fill = job))+
  geom_bar(stat = "identity", position = position_dodge(.9))


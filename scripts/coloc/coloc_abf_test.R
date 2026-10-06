
library(dplyr)
library(glue)
library(data.table)
library(coloc)



path_coloc <- "/scratch/dariush.ghasemi/projects/pqtl_coloc/results/"
#path_dicot <- "/scratch/dariush.ghasemi/projects/pqtl_coloc/tiledb/results/dicot_seq.21547.6_7_75024_284666/_genesandhealth_v010_binary_traits_3digitICD10_b934bf16f4.csv.gz"
#path_cases <- "/scratch/dariush.ghasemi/projects/pqtl_coloc/gnh_cases.txt"
path_dicot <- glue(path_coloc, "believe_dic/tmp/gwas/seq.16300.4_22_43928847_43998522/seq.16300.4_22_43928847_43998522_genesandhealth_v010_binary_traits_3digitICD10_ff8ea4b2da.csv.gz")
path_quant <- glue(path_coloc, "believe_test_quant/tmp/gwas/seq.16300.4_22_43928847_43998522/seq.16300.4_22_43928847_43998522_genesandhealth_v010_quantitative_traits_median_values_0e55cdb425.csv.gz")
#path_pwas  <- glue(path_coloc, "believe_dic/tmp/pwas/seq.15534.26_22_43928847_43985602_sumstat.csv.gz")
path_pwas  <- glue(path_coloc, "believe_test_quant/tmp/pwas/seq.16300.4_22_43928847_43998522_sumstat.csv.gz")


# Read protein GWAS
pwas <- fread(path_pwas)

# Read G&H GWAS
dicot_gwas <- fread(path_dicot)
quant_gwas <- fread(path_quant)

# fread(path_cases) %>%
#   filter(!count_v010 %in% c("", "<10"),
#          `3-digit ICD10 term` == "Congenital malformations of breast"
#          ) %>%
#   distinct(count_v010)  


# files from tileDB
sums_lists <- list.files(
  path = dirname(path_quant),
  pattern = ".csv.gz",
  full.names = TRUE
)


# all possible protein-phenotype pairs
#traits2test <- expand.grid(path_pwas, sums_lists, stringsAsFactors = FALSE)

#-------------------------------#
# -----      FUNCTIONS     -----
#-------------------------------#

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
  
  # compute z-score and take absolute, 
  # raise digits with mpfr, 
  # apply pnorm for non-missing values
  z_score <- b[i] / se[i]
  z_mpfr <- Rmpfr::mpfr(- abs(z_score), 120)
  p_mpfr <- 2 * Rmpfr::pnorm(z_mpfr) # dispatch the correct pnorm()  
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


prepare4coloc <- function(data, dichotomous = FALSE){
  
  temp  <- data |>
    # group_by(position) |>
    # dplyr::slice_max(MLOG10P, n = 1) |> # handle multi-allelic variants
    # ungroup() |>
    dplyr::mutate(
      #snp = paste0(CHR, ":", POS),
      varbeta = SE^2,
      pvalues = safe_pnorm(BETA, SE, p = TRUE),
      MAF = ifelse(EAF < 0.5, EAF, 1- EAF)
    ) |>
    dplyr::rename_with(
      ~gsub("meta_total_samples", "N", .x)
    ) |>
    dplyr::rename(
      snp = SNPID,
      beta = BETA,
      position = POS
    )
  
  if(dichotomous){
    
    temp <- temp |>
      dplyr::mutate(s = meta_total_cases/N) %>%
      dplyr::select(position, snp, beta, varbeta, MAF, pvalues, s)
  
    } else {
      
      temp <- temp %>%
        dplyr::mutate(sdY = coloc:::sdY.est(varbeta, MAF, N)) %>%
        dplyr::select(position, snp, beta, varbeta, MAF, pvalues, sdY)
    }
  
  
  odata <- as.list(na.omit(temp))
  
  if(dichotomous) {
    
    odata$type <- "cc"
    odata$s <- unique(odata$s)
    
    } else {
      
      odata$type <- "quant"
      odata$sdY <- unique(odata$sdY)
  }
  
  return(odata)
}


run_coloc <- function(gfile){
  
  gwas <- fread(gfile)
  
  annot_gwas <- prepare4coloc(gwas, dichotomous = bin_gwas)
  
  # run coloc standard
  res <- coloc::coloc.abf(annot_pwas, annot_gwas)
  
  res_h4 <- res$summary %>% t() %>% as.data.frame()
  
  locus <- path_pwas %>% basename() %>%
    stringr::str_remove("_sumstat.csv.gz") %>%
    stringr::str_remove("seq.(\\d)+.(\\d)+_")
  
  seqid   <- unique(pwas$meta_notes_source_id)
  protein <- unique(pwas$meta_trait_desc)
  pheno <- unique(gwas$meta_trait_desc)
  
  res_final <- data.frame(
    "seqid" = seqid,
    "locus" = locus,
    "protein" = protein,
    "phenotype" = pheno
  ) %>%
    cbind(res_h4)
  
  return(res_final)
}


plot_datasets <- function(d1, d2, color = "dodgerblue3") {
  
  if(!("position" %in% names(d1))) stop("no position element given")
  if(!("position" %in% names(d2))) stop("no position element given")
  
  xlim=c(min(c(d1$position, d2$position)),
         max(c(d1$position, d2$position)))
  
  intsnps=intersect(d1$snp, d2$snp)
  
  par(mfrow=c(2,1))
  
  plot_dataset(
    d1,
    highlight_list=list(intsnps),
    color=color[1],
    main="Dataset 1",
    xlim=xlim,
    legend = FALSE
    )
  
  plot_dataset(
    d2,
    highlight_list=list(intsnps),
    color=color[1],
    main="Dataset 2",
    xlim=xlim,
    legend = FALSE
    )
}


#-------------------------------#
# -----        CHECKS      -----
#-------------------------------#

# check removal of multi-allelic variants
sums_lists[1] %>% fread() %>% arrange(POS)
sums_lists[1] %>% fread() %>% nrow()
sums_lists[1] %>% fread() %>% distinct(POS) %>% nrow()


# quant_gwas |>
#   dplyr::rename(
#     snp = SNPID,
#     position = POS,
#     beta = BETA
#     ) |>
#   group_by(position) |>
#   dplyr::slice_max(MLOG10P, n = 1) |> # handle multi-allelic variants
#   ungroup() |>
#   dplyr::rename_with(~gsub("meta_total_samples", "N", .x)) |>
#   dplyr::mutate(
#     snp = paste0(CHR, ":", position),
#     varbeta = SE^2,
#     pvalues = safe_pnorm(beta, SE, p = TRUE),
#     MAF = ifelse(EAF < 0.5, EAF, 1- EAF)
#   )

bin_gwas <- FALSE
bin_pwas <- FALSE

prepare4coloc(quant_gwas, FALSE)
prepare4coloc(pwas, FALSE) |> coloc::plot_dataset()
prepare4coloc(quant_gwas, FALSE) |> coloc::plot_dataset()
prepare4coloc(dicot_gwas, TRUE) |> plot_dataset()
prepare4coloc(dicot_gwas, TRUE) |> check_dataset(suffix = "")
prepare4coloc(dicot_gwas, TRUE) |> process.dataset(suffix="")


# find duplicate variants
n_distinct(quant_gwas$SNPID)
n_distinct(quant_gwas$POS)

# find duplicate SNPs -- mirrored variants
quant_gwas %>% add_count(SNPID) %>% filter(n > 1)

quant_gwas %>%
  mutate(npos = n_distinct(SNPID), .by = POS) %>% 
  filter(npos >1) %>% dim()


#-------------------------------#
# -----      PLOTTING      -----
#-------------------------------#

pwas %>% 
  dplyr::select(POS, SNPID, MLOG10P) %>%
  inner_join(
    quant_gwas %>% dplyr::select(POS, SNPID, MLOG10P),
    join_by(POS, SNPID), suffix = c("_pwas", "_gwas")
    ) %>%
  ggplot(aes(x = POS)) +
  geom_point(aes(y = MLOG10P_pwas), size = 2, color = "#618685", shape = 21) +
  geom_point(aes(y = MLOG10P_gwas), size = 2, color = "#eca1a6", shape = 21) +
  labs(
    x = "Genomic position",
    y = "-log10(P)",
    title = "Regional association plot: seq.16300.4 (green) vs ALT (pink)"
  ) +
  theme_light()
  



# Reshape PWAS: adding sdY and type = 'quant'
annot_pwas <- prepare4coloc(pwas, dichotomous = FALSE)
#annot_gwas <- prepare4coloc(dicot_gwas, dichotomous = bin_gwas)
annot_gwas <- prepare4coloc(quant_gwas, dichotomous = FALSE)

# plot region
plot_datasets(annot_pwas, annot_gwas)

# test coloc with annotated data
coloc::coloc.abf(annot_pwas, annot_gwas)

# GWASs with i=1 or i=4 include multi-allelic sites 

# iterate run_coloc()
res <- purrr::map(sums_lists[1:10], run_coloc)


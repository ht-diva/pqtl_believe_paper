
library(pgenlibr)

my_seqid <- "seq.9253.52_res"
my_locuseq <- "seq.9253.52_9_132740402_133566789"

# Actual vs. residual proteins
#path_pheno <- "/exchange/healthds/pQTL/BELIEVE/Proteomics_QC_files/Residuals_for_regenie_pQTL/BELIEVE_proteins_NonImp_full_14062024.txt"
path_pheno <- "/exchange/healthds/pQTL/BELIEVE/Proteomics_QC_files/Residuals_for_regenie_pQTL/BELIEVE_NonImp_residuals_new_14_062024.txt"


# Read protein residulas
pheno <- read.table(path_pheno, header = TRUE, sep = "\t")

hist(pheno$seq.9253.52)

path_susie <- "/scratch/dariush.ghasemi/projects/pqtl_susie/"

path_pgen <- glue(path_susie, "results/believe_nomir/tmp/", my_locuseq, "_dosage.pgen")
path_psam <- gsub(".pgen", ".psam", path_pgen)
path_pvar <- gsub(".pgen", ".pvar", path_pgen)

# Read genotype files
pvar <- pgenlibr::NewPvar(path_pvar)
pgen    <- pgenlibr::NewPgen(path_pgen, pvar = pvar)
psam_df <- read.delim(path_psam, header = TRUE, comment.char = "")
pvar_df <- read.delim(path_pvar, header = TRUE, comment.char = "")


# 2. Get the indices for specific variants (e.g., rsIDs or positions)
variant_indices <- GetVariantsById(pvar, c("rs123456", "rs7891011"))

# 4. Clean up connections
pgenlibr::ClosePgen(pgen)
pgenlibr::ClosePvar(pvar)


# Count number of variants & samples in raw data
n_variants <- nrow(pvar_df) # Or → pgenlibr::GetVariantCt(pgen)
n_samples  <- nrow(psam_df) # Or → pgenlibr::GetRawSampleCt(pgen)

# Extract dosages for all of variants
# meanimpute: if true, missing values are mean-imputed instead of NA.
dosage <- pgenlibr::ReadList(pgen, 1:n_variants, meanimpute = FALSE)

# Add variant IDs as column names
colnames(dosage) <- pvar_df$ID

# Add sample IDs as row names
#rownames(dosage) <- psam_df$IID


# Define sample IDs colomn
dosage <- as.data.frame(dosage)
dosage$IID <- psam_df$IID



geno_seqid <- pheno %>%
  select(all_of(c("IID", my_seqid))) %>%
  left_join(
    dosage, join_by(IID)
    )


X <- geno_seqid %>% select(- IID, - seq.9253.52_res) %>% as.matrix()
Y <- geno_seqid$seq.9253.52_res

# Run SuSiE 
res_susie_ild <- tryCatch(
  susie(X, y = Y, max_iter = 1000, min_abs_corr = 0.05),
  error = err_handling
)





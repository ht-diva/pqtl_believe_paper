#========================#
# This script is written to find missing in
# the genotype data of INTERVAL hg38.
# Author: Dariush Ghasemi
# Date: August 19, 2026
#========================#

library(vcfR)
library(pgenlibr)

# Full VCF file CHR21
vcf_path <- "/processing_data/shared_datasets/plasma_proteome/interval/genotypes/b38_TOPMED_imputed/filtered/chr_21_interval.vcf.gz"

# Reading entire VCF took >7 hours -- not worth it!
#int_chr21 <- vcfR::read.vcfR(vcf_path)

# PGEN files for subset region in CHR21
path_pvar <- "/scratch/dariush.ghasemi/projects/pqtl_believe_paper/chr21_35467738_36447522.pvar"
path_pgen <- str_replace_all(path_pvar, ".pvar", ".pgen")
path_psam <- str_replace_all(path_pvar, ".pvar", ".psam")

# Import genotype
pvar <- pgenlibr::NewPvar(path_pvar)
pgen <- pgenlibr::NewPgen(path_pgen, pvar=pvar)

pvar_df <- fread(path_pvar)
psam_df <- fread(path_psam)

n_variants <- pgenlibr::GetVariantCt(pgen)
n_samples  <- pgenlibr::GetRawSampleCt(pgen)

# Extract dosages for all of variants
dosage <- pgenlibr::ReadList(pgen, 1:n_variants, meanimpute = FALSE)

# Add variant IDs as column names and sample IDs as row names
colnames(dosage) <- pvar_df$ID
rownames(dosage) <- psam_df$IID


# quantify missing rate for each variant in the PGEN
na_rate <- sapply(
  colnames(dosage), 
  function (col) sum(is.na(dosage[, col]))/ nrow(dosage)
  )

# Prepare results for merge
df_missing <- as.data.frame(na_rate) %>%
  rownames_to_column("ID") %>%
  mutate(has_na = na_rate != 0)

# Merge PVAR with missing rates to compare MAF
pvar_df %>%
  mutate(
    MAF = str_split_fixed(INFO, ";", 4)[,2],
    MAF = str_remove(MAF, "MAF=") %>% as.numeric()
    ) %>%
  left_join(df_missing, join_by(ID)) %>%
  #count(has_na)
  summarize(
    MAF_min = min(MAF),
    MAF_median = median(MAF),
    MAF_max = max(MAF),
    .by = has_na
    ) #%>%
  # ggplot(aes(x = MAF, y = na_rate)) +
  # geom_point(size = 3, shape = 21)


#=============================#
# June 4, 2026
# Dariush Ghasemi

# This script compares min MAF and 
# lambda index for BELIEVE's pQTL.
#=============================#

library(tidyverse)


path_susie <- "/scratch/dariush.ghasemi/projects/pqtl_susie/"

# Believe
headers <- c("CHR", "POS", "SNPID", "EA", "NEA", "EAF", "BETA", "SE", "P", "MLOG10P", "Z")



#------------------------------#

# Check MAF in each GWAS subset
find_maf <- function(seqid, locus){
    
  locuseq <- paste0(seqid, "_", gsub("chr", "", locus))
  
  gwas_file <- glue(path_susie, "results/believe_estVarF/tmp/", locuseq, "_sumstat.csv")
  gwas  <- fread(gwas_file, header = FALSE, col.names = headers, sep = "\t", data.table = FALSE)
  
  #gwas$seqid <- str_split_fixed(locuseq, "_", 2)[1]
  #gwas$locus <- str_split_fixed(locuseq, "_", 2)[2] %>% paste0("chr", .)
  
  # Flag if any variant passes threshold
  gwas$MAF <- pmin(gwas$EAF, 1 - gwas$EAF)
  
  # Maximum value
  min_maf <- min(gwas$MAF, na.rm = TRUE)
  
  # Accounts for multiple rows
  gwas[gwas$MAF == min_maf, ]
  
}

# test function if works
find_maf("seq.23690.228", "chr5_177388765_177415473")
find_maf("seq.6364.7_12_6288176_7250542")

#------------------------------#
# Append min MAF to each locus in the report file
cs_report_maf <- cs_report %>%
  #slice_min(lambda, n=5) %>%
  filter(lambda > 0.09 | lambda < 5.4e-4) %>%  # 203 loci with extreme lambda
  mutate(locuseq = paste0(seqid, "_", gsub("chr", "", locus))) %>%
  nest_by(locuseq) %>%
  mutate(
    rare_snps = list(find_maf(data$seqid, data$locus)),
    maf = min(rare_snps$MAF)
    ) %>%
  ungroup() %>%
  unnest(data)

# Correlation: r = -0.55, p=2.2e-16
cor.test(cs_report_maf$lambda, cs_report_maf$maf)

#------------------------------#
# Scatter plot: λ vs. MAF
cs_report_maf %>%
  ggplot(aes(x = lambda, y = maf)) +
  geom_point(size = 3, shape = 21, fill = "red", color = "grey20", alpha = .6) +
  scale_x_continuous(breaks = seq(0, 0.41, 0.05)) +
  scale_y_continuous(breaks = seq(0.0001, 0.0005, 0.00005)) +
  theme_light() + 
  labs(
    y = "Minimum of Region MAF",
    x = "Extreme tails for λ at 203 loci"
    )


ggsave(filename = "04-Jun-26_extreme_lambda_values_200loci.png",
       width = 8, height = 5.5, dpi = 300, units = "in")

#------------------------------#
# Sent report/plot to Claudia G at the end of the day.

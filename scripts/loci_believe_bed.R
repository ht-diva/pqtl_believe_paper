
library(tidyverse)
library(data.table)

# 8904 loci excluding HLA
path_lb_annot <- "/exchange/healthds/pQTL/BELIEVE/Working_shared/LB/mapped_LB_gp_ann_va_ann_bl_ann_collapsed_hf_ann_wo_HLA.csv"

# Started with all 8,904 Believe pQTLs (seqid_locus pairs), then:
#   - omitted the duplicate regions (likely for aptamers with similar biological function)
#   - converted staring position of the region to 0-based (start - 1)
#   - sorted the regions based on both chromosomal and 0-based starting position of the regions
#   - combined the 4,958 pQTL regions into 768 non-overlapping windows using bedtools merge -i <myfile.bed>
#   - extracted multiple regions in a chromosome from VCF file via bcftools view -R <merged.bed>. 


#----------------------------# 
#         BELIEVE pQTLs
#----------------------------# 

lb_believe_annot %>%
  mutate(beg = start - 1) %>% # 0-based to satisfy bed format
  distinct(chr, beg, end) %>% # unique regions --> ended up with 4,958 regions
  arrange(chr, beg) %>%
  #filter(chr%in% 21:22) %>% # for test
  mutate(chr = str_c("chr", chr)) %>%
  write.table(
    "believe/believe_loci_unique_100k_extended.bed", 
    quote = F, row.names = F, sep = "\t", col.names = FALSE
    )


#----------------------------# 
#         Enlarge pQTLs
#----------------------------# 

# Will be added to each tail
buffer <- 100000

lb_believe_annot %>%
  # extend locus boundaries: +/- buffer size
  dplyr::mutate(
    zero_proximate = start - buffer <= 1, # 12 loci being affected
    beg_ext = ifelse(start - buffer - 1 < 1, 0, start - buffer), # (-1) 0-based it to satisfy bed format
    end_ext = end + buffer
    ) %>% 
  # Count any affected loci
  #filter(zero_proximate) %>% select(locus:loci_cat) %>% distinct(locus)
  distinct(chr, beg_ext, end_ext) %>% # 4,958 unique regions
  arrange(chr, beg_ext) %>%
  mutate(chr = str_c("chr", chr)) %>%
  write.table(
    "believe/21-May-26_believe_loci_unique_100kb_extended.bed", 
    quote = F, row.names = F, sep = "\t", col.names = FALSE
    )


#----------------------------# 
#         Merge pQTLs
#----------------------------# 

# ***** Bash commands to merge the overlapping loci  *****

# Ended up with below non-overlapping, merged pQTL regions:
#  - without extension: 768
#  - with buffer=100k: 653

#bedtools merge -i  ~/believe/believe_loci_unique.bed > ~/believe/believe_loci_merged.bed
#bedtools merge -i  ~/believe/21-May-26_believe_loci_unique_100kb_extended.bed > ~/believe/believe_loci_merged_extended.bed


# Explore merged region
lb_believe_merged <- fread("~/believe/believe_loci_merged_extended.bed",
      col.names = c("chr", "start", "end")) %>%
  chop_locus() #%>% count(loci_cat)

#      loci_cat     n
# 1:     1-2Mbp   147
# 2: 100-500Kbp   178
# 3:     2-5Mbp   115
# 4:  500K-1Mbp   170
# 5:      >5Mbp    43


summary(lb_believe_merged$loci_width)

# Min: 200 kb
# Max: 12 Mb
# Q1: 465 kb
# Q2: 919 kb
# Q3: 1.9 Mb


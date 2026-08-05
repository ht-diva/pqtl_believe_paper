
library(tidyverse)
library(data.table)

path_lb_pipe <- "/scratch/dariush.ghasemi/projects/pqtl_conditional/results/believe_hla/combined_loci.csv"
#path_lb_annot <- "/exchange/healthds/pQTL/BELIEVE/Working_shared/LB/mapped_LB_gp_ann_va_ann_bl_ann_collapsed_hf_ann.csv"
path_lb_woHLA <- "/exchange/healthds/pQTL/BELIEVE/Working_shared/LB/mapped_LB_gp_ann_va_ann_bl_ann_collapsed_hf_ann_wo_HLA.csv"

path_lb_out <- "/scratch/dariush.ghasemi/projects/pqtl_coloc/config/19-Nov-25_believe_loci.csv"
path_lb_susie <- "/scratch/dariush.ghasemi/projects/pqtl_susie/config/believe_loci.csv"
path_lb_out_cis <- "/scratch/dariush.ghasemi/projects/pqtl_coloc/config/believe_loci_cis.csv"
path_lb_out_large <- "/scratch/dariush.ghasemi/projects/pqtl_locuszoom/conf/believe_loci_large.csv"
path_lb_test_gnh <- "believe_loci_test_gnh.csv"

plt_hist <- "26-Feb-26_histogram_loci_width_believe.png"
plt_hla <- "26-Feb-26_histogram_loci_width_believe_HLA_comparison.png"
plt_vioin <- "19-Nov-25_density_loci_width_believe.png"

#-------------------------------#
# -----    Locus Breaker   -----
#-------------------------------#

lb_believe <- fread(path_lb_pipe)
lb_believe_annot <- fread(path_lb_woHLA) %>% arrange(chr, phenotype_id)


chop_locus <- function(df){
  
  df %>%
    dplyr::rename(seqid = phenotype_id) %>%
    mutate(
      locus = str_c("chr", chr, "_", start, "_", end),
      loci_width = end - start,
      loci_cat = case_when(
        loci_width == 0 ~ "1-SNP",
        loci_width > 0       & loci_width <= 100000  ~ "1bp-100Kbp",
        loci_width > 100000  & loci_width <= 500000  ~ "100-500Kbp",
        loci_width > 500000  & loci_width <= 1000000 ~ "500K-1Mbp",
        loci_width > 1000000 & loci_width <= 2000000 ~ "1-2Mbp",
        loci_width > 2000000 & loci_width <= 5000000 ~ "2-5Mbp",
        TRUE ~ ">5Mbp"
      )
    )
}


# Categorize loci
lb_believe <- lb_believe %>% chop_locus()
lb_believe_annot <- lb_believe_annot %>% chop_locus()


write.csv(lb_believe, file = path_lb_out, quote = F, row.names = F)

# to run SuSiE
lb_believe_annot %>% 
  select(chr:cis_or_trans, loci_width, loci_cat) %>%
  write.csv(file = path_lb_susie, quote = F, row.names = F)
  

# cis vs. trans
lb_believe %>%
  count(loci_cat, cis_or_trans) %>%
  spread(cis_or_trans, n) %>% DT::datatable()

# save cis loci -- last save on 20-Apr-26
lb_believe_annot %>%
  dplyr::filter(cis_or_trans == "cis") %>%
  dplyr::select(chr:cis_or_trans) %>%
  chop_locus() %>% 
  write.csv(file = path_lb_out_cis, quote = F, row.names = F)

# to test in Genes&Health TRE
lb_believe %>%
  filter(
    chr %in% 21:22, 
    !loci_cat %in% c("1-SNP", "1bp - 100Kbp")
    ) %>% #count(loci_cat)
  select(chr:end, SNPID, seqid, loci_cat) %>%
  write.csv(path_lb_test_gnh, row.names = F)


# Unique cis pQTLs to compute LD in G&H TRE
lb_believe_annot %>% 
  select(chr:end, cis_or_trans, loci_cat) %>%
  distinct() %>%  # count(loci_cat)  # Ended with 4958 unique regions
  filter(cis_or_trans == "cis") %>%
  # filter(loci_cat=="1-SNP")        # No cis 1-SNP region exist
  fwrite("believe/believe_loci_uniq_cis_4TRE.tsv", row.names = F, sep = "\t")


lb_believe_4TRE <- list.files(
  "/scratch/dariush.ghasemi/projects/pqtl_susie/ld/ld_matrix/",
  pattern = "_ld.matrix"
  ) %>%
  data_frame("ldfile" = .) %>%
  mutate(locus = str_c("chr", str_remove(ldfile, "_ld.matrix")))


lb_believe_annot %>% 
  select(chr:cis_or_trans, locus, loci_cat) %>%
  filter(cis_or_trans == "cis",
         locus %in% lb_believe_4TRE$locus) %>%
  write.csv("/scratch/dariush.ghasemi/projects/pqtl_susie/config/gnh_93.csv",
            quote = F, row.names = F)


#-------------------------------#
# -----   Loci Histogram   -----
#-------------------------------#

hist(lb_believe$loci_width, nclass = 50)
lb_believe %>% count(loci_cat)

# Size of the largest locus in Believe
mx_lb <- round(max(lb_believe$loci_width)/10^6, 1)

level_order <- c('1-SNP', '1bp-100Kbp', '100-500Kbp', '500K-1Mbp', '1-2Mbp', '2-5Mbp', '>5Mbp')

lb_believe %>%
  ggplot(aes(factor(loci_cat, levels = level_order))) +
  geom_histogram(stat = "count", color = "#003421", fill = "#095F54") +
  stat_count(
    binwidth = 1, 
    geom = 'text', 
    color = '#E69F00', 
    #aes(label = after_stat(count)),
    #position = position_nudge(vjust = 0),
    aes(label=..count.., y=..count.. + 130)
  ) +
  scale_y_continuous(breaks = seq(0, 5000, 500), limits = c(0, 5000)) +
  labs(
    x = paste0("\nLocus Size (max: ", mx_lb, " Mbp)"),
    y = "#Loci"
  ) +
  ggtitle("8,910 Loci in Believe Study (excluding HLA)") +
  theme_light() +
  theme(
    axis.ticks.length = unit(2, "mm"),
    axis.title  = element_text(size = 13, face = 1),
    axis.text.x = element_text(size = 11, face = 1)
  )


ggsave(filename = plt_hist, width = 8, height = 5.5, dpi = 300, units = "in")

#-------------------------------#
# -----   HLA Comparison   -----
#-------------------------------#

# Bar Plot: Annotated with frequency and percentage
lb_believe %>%
  dplyr::mutate(version = "With HLA (n=9,679)") %>%
  select(version, loci_cat) %>%
  rbind(
    lb_believe_annot %>%
      dplyr::mutate(version = "Without HLA (n=8,910)") %>%
      select(version, loci_cat)
  ) %>%
  summarize(n = n(), .by = c("version", "loci_cat")) %>%
  dplyr::mutate(
    prop = n/sum(n),
    percent = str_c(n, " (", round(prop*100, 2), "%)"),
    .by = "version"
    ) %>%
  ggplot(aes(x = factor(loci_cat, levels = level_order),
             y = n, fill = version)) +
  #geom_bar(position = position_dodge(0.9), stat = "count") +
  geom_bar(stat = "identity", position = "dodge") +
  geom_text(
    aes(label = percent),
    position = position_dodge(1),
    size = 3.5, angle = 30,
    vjust = -1, hjust = .25
    ) +
  scale_fill_manual(values = c("#E69F00", "#095F54")) +
  scale_y_continuous(breaks = seq(0, 5000, 500), limits = c(0, 5000)) +
  labs(x = paste0("\nLocus Size"),
       y = "#Loci",
       fill = "LB version"
       ) +
  theme_light() +
  theme(
    legend.position = c(.85, .85),
    legend.background = element_blank(),
    legend.title = element_text(size = 13, face = 2),
    legend.text = element_text(size = 12, face = 3),
    axis.ticks.length = unit(2.5,"mm")
  )


ggsave(filename = plt_hla, width = 8, height = 6, dpi = 300, units = "in")


#-------------------------------#
# -----   MAF Comparison   -----
#-------------------------------#

# violin plot
lb_believe %>%
  mutate(MAF = ifelse(EAF < 0.5, EAF, 1-EAF),
         locus = factor(loci_cat, levels = level_order)
  ) %>%
  ggplot() +
  #geom_density(aes(x = EAF, fill = loci_cat))
  geom_violin(aes(y = MAF, x = locus, color = locus), fill = NA) +
  labs(x = "Width of loci", y = "MAF of index SNP", color = "Locus") +
  theme_dark() +
  theme(
    axis.title  = element_text(size = 12, face = 2),
    axis.text.x = element_text(size = 8, face = 2)
  )

ggsave(filename = plt_violin, last_plot(), width = 8.5, height = 5.5, dpi = 300, units = "in")

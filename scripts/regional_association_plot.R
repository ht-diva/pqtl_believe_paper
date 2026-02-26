#===============================#
# Drawing regional association
# plot for large loci in Believe
# February 26, 2026
# Dariush Ghasemi
#===============================#


path_gwas <- "/exchange/healthds/pQTL/BELIEVE/harmonized_gwas/seq.14151.4/seq.14151.4.gwaslab.tsv.gz"
path_sumstat <- "/scratch/dariush.ghasemi/projects/pqtl_susie/results/believe_with_multiallelic/tmp/"


#outputs
regional_plots <- "/home/dariush.ghasemi/26-Feb-26_Believe_regional_association_plots_large_loci.pdf"

#-------------------------------#
# -----      LocusZoom     -----
#-------------------------------#

# large loci for LocusZoom pipe
lb_believe_large <- lb_believe %>%
  filter(loci_cat %in% c("1-2Mbp", "2-5Mbp", ">5Mbp")) %>%
  group_by(loci_cat) %>%
  slice_max(loci_width, n = 5) %>% 
  ungroup() #distinct(loci_width)


write.csv(x = lb_believe_large, file = path_lb_out_large, quote = F, row.names = F)


# Regional Association Plots
headers <- c("CHR", "POS", "SNPID", "EA", "NEA", "EAF", "BETA", "SE", "P", "MLOG10P", "Z")


pdf(regional_plots, height = 5.5, width = 7.5)

lapply(
  lb_believe_large$seqid, function(seqname) {
    
    # Extract loci characteristics from LB file
    my_locus <- lb_believe_large %>% filter(seqid == seqname)
    # Define locuseq
    seqid_locus <- str_c(my_locus$seqid, "_", my_locus$chr, "_", my_locus$start, "_", my_locus$end)
    # Create approperiate path to merginal GWAS exist in SuSiE pipe
    path_df <- glue(path_sumstat, seqid_locus, "_sumstat.csv")
    # Read marginal GWAS sumstat for Believe study
    df_sumstat <- fread(path_df, header = F, col.names = headers)
    
    # Regional plot
    df_sumstat %>%
      ggplot(aes(x = POS, y = MLOG10P)) +
      geom_point(size = 3, fill = "brown4", shape = 21) +
      labs(x = "Genomic Position (hg38)") + 
      ggtitle(paste0(seqid_locus, "(size: ", my_locus$loci_cat,")")) + 
      theme_light() +
      theme(
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        axis.ticks.length = unit(2.5,"mm")
      )
    
  })

dev.off()


# example regional plot
path_df <- glue(path_sumstat, "seq.8398.277_3_157527732_159526583_sumstat.csv")
df_sumstat <- fread(path_df, header = F, col.names = headers)

df_sumstat %>%
  ggplot(aes(x = POS, y = MLOG10P)) +
  geom_point(size = 3, fill = "brown4", shape = 21) +
  labs(x = "Genomic Position (hg38)") + 
  theme_light() +
  theme(
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    axis.ticks.length = unit(2.5,"mm")
  )




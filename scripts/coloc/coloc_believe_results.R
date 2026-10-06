
library(glue)
library(purrr)

path_coloc <- "/scratch/dariush.ghasemi/projects/pqtl_coloc/results"
path_res_binary <- glue(path_coloc, "/believe_cis_binary/combined_coloc_results.tsv")
path_res_quant  <- glue(path_coloc, "/believe_cis_quant/combined_coloc_results.tsv")

#-------------------------------#
# -----      Read data     -----
#-------------------------------#

# # Files from coloc pipe
# files_list <- list.files(
#   path = path_coloc,
#   pattern = ".csv",
#   full.names = TRUE
# )
# 
# # Merge files
# coloc_res_dic <- map_dfr(files_list, fread)


# Read results
coloc_res_dic   <- fread(path_res_binary)
coloc_res_quant <- fread(path_res_quant)


# New cis proteins from Meta-analysis
new_seqids <- fread("/exchange/healthds/pQTL/BELIEVE/Downstream_LB_analysis/subset_new_cis_loci.txt") %>%
  pull(phenotype_id)



#-------------------------------#
# -----    Check results   -----
#-------------------------------#

# coloc_signif <- combined_results %>%
#   inner_join(
#     lb_believe %>%
#       filter(cis_or_trans == "cis") %>%
#       mutate(locus = str_c(seqid, chr, start, end, sep = "_")),
#     join_by(locus)
#   )

# Sanity checks

# No. unique phenotypes
n_distinct(coloc_res_quant$phenotype)
n_distinct(coloc_res_dic$phenotype)


# No. unique seqid-locus
n_distinct(coloc_res_quant$seqid, coloc_res_quant$locus)
n_distinct(coloc_res_dic$seqid, coloc_res_dic$locus)

# No. of Tests
nrow(coloc_res_quant)
nrow(coloc_res_dic)

# No. of Tests with PP.H4 > 0.8
coloc_res_quant[coloc_res_quant$PP.H4.abf > 0.8] %>% nrow()
coloc_res_dic[coloc_res_dic$PP.H4.abf > 0.8] %>% nrow()

# Percentage
prop.table(table(coloc_res_quant$PP.H4.abf > 0.8)) %>% round(4)
prop.table(table(coloc_res_dic$PP.H4.abf > 0.8)) %>% round(4)

# Colocalized traits
#coloc_res_dic %>%
coloc_res_quant %>%
  dplyr::filter(PP.H4.abf > 0.8) %>%
  #count(seqid, locus) %>% arrange(n)
  #count(phenotype) %>% arrange(n)
  #distinct(seqid, locus)
  #add_count(phenotype) %>%
  summarize(npheno = n(), .by = phenotype) %>%
  dplyr::filter(npheno > 2)


coloc_res_quant %>%
  dplyr::filter(
    PP.H4.abf > 0.8,
    seqid %in% new_seqids
    ) %>%
  dplyr::select(seqid:nsnps, PP.H4.abf) %>%
  write.csv(quote = F)


coloc_res_quant %>%
  #dplyr::filter(PP.H4.abf > 0.8) %>% 
  View()


#-------------------------------#
# -----  Visualize results -----
#-------------------------------#

coloc_res_quant %>%
  ggplot(aes(y = nsnps))+ geom_violin()
  geom_boxplot() +
  #geom_histogram(fill = "steelblue", color = "white", binwidth = 1000) +
  scale_y_continuous(breaks = seq(0,70000,10000)) +
  theme_light() +
  theme(axis.title = element_text(size =12, face = 2),
        axis.text = element_text(size =12))


coloc_res_quant %>%
  ggplot(aes(PP.H4.abf))+
  geom_histogram(fill = "steelblue", color = "white", binwidth = 0.015) +
  scale_x_continuous(breaks = seq(0,1,0.2)) +
  geom_vline(xintercept = 0.8, linetype = 2) +
  theme_light() +
  theme(axis.title = element_text(size =12, face = 2),
        axis.text = element_text(size =12))

ggsave("08-Jan-26_believe_coloc_cis_pp4.png", width = 8, height = 5.5, dpi = 300)

# histogram of PP.H0-4
coloc_res_quant %>%
  pivot_longer(
    cols = PP.H0.abf:PP.H4.abf,
    names_to = "Hypothesis",
    values_to = "value"
  ) %>%
  ggplot(aes(value))+
  geom_histogram(fill = "steelblue", color = "white", binwidth = 0.035) +
  facet_wrap(~Hypothesis, scales = "free") +
  theme_light() +
  theme(
    axis.title.x = element_blank(),
    strip.placement = "outside",
    strip.background = element_blank(),
    strip.text.x = element_text(size = 12, color = "Black", face = 4),
    panel.border = element_rect(colour = "black", fill = NA),
  )

ggsave("08-Jan-26_believe_coloc_cis_pps.png", width = 12, height = 8, dpi = 300)

# scatter
coloc_res_quant %>% #count(phenotype)
  dplyr::filter(PP.H4.abf > 0.8) %>%
  mutate(pheno_sig = ifelse(PP.H4.abf > 0.8, phenotype, "")) %>%
  ggplot(aes(x = phenotype,  y = PP.H4.abf, color = protein))+
  geom_point(show.legend = F) +
  theme(axis.text.x = element_text(angle =-90, vjust = .5, hjust = 0))

# scatter
coloc_res_quant %>%
  #ggplot(aes(x = nsnps, y = PP.H4.abf)) +
  ggplot(aes(x = PP.H3.abf, y = PP.H4.abf)) +
  geom_point(alpha = 0.6) +
  theme_light() +
  theme(axis.title = element_text(size =12, face = 2),
        axis.text = element_text(size =12))

ggsave("08-Jan-26_believe_coloc_cis_pp3n4.png", width = 8, height = 5.5, dpi = 300)

# Colocalized traits
coloc_res_dic %>%
  #coloc_res_quant %>%
  dplyr::filter(PP.H4.abf > 0.8) %>%
  #filter(seqid == "seq.20091.138")
  #count(seqid, locus) %>% arrange(n)
  #count(phenotype) %>% arrange(n)
  #distinct(seqid, locus)
  #add_count(phenotype) %>%
  summarize(npheno = n(), .by = phenotype) %>% 
  #filter(npheno > 2) %>%
  ggplot(aes(y = fct_reorder(phenotype, npheno), x = npheno))+
  geom_col() +
  labs(
    x = "#Unique seqid-locus colocalized with each trait",
    y = "Phenotype"
  ) #+ theme_minimal()

ggsave("30-Apr-26_believe_coloc_cis_signif_quant.png",
       width = 8, height = 13, dpi = 300)




library(forcats)

# Heatmap
coloc_res_quant %>%
  dplyr::filter(PP.H4.abf > 0.8) %>%
  mutate(ncoloc_pheno = n_distinct(seqid),     .by = "phenotype") %>%
  mutate(ncoloc_prot  = n_distinct(phenotype), .by = "seqid") %>%
  mutate(seqid_prot = str_c(protein," (", seqid, ")")) %>% 
  ggplot(aes(x = fct_reorder(phenotype, ncoloc_pheno),
             y = fct_reorder(seqid_prot, ncoloc_prot),
             fill = PP.H4.abf)) +
  geom_tile(color = "white") +
  scale_x_discrete(position = "top") +
  scale_y_discrete(position = "right") +
  #scale_fill_viridis_c(name = "PP.H4", option = "H") +
  scale_fill_gradient(low = "gold2", high = "red") +
  labs(
    title = "Colocalized pQTLs with Blood Biomarkers",
    x = "Blood biomarker",
    y = "Protein"
  )+
  theme_light() +
  theme(
    #plot.margin = margin(r = 10),
    legend.position = c(1.35, 1.07),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_text(size=9, face=1, angle = 90, hjust = 0),
    axis.text.y = element_text(size=8, face=1, hjust = 0)
  )

ggsave("27-Apr-26_believe_coloc_signif_cis_quant.png", width = 16, height = 21, dpi = 300)

# Dot plot
combined_results %>%
  dplyr::filter(PP.H4.abf > 0.8) %>%
  ggplot(aes(x = phenotype,
             y = protein,
             color = PP.H4.abf,
             size = nsnps)) +
  geom_point(alpha = 0.8) +
  scale_color_viridis_c(name = "PP.H4") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Protein–Biomarker Colocalization",
       x = "Blood biomarker",
       y = "Protein")


# pretty heatmap
coloc_signif_long <- coloc_res_dic %>%
  dplyr::filter(PP.H4.abf > 0.8) %>%
  # pivot_wider(
  #   id_cols = locus, 
  #   names_from = phenotype,
  #   values_from = PP.H4.abf
  #   ) %>%
  # mutate(across(is.na(), ~ replace_na(0))) %>%
  select(phenotype, PP.H4.abf, locus) %>%
  spread(phenotype, PP.H4.abf, fill = 0) %>%
  tibble::column_to_rownames(var = "locus")



library(dplyr)
library(tidyr)

# Filter proteins & phenotypes with ≥1 significant hit
sig_pairs <- coloc_res_dic %>% filter(PP.H4.abf > 0.8)

coloc_signif_long <- coloc_res_dic %>%
  filter(
    protein %in% sig_pairs$protein,
    phenotype %in% sig_pairs$phenotype
  ) %>%
  group_by(protein, phenotype) %>%
  summarise(PP.H4.abf = max(PP.H4.abf), .groups = "drop") %>%
  pivot_wider(
    names_from = phenotype,
    values_from = PP.H4.abf,
    #values_fill = 0
  ) %>%
  tibble::column_to_rownames("protein") %>%
  as.matrix()

# Emphasize strong signals
coloc_signif_plot <- coloc_signif_long
coloc_signif_plot[coloc_signif_plot < 0.7] <- 0

# Or binary map
coloc_signif_bin <- ifelse(coloc_signif_long > 0.8, 1, 0)

png("23-Apr-26_believe_coloc_signif_cis_bin_full.png", 
    width = 21, height = 18, units = "in", res = 300)


library(pheatmap)

pheatmap(
  coloc_signif_long[,-1],
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  #show_rownames = FALSE,
  #show_colnames = FALSE,
  color = colorRampPalette(c("white", "red"))(100)
)

dev.off()

coloc_res %>%
  dplyr::filter(PP.H4.abf > 0.8) %>% #distinct(protein, phenotype)
  group_by(protein) %>%
  arrange(phenotype) %>% 
  summarise(
    n_traits = n(),
    traits = paste0(unique(phenotype), collapse = ", ")
  ) %>% View()




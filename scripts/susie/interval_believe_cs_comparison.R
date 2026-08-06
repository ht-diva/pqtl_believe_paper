
path_lifted_intl_chr22 <- "/scratch/dariush.ghasemi/projects/pqtl_liftover/results/chr22_harm_renormalized/VCF_lifted/lifted_variants.txt"

lifted_intl <- fread(path_lifted_intl_chr22, header = F,
                     col.names = c("chr", "pos_b38", "id_b37"))

# Note:
# In the lifted file, we may have a genomic cite mapped to 
# a different chromosome --> table(df$chr) shows this.

genomic_key <- lifted_intl %>%
  dplyr::filter(chr == "chr22") %>%
  mutate(
    pos_b37 = str_split_fixed(id_b37, ":", 4)[,2],
    chr_pos = str_c(str_remove(chr, "chr"), pos_b37, sep = ":") 
  )

# sanity checks
genomic_key %>% 
  mutate(nb37 = n_distinct(pos_b37), .by = pos_b38) %>% 
  filter(nb37 > 1)


n_distinct(genomic_key$pos_b38)



lb_int_annot %>%
  dplyr::filter(chr == 22) %>%
  mutate(
    chr_beg = str_c(chr, start, sep = ":"),
    chr_end = str_c(chr, end,   sep = ":"),
  ) %>%
  left_join(genomic_key[, c(1,3)], join_by(chr_end == chr_pos)) %>%
  dplyr::rename(start_b38 = pos_b38) %>%
  left_join(genomic_key[, c(1,3)], join_by(chr_beg == chr_pos), suffix = c("_beg", "_end"))
            





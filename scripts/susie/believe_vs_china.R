
path_pip_group  <- "/group/diangelantonio/users/dariush/results/pqtl_susie"
path_ckb_cslist <- "/scratch/dariush.ghasemi/projects/china_kadoori_suppl.tabels.xlsx"
path_pip_cslist <- glue(path_pip_group, "/believe_bgen/combined_credibles.tsv")


cssums_believe <- fread(path_pips_combined)

cs_believe <- fread(path_pip_cslist)

cs_china <- readxl::read_xlsx(path_ckb_cslist, skip = 1, sheet = 7) #%>% janitor::clean_names()



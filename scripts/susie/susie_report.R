#=============================#
# June 4, 2026
# Dariush Ghasemi

# Here I am modeling locus characteristics
# that influence lambda in BELIEVE study.
#=============================#

library(tidyverse)

# Base directory
path_susie <- "/scratch/dariush.ghasemi/projects/pqtl_susie/"
path_group <- "/group/diangelantonio/users/dariush/results/"

# Latest susie run on Believe's pQTLs
path_pip_report <- glue(path_group, "pqtl_susie/believe_estVarF/susierss/combined_reports.tsv")

# Report of loci with Negative Variance
path_report_negvar <-  glue(path_susie, "results/believe_negvar_89loci/susierss/combined_reports.tsv")
path_repott_issu11 <-  glue(path_susie, "results/issue_11/susierss/combined_reports.tsv")

# Read report file
cs_report        <- fread(path_pip_report)
cs_report_negvar <- fread(path_report_negvar)
cs_report_issu11 <- fread(path_repott_issu11)

list_report <- list.files(
  path = glue(path_susie, "results/believe_estVarF/susierss/cs_report"),
  pattern = "seq.(\\d)+.(\\d)+_(\\d)+_(\\d)+_(\\d)+.report$",
  full.names = T
)

# combine reports
cs_report <- map_dfr(list_report, ~ fread(.x, colClasses = "character"))


#----------------------------#
# ----       Plots       ----
#----------------------------#

# Histogram of Lambda
hist(cs_report$lambda, nclass = 100)


# Box plot + histogram 
cs_report %>%
  left_join(
    lb_believe_annot %>% select(seqid, locus, cis_or_trans, loci_width, loci_cat),
    join_by(seqid, locus)
  ) %>%
  mutate(loci_order = factor(loci_cat, levels = level_order)) %>% 
  ggplot() +
  #geom_boxplot(aes(x = lambda)) +
  geom_histogram(aes(x = lambda), color="white", bins = 50) +
  #\scale_x_continuous(breaks = seq(0, 0.41, 0.05))+
  #geom_boxplot(aes(x = loci_order, y = lambda)) +
  theme_light() +
  labs(x = "λ", y = "#Loci")

ggsave("10-Jun-26_histogram_lambda.png", width = 8, height = 7, dpi = 200)



cs_report %>% #count(nmulti_allelic)
  ggplot(aes(x = nindels_shared, y = lambda)) +
  geom_point() +
  labs(x = "# Insertion/deletion variants in the region")


# If higher lambda is due to lower P-value of index SNP

cs_report_pval <- cs_report %>%
  select(- lambda_warning) %>%
  # Add locus annotation + index SNP info
  left_join(
    lb_believe_annot %>% select(seqid, locus:loci_cat, SNPID, EAF:MLOG10P),
    join_by(seqid, locus)
  ) %>%
  mutate(
    loci_order = factor(loci_cat, levels = level_order),
    MAF = ifelse(EAF <= 0.5, EAF, 1 - EAF),
    Z = BETA/SE,
    log_lambda = log10(lambda),
    prop_indels = nindels_shared/nvar_shared
    ) %>%
  # Add CS size for each locus
  left_join(cs_believe_sum, join_by(seqid, locus))



# Scatter plot
cs_report_pval %>% 
  #ggplot(aes(x = log_lambda, y = abs(Z), color = loci_order)) +
  ggplot(aes(x = lambda, y = cs_avgr2, color = loci_order)) +
  geom_point(alpha = .85, size = 3, shape = 19) +
  scale_x_continuous(breaks = seq(0, 0.41, 0.05))+
  #scale_y_continuous(breaks = seq(0, 1500, 100))+
  scale_color_brewer(palette = "Set1") +
  labs(color = "Locus", x = "λ") +  #x=paste("log10(λ)")
  theme_light() +
  theme(
    legend.position = c(0.85, .40),
    legend.background = element_blank(),
    legend.box.background = element_rect(fill = NA, linetype = 2)
    )


ggsave("10-Jun-26_lambda_vs_csavgr2.png", width = 8, height = 5.5, dpi = 200)


#----------------------------#
# ----    Model Lambda   ----
#----------------------------#

# Pearson's correlation: 0.717
cor.test(cs_report_pval$lambda, cs_report_pval$MLOG10)

# Regression model
m_formula <- "log_lambda ~ MAF"
m_formula <- "log_lambda ~ nvar_shared"
m_formula <- "log_lambda ~ prop_indels"
m_formula <- "log_lambda ~ n_cs"
m_formula <- "log_lambda ~ n_snps"
m_formula <- "log_lambda ~ cs_power"
m_formula <- "log_lambda ~ cs_avgr2"
m_formula <- "log_lambda ~ cs_minr2"
m_formula <- "log_lambda ~ Z + loci_width"
m_formula <- "log_lambda ~ Z + loci_width + cs_avgr2 + MAF"

lm(m_formula, data = cs_report_pval) %>%
  summary()# %>% broom::tidy() %>% knitr::kable()


#----------------------------#
# ----   Report compare  ----
#----------------------------#

cs_report %>%
  mutate(est_var = "TRUE") %>%
  rbind(
    cs_report_negvar %>% mutate(est_var = "FALSE")
  ) %>%
  ggplot(aes(y = lambda, x = est_var, color = est_var)) +
  geom_jitter(show.legend = F, color = "grey30", shape = 21)+
  geom_boxplot(width = .2, show.legend = F) +
  #geom_density() +
  #guides()
  labs(x = "Estimate residual variance") +
  theme_classic()


ggsave("18-May-26_est_var_comparison.png", width = 7, height = 7, dpi = 200)


# Assess effectiveness of score option in reducing Lambda
cs_report_negvar %>%
  select(seqid, locus, lambda) %>%
  inner_join(
    cs_report_issu11 %>% select(seqid, locus, lambda),
    join_by(seqid, locus),
    suffix = c("_default", "_score")
    ) %>%
  # pivot_longer(
  #   cols = starts_with("lambda"),
  #   names_to = "method",
  #   values_to = "lambda"
  #   ) %>%
  #slice_max(lambda_score, n = 10)
  ggplot(aes(lambda_default, lambda_score)) +
  geom_point(size = 3.5, fill = "#FDC700", shape = 21) +
  geom_abline(slope = 1, lty = 2) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  labs(x = "λ with default params",
       y = "λ with score option",
       title = "Lambda accounted for Z-scores from LLM",
       subtitle = "(37 highly significant pQTLs)") +
  theme_light()

ggsave("02-Sep-26_lambda_comparison_score_option.png",
       width = 7.5, height = 6.5, dpi = 200)




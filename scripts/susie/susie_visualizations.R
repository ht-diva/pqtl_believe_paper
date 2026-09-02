
#=============================#
# September 1, 2026
# Dariush Ghasemi

# This script put together all
# the plots depicting susie results.
#=============================#

library(ggcorrplot)
library(ggplotify)
library(ggpubr)

#----------------------------#
# ----   Regional plot   ----
#----------------------------#

# Regional plot
plt_lz <- sumstat %>%
  ggplot(aes(x = POS, y = MLOG10P)) +
  geom_point(size = 3, color = "#FDC700", shape = 16, alpha = .6) +
  geom_point(data = cs_sum_plot, aes(color = cs_id), size = 4, shape = 21, stroke = 1.3) +
  scale_x_continuous(labels = function(x) round(x/1e6, 2)) +
  scale_color_discrete() +
  labs(
    title = paste("Region:", my_locuseq),
    x = "Genomic Position (Mb)"
    ) + 
  theme_light() +
  theme(
    legend.position = c(.85, .55),
    legend.background = element_blank(),
    plot.title = element_text(size = 10, face = 2, hjust = 0.5),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 9),
    axis.ticks.length = unit(2.5, "mm")
  )


#----------------------------#
# ----    Kriging plot   ----
#----------------------------#

# Compute λ
lambda <- susieR::estimate_s_rss(z = z_scores, R = R)

# Consistency between alleles in LD and sumstat
condz <- susieR::kriging_rss(z = z_scores, R = R)

# Plot + lambda
plt_kriging <- condz$plot + 
  labs(
    #title = paste("SeqID-Locus= ", my_locuseq),
    title = paste0("LD matrix consistency with Z-scores (λ=", signif(lambda, 4), ")")
  ) +
  theme_light() +
  theme(
    plot.title = element_text(size = 10, face = 2, hjust = 0.5),
    axis.title = element_text(size = 12),
    panel.background = element_blank()
  )


#----------------------------#
# ----   SuSiE plot   ----
#----------------------------#

# SuSiE Plot
plt_cs <- ggplotify::as.ggplot(function() {
  
  par(
    mar = c(7,5,3,3),
    bg = "white",
    fg = "black"
  )
  
  susieR::susie_plot(
    res_rss_score,
    y = "PIP",
    b = betas,
    xlab = "Variants",
    add_bar = FALSE,
    add_legend = TRUE,
    cex.axis = .75
  )
  
})


plt_cs <- plt_cs +
  labs(
    title = "Posterior Inclusion Probablity of variants",
    #subtitle = paste("SeqID_locus:", my_locuseq)
  ) +
  theme(
    plot.title = element_text(size = 10, face = 2, hjust = 0.5),
    plot.background  = element_rect(fill = "white", colour = NA),
    panel.background = element_rect(fill = "white", colour = NA)
  )


#----------------------------#
# ----     CS LD plot    ----
#----------------------------#

# CS LD Correlation
cs_ld <- susieR::get_cs_correlation(res_rss_score, Xcorr = R)

# CS LD Corr. plot
plt_cs_ld <- ggcorrplot(
  cs_ld,
  method = "square",
  type = "lower",
  lab = TRUE,
  tl.cex = 10,
  lab_size = 2,
  digits = 2,
  colors = c("#3B4CC0", "white", "#B40426"),
  outline.color = "white",
  ggtheme = "theme_minimal"
) +
  labs(
    title = "LD correlation between credible sets",
    x = NULL,
    y = NULL
  ) +
  theme(
    plot.title = element_text(size = 10, face = 2, hjust = 0.5),
    panel.background = element_blank()
  )


#----------------------------#
# ----     Joint plot    ----
#----------------------------#

# Combine plots
ggpubr::ggarrange(
  plt_lz,
  plt_kriging,
  plt_cs,
  plt_cs_ld,
  nrow = 2,
  ncol = 2,
  #align = "hv",
  heights = c(1, 1)
)


ggsave("31-Aug-26_combined_susie_plots.png", 
       width = 9, height = 7.75, dpi = 200)


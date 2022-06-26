# Make figure 4

# RIL effect on all polymorphs
fig4 <- caco3_comp %>%
  modelr::data_grid(log_weight = mean(log_weight),
                    log_ril =  modelr::seq_range(log_ril, n = 101),
                    mean_T = mean(mean_T),
                    sqrt_asp_ratio = mean(sqrt_asp_ratio),
                    method = "double") %>%
  tidybayes::add_epred_draws(m3_caco3_comp, resp = resp, re_formula = NA) %>%
  mutate(mineral = factor(.category,
                          levels = paste0(c("L", "AR", "H", "M", "AC"), "umolh"),
                          labels = c("LMC", "Aragonite", "HMC", "MHC", "ACMC")
                          )
         ) %>%
  ggplot(aes(x = exp(log_ril))) +
  tidybayes::stat_lineribbon(aes(y = .epred, color = mineral), size = 0.5, .width = c(0.5, 0.8, 0.95)) +
  scale_color_viridis_d(option = "C", end = 0.95, guide = "none") +
  scale_fill_brewer("CI", palette = "Greys") +
  facet_grid(rows = "mineral", scales = "free_y") +
  ylab("Excretion rate (&mu;mol h<sup>-1</sup>)") +
  xlab("Relative intestinal length") +
  theme(axis.title.y = ggtext::element_markdown(),
        legend.position = "bottom",
        legend.title = element_text(size = 11),
        legend.text = element_text(size = 10),
        legend.key.size = unit(4, "mm"),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold", size = 12))

ggsave(here::here("outputs", "figures", "fig4.png"), fig4,
       width = 10, height = 24, units = "cm", dpi = 600, type = "cairo")
ggsave(here::here("outputs", "figures", "fig4.pdf"), fig4,
         width = 10, height = 24, units = "cm", device = cairo_pdf)

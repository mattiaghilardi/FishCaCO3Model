# Make figure 1

fig1 <- ggpubr::ggarrange(
  # Fixed effects
  m_caco3_sigma %>%
    tidybayes::gather_draws(b_scalelog_weight, b_scalelog_ril, b_scalemean_T, b_scalesqrt_asp_ratio) %>%
    tidybayes::median_qi(.width = c(.5, .95)) %>%
    mutate(.variable = case_when(.variable == "b_scalelog_weight" ~ "*ln* Body mass",
                                 .variable == "b_scalelog_ril" ~ "*ln* Relative<br>intestinal length",
                                 .variable == "b_scalemean_T" ~ "Temperature",
                                 .variable == "b_scalesqrt_asp_ratio" ~ "*sqrt* Aspect ratio")) %>%
    ggplot(aes(y = reorder(as.factor(.variable), .value), x = .value, xmin = .lower, xmax = .upper)) +
    tidybayes::geom_pointinterval(interval_size_range = c(0.5, 1), point_size = 1.2) +
    ylab("") + xlab("Effect size") +
    geom_vline(xintercept = 0, linetype = 2) +
    theme(axis.text.y = ggtext::element_markdown(size = 10)),
  # Family effect
  m_caco3_sigma %>%
    tidybayes::spread_draws(b_Intercept, r_family[family,]) %>%
    mutate(family = recode(family,
                           "Scaridae" = "Labridae-S",
                           "Labridae" = "Labridae-NS")) %>%
    ggplot(aes(y = family, x = r_family)) +
    geom_vline(aes(xintercept = 0), linetype = 2, color = "black") +
    tidybayes::stat_pointinterval(point_size = 1.2, .width = c(0.5, 0.95),
                                  interval_size_range = c(0.5, 1)) +
    ylab("") + xlab("Effect size") +
    theme(axis.title = element_text(size = 12),
          axis.text = element_text(size = 10)),
  ncol = 2, align = "hv", labels = letters[1:2])

ggsave(here::here("outputs", "figures", "fig1.png"), fig1,
       width = 17, height = 10, units = "cm", dpi = 600, type = "cairo")
ggsave(here::here("outputs", "figures", "fig1.pdf"), fig1,
       width = 17, height = 10, units = "cm", device = cairo_pdf)

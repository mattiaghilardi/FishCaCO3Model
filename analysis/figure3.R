# Make figure 3

fig3 <- ggpubr::ggarrange(
  # Fixed effects lognormal
  m3_caco3_comp %>%
    tidybayes::gather_draws(b_Humolh_scalelog_weight,
                            b_Lumolh_scalelog_weight,
                            b_Mumolh_scalelog_weight,
                            b_ARumolh_scalelog_weight,
                            b_ACumolh_scalelog_weight,
                            b_Humolh_scalelog_ril,
                            b_Mumolh_scalelog_ril,
                            b_ARumolh_scalelog_ril,
                            b_ACumolh_scalelog_ril,
                            b_Humolh_scalemean_T,
                            b_ARumolh_scalemean_T,
                            b_ACumolh_scalemean_T,
                            b_Humolh_scalesqrt_asp_ratio,
                            b_Mumolh_scalesqrt_asp_ratio,
                            b_ARumolh_scalesqrt_asp_ratio,
                            b_ACumolh_scalesqrt_asp_ratio) %>%
    tidybayes::median_qi(.width = c(0.5, 0.95)) %>%
    mutate(.variable = gsub("b_", "", .variable),
           var = gsub(".*\\umolh_", "", .variable),
           var = factor(var, levels = c("scalelog_ril",
                                        "scalesqrt_asp_ratio",
                                        "scalemean_T",
                                        "scalelog_weight")
                        ),
           caco3_phase = gsub("\\umolh.*", "", .variable),
           caco3_phase = factor(caco3_phase,
                                levels = c("L", "AR", "H", "M", "AC"),
                                labels = c("LMC", "Aragonite", "HMC", "MHC", "ACMC")
                                )
           ) %>%
    ggplot(aes(y = var, x = .value, xmin = .lower, xmax = .upper, color = caco3_phase)) +
    geom_vline(xintercept = 0, linetype = 2) +
    tidybayes::geom_pointinterval(position = position_dodge(width = 0.7),
                                  interval_size_range = c(0.5, 1), point_size = 1.2) +
    scale_color_viridis_d(option = "C", end = 0.95) +
    labs(title = "Excretion rate", x = "Effect size", y = "") +
    scale_y_discrete(labels = c("*ln* Relative<br>intestinal length",
                                "*sqrt* Aspect ratio",
                                "Temperature",
                                "*ln* Body mass")
                     ) +
    theme(axis.text.y = ggtext::element_markdown(size = 10),
          plot.title = element_text(face = "bold", size = 12, hjust = 0.5),
          legend.title = element_blank()),
  # Fixed effects hurdle probability
  m3_caco3_comp %>%
    tidybayes::gather_draws(b_hu_Humolh_scalelog_ril,
                            b_hu_Lumolh_scalelog_ril,
                            b_hu_Mumolh_scalelog_ril,
                            b_hu_ARumolh_scalelog_ril,
                            b_hu_ACumolh_scalelog_ril,
                            b_hu_Humolh_scalemean_T,
                            b_hu_Lumolh_scalemean_T,
                            b_hu_Mumolh_scalemean_T,
                            b_hu_ARumolh_scalemean_T,
                            b_hu_ACumolh_scalemean_T) %>%
    mutate(.value = -.value) %>%
    tidybayes::median_qi(.width = c(0.5, 0.95)) %>%
    mutate(.variable = gsub("b_hu_", "", .variable),
           var = gsub(".*\\umolh_", "", .variable),
           caco3_phase = gsub("\\umolh.*", "", .variable),
           caco3_phase = factor(caco3_phase,
                                levels = c("L", "AR", "H", "M", "AC"),
                                labels = c("LMC", "Aragonite", "HMC", "MHC", "ACMC")
           )
    ) %>%
    ggplot(aes(y = var, x = .value, xmin = .lower, xmax = .upper, color = caco3_phase)) +
    geom_vline(xintercept = 0, linetype = 2) +
    tidybayes::geom_pointinterval(position = position_dodge(width = 0.7),
                                  interval_size_range = c(0.5, 1), point_size = 1.2) +
    scale_color_viridis_d(option = "C", end = 0.95) +
    labs(title = "Probability of excretion", x = "Effect size", y = "") +
    scale_y_discrete(labels = c("*ln* Relative<br>intestinal length",
                                "Temperature")
    ) +
    theme(axis.text.y = ggtext::element_markdown(size = 10),
          plot.title = element_text(face = "bold", size = 12, hjust = 0.5),
          legend.title = element_blank()),
  ncol = 2, align = "hv", common.legend = TRUE, labels = letters[1:2], label.x = 0.2)

ggsave(here::here("outputs", "figures", "fig3.png"), fig3,
       width = 17, height = 10, units = "cm", dpi = 600, type = "cairo")
ggsave(here::here("outputs", "figures", "fig3.pdf"), fig3,
       width = 17, height = 10, units = "cm", device = cairo_pdf)

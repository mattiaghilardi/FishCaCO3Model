# Make figure 6

nd <- expand.grid(t = 25,
                  w = seq(1, 5000, 100),
                  ril = c(seq(0.25, 1.5, 0.25), 2, 4, 6),
                  ar = seq(0.5, 4, 0.5))

nd <- cbind(nd, pred_C(t = nd$t, w = nd$w, ril = nd$ril, ar = nd$ar))

fig6 <- nd %>%
  mutate(w = w/1000) %>%
  tidyr::pivot_longer(cols = c(C_wilson_p1, C_wilson_p2.4), names_to = "p", values_to = "C_wilson") %>%
  ggplot(aes(x = w, y = C_ghilardi, color = as.factor(ril))) +
  geom_line(aes(y = C_wilson, linetype = p), color = "black", size = 1) +
  geom_line() +
  facet_wrap(facets = vars(paste("AR = ", ar)), nrow = 2) +
  scale_color_viridis_d("RIL", option = "D", direction = -1, end = 0.9) +
  scale_linetype_manual("Wilson et al.<br>(2009)",
                        values = c(1, 3),
                        labels = paste("<i>&rho;</i> =", c(1, 2.4))) +
  xlab("Body mass (kg)") +
  ylab("Ca(Mg)CO<sub>3</sub> excretion (&mu;mol h<sup>-1</sup>)") +
  theme(axis.title.y = ggtext::element_markdown(),
        legend.text = ggtext::element_markdown(),
        legend.title = ggtext::element_markdown())

ggsave(here::here("outputs", "figures", "fig6.png"), fig6,
       width = 20, height = 12, unit = "cm", dpi = 600, type = "cairo")
ggsave(here::here("outputs", "figures", "fig6.pdf"), fig6,
       width = 20, height = 12, unit = "cm", device = cairo_pdf)

rm(nd)

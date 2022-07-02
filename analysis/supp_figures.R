# Make supplementary figures

# Fig. S1
# Plot observed vs predicted for total excretion rate
predict(m_caco3_sigma) %>%
  as_tibble() %>%
  bind_cols(data_caco3_f_weight) %>%
  plot_obs_vs_pred(pred = "Estimate", obs = "log_exc_rate",
                   xlab = "Predicted *ln* Ca(Mg)CO<sub>3</sub> excretion (&mu;mol h<sup>-1</sup>)",
                   ylab = "Observed *ln* Ca(Mg)CO<sub>3</sub> excretion (&mu;mol h<sup>-1</sup>)",
                   point_size = 1.5,
                   text_size = 4) +
  theme(axis.title.x = ggtext::element_markdown(),
        axis.title.y = ggtext::element_markdown())

ggsave(here::here("outputs", "figures", "obs_vs_pred_caco3.pdf"),
       width = 14, height = 10, units = "cm", device = cairo_pdf)

# Fig. S2
# PP check proportion of zeroes
ggpubr::ggarrange(plotlist = Map(function(x, y) {
  pp_check(m3_caco3_comp, resp = x, nsamples = 1000,
           type = "stat", stat = "prop_zero", binwidth = 0.005) +
    ggtitle(y) +
    theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 12))
}, x = c("Lumolh", "ARumolh", "Humolh", "Mumolh", "ACumolh"),
y = c("LMC", "Aragonite", "HMC", "MHC", "ACMC")),
ncol = 3, nrow = 2, legend = "none")

ggsave(here::here("outputs", "figures", "proportion_of_zeros.pdf"),
       width = 16, height = 8, units = "cm", device = cairo_pdf)

# Fig. S3
# Plot observed vs predicted for carbonate composition
set.seed(98765)
caco3_comp_pred <- as_tibble(predict(m3_caco3_comp, newdata = caco3_comp, robust = TRUE))

ggpubr::ggarrange(plotlist = Map(function(x, y) {
  plot_obs_vs_pred(obs = log(caco3_comp[, paste0(x, "_umol_h"), drop = TRUE] + 1),
                   pred = log(caco3_comp_pred[, paste0("Estimate.", x, "umolh"), drop = TRUE] + 1),
                   xlab = "*ln* (x+1) Predicted excretion (&mu;mol h<sup>-1</sup>)",
                   ylab = "*ln* (x+1) Observed excretion (&mu;mol h<sup>-1</sup>)",
                   point_size = 1.5) +
    ggtitle(y) +
    theme(axis.title.x = ggtext::element_markdown(size = 9),
          axis.title.y = ggtext::element_markdown(size = 9),
          plot.title = element_text(face = "bold", hjust = 0.5, size = 12))
}, x = c("L", "AR", "H", "M", "AC"),
y = c("LMC", "Aragonite", "HMC", "MHC", "ACMC")),
ncol = 3, nrow = 2)

ggsave(here::here("outputs", "figures", "obs_vs_pred_caco3_comp.pdf"),
       width = 22, height = 14, units = "cm", device = cairo_pdf)

# Fig. S4
# Plot family effect on carbonate composition

# Prepare data
df <- bind_rows(m3_caco3_comp %>%
                  tidybayes::gather_draws(r_family__hu_Humolh[family,],
                                          r_family__hu_Lumolh[family,],
                                          r_family__hu_Mumolh[family,],
                                          r_family__hu_ACumolh[family,],
                                          r_family__hu_ARumolh[family,]) %>%
                  mutate(.value = -.value,
                         pars = "hu"),
                m3_caco3_comp %>%
                  tidybayes::gather_draws(r_family__Humolh[family,],
                                          r_family__Lumolh[family,],
                                          r_family__Mumolh[family,],
                                          r_family__ACumolh[family,],
                                          r_family__ARumolh[family,]) %>%
                  mutate(pars = "mu")) %>%
  group_by(family, .variable, pars) %>%
  tidybayes::median_qi(.width = c(0.5, 0.95)) %>%
  mutate(var = if_else(pars == "hu",
                       gsub("r_family__hu_", "", .variable),
                       gsub("r_family__", "", .variable)),
         caco3_phase = gsub("\\umolh.*", "", var),
         caco3_phase = factor(caco3_phase,
                              levels = c("L", "AR", "H", "M", "AC"),
                              labels = c("LMC", "Aragonite", "HMC", "MHC", "ACMC")),
         pars = factor(pars,
                       levels = c("hu", "mu"),
                       labels = c("Probability of excretion", "Excretion rate")),
         effect = if_else(.lower > 0, "pos",
                          if_else(.upper < 0, "neg", "null")),
         family = recode(family,
                         "Scaridae" = "Labridae-S",
                         "Labridae" = "Labridae-NS")
  )

# Plot
figS4 <- lapply(levels(df$pars), function(i) {
  ggplot(filter(df, pars == i), aes(y = family, x = .value, xmin = .lower, xmax = .upper)) +
    geom_vline(xintercept = 0, linetype = 2) +
    tidybayes::geom_pointinterval(aes(color = effect), point_size = 1.2,
                                  interval_size_range = c(0.5, 1)) +
    scale_color_manual(breaks = c("neg", "null", "pos"),
                       values = c("#B2182B", "#969696", "#2166AC"),
                       guide = "none") +
    facet_grid(cols = vars(caco3_phase),
               scales = "free_x") +
    xlab("Effect size") + ylab("") +
    theme(axis.text = element_text(size = 10),
          strip.background = element_blank(),
          strip.text = element_text(face = "bold", size = 12))
  })

# Arrange plots
ggpubr::ggarrange(ggpubr::annotate_figure(figS4[[1]] + ggpubr::rremove("x.title"),
                                          right = ggpubr::text_grob("Probability of excretion",
                                                                    rot = 270, family = "serif", face = "bold", size = 12)) ,
                  ggpubr::annotate_figure(figS4[[2]] + theme(strip.text = element_blank()),
                                          right = ggpubr::text_grob("Excretion rate",
                                                                    rot = 270, family = "serif", face = "bold")),
                  nrow = 2, labels = letters[1:2], label.x = 0.01)

ggsave(here::here("outputs", "figures", "family_effect_composition.pdf"),
       width = 22, height = 20, units = "cm", device = cairo_pdf)

rm(df, figS4)

# Fig. S5
# Plot average carbonate composition by family

# Get family-level traits
fam_traits <- family_traits(unique(caco3_comp$family), type = "food items")

# Average adult (1/2 max length) biomass for each family
taxo <- rfishbase::load_taxa() %>%
  as.data.frame() %>%
  filter(!is.na(Species) & Family %in% fam_traits$family)

biom_fam <- rfishbase::estimate(taxo$Species, fields = c("Species", "MaxLengthTL")) %>%
  left_join(taxo[, c("Family", "Species")]) %>%
  rename(family = Family, species = Species) %>%
  mutate(length = MaxLengthTL/2) %>%
  biomass_FB() %>%
  group_by(family) %>%
  summarise(biomass = round(mean(biomass, na.rm = TRUE), 2)) %>%
  ungroup()

# Prepare dataset (two T levels: 25 and 30°C)
nd <- left_join(fam_traits, biom_fam) %>%
  mutate(log_weight = log(biomass/1000),
         log_ril = log(ril),
         sqrt_asp_ratio = sqrt(aspect_ratio),
         method = "double") %>%
  tidyr::expand_grid(mean_T = c(25, 30))

# Prediction and computation of relative composition
nd <- predict(m3_caco3_comp, newdata = nd, robust = TRUE) %>%
  as_tibble() %>%
  select(starts_with("Estimate")) %>%
  mutate(total = rowSums(across()),
         propL = Estimate.Lumolh/total,
         propH = Estimate.Humolh/total,
         propAR = Estimate.ARumolh/total,
         propM = Estimate.Mumolh/total,
         propAC = Estimate.ACumolh/total) %>%
  bind_cols(nd) %>%
  tidyr::pivot_longer(cols = starts_with("prop"), names_to = "caco3_phase", values_to = "prop") %>%
  mutate(caco3_phase = factor(caco3_phase,
                              levels = paste0("prop", c("L", "AR", "H", "M", "AC")),
                              labels = c("LMC", "Aragonite", "HMC", "MHC", "ACMC")
  ),
  family = recode(family,
                  "Scaridae" = "Labridae-S",
                  "Labridae" = "Labridae-NS")
  )

# Plot
ggplot(nd, aes(x = prop, y = as.factor(mean_T), fill = caco3_phase)) +
  geom_col(orientation = "y", position = "stack", col = "black", size = 0.2, width = 1) +
  facet_grid(rows = vars(forcats::fct_rev(family))) +
  scale_fill_viridis_d("", option = "C") +
  xlab("Proportion of excreted Ca(Mg)CO<sub>3</sub>") + ylab("Temperature (°C)") +
  scale_x_continuous(limits = c(0, 1), expand = expansion(mult = 0)) +
  theme(legend.position = "top",
        strip.background = element_blank(),
        strip.text.y = element_text(angle = 0, hjust = 0, size = 10),
        axis.text = element_text(size = 10),
        panel.spacing.y = unit(0.2, "lines"),
        axis.title.x = ggtext::element_markdown())

ggsave(here::here("outputs", "figures", "family_composition.pdf"),
       width = 18, height = 18, units = "cm", device = cairo_pdf)

rm(nd, fam_traits, biom_fam, taxo)

# Fig. S6
# Compare predictions of Wilson's model with ours
nd <- expand.grid(t = c(20, 25, 30),
                  w = seq(1, 5000, 100),
                  ril = 0.5,
                  ar = 1.5)

nd <- cbind(nd, pred_C(t = nd$t, w = nd$w, ril = nd$ril, ar = nd$ar))

ggplot(nd %>%
         tidyr::pivot_longer(cols = c(C_wilson_p1, C_ghilardi),
                             names_to = "var",
                             values_to = "value") %>%
         mutate(w = w/1000),
       aes(x = w, y = value, color = var)) +
  geom_line() +
  facet_wrap(facets = vars(paste("T =", t, "°C"))) +
  scale_color_manual("",
                     values = c("firebrick", "deepskyblue2"),
                     labels = c("This study",
                                "Wilson et al. (2009)")) +
  xlab("Body mass (kg)") +
  ylab("Ca(Mg)CO<sub>3</sub> excretion (&mu;mol h<sup>-1</sup>)") +
  theme(axis.title.y = ggtext::element_markdown(),
        legend.position = c(0.15, 0.9),
        legend.background = element_blank())

ggsave(here::here("outputs", "figures", "comparison_with_wilson.pdf"),
       width = 18, height = 8, unit = "cm", device = cairo_pdf)

rm(nd)

# Fig. S7
# Sensitivity - Fig. 1
ggpubr::ggarrange(
  # Fixed effects
  m_caco3_sigma_corrected %>%
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
  m_caco3_sigma_corrected %>%
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

ggsave(here::here("outputs", "figures", "sensitivity_fig1.pdf"),
       width = 17, height = 10, units = "cm", device = cairo_pdf)

# Fig. S8
# Sensitivity - Fig. 3
ggpubr::ggarrange(
  # Fixed effects lognormal
  m3_caco3_comp_corrected %>%
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
  m3_caco3_comp_corrected %>%
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

ggsave(here::here("outputs", "figures", "sensitivity_fig3.pdf"),
       width = 17, height = 10, units = "cm", device = cairo_pdf)

# Fig. S9
# Plot observed vs predicted intestinal length by observed/unobserved species
left_join(int_man_tet, int_man_tet_pred) %>%
  group_by(obs) %>%
  mutate(nobs = n(),
         label = paste0(obs, " (", nobs, ")")) %>%
  ggplot(aes(x = int_length, y = il_log, color = label, fill = label)) +
  geom_point(size = 1.5, alpha = 0.5, stroke = 0.1, show.legend = FALSE) +
  geom_abline(linetype = 2, color = "red", size = 0.7) +
  geom_smooth(method = "lm", size = 0.7, alpha = 0.3) +
  scale_color_viridis_d(option = "C", end = 0.8) +
  scale_fill_viridis_d(option = "C", end = 0.8) +
  ggpubr::stat_regline_equation(aes(label =  paste(..eq.label.., ..rr.label.., sep = "~~~~")),
                                family = "serif", show.legend = FALSE, size = 4) +
  labs(x = "Predicted"~italic(ln)~"int. length (mm)",
       y = "Observed"~italic(ln)~"int. length (mm)") +
  theme(legend.title = element_blank(),
        legend.position = c(0.8, 0.2))

ggsave(here::here("outputs", "figures", "int_validation_observed.pdf"),
       width = 14, height = 10, units = "cm", device = cairo_pdf)

# Fig. S10
# Plot observed vs predicted intestinal length by location
left_join(int_man_tet, int_man_tet_pred) %>%
  filter(obs == "observed") %>%
  group_by(location) %>%
  mutate(nobs = n(),
         label = paste0(location, " (", nobs, ")")) %>%
  ggplot(aes(x = int_length, y = il_log, color = label, fill = label)) +
  geom_point(size = 1.5, alpha = 0.5, stroke = 0.1, show.legend = FALSE) +
  geom_abline(linetype = 2, color = "red", size = 0.7) +
  geom_smooth(method = "lm", size = 0.7, alpha = 0.3) +
  scale_color_viridis_d(option = "C", end = 0.8) +
  scale_fill_viridis_d(option = "C", end = 0.8) +
  ggpubr::stat_regline_equation(aes(label =  paste(..eq.label.., ..rr.label.., sep = "~~~~")),
                                family = "serif", show.legend = FALSE, size = 4) +
  labs(x = "Predicted"~italic(ln)~"int. length (mm)",
       y = "Observed"~italic(ln)~"int. length (mm)") +
  theme(legend.title = element_blank(),
        legend.position = c(0.8, 0.2))

ggsave(here::here("outputs", "figures", "int_validation_location.pdf"),
       width = 14, height = 10, units = "cm", device = cairo_pdf)

# Fig. S11
# Plot observed vs predicted intestinal length at genus-level
left_join(int_sp, int_sp_pred) %>%
  plot_obs_vs_pred(pred = "int_length", obs = "il_log",
                   xlab = "Predicted"~italic(ln)~"int. length (mm)",
                   ylab = "Observed"~italic(ln)~"int. length (mm)",
                   point_size = 1.5,
                   text_size = 4)

ggsave(here::here("outputs", "figures", "int_validation_genus.pdf"),
       width = 14, height = 10, units = "cm", device = cairo_pdf)

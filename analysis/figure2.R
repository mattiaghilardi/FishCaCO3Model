# Make figure 2

# Tibble with specific arguments for each plot
nd <- tibble(x = c("log_weight", "mean_T", "sqrt_asp_ratio", "log_ril"),
             xlab = c("Body mass (kg)", "Temperature (°C)",
                      "Caudal fin aspect ratio", "Relative intestinal length"),
             trans_x = c("exp", "", "square", "exp"),
             newdata = list(data_caco3_f_weight %>%
                              modelr::data_grid(log_weight = modelr::seq_range(log_weight, n = 101),
                                                log_ril = mean(log_ril),
                                                mean_T = mean(mean_T),
                                                sqrt_asp_ratio = mean(sqrt_asp_ratio),
                                                method = "double"),
                            data_caco3_f_weight %>%
                              modelr::data_grid(log_weight = mean(log_weight),
                                                log_ril = mean(log_ril),
                                                mean_T = modelr::seq_range(mean_T, n = 101),
                                                sqrt_asp_ratio = mean(sqrt_asp_ratio),
                                                method = "double"),
                            data_caco3_f_weight %>%
                              modelr::data_grid(log_weight = mean(log_weight),
                                                log_ril = mean(log_ril),
                                                mean_T = mean(mean_T),
                                                sqrt_asp_ratio = modelr::seq_range(sqrt_asp_ratio, n = 101),
                                                method = "double"),
                            data_caco3_f_weight %>%
                              modelr::data_grid(log_weight = mean(log_weight),
                                                log_ril = modelr::seq_range(log_ril, n = 101),
                                                mean_T = mean(mean_T),
                                                sqrt_asp_ratio = mean(sqrt_asp_ratio),
                                                method = "double")
             )
)

# List of tibbles for dots
dots <- list(tibble(x = c(0.5, 3),
                    log_exc_rate = fitted(m_caco3_sigma,
                                          newdata = tibble(log_weight = log(c(0.5, 3))) %>%
                                            bind_cols(nd$newdata[[1]] %>%
                                                        select(-log_weight) %>%
                                                        unique()
                                            ),
                                          re_formula = NA)[, 1],
                    color = c("black", "black"),
                    fill = c("white", "white"),
                    label = c("", "")
             ),
             tibble(x = c(24, 29.3),
                    log_exc_rate = fitted(m_caco3_sigma,
                                          newdata = tibble(mean_T = c(24, 29.3)) %>%
                                            bind_cols(nd$newdata[[2]] %>%
                                                        select(-mean_T) %>%
                                                        unique()
                                            ),
                                          re_formula = NA)[, 1],
                    color = c("black", "black"),
                    fill = c("#FF8C00", "#CC0018"),
                    label = paste(x, "°C")
             ),
             tibble(x = c(min(data_caco3_f_weight$aspect_ratio),
                          max(data_caco3_f_weight$aspect_ratio)),
                    log_exc_rate = fitted(m_caco3_sigma,
                                          newdata = tibble(sqrt_asp_ratio = c(min(data_caco3_f_weight$sqrt_asp_ratio),
                                                                              max(data_caco3_f_weight$sqrt_asp_ratio))) %>%
                                            bind_cols(nd$newdata[[3]] %>%
                                                        select(-sqrt_asp_ratio) %>%
                                                        unique()
                                            ),
                                          re_formula = NA)[, 1],
                    color = c("black", "black"),
                    fill = c("white", "white"),
                    label = x
             ),
             tibble(x = c(1, 5),
                    log_exc_rate = fitted(m_caco3_sigma,
                                          newdata = tibble(log_ril = log(c(1, 5))) %>%
                                            bind_cols(nd$newdata[[4]] %>%
                                                        select(-log_ril) %>%
                                                        unique()
                                            ),
                                          re_formula = NA)[, 1],
                    color = c("black", "black"),
                    fill = c("white", "white"),
                    label = x
             )
)

y_lab <- "Ca(Mg)CO<sub>3</sub> excretion (&mu;mol h<sup>-1</sup>)"
seed <- 98765

# Plot
fig2 <- lapply(1:4, function(i) {
  plot_partial_eff_brms(model = m_caco3_sigma,
                        newdata = nd$newdata[[i]],
                        data = data_caco3_f_weight,
                        x = nd$x[i],
                        xlab = nd$xlab[i],
                        ylab = y_lab,
                        seed = seed,
                        backtrans_y = "exp",
                        backtrans_x = nd$trans_x[i]) +
    theme(text = element_text(size = 10),
          axis.title.y = ggtext::element_markdown()) +
    geom_point(data = dots[[i]],
               aes(x = x,
                   y = exp(log_exc_rate),
                   color = levels(as.factor(color)),
                   fill = levels(as.factor(fill))),
               shape = c(21, 24),
               size = 2) +
    scale_color_manual(values = dots[[i]]$color, guide = "none") +
    scale_fill_manual(values = dots[[i]]$fill, guide = "none") +
    annotate(geom = "text", x = dots[[i]]$x, y = exp(dots[[i]]$log_exc_rate),
             label = dots[[i]]$label, vjust = -2,
             size = 3.5, family = "serif", fill = NA)
    })

# Equation parameters
# Body mass
# Unscaled slope
b1w <- fixef(m_caco3_sigma)["scalelog_weight", 1] / sd(data_caco3_f_weight$log_weight)
# Intercept W = 1 kg
lnb0w <- fixef(m_caco3_sigma)["Intercept", 1] +
  b1w * (log(1) - mean(data_caco3_f_weight$log_weight))
b0w <- round(exp(lnb0w), 2)
b1w <- round(b1w, 2)

# Temperature
# Unscaled slope
b1t <- fixef(m_caco3_sigma)["scalemean_T", 1] / sd(data_caco3_f_weight$mean_T)
# Intercept at T = 0°C
b0t <- exp(fixef(m_caco3_sigma)["Intercept", 1] +
             b1t * (0 - mean(data_caco3_f_weight$mean_T)))
b0t <- round(b0t, 2)
b1t <- round(b1t, 2)

b0t * exp(b1t * T)

# Aspect ratio
# Unscaled slope
b1ar <- fixef(m_caco3_sigma)["scalesqrt_asp_ratio", 1] / sd(data_caco3_f_weight$sqrt_asp_ratio)
# Intercept at AR = 0
b0ar <- exp(fixef(m_caco3_sigma)["Intercept", 1] +
              b1ar * (sqrt(0) - mean(data_caco3_f_weight$sqrt_asp_ratio)))
b0ar <- round(b0ar, 2)
b1ar <- round(b1ar, 2)

# RIL
# Unscaled slope
b1ril <- fixef(m_caco3_sigma)["scalelog_ril", 1] / sd(data_caco3_f_weight$log_ril)
# Intercept RIL = 1
b0ril <- exp(fixef(m_caco3_sigma)["Intercept", 1] +
               b1ril * (log(1) - mean(data_caco3_f_weight$log_ril)))
b0ril <- round(b0ril, 2)
b1ril <- round(b1ril, 2)

# Isometric scaling
isometry <- tibble(log_weight = modelr::seq_range(data_caco3_f_weight$log_weight, n = 101),
                   log_exc_rate = lnb0w + 1 * log_weight
                   )

dots_iso <- tibble(log_weight = log(c(0.5, 3)),
                   log_exc_rate = lnb0w + 1 * log_weight
                   )

# Function to add equations
add_eq <- function(label) {
  annotate(geom = "richtext", label = label, size = 3.5,
           x = -Inf, y = Inf, hjust = -0.2, vjust = 1.2,
           family = "serif", label.color = NA, fill = NA)
}

# Add isometry and equations
fig2[[1]] <- fig2[[1]] +
  geom_line(data = isometry,
            mapping = aes(x = exp(log_weight), y = exp(log_exc_rate)),
            linetype = "dashed", color = "#CC0018") +
  geom_point(data = dots_iso, aes(x = exp(log_weight), y = exp(log_exc_rate)),
             fill = "#CC0018", color = "black", shape = c(21, 24), size = 2) +
  annotate(geom = "text", x = exp(dots_iso$log_weight), y = exp(dots_iso$log_exc_rate),
           label = paste(c("0.5", "3"), "kg"), vjust = -2,
           size = 3.5, family = "serif", fill = NA) +
  annotate(geom = "richtext", x = 10, y = 145, label = "*&beta;* = 1",
           color = "#CC0018", size = 3.5, family = "serif", label.color = NA, fill = NA) +
  annotate(geom = "richtext", x = 10, y = 30, label = glue::glue("*&beta;* = {b1w}"),
           size = 3.5, family = "serif", label.color = NA, fill = NA) +
  add_eq(glue::glue("*y* = {b0w} · *x^&beta;*"))

fig2[[2]] <- fig2[[2]] +
  add_eq(glue::glue("*y* = {b0t} · e^{b1t}*x*"))

fig2[[3]] <- fig2[[3]] +
  add_eq(glue::glue("*y* = {b0ar} · e^{b1ar}&radic;*x*"))

fig2[[4]] <- fig2[[4]] +
  add_eq(glue::glue("*y* = {b0ril} · *x*^{b1ril}"))

# Final plot
fig2 <- ggpubr::ggarrange(plotlist = fig2, ncol = 2, nrow = 2, align = "hv",
                          labels = letters[1:4], label.x = 0.02)

ggsave(here::here("outputs", "figures", "fig2.png"), fig2,
       width = 17, height = 15, units = "cm", dpi = 600, type = "cairo")
ggsave(here::here("outputs", "figures", "fig2.pdf"), fig2,
       width = 17, height = 15, units = "cm", device = cairo_pdf)

rm(nd, y_lab, seed, isometry, dots, dots_iso, add_eq,
   lnb0w, b0w, b1w, b0t, b1t, b0ar, b1ar, b0ril, b1ril)

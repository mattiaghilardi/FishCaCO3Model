# Make figure 5

# Extract and summarise posterior samples of correlation parameters between hurdle probabilities
hu_cor <- posterior_samples(m3_caco3_comp, pars = "^cor_") %>%
  as_tibble() %>%
  tidyr::pivot_longer(names_to = "var", values_to = "value", cols = 1:45) %>%
  mutate(var = gsub(".*\\cor_family__", "", var),
         var1 = gsub("\\_Intercept.*", "", var),
         var2 = gsub(".*\\_Intercept__", "", var),
         var2 = gsub("\\_Intercept.*", "", var2),
         var1 = ifelse(startsWith(var1, "hu_"), var1, paste0("mu_", var1)),
         var2 = ifelse(startsWith(var2, "hu_"), var2, paste0("mu_", var2))
         ) %>%
  group_by(var1, var2) %>%
  summarise(median = median(value),
            q2.5 = quantile(value, probs = 0.025),
            q97.5 = quantile(value, probs = 0.975)) %>%
  ungroup() %>%
  filter(startsWith(var1, "hu_") & startsWith(var2, "hu_"))

# Duplicate tibble by inverting variable names
hu_cor <- bind_rows(hu_cor,
          hu_cor %>%
            rename(var1 = var2, var2 = var1)) %>%
  mutate(var1 = factor(var1,
                       levels = paste0("hu_", c("L", "AR", "H", "M", "AC"), "umolh"),
                       labels = c("LMC", "Aragonite", "HMC", "MHC", "ACMC")),
         var2 = factor(var2,
                       levels = paste0("hu_", c("L", "AR", "H", "M", "AC"), "umolh"),
                       labels = c("LMC", "Aragonite", "HMC", "MHC", "ACMC"))
         )

# Add terms for diagonal of correlation matrix
hu_cor <- bind_rows(hu_cor,
                    tibble(var1 = factor(levels(hu_cor$var1)),
                           var2 = factor(levels(hu_cor$var2)),
                           median = 1))

# Transform in correlation matrix
hu_cor_mat <- hu_cor %>%
  select(1:3) %>%
  arrange(var1, var2) %>%
  tidyr::pivot_wider(values_from = median, names_from = var2) %>%
  tibble::column_to_rownames("var1") %>%
  as.matrix()

# Helper functions to get lower and upper triangle of a correlation matrix
get_lower_tri <- function(m, diag = FALSE){
  m[upper.tri(m, diag = diag)] <- NA
  m
}
get_upper_tri <- function(m, diag = FALSE){
  m[lower.tri(m, diag = diag)] <- NA
  m
}

# Tibble with median correlations on the lower triangular part
hu_cor_low <- as_tibble(get_lower_tri(hu_cor_mat, diag = TRUE)) %>%
  tidyr::pivot_longer(cols = 1:5, names_to = "var2", values_to = "median") %>%
  mutate(var2 = factor(var2, levels = levels(hu_cor$var2)),
         var1 = factor(rep(levels(var2), each = nlevels(var2)), levels = levels(var2)))

# Tibble with median and 95% CI on the upper triangular part
hu_cor_up <- as_tibble(get_upper_tri(hu_cor_mat, diag = TRUE)) %>%
  tidyr::pivot_longer(cols = 1:5, names_to = "var2", values_to = "median") %>%
  mutate(var2 = factor(var2, levels = levels(hu_cor$var2)),
         var1 = factor(rep(levels(var2), each = nlevels(var2)), levels = levels(var2))) %>%
  left_join(hu_cor) %>%
  mutate(label = ifelse(!is.na(median),
                        paste0(round(median, 2), "\n[", round(q2.5, 2), ", ", round(q97.5, 2), "]"),
                        NA)
  )

# Tibble with polymorphs' names on the diagonal
hu_cor_diag <- as_tibble(get_upper_tri(hu_cor_mat)) %>%
  tidyr::pivot_longer(cols = 1:5, names_to = "var2", values_to = "median") %>%
  mutate(var2 = factor(var2, levels = levels(hu_cor$var2)),
         var1 = factor(rep(levels(var2), each = nlevels(var2)), levels = levels(var2)),
         name = ifelse(median == 1, as.character(var2), NA))

# Plot
fig5 <- ggplot(mapping = aes(x = var1, y = var2)) +
  geom_tile(data = hu_cor_low, aes(fill = median)) +
  scale_fill_gradientn(colors = colorRampPalette(RColorBrewer::brewer.pal(11, "RdBu"))(100),
                       limits = c(-1, 1), na.value = NA) +
  geom_text(data = hu_cor_up,
            aes(label = label), size = 3.5, color = "black") +
  geom_text(data = hu_cor_diag,
            aes(label = name), size = 4, color = "black") +
  theme_void() +
  theme(legend.title = element_blank(),
        legend.key.height = unit(3.5, units = "lines"))

ggsave(here::here("outputs", "figures", "fig5.png"), fig5,
       units = "cm", width = 13, height = 10, dpi = 600, type = "cairo")
ggsave(here::here("outputs", "figures", "fig5.pdf"), fig5,
       units = "cm", width = 13, height = 10, device = cairo_pdf)

rm(hu_cor, hu_cor_diag, hu_cor_low, hu_cor_up, hu_cor_mat, get_lower_tri, get_upper_tri)

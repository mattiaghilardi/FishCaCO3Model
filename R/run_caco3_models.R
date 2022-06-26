#' Run a set of 20 models on species-level carbonate excretion rates
#'
#' @param data A data frame containing 9 columns:
#'             "log_exc_rate": natural-log transformed carbonate excretion rate in "µmol/h";
#'             "log_weight": natural-log transformed fish weight in "kg";
#'             "log_ril": natural-log transformed relative intestinal length;
#'             "sqrt_asp_ratio": square-root transformed aspect ratio of caudal fin;
#'             "mean_S": salinity;
#'             "mean_T": water temperature in "°C";
#'             "family": taxonomic family;
#'             "sqrt_sampling": square-root transformed sampling interval in hours
#'             "weights": the weight of each observation.
#' @param ... Additional parameters passed to \code{\link[brms:brm]{brms::brm()}}
#'
#' @return A tibble with four columns:
#'         "id: the model names;
#'         "bf": the model formulas;
#'         "prior": the priors;
#'         "fit": the model objects.
#'
#' @importFrom brms brm prior student
#' @importFrom dplyr tibble group_split mutate bind_rows
#'
run_caco3_models <- function(data, ...) {

  # model formulas
  bf <- list(
    # linear
    log_exc_rate|weights(weights) ~ 1,
    log_exc_rate|weights(weights) ~ scale(log_weight),
    log_exc_rate|weights(weights) ~ scale(log_weight) + scale(log_ril),
    log_exc_rate|weights(weights) ~ scale(log_weight) + scale(mean_S),
    log_exc_rate|weights(weights) ~ scale(log_weight) + scale(log_ril) + scale(mean_T),
    log_exc_rate|weights(weights) ~ scale(log_weight) + scale(log_ril) + scale(sqrt_asp_ratio),
    log_exc_rate|weights(weights) ~ scale(log_weight) + scale(mean_S) + scale(mean_T),
    log_exc_rate|weights(weights) ~ scale(log_weight) + scale(mean_S) + scale(sqrt_asp_ratio),
    log_exc_rate|weights(weights) ~ scale(log_weight) + scale(log_ril) + scale(mean_T) + scale(sqrt_asp_ratio),
    log_exc_rate|weights(weights) ~ scale(log_weight) + scale(mean_S) + scale(mean_T) + scale(sqrt_asp_ratio),
    log_exc_rate|weights(weights) ~ scale(log_weight) + scale(log_ril) + scale(sqrt_sampling),
    log_exc_rate|weights(weights) ~ scale(log_weight) + scale(mean_S) + scale(sqrt_sampling),
    log_exc_rate|weights(weights) ~ scale(log_weight) + scale(log_ril) + scale(mean_T) + scale(sqrt_sampling),
    log_exc_rate|weights(weights) ~ scale(log_weight) + scale(mean_S) + scale(mean_T) + scale(sqrt_sampling),
    log_exc_rate|weights(weights) ~ scale(log_weight) + scale(log_ril) + scale(sqrt_asp_ratio) + scale(sqrt_sampling),
    log_exc_rate|weights(weights) ~ scale(log_weight) + scale(mean_S) + scale(sqrt_asp_ratio) + scale(sqrt_sampling),
    log_exc_rate|weights(weights) ~ scale(log_weight) + scale(log_ril) + scale(mean_T) + scale(sqrt_asp_ratio) + scale(sqrt_sampling),
    log_exc_rate|weights(weights) ~ scale(log_weight) + scale(mean_S) + scale(mean_T) + scale(sqrt_asp_ratio) + scale(sqrt_sampling),

    # multilevel - add family as group level effect on the intercept
    log_exc_rate|weights(weights) ~ (1| family),
    log_exc_rate|weights(weights) ~ scale(log_weight) + (1| family),
    log_exc_rate|weights(weights) ~ scale(log_weight) + scale(log_ril) + (1| family),
    log_exc_rate|weights(weights) ~ scale(log_weight) + scale(mean_S) + (1| family),
    log_exc_rate|weights(weights) ~ scale(log_weight) + scale(log_ril) + scale(mean_T) + (1| family),
    log_exc_rate|weights(weights) ~ scale(log_weight) + scale(log_ril) + scale(sqrt_asp_ratio) + (1| family),
    log_exc_rate|weights(weights) ~ scale(log_weight) + scale(mean_S) + scale(mean_T) + (1| family),
    log_exc_rate|weights(weights) ~ scale(log_weight) + scale(mean_S) + scale(sqrt_asp_ratio) + (1| family),
    log_exc_rate|weights(weights) ~ scale(log_weight) + scale(log_ril) + scale(mean_T) + scale(sqrt_asp_ratio) + (1| family),
    log_exc_rate|weights(weights) ~ scale(log_weight) + scale(mean_S) + scale(mean_T) + scale(sqrt_asp_ratio) + (1| family),
    log_exc_rate|weights(weights) ~ scale(log_weight) + scale(log_ril) + scale(sqrt_sampling) + (1| family),
    log_exc_rate|weights(weights) ~ scale(log_weight) + scale(mean_S) + scale(sqrt_sampling) + (1| family),
    log_exc_rate|weights(weights) ~ scale(log_weight) + scale(log_ril) + scale(mean_T) + scale(sqrt_sampling) + (1| family),
    log_exc_rate|weights(weights) ~ scale(log_weight) + scale(mean_S) + scale(mean_T) + scale(sqrt_sampling) + (1| family),
    log_exc_rate|weights(weights) ~ scale(log_weight) + scale(log_ril) + scale(sqrt_asp_ratio) + scale(sqrt_sampling) + (1| family),
    log_exc_rate|weights(weights) ~ scale(log_weight) + scale(mean_S) + scale(sqrt_asp_ratio) + scale(sqrt_sampling) + (1| family),
    log_exc_rate|weights(weights) ~ scale(log_weight) + scale(log_ril) + scale(mean_T) + scale(sqrt_asp_ratio) + scale(sqrt_sampling) + (1| family),
    log_exc_rate|weights(weights) ~ scale(log_weight) + scale(mean_S) + scale(mean_T) + scale(sqrt_asp_ratio) + scale(sqrt_sampling) + (1| family)
  )

  # priors
  priors_I  <- c(brms::prior(normal(0, 5), class = Intercept),
                 brms::prior(student_t(3, 0, 2.5), class = sigma),
                 brms::prior(gamma(2, 0.1), class = nu)) #default in brms
  prior_b  <- brms::prior(normal(0, 5), class = b)
  prior_sd <- brms::prior(normal(0, 5), class = sd)
  priors_mI <- c(priors_I, prior_sd)
  priors_l <- c(priors_I, prior_b)
  priors_ml <- c(priors_I, prior_b, prior_sd)

  # Tibble with model names, formulas and priors
  model_tbl <- dplyr::tibble(id = c(paste0("m", 1:18, "_caco3"), paste0("m", 1:18, "_fam_caco3")),
                             bf = bf,
                             prior = c(list(priors_I), rep(list(priors_l), 17), list(priors_mI), rep(list(priors_ml), 17)))

  # Convert tibble to list
  model_list <- dplyr::group_split(model_tbl, id, .keep = TRUE)

  # Run models
  models <- lapply(seq_along(model_list), function(x) {
    message(paste("Running model", x, "of", length(model_list)))
    fit <- brms::brm(model_list[[x]]$bf[[1]],
                     data = data,
                     family = brms::student(),
                     prior = model_list[[x]]$prior[[1]],
                     ...,
                     file = paste0("outputs/models/", model_list[[x]]$id)
    )
    x <- dplyr::mutate(model_list[[x]], fit = list(fit))
  })

  # Convert list to tibble
  models <- dplyr::bind_rows(models)
  names(models$fit) <- models$id

  models
}

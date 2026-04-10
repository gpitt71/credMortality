#' @import data.table
#' @import rpart
NULL

# Preprocessing

wide_mat <- function(DT, value_col) {
  w <- dcast(DT, Age ~ Period, value.var = value_col)
  m <- as.matrix(w)
  rownames(m) <- m[, "Age"]
  m <- m[, -1, drop = FALSE]
  m
}

#' @export
fit_mortality_models <- function(data_pp,
                                 years_fit,
                                 ages_fit,
                                 global_mortality_model_candidates = c("lc", "apc", "rh")) {
  l <- list()

  mortality_model_lc <- lc(link = "log")
  mortality_model_apc <- apc(link = "log")
  mortality_model_rh <- rh(link = "log")

  # Lee-Carter

  if ("lc" %in% global_mortality_model_candidates) {
    model_option <- "lc"

    suppressWarnings(
      mortality_model_fit <- fit(
        mortality_model_lc,
        data = data_pp$superpopulation,
        years.fit = years_fit,
        ages.fit = ages_fit,
        verbose = TRUE
      )
    )



    l[[model_option]] <- mortality_model_fit

    start.ax = mortality_model_fit$ax
    start.bx = mortality_model_fit$bx
    start.kt = mortality_model_fit$kt
  }
  # Age-Period-Cohort

  if ("apc" %in% global_mortality_model_candidates) {
    model_option <- "apc"

    suppressWarnings(
      mortality_model_fit <- fit(
        mortality_model_apc,
        data = data_pp$superpopulation,
        years.fit = years_fit,
        ages.fit = ages_fit,
        verbose = FALSE
      )
    )

    l[[model_option]] <- mortality_model_fit
  }

  # Renshaw- Habermann
  if ("rh" %in% global_mortality_model_candidates) {
    model_option <- "rh"

    datahat <- data_pp$superpopulation

    suppressWarnings(
      mortality_model_fit <- fit(
        mortality_model_rh,
        data = data_pp$superpopulation,
        years.fit = years_fit,
        ages.fit = ages_fit,
        verbose = FALSE
      )
    )


    l[[model_option]] <- mortality_model_fit
  }


  return(l)


}

# Fit credibility models ----
#' @export
fit_and_predict_total_model <- function(data,
                                        N_groups,
                                        mortality_models_fit,
                                        years_fit,
                                        ages_fit,
                                        bias = 70,
                                        forecasting_horizon) {
  out <- list()
  observed_rates <- list()

  ## Pre-process the data to get out the tprime plus h year.
  ## The fitting pre-processing only considers up to year tprime.
  data_pp_2 <- data_preprocessing(
    data = data,
    N_groups = N_groups,
    bias = bias,
    ages_fit = ages_fit,
    years_fit = c(years_fit, (last(years_fit) + 1):(last(years_fit) + forecasting_horizon))
    # scenario = scenario
  )


  ## For each model fit, we predict the future central mortality rates
  for (model_option in names(mortality_models_fit)) {
    l <- list()

    mortality_model_fit <- mortality_models_fit[[model_option]]

    year.predict <- max(years_fit) + forecasting_horizon

    cv.arima.kt <- auto.arima(as.numeric(mortality_model_fit$kt), ic = "bic")

    if (model_option != "lc") {
      cv.arima.gc <- auto.arima(as.numeric(mortality_model_fit$gc), ic = "bic")
      gc.order <- unname(arimaorder(cv.arima.gc))
    } else{
      gc.order <- c(1, 1, 0)
    }

    ## actual forecasts
    mortality_model_forecast <- forecast(
      mortality_model_fit,
      kt.method = "iarima",
      gc.order = gc.order,
      kt.order = unname(arimaorder(cv.arima.kt)),
      h = forecasting_horizon
    )

    muxt_hat <- mortality_model_forecast$fitted

    ### different handling if projecting more than one year ahead
    if (forecasting_horizon > 1) {
      muxt_hat_predicted <- mortality_model_forecast$rates[, as.character(last(years_fit) + forecasting_horizon)]


    } else{
      muxt_hat_predicted <- mortality_model_forecast$rates
    }


    ## Predictions for the superpopulation
    l[['muhatsuperpop']] <- muxt_hat_predicted
    observed_rates[['muxt_actual_superpop']] <- data_pp_2$superpopulation$Dxt[, (length(years_fit) + forecasting_horizon)] / data_pp_2$superpopulation$Ext[, (length(years_fit) + forecasting_horizon)]


    ## Predictions for the subpopulations
    for (i in 1:length(data_pp_2$subpopulations)) {
      l[[paste0("muhat", i)]] <-  muxt_hat_predicted

      subgroups_forecast <- data_pp_2$subpopulations[[i]]

      observed_rates[[paste0('muxt_actual_', i)]] <- subgroups_forecast$Dxt[, (length(years_fit) + forecasting_horizon)] / subgroups_forecast$Ext[, (length(years_fit) + forecasting_horizon)]


    }

    # observed_rates[['full_pp_data']] <- data_pp_2
    # out[['actual_data']] <- observed_rates
    # out[[model_option]] <- l

    observed_rates[['pre_processed_data_forecasting']] <- data_pp_2
    out[['actual_data']] <- observed_rates
    out[[model_option]] <- l

  }

  return(out)
}

#' @export
fit_and_predict_credibility_models <- function(data,
                                               N_groups,
                                               mortality_models_fit,
                                               years_fit,
                                               ages_fit,
                                               bias = 70,
                                               forecasting_horizon) {
  out <- list()
  observed_rates <- list()

  data_pp_2 <- data_preprocessing(
    data = data,
    N_groups = N_groups,
    bias = bias,
    ages_fit = ages_fit,
    years_fit = c(years_fit, (last(years_fit) + 1):(last(years_fit) + forecasting_horizon))
  )

  for (model_option in names(mortality_models_fit)) {
    l <- list()

    mortality_model_fit <- mortality_models_fit[[model_option]]

    year.predict <- max(years_fit) + forecasting_horizon

    cv.arima.kt <- auto.arima(as.numeric(mortality_model_fit$kt), ic = "bic")

    if (model_option != "lc") {
      cv.arima.gc <- auto.arima(as.numeric(mortality_model_fit$gc), ic = "bic")
      gc.order <- unname(arimaorder(cv.arima.gc))
    } else{
      gc.order <- c(1, 1, 0)
    }

    #

    mortality_model_forecast <- forecast(
      mortality_model_fit,
      kt.method = "iarima",
      gc.order = gc.order,
      kt.order = unname(arimaorder(cv.arima.kt)),
      h = forecasting_horizon
    )

    muxt_hat <- mortality_model_forecast$fitted

    if (forecasting_horizon > 1) {
      muxt_hat_predicted <- mortality_model_forecast$rates[, as.character(last(years_fit) + forecasting_horizon)]


    } else{
      muxt_hat_predicted <- mortality_model_forecast$rates


    }

    for (i in 1:N_groups) {
      ## subpopulation specific effect
      tmp_occ <- data_pp_2$subpopulations[[i]]$Dxt[as.character(ages_fit), as.character(years_fit)]
      tmp_exp <- data_pp_2$subpopulations[[i]]$Ext[as.character(ages_fit), as.character(years_fit)]

      assign(
        paste0("C", i),
        apply(tmp_occ, 1, sum, na.rm = T) / apply(tmp_exp *
                                                    muxt_hat, 1, sum, na.rm = T)
      )

      fit1 <- rpart(get(paste0("C", i)) ~ ages_fit)
      assign(paste0("C", i), predict(fit1, data.frame(ages_fit)))

      assign(paste0("Fxt_", i), tmp_occ / tmp_exp)


      ## Variance calculations
      assign(
        paste0("tmpvar", i),
        #((

        (
          apply(get(paste0("Fxt_", i)), 1, sum, na.rm = T) - apply(muxt_hat, 1, sum, na.rm =
                                                                     T)
        )^2 - apply(muxt_hat / tmp_exp, 1, sum, na.rm = T)
      )

      denominator <- apply(muxt_hat, 1, sum, na.rm = TRUE)^2

      assign(paste0("varthetax_", i),
             pmax(0, get(paste0("tmpvar", i)) / denominator, na.rm = T))

      ### Variance smoothing

      fit1 <- rpart(get(paste0("varthetax_", i)) ~ ages_fit)
      assign(paste0("varthetax_", i), predict(fit1, data.frame(ages_fit)))

      assign(paste0("Z_", i), 1 / (1 + get(paste0(
        "varthetax_", i
      )) *
        apply(tmp_exp * muxt_hat, 1, sum, na.rm = T)))




      l[[paste0("C", i)]] <- get(paste0("C", i))
      l[[paste0("Fxt_", i)]] <- get(paste0("Fxt_", i))
      l[[paste0("Z_", i)]] <- get(paste0("Z_", i))
      l[[paste0("varthetax_", i)]] <- get(paste0("varthetax_", i))

      l[[paste0("muhat", i)]] <- get(paste0("Z_", i)) * muxt_hat_predicted + (1 - get(paste0("Z_", i))) * get(paste0("C", i)) * muxt_hat_predicted
      l[[paste0("muhat", i, "_MLE")]] <- get(paste0("C", i)) * muxt_hat_predicted


      subgroups_forecast <- data_pp_2$subpopulations[[i]]


      observed_rates[[paste0('muxt_actual_', i)]] <- subgroups_forecast$Dxt[, (length(years_fit) + forecasting_horizon)] / subgroups_forecast$Ext[, (length(years_fit) + forecasting_horizon)]

    }

    # observed_rates[['full_pp_data']] <- data_pp_2
    # out[['actual_data']] <- observed_rates
    # out[[model_option]] <- l

    observed_rates[['pre_processed_data_forecasting']] <- data_pp_2
    out[['actual_data']] <- observed_rates
    out[[model_option]] <- l

  }



  return(out)
}


is_implausible_muhat <- function(muhat,
                                 global_muhat,
                                 max_ratio = 50,
                                 max_log_jump = 3.5,
                                 eps = 1e-12) {
  muhat <- as.numeric(muhat)
  global_muhat <- as.numeric(global_muhat)

  if (length(muhat) != length(global_muhat)) {
    return(TRUE)
  }

  if (any(!is.finite(muhat)) || any(muhat <= 0)) {
    return(TRUE)
  }

  log_ratio <- abs(log(pmax(muhat, eps) / pmax(global_muhat, eps)))
  if (any(!is.finite(log_ratio)) || any(log_ratio > log(max_ratio))) {
    return(TRUE)
  }

  log_muhat <- log(pmax(muhat, eps))
  if (length(log_muhat) > 1 && any(abs(diff(log_muhat)) > max_log_jump, na.rm = TRUE)) {
    return(TRUE)
  }

  FALSE
}

#' @export
fit_and_predict_separate_models <- function(data,
                                            N_groups,
                                            mortality_models_fit,
                                            years_fit,
                                            ages_fit,
                                            bias = 70,
                                            forecasting_horizon,
                                            max_ratio_implausible = 50,
                                            max_log_jump_implausible = 3.5) {
  out <- list()
  observed_rates <- list()
  failed_fits <- data.frame(group = 9999,
                            failed_fit = 99999,
                            model = "bla")


  mortality_model_list = list(
    lc = lc(link = "log"),
    apc = apc(link = "log"),
    rh = rh(link = "log")
  )


  data_pp_2 <- data_preprocessing(
    data = data,
    N_groups = N_groups,
    bias = bias,
    ages_fit = ages_fit,
    years_fit = c(years_fit, (last(years_fit) + 1):(last(years_fit) + forecasting_horizon))
    # scenario = scenario
  )

  for (model_option in names(mortality_models_fit)) {
    tmp_failed = FALSE
    l <- list()
    i = 1
    mortality_model = mortality_model_list[[model_option]]

    ## -------- added: global forecast for fallback --------
    global_fit <- mortality_models_fit[[model_option]]

    global_cv.arima.kt <- auto.arima(as.numeric(global_fit$kt), ic = "bic")

    if (model_option != "lc") {
      global_cv.arima.gc <- auto.arima(as.numeric(global_fit$gc), ic = "bic")
      global_gc.order <- unname(arimaorder(global_cv.arima.gc))
    } else {
      global_gc.order <- c(1, 1, 0)
    }

    global_forecast_try <- tryCatch(
      forecast(
        global_fit,
        gc.order = global_gc.order,
        kt.method = "iarima",
        kt.order = unname(arimaorder(global_cv.arima.kt)),
        h = forecasting_horizon
      ),
      error = function(e) {
        return("mrw")
      }
    )

    if (is.character(global_forecast_try)) {
      global_forecast <- forecast(global_fit, h = forecasting_horizon)
    } else {
      global_forecast <- forecast(
        global_fit,
        gc.order = global_gc.order,
        kt.method = "iarima",
        kt.order = unname(arimaorder(global_cv.arima.kt)),
        h = forecasting_horizon
      )
    }

    if (forecasting_horizon > 1) {
      global_muhat <- global_forecast$rates[, as.character(last(years_fit) + forecasting_horizon)]
    } else {
      global_muhat <- global_forecast$rates
    }
    ## -------- end added block --------

    for (i in 1:N_groups) {
      failure_logged <- FALSE

      popspecdata <- structure(
        list(
          Dxt = data_pp_2$subpopulations[[i]]$Dxt[as.character(ages_fit), as.character(years_fit)],
          Ext = data_pp_2$subpopulations[[i]]$Ext[as.character(ages_fit), as.character(years_fit)],
          ages = ages_fit,
          years = years_fit,
          type = 'central',
          series = 'male',
          label =  paste0("superpop", i)
        ),
        class = "StMoMoData"
      )

      error_check <- tryCatch(
        assign(
          paste0("mortality_model_fit", i),
          fit(
            mortality_model,
            data = popspecdata,
            start.ax = mortality_models_fit[[model_option]]$ax,
            start.bx = mortality_models_fit[[model_option]]$bx,
            start.kt = mortality_models_fit[[model_option]]$kt,
            start.gc = mortality_models_fit[[model_option]]$gc,
            years.fit = years_fit,
            ages.fit = ages_fit,
            verbose = FALSE
          )
        ),
        error = function(e) {
          return("failed")
        }
      )

      tmp_failed = FALSE

      if (is.character(error_check)) {
        assign(paste0("mortality_model_fit", i),
               mortality_models_fit[[model_option]])

        tmp_failed = TRUE

      } else{
        assign(paste0("mortality_model_fit", i), error_check)

        if (!get(paste0("mortality_model_fit", i))$conv) {
          assign(paste0("mortality_model_fit", i),
                 mortality_models_fit[[model_option]])

          tmp_failed = TRUE
        }
      }

      if (!get(paste0("mortality_model_fit", i))$conv) {
        assign(paste0("mortality_model_fit", i),
               mortality_models_fit[[model_option]])

        tmp_failed = TRUE
      }

      if (tmp_failed) {
        failed_fits = rbind(failed_fits,
                            data.frame(
                              group = i,
                              failed_fit = 1,
                              model = model_option
                            ))
        failure_logged <- TRUE
      }

      assign(paste0("cv.arima.kt_", i),
             auto.arima(as.numeric(get(
               paste0("mortality_model_fit", i)
             )$kt), ic = "bic"))

      if (model_option != "lc") {
        assign(paste0("cv.arima.gc_", i),
               auto.arima(as.numeric(get(
                 paste0("mortality_model_fit", i)
               )$gc), ic = "bic"))
        assign(paste0("gc.order_", i), unname(arimaorder(get(
          paste0("cv.arima.gc_", i)
        ))))

      } else{
        assign(paste0("gc.order_", i), c(1, 1, 0))
      }

      error_check <- tryCatch(
        assign(
          paste0("mortality_model_forecast", i) ,
          forecast(
            get(paste0("mortality_model_fit", i)),
            gc.order = get(paste0("gc.order_", i)),
            kt.method = "iarima",
            kt.order =
              unname(arimaorder(get(
                paste0("cv.arima.kt_", i)
              ))),
            h = forecasting_horizon
          )
        ),
        error = function(e) {
          return("mrw")
        }
      )

      if (is.character(error_check)) {
        assign(
          paste0("mortality_model_forecast", i) ,
          forecast(get(
            paste0("mortality_model_fit", i)
          ), h = forecasting_horizon)
        )
      } else{
        assign(
          paste0("mortality_model_forecast", i) ,
          forecast(
            get(paste0("mortality_model_fit", i)),
            gc.order = get(paste0("gc.order_", i)),
            kt.method = "iarima",
            kt.order =
              unname(arimaorder(get(
                paste0("cv.arima.kt_", i)
              ))),
            h = forecasting_horizon
          )
        )
      }

      if (forecasting_horizon > 1) {
        assign(paste0("muhat", i) ,
               get(paste0("mortality_model_forecast", i))$rates[, as.character(last(years_fit) + forecasting_horizon)])
      } else{
        assign(paste0("muhat", i) ,
               get(paste0("mortality_model_forecast", i))$rates)
      }

      ## -------- added: implausible forecast correction --------
      if (is_implausible_muhat(
        muhat = get(paste0("muhat", i)),
        global_muhat = global_muhat,
        max_ratio = max_ratio_implausible,
        max_log_jump = max_log_jump_implausible
      )) {
        assign(paste0("muhat", i), global_muhat)
        assign(paste0("mortality_model_fit", i), mortality_models_fit[[model_option]])

        if (!failure_logged) {
          failed_fits = rbind(failed_fits,
                              data.frame(
                                group = i,
                                failed_fit = 1,
                                model = model_option
                              ))
        }
      }
      ## -------- end added block --------

      l[[paste0("muhat", i)]] <- get(paste0("muhat", i))
      l[[paste0("model", i)]] <- get(paste0("mortality_model_fit", i))

      subgroups_forecast <- data_pp_2$subpopulations[[i]]

      observed_rates[[paste0('muxt_actual_', i)]] <- subgroups_forecast$Dxt[, (length(years_fit) + forecasting_horizon)] / subgroups_forecast$Ext[, (length(years_fit) + forecasting_horizon)]
    }

    observed_rates[['pre_processed_data_forecasting']] <- data_pp_2
    out[['actual_data']] <- observed_rates
    out[[model_option]] <- l
  }

  attr(out, "failed_fit") <-  failed_fits[-1, ]

  return(out)
}
# Sample size study ----

#' @export
sample_size_effect_on_predictions <- function(data,
                                              years_fit_basic,
                                              number_of_rolls,
                                              ages_fit,
                                              N_groups = 3,
                                              bias = 70,
                                              global_mortality_model = "lc") {
  out <- NULL


  for (years_ix in 1:number_of_rolls) {
    dt <- NULL
    setDT(dt)

    years_fit <- c(years_fit_basic, (last(years_fit_basic) + 1):(last(years_fit_basic) +
                                                                   years_ix))


    data_pp <- data_preprocessing(data = data,
                                  ages_fit = ages_fit,
                                  years_fit = years_fit)

    mortality_models_fit <-     fit_mortality_models(
      data_pp,
      years_fit - min(years_fit),
      ages_fit,
      global_mortality_model_candidates = global_mortality_model
    )

    credibility_model <- fit_and_predict_credibility_models(
      data = data,
      N_groups = N_groups,
      mortality_models_fit = mortality_models_fit,
      years_fit = years_fit - min(years_fit),
      ages_fit = ages_fit,
      forecasting_horizon = 1,
      bias = min(years_fit)
    )

    tmp <- credibility_model[[global_mortality_model]]

    dt[, `:=`(
      ages_fit = ages_fit,
      muhat_c = tmp$muhat1,
      theta = tmp$C1,
      muhat_mle = tmp$muhat1_MLE,
      label = "1"
    )]

    dt <- rbind(
      dt,
      data.table(
        ages_fit = ages_fit,
        muhat_c = tmp$muhat2,
        theta = tmp$C2,
        muhat_mle = tmp$muhat2_MLE,
        label = "2"
      )

    )

    dt <- rbind(
      dt,
      data.table(
        ages_fit = ages_fit,
        muhat_c = tmp$muhat3,
        theta = tmp$C3,
        muhat_mle = tmp$muhat3_MLE,
        label = "3"
      )

    )


    dt[, number_of_sample_years := length(years_fit)]

    out <- rbind(out, dt)

  }

  return(out)



}

#' @export
sample_size_effect_plotter <- function(ss_effect,
                                       color_group,
                                       ages_to_plot,
                                       label_to_plot) {
  plot_dt <- data.table::copy(ss_effect)[label == label_to_plot &
                                           ages_fit %in% ages_to_plot][order(ages_fit, number_of_sample_years)][, `:=`(
                                             diff_muhat_c   = c(NA_real_, diff(muhat_c)),
                                             diff_muhat_mle = c(NA_real_, diff(muhat_mle)),
                                             diff_global    = c(NA_real_, diff(muhat_mle / theta)),
                                             global =      muhat_mle / theta
                                           ), by = ages_fit]

  # choose as many distinct shapes as needed
  available_shapes <- c(16, 17, 15, 18, 8, 3, 7, 9, 10, 12, 13)
  shape_values <- setNames(available_shapes[seq_along(ages_to_plot)], as.character(sort(ages_to_plot)))

  ggplot(plot_dt, aes(x = number_of_sample_years, group = factor(ages_fit))) +
    geom_point(
      aes(y = muhat_c, shape = factor(ages_fit)),
      color = color_group,
      size = 3,
      na.rm = TRUE
    ) +
    geom_line(
      aes(y = muhat_mle),
      color = color_group,
      alpha = 0.35,
      linewidth = 1,
      na.rm = TRUE
    ) +
    geom_line(
      aes(y = global),
      color = "gray40",
      linewidth = 1,
      na.rm = TRUE
    ) +
    scale_shape_manual(values = shape_values, name = "Age") +
    labs(x = "", y = "", title = NULL) +
    theme_minimal(base_size = 14) + theme(legend.position = "none")
}

# Assess mortality model performance ----
#' @export
poisson_nll <- function(occurrence, exposure, muxt_hat) {
  out <- -sum(occurrence * log(exposure * muxt_hat) - exposure * muxt_hat -
                lfactorial(occurrence))

  return(out)

}


# -------------------------------------------------------------------
# Helpers for the Section 4.5 simulation study (main.R line 764 ff.)
# -------------------------------------------------------------------

#' @export
`%||%` <- function(x, y) {
  if (is.null(x))
    y
  else
    x
}
#' @export
.get_forecast_pp <- function(model_fit_and_prediction) {
  model_fit_and_prediction$actual_data$pre_processed_data_forecasting %||%
    model_fit_and_prediction$actual_data$full_pp_data
}
#' @export
.get_superpopulation <- function(pp) {
  pp$superpopulation %||% pp$superpopulation
}
#' @export
.get_subpopulation <- function(pp, i) {
  if (!is.null(pp$subpopulations)) {
    return(pp$subpopulations[[i]])
  }

  # legacy layout: group 1 lives in datahat, groups 2:3 in list_of_extra_exposures
  if (i == 1 && !is.null(pp$datahat)) {
    return(list(Dxt = pp$datahat$Dxt, Ext = pp$datahat$Ext))
  }

  pp$list_of_extra_exposures[[i - 1]]
}


#' @export
.make_selection_mask <- function(n,
                                 subset_rates = NULL,
                                 bias = 15) {
  out <- rep(TRUE, n)

  if (!is.null(subset_rates)) {
    out <- rep(FALSE, n)
    idx <- subset_rates - bias + 1L
    idx <- idx[idx >= 1L & idx <= n]
    out[idx] <- TRUE
  }

  out
}

#' @export
.poisson_deviance_cells <- function(dxt, ext, mu_hat) {
  mu_hat <- as.numeric(mu_hat)
  dxt <- as.numeric(dxt)
  ext <- as.numeric(ext)

  out <- rep(NA_real_, length(mu_hat))
  ok <- is.finite(dxt) & dxt >= 0 &
    is.finite(ext) & ext > 0 &
    is.finite(mu_hat) & mu_hat >= 0

  if (!any(ok))
    return(out)

  fxt <- dxt / ext

  zero_idx <- ok & dxt == 0
  pos_idx  <- ok & dxt > 0 & mu_hat > 0
  inf_idx  <- ok & dxt > 0 & mu_hat == 0

  out[zero_idx] <- 2 * ext[zero_idx] * mu_hat[zero_idx]
  out[pos_idx]  <- 2 * ext[pos_idx] * (mu_hat[pos_idx] -
                                         fxt[pos_idx] +
                                         fxt[pos_idx] * log(fxt[pos_idx] / mu_hat[pos_idx]))
  out[inf_idx] <- Inf

  out
}

#' @export
.assess_one_group <- function(mu_hat, actual_rate, dxt, ext, selected) {
  mu_hat <- as.numeric(mu_hat)
  actual_rate <- as.numeric(actual_rate)
  dxt <- as.numeric(dxt)
  ext <- as.numeric(ext)

  valid <- is.finite(actual_rate) &
    is.finite(mu_hat) & is.finite(ext) & ext > 0 & selected

  if (!any(valid)) {
    return(list(
      error = NA_real_,
      oos_deviance = NA_real_,
      n = 0L
    ))
  }

  obs <- actual_rate[valid]
  pred <- pmax(mu_hat[valid], .Machine$double.eps)
  d_use <- dxt[valid]
  e_use <- ext[valid]

  se <- (pred - obs)^2
  dev <- .poisson_deviance_cells(d_use, e_use, pred)

  list(
    error = mean(se, na.rm = TRUE),
    oos_deviance = mean(dev, na.rm = TRUE),
    n = sum(valid)
  )
}

#' Updated model assessment compatible with both old and new fit_and_predict_* outputs.
#' MARE follows the manuscript definition, including ARE = 0 when the observed rate is 0.
#' Poisson deviance is returned as an average over the selected cells, not a sum.
#' @export
model_assessment <- function(model_fit_and_prediction,
                             N_groups,
                             years_fit,
                             forecasting_horizon,
                             bias = 15,
                             subset_rates = NULL) {
  out <- list()
  tmp_actual_data <- model_fit_and_prediction[["actual_data"]]
  pp_forecast <- .get_forecast_pp(model_fit_and_prediction)
  forecast_col <- length(years_fit) + forecasting_horizon

  model_names <- intersect(c("lc", "apc", "rh"), names(model_fit_and_prediction))

  for (model_ix in model_names) {
    tmp_model <- model_fit_and_prediction[[model_ix]]
    l <- list()
    group_metrics <- vector("list", N_groups)

    for (i in seq_len(N_groups)) {
      actual_name <- paste0("muxt_actual_", i)
      fitted_name <- paste0("muhat", i)

      selected <- .make_selection_mask(
        n = length(tmp_actual_data[[actual_name]]),
        subset_rates = subset_rates,
        bias = bias
      )

      subpop <- .get_subpopulation(pp_forecast, i)

      group_metrics[[i]] <- .assess_one_group(
        mu_hat = tmp_model[[fitted_name]],
        actual_rate = tmp_actual_data[[actual_name]],
        dxt = subpop$Dxt[, forecast_col],
        ext = subpop$Ext[, forecast_col],
        selected = selected
      )

      l[[paste0("error_", i)]] <- group_metrics[[i]]$error
      l[[paste0("oos_deviance_", i)]] <- group_metrics[[i]]$oos_deviance
    }



    n_vec <- vapply(group_metrics, `[[`, numeric(1), "n")
    err_vec <- vapply(group_metrics, `[[`, numeric(1), "error")
    dev_vec <- vapply(group_metrics, `[[`, numeric(1), "oos_deviance")

    keep <- n_vec > 0 & is.finite(err_vec) & is.finite(dev_vec)

    l[["error"]] <- if (any(keep))
      stats::weighted.mean(err_vec[keep], n_vec[keep])
    else
      NA_real_
    l[["oos_deviance"]] <- if (any(keep))
      stats::weighted.mean(dev_vec[keep], n_vec[keep])
    else
      NA_real_

    out[[model_ix]] <- l
  }

  out
}

#' @export
model_assessment_superpop <- function(model_fit_and_prediction,
                                      N_groups,
                                      years_fit,
                                      forecasting_horizon,
                                      bias = 15,
                                      subset_rates = NULL) {
  out <- list()
  tmp_actual_data <- model_fit_and_prediction[["actual_data"]]
  pp_forecast <- .get_forecast_pp(model_fit_and_prediction)
  forecast_col <- length(years_fit) + forecasting_horizon
  superpop <- .get_superpopulation(pp_forecast)

  model_names <- intersect(c("lc", "apc", "rh"), names(model_fit_and_prediction))

  for (model_ix in model_names) {
    tmp_model <- model_fit_and_prediction[[model_ix]]
    selected <- .make_selection_mask(
      n = length(tmp_actual_data[["muxt_actual_superpop"]]),
      subset_rates = subset_rates,
      bias = bias
    )

    met <- .assess_one_group(
      mu_hat = tmp_model[["muhatsuperpop"]],
      actual_rate = tmp_actual_data[["muxt_actual_superpop"]],
      dxt = superpop$Dxt[, forecast_col],
      ext = superpop$Ext[, forecast_col],
      selected = selected
    )

    l <- list(error = met$error,
              oos_deviance = met$oos_deviance)

    for (i in seq_len(N_groups)) {
      l[[paste0("error_", i)]] <- NA_real_
      l[[paste0("oos_deviance_", i)]] <- NA_real_
    }

    out[[model_ix]] <- l
  }

  out
}

#' @export
model_assessment_full_mle <- function(model_fit_and_prediction,
                                      N_groups,
                                      years_fit,
                                      forecasting_horizon,
                                      bias = 15,
                                      subset_rates = NULL) {
  out <- list()
  tmp_actual_data <- model_fit_and_prediction[["actual_data"]]
  pp_forecast <- .get_forecast_pp(model_fit_and_prediction)
  forecast_col <- length(years_fit) + forecasting_horizon

  model_names <- intersect(c("lc", "apc", "rh"), names(model_fit_and_prediction))

  for (model_ix in model_names) {
    tmp_model <- model_fit_and_prediction[[model_ix]]
    l <- list()
    group_metrics <- vector("list", N_groups)

    for (i in seq_len(N_groups)) {
      actual_name <- paste0("muxt_actual_", i)
      fitted_name <- paste0("muhat", i, "_MLE")

      selected <- .make_selection_mask(
        n = length(tmp_actual_data[[actual_name]]),
        subset_rates = subset_rates,
        bias = bias
      )

      subpop <- .get_subpopulation(pp_forecast, i)

      group_metrics[[i]] <- .assess_one_group(
        mu_hat = tmp_model[[fitted_name]],
        actual_rate = tmp_actual_data[[actual_name]],
        dxt = subpop$Dxt[, forecast_col],
        ext = subpop$Ext[, forecast_col],
        selected = selected
      )

      l[[paste0("error_", i)]] <- group_metrics[[i]]$error
      l[[paste0("oos_deviance_", i)]] <- group_metrics[[i]]$oos_deviance
    }

    n_vec <- vapply(group_metrics, `[[`, numeric(1), "n")
    err_vec <- vapply(group_metrics, `[[`, numeric(1), "error")
    dev_vec <- vapply(group_metrics, `[[`, numeric(1), "oos_deviance")

    keep <- n_vec > 0 & is.finite(err_vec) & is.finite(dev_vec)

    l[["error"]] <- if (any(keep))
      stats::weighted.mean(err_vec[keep], n_vec[keep])
    else
      NA_real_
    l[["oos_deviance"]] <- if (any(keep))
      stats::weighted.mean(dev_vec[keep], n_vec[keep])
    else
      NA_real_

    out[[model_ix]] <- l
  }

  out
}

#' @export
.call_data_preprocessing <- function(data, N_groups, ages_fit, years_fit) {
  fn_formals <- names(formals(data_preprocessing))
  args <- list(
    data = data,
    N_groups = N_groups,
    ages_fit = ages_fit,
    years_fit = years_fit
  )

  if ("scenario" %in% fn_formals) {
    args$scenario <- 2
  }

  do.call(data_preprocessing, args)
}

#' @export
.call_fit_mortality_models <- function(data_pp,
                                       years_fit,
                                       ages_fit,
                                       global_mortality_model_candidates) {
  fn_formals <- names(formals(fit_mortality_models))
  args <- list(
    data_pp = data_pp,
    years_fit = years_fit,
    ages_fit = ages_fit
  )

  if ("global_mortality_model_candidates" %in% fn_formals) {
    args$global_mortality_model_candidates <- global_mortality_model_candidates
  } else if ("global_mortality_model" %in% fn_formals) {
    args$global_mortality_model <- global_mortality_model_candidates
  }

  if ("separate_exposures" %in% fn_formals) {
    args$separate_exposures <- FALSE
  }

  do.call(fit_mortality_models, args)
}

#' @export
.call_fit_predict <- function(fun,
                              data,
                              data_pp,
                              N_groups,
                              mortality_models_fit,
                              years_fit,
                              ages_fit,
                              forecasting_horizon,
                              bias) {
  fn_formals <- names(formals(fun))
  args <- list(
    data = data,
    N_groups = N_groups,
    mortality_models_fit = mortality_models_fit,
    years_fit = years_fit,
    ages_fit = ages_fit,
    forecasting_horizon = forecasting_horizon,
    bias = bias
  )

  if ("data_pp" %in% fn_formals) {
    args$data_pp <- data_pp
  }

  if ("scenario" %in% fn_formals) {
    args$scenario <- 2
  }

  do.call(fun, args)
}

#' @export
age_breaks_tbl <- function(age_min = 16,
                           age_max = 85,
                           width = 5) {
  starts <- seq(age_min, age_max, by = width)
  data.frame(
    interval_start = starts,
    interval_end = pmin(starts + width - 1L, age_max),
    stringsAsFactors = FALSE
  ) |>
    dplyr::mutate(age_breakets = paste(interval_start, interval_end, sep = "-"))
}

#' @export
.as_assessment_df <- function(perf,
                              seed_ix,
                              years_ix,
                              model_label,
                              forecasting_horizon,
                              age_breaket) {
  df <- do.call(rbind, lapply(perf, data.frame, stringsAsFactors = FALSE))
  df[["seed_ix"]] <- seed_ix
  df[["seed"]] <- seed_ix
  df[["model"]] <- model_label
  df[["predictor"]] <- rownames(df)
  df[["years_ix"]] <- years_ix
  df[["forecasting_horizon"]] <- forecasting_horizon
  df[["age_breakets"]] <- age_breaket
  rownames(df) <- NULL
  df
}

#' @export
.simulation_cell_setup <- function(seed_ix,
                                   years_ix,
                                   ages_fit = 15:85,
                                   years_fit_basic = 118:140,
                                   N_groups = 3,
                                   exposure_sup = 94500,
                                   exposure_sub = c(5000, 500),
                                   global_mortality_model_candidates = c("lc", "apc", "rh"),
                                   hmd_user = "",
                                   hmd_pass = "") {
  dt_0 <- data_generator_hmd_lt(
    seed_input = seed_ix,
    hmd_username = hmd_user,
    hmd_password = hmd_pass,
    exposure_sup = exposure_sup,
    exposure_sub = exposure_sub
  )

  years_fit <- c(
    years_fit_basic,
    seq.int(max(years_fit_basic) + 1L, max(years_fit_basic) + years_ix)
  )

  data_pp <- .call_data_preprocessing(
    data = dt_0,
    N_groups = N_groups,
    ages_fit = ages_fit,
    years_fit = years_fit
  )

  mortality_models_fit <- .call_fit_mortality_models(
    data_pp = data_pp,
    years_fit = years_fit - min(years_fit),
    ages_fit = ages_fit,
    global_mortality_model_candidates = global_mortality_model_candidates
  )

  list(
    dt_0 = dt_0,
    years_fit = years_fit,
    data_pp = data_pp,
    mortality_models_fit = mortality_models_fit
  )
}

#' @export
.summarise_separate_failures <- function(separate_model,
                                         predictors = c("lc", "apc", "rh")) {
  failed_fits <- attr(separate_model, "failed_fit", exact = TRUE)

  out <- data.table::data.table(
    predictor = predictors,
    n_failed_separate = 0L
  )

  if (is.null(failed_fits) || NROW(failed_fits) == 0L) {
    return(out[])
  }

  failed_fits <- data.table::as.data.table(failed_fits)
  failed_fits <- failed_fits[
    !is.na(model) & !is.na(failed_fit) & failed_fit > 0,
    .(n_failed_separate = as.integer(sum(failed_fit > 0))),
    by = .(predictor = as.character(model))
  ]

  out[failed_fits, on = "predictor", n_failed_separate := i.n_failed_separate]
  out[is.na(n_failed_separate), n_failed_separate := 0L]
  out[]
}

#' @export
.summarise_insample_bic <- function(mortality_models_fit,
                                    data_pp) {
  out <- lapply(names(mortality_models_fit), function(model_ix) {
    tmp <- compute_insample_bic_scratch(
      model_fit = mortality_models_fit[[model_ix]],
      data_stmomo = data_pp$superpopulation,
      model_label = model_ix
    )

    data.table::setDT(tmp)
    data.table::setnames(tmp, "model", "predictor")
    tmp[]
  })

  data.table::rbindlist(out, use.names = TRUE)
}

#' @export
run_simulation_study_cell <- function(seed_ix,
                                      years_ix,
                                      ages_fit = 15:85,
                                      years_fit_basic = 118:140,
                                      N_groups = 3,
                                      forecasting_horizon = 1,
                                      ageb_width = 5,
                                      exposure_sup = 94500,
                                      exposure_sub = c(5000, 500),
                                      global_mortality_model = c("lc", "apc", "rh"),
                                      hmd_user = "",
                                      hmd_pass = "") {

  dt_0 <- data_generator_hmd_lt(
    seed_input = seed_ix,
    hmd_username = hmd_user,
    hmd_password = hmd_pass,
    exposure_sup = exposure_sup,
    exposure_sub = exposure_sub
  )

  years_fit <- c(
    years_fit_basic,
    seq.int(max(years_fit_basic) + 1L, max(years_fit_basic) + years_ix)
  )

  data_pp <- data_preprocessing(
    data = dt_0,
    N_groups = N_groups,
    ages_fit = ages_fit,
    years_fit = years_fit
  )

  mortality_models_fit <- fit_mortality_models(
    data_pp = data_pp,
    years_fit = years_fit - min(years_fit),
    ages_fit = ages_fit,
    global_mortality_model_candidates = global_mortality_model
  )

  years_fit_shifted <- years_fit - min(years_fit)

  total_model <- fit_and_predict_total_model(
    data = dt_0,
    N_groups = N_groups,
    mortality_models_fit = mortality_models_fit,
    years_fit = years_fit_shifted,
    ages_fit = ages_fit,
    forecasting_horizon = forecasting_horizon,
    bias = min(years_fit)
  )

  credibility_model <- fit_and_predict_credibility_models(
    data = dt_0,
    N_groups = N_groups,
    mortality_models_fit = mortality_models_fit,
    years_fit = years_fit_shifted,
    ages_fit = ages_fit,
    forecasting_horizon = forecasting_horizon,
    bias = min(years_fit)
  )

  separate_model <- fit_and_predict_separate_models(
    data = dt_0,
    N_groups = N_groups,
    mortality_models_fit = mortality_models_fit,
    years_fit = years_fit_shifted,
    ages_fit = ages_fit,
    forecasting_horizon = forecasting_horizon,
    bias = min(years_fit)
  )

  bic_dt <- .summarise_insample_bic(
    mortality_models_fit = mortality_models_fit,
    data_pp = data_pp
  )
  bic_dt[, `:=`(
    seed_ix = as.integer(seed_ix),
    years_ix = as.integer(years_ix),
    n_years_fit = length(years_fit)
  )]

  failed_dt <- .summarise_separate_failures(
    separate_model = separate_model,
    predictors = names(mortality_models_fit)
  )

  extra_cols <- merge(
    bic_dt,
    failed_dt,
    by = "predictor",
    all.x = TRUE
  )

  brks <- age_breaks_tbl(age_min = min(ages_fit),
                         age_max = max(ages_fit),
                         width = ageb_width)

  out <- lapply(seq_len(nrow(brks)), function(j) {
    subset_rates <- brks$interval_start[j]:brks$interval_end[j]
    interval_name <- brks$age_breakets[j]

    rbind(
      .as_assessment_df(
        perf = model_assessment(
          model_fit_and_prediction = total_model,
          N_groups = N_groups,
          years_fit = years_fit,
          forecasting_horizon = forecasting_horizon,
          subset_rates = subset_rates,
          bias = min(ages_fit)
        ),
        seed_ix = seed_ix,
        years_ix = years_ix,
        model_label = "total",
        forecasting_horizon = forecasting_horizon,
        age_breaket = interval_name
      ),
      .as_assessment_df(
        perf = model_assessment(
          model_fit_and_prediction = credibility_model,
          N_groups = N_groups,
          years_fit = years_fit,
          forecasting_horizon = forecasting_horizon,
          subset_rates = subset_rates,
          bias = min(ages_fit)
        ),
        seed_ix = seed_ix,
        years_ix = years_ix,
        model_label = "credibility",
        forecasting_horizon = forecasting_horizon,
        age_breaket = interval_name
      ),
      .as_assessment_df(
        perf = model_assessment_full_mle(
          model_fit_and_prediction = credibility_model,
          N_groups = N_groups,
          years_fit = years_fit,
          forecasting_horizon = forecasting_horizon,
          subset_rates = subset_rates,
          bias = min(ages_fit)
        ),
        seed_ix = seed_ix,
        years_ix = years_ix,
        model_label = "full_mle",
        forecasting_horizon = forecasting_horizon,
        age_breaket = interval_name
      ),
      .as_assessment_df(
        perf = model_assessment_superpop(
          model_fit_and_prediction = total_model,
          N_groups = N_groups,
          years_fit = years_fit,
          forecasting_horizon = forecasting_horizon,
          subset_rates = subset_rates,
          bias = min(ages_fit)
        ),
        seed_ix = seed_ix,
        years_ix = years_ix,
        model_label = "superpopulation",
        forecasting_horizon = forecasting_horizon,
        age_breaket = interval_name
      ),
      .as_assessment_df(
        perf = model_assessment(
          model_fit_and_prediction = separate_model,
          N_groups = N_groups,
          years_fit = years_fit,
          forecasting_horizon = forecasting_horizon,
          subset_rates = subset_rates,
          bias = min(ages_fit)
        ),
        seed_ix = seed_ix,
        years_ix = years_ix,
        model_label = "separate",
        forecasting_horizon = forecasting_horizon,
        age_breaket = interval_name
      )
    )
  })

  out <- data.table::as.data.table(dplyr::bind_rows(out))
  out <- merge(out, extra_cols, by = c("predictor", "seed_ix", "years_ix"), all.x = TRUE)
  data.table::setorder(out, seed_ix, years_ix, forecasting_horizon, age_breakets, predictor, model)
  out[]
}


# Convergence and BIC ----
#' @export
empty_failed_fits_dt <- function() {
  data.table::data.table(
    seed_ix = integer(),
    years_ix = integer(),
    predictor = character(),
    group = integer(),
    failed_fit = integer()
  )
}
#' @export
run_separate_failures_cell <- function(seed_ix,
                                       years_ix,
                                       ages_fit = 15:85,
                                       years_fit_basic = 118:140,
                                       N_groups = 3,
                                       forecasting_horizon = 1,
                                       exposure_sup = 94500,
                                       exposure_sub = c(5000, 500),
                                       global_mortality_model = c("lc", "apc", "rh"),
                                       hmd_user = "",
                                       hmd_pass = "") {
  dt_0 <- data_generator_hmd_lt(
    seed_input = seed_ix,
    hmd_username = hmd_user,
    hmd_password = hmd_pass,
    exposure_sup = exposure_sup,
    exposure_sub = exposure_sub
  )

  years_fit <- c(
    years_fit_basic,
    seq.int(max(years_fit_basic) + 1L, max(years_fit_basic) + years_ix)
  )

  data_pp <- data_preprocessing(
    data = dt_0,
    N_groups = N_groups,
    ages_fit = ages_fit,
    years_fit = years_fit
  )

  mortality_models_fit <- fit_mortality_models(
    data_pp = data_pp,
    years_fit = years_fit - min(years_fit),
    ages_fit = ages_fit,
    global_mortality_model_candidates = global_mortality_model
  )

  separate_model <- fit_and_predict_separate_models(
    data = dt_0,
    N_groups = N_groups,
    mortality_models_fit = mortality_models_fit,
    years_fit = years_fit - min(years_fit),
    ages_fit = ages_fit,
    forecasting_horizon = forecasting_horizon,
    bias = min(years_fit)
  )

  failed_fits <- attr(separate_model, "failed_fit", exact = TRUE)

  if (is.null(failed_fits) || NROW(failed_fits) == 0L) {
    return(empty_failed_fits_dt())
  }

  failed_fits <- data.table::as.data.table(failed_fits)
  failed_fits <- failed_fits[!is.na(model)]

  if (!nrow(failed_fits)) {
    return(empty_failed_fits_dt())
  }

  failed_fits[, `:=`(
    seed_ix = as.integer(seed_ix),
    years_ix = as.integer(years_ix),
    predictor = as.character(model)
  )]
  failed_fits[, model := NULL]

  data.table::setcolorder(failed_fits,
                          c("seed_ix", "years_ix", "predictor", "group", "failed_fit"))
  data.table::setorder(failed_fits, seed_ix, years_ix, predictor, group)
  failed_fits[]
}

#' @export
bind_separate_failures <- function(x) {
  if (length(x) == 0L) {
    return(empty_failed_fits_dt())
  }

  out <- data.table::rbindlist(x, use.names = TRUE, fill = TRUE)

  if (nrow(out) == 0L) {
    return(empty_failed_fits_dt())
  }

  data.table::setorder(out, seed_ix, years_ix, predictor, group)
  out[]
}


# Plots ----
#' @export
weights_plotter <- function(credibility_model,
                            N_groups,
                            ages_fit,
                            predictor = "lc",
                            ages_breaks = c(50, 60, 70, 80, 90)) {
  tmp <- credibility_model[[predictor]]



  # Build long table of all Z_i columns stacked in i = 1..N_groups order,
  # and within each i keep the original vector order.
  dt_Zs <- rbindlist(lapply(1:N_groups, function(i) {
    wname <- paste0("Z_", i)
    data.table(Zs = 1 - tmp[[wname]],
               label = wname,
               ages.code = ages_fit)
  }), use.names = TRUE)

  # Match original factor(levels=...) exactly
  dt_Zs[, label := factor(label, levels = c("Z_2", "Z_3", "Z_1"))]

  text_size <- 28

  # Keep ggplot part unchanged (ggplot works fine with data.table)
  ggplot(data = dt_Zs, aes(x = ages.code, y = Zs)) +
    geom_point(aes(colour = label), size = 3, alpha = .7) +
    theme_bw() +
    scale_x_continuous(breaks = ages_breaks) +
    theme(
      text = element_text(size = text_size),
      legend.position = "inside",
      legend.justification = c("right", "bottom"),
      legend.background = element_blank(),
      legend.key = element_blank()
    ) +
    ylab("") +
    xlab("") +
    labs(color = "")
}
#' @export
all_thetas_plotter <- function(credibility_model,
                               N_groups,
                               ages_fit,
                               predictor = "lc") {
  tmp <- credibility_model[[predictor]]


  cis <- c()
  names_to_add <- c()


  for (i in 1:N_groups) {
    cs_to_subset <- paste0("C", i)
    names_to_add <- c(names_to_add, cs_to_subset)

    cis <- c(cis, tmp[[cs_to_subset]])

  }

  dt_cis <- data.frame(
    cis = cis,
    label = factor(rep(names_to_add, each = length(ages_fit)), levels =
                     c('C2', 'C3', 'C1')),
    ages.code = rep(ages_fit, N_groups)
  )



  text_size = 28


  dt_cis %>%
    ggplot(aes(x = ages.code, y = cis)) +
    geom_point(aes(colour = label), size = 3, alpha = .7) +
    ylim(0.5, 2.1) +
    theme_bw() +
    scale_x_continuous(breaks = (1:9) * 10) +
    theme(
      text = element_text(size = 28),
      legend.position = "inside",
      # must now be "inside"
      # legend.position.inside = c(0.85, 0.85),              # numeric position
      legend.justification = c("right", "top"),
      legend.key.height = unit(1.5, "lines"),
      legend.background = element_blank(),
      # removes legend box
      legend.key = element_blank()
    ) +
    ylab("") +
    xlab("") +
    labs(color = "")




}


#' @export
exposure_plotter <- function(data_pp, predictor = "lc") {
  E1 <- apply(data_pp$datahat$Ext, 1, sum, na.rm = T)
  E2 <-  apply(data_pp$list_of_extra_exposures[[1]]$Ext, 1, sum, na.rm =
                 T)
  E3 <-  apply(data_pp$list_of_extra_exposures[[2]]$Ext, 1, sum, na.rm =
                 T)


  dt <- data.frame(
    exposure = c(E1, E2, E3),
    Age = rep(names(E1), 3),
    group = rep(c("p", "p1", "p2"), each = length(names(E1)))
  )

  p <- dt %>%
    mutate(group = recode(group, p = "Group 0", p1 = "Group 1", p2 = "Group 2")) %>%
    ggplot(aes(x = as.integer(Age), y = exposure)) +
    geom_line(aes(colour = group), size = 1.5) +
    theme_bw() +
    theme(text = element_text(size = 28), legend.position = "top") +
    scale_color_manual(values = c(
      "Group 0" = "#a71429",
      "Group 1" = "#4169E1",
      "Group 2" = "#2E8B57"
    )) +
    ylab("") +
    xlab("") +
    labs(color = "")

  return(p)


}


#' @export
variance_plotter <- function(credibility_model, predictor = "lc", ages_fit) {
  #
  v1 <- credibility_model[[predictor]]$varthetax_1
  v2 <-  credibility_model[[predictor]]$varthetax_2
  v3 <-  credibility_model[[predictor]]$varthetax_3


  dt <- data.frame(
    variance = c(v1, v2, v3),
    Age = rep(as.character(ages_fit), 3),
    group = rep(c("p", "p1", "p2"), each = length(as.character(ages_fit)))
  )


  p <- dt %>%
    mutate(group = recode(group, p1 = "Group 1", p2 = "Group 2", p = "Group 3")) %>%
    ggplot(aes(
      x = as.integer(Age),
      y = (variance),
      color = interaction(group),
      shape = interaction(group)
    )) +
    geom_point(size = 3) +
    theme_bw() +
    theme(
      text = element_text(size = 28),
      legend.position = "inside",
      # must now be "inside"
      # legend.position.inside = c(0.85, 0.85),              # numeric position
      legend.justification = c("right", "top"),
      legend.background = element_blank(),
      # removes legend box
      legend.key = element_blank()
    ) +
    scale_color_manual(
      values = c(
        "Group 1" = "#4169E1",
        "Group 2" = "#2E8B57",
        "Group 3" = "#a71429"
      ),
      labels = c(expression(hat(Var)(Theta[x]^1)), expression(hat(Var)(Theta[x]^2)), expression(hat(Var)(Theta[x]^3)))
    ) +
    scale_shape_manual(
      values = c(
        "Group 1" = 16,
        "Group 2" = 17,
        "Group 3" = 18
      ),
      labels = c(expression(hat(Var)(Theta[x]^1)), expression(hat(Var)(Theta[x]^2)), expression(hat(Var)(Theta[x]^3)))
    ) +
    labs(color = "Group", shape = "Group") +  # Customize the legend titles
    ylab("") +
    xlab("") +
    # ylim(-0.0001, 1)+
    labs(color = "", shape = "")
  guides(color = guide_legend(override.aes = list(shape = c(16, 17, 18))))  # Ensures the legend reflects both shape and color


  return(p)


}

#' @export
thetas_plotter <- function(credibility_model,
                           subgroup,
                           ages_fit,
                           years_fit,
                           mortality_models_fit,
                           predictor = "lc",
                           chosen_age = 70) {
  tmp <- credibility_model[[predictor]]

  tmp_fit <- mortality_models_fit[[predictor]]


  mortality_model_forecast <- forecast(
    tmp_fit,
    h = 1 #doesnt really matter here. I want the fitted values
  )

  muxt_hat <- mortality_model_forecast$fitted


  chosen_age = as.character(chosen_age)

  selection_ratios <- paste0("Fxt_", subgroup)



  Cs_frame <- data.frame(y = (tmp[[selection_ratios]] * 1 / muxt_hat)[chosen_age, ], x =
                           years_fit) %>%
    filter(y > 0)
  selection_ratios <- paste0("C", subgroup)
  C1 <- tmp[[selection_ratios]]
  selection_ratios <- paste0("varthetax_", subgroup)
  varthetax_1 <- tmp[[selection_ratios]]

  names(varthetax_1) <- names(C1)

  text_size = 28

  Cs_frame %>%
    ggplot(aes(x = x, y = y)) +
    geom_point(size = 3) +
    geom_hline(
      yintercept = C1[chosen_age],
      linetype = "dotted",
      linewidth = 1.2,
      color = "#a71429"
    ) +
    geom_hline(
      yintercept = C1[chosen_age] + sqrt(varthetax_1[chosen_age]),
      linewidth = 1.2,
      color = "#4169E1",
      linetype = "dotted"
    ) +
    geom_hline(
      yintercept = C1[chosen_age] - sqrt(varthetax_1[chosen_age]),
      linewidth = 1.2,
      color = "#4169E1",
      linetype = "dotted"
    ) +
    theme_bw() +
    scale_x_continuous(breaks = c(10, 20, 30, 40, 50)) +
    theme(text = element_text(size = text_size),
          legend.position = "top") +
    ylab("") +
    xlab("") +
    labs(color = "")

}



#' @export
plot_credibility_mle_global <- function(credibility_model,
                                        N_groups,
                                        ages_fit,
                                        predictor = "lc",
                                        ages_breaks = c(50, 60, 70, 80, 90)) {
  tmp <- credibility_model[[predictor]]



  # Build long table of all Z_i columns stacked in i = 1..N_groups order,
  # and within each i keep the original vector order.
  dt_mles <- rbindlist(lapply(1:N_groups, function(i) {
    wname <- paste0("muhat", i, "_MLE")
    data.table(mortality = tmp[[wname]],
               label = wname,
               ages.code = ages_fit)
  }), use.names = TRUE)

  # any of them is good
  dt_global <- data.table(
    mortality = tmp$muhat1_MLE / tmp$C1,
    label = "global_mortality_model",
    ages.code = ages_fit
  )



  dt_credibility <- rbindlist(lapply(1:N_groups, function(i) {
    wname <- paste0("muhat", i)
    data.table(mortality = tmp[[wname]],
               label = wname,
               ages.code = ages_fit)
  }), use.names = TRUE)

  dt <- rbind(dt_mles, dt_global, dt_credibility)

  dt_plot <- copy(dt)

  dt_plot[, model :=
            fifelse(grepl("_MLE$", label),
                    "MLE",
                    fifelse(
                      grepl("^muhat[0-9]+$", label),
                      "credibility",
                      fifelse(label == "global_mortality_model", "global", NA_character_)
                    ))]

  dt_plot[, group :=
            fifelse(grepl("1", label), "1", fifelse(grepl("2", label), "2", fifelse(grepl("3", label), "3", NA_character_)))]
  dt_plot[, log_mortality := log(mortality)]
  cols <- c(
    "1" = "#a71429",
    "2" = "#4169E1",
    "3" = "#2E8B57",
    "global" = "darkgray"
  )

  p <- ggplot() +

    geom_line(
      data = dt_plot[model == "MLE"],
      aes(
        x = ages.code,
        y = log_mortality,
        color = group,
        group = label
      ),
      alpha = 0.9
    ) +

    geom_line(
      data = dt_plot[model == "global"],
      aes(x = ages.code, y = log_mortality),
      color = cols["global"],
      linewidth = 1.5
    ) +

    geom_point(
      data = dt_plot[model == "credibility"],
      aes(x = ages.code, y = log_mortality, color = group),
      size = 1.5
    ) +

    scale_color_manual(values = cols) +
    guides(color = "none") +
    labs(x = "", y = "") +
    theme_bw() +
    theme(text = element_text(size = 28))

  return(p)

}

#' @export
binomial_simulator <- function(starting_exposure, probabilities) {
  starting_period_lives <- starting_exposure
  starting_period_events <- NULL

  for (ix in 1:length(probabilities)) {
    events <- rbinom(1, size = starting_period_lives[ix], p = probabilities[ix])

    starting_period_lives <- c(starting_period_lives, starting_period_lives[ix] -
                                 events)

    starting_period_events <- c(starting_period_events, events)

  }


  l <- list(Lx = starting_period_lives[-length(starting_period_lives)], Dx =
              starting_period_events)

  return(l)


}

#' @export
average_exposure_and_deaths <- function(raw_exposure) {
  out <- raw_exposure
  events <- NULL



  # One would need some adjustments for older ages and age 0. Since we
  # do not care about those we simply disgregard them.
  for (ix in 1:(length(raw_exposure)) - 1) {
    out[ix] <- raw_exposure[ix] / 3 + raw_exposure[ix + 1] / 6 + raw_exposure[ix] /
      6 + raw_exposure[ix + 1] / 3

    events[ix] <- (raw_exposure[ix] - raw_exposure[ix + 1]) / 2
  }



  events <- c(events, last(raw_exposure))
  l <- list(Ext = as.numeric(out), Dxt = as.numeric(events))
  return(l)

}

#' @export
plot_predicted_event_trends <- function(mortality_models_fit,
                                        data_pp,
                                        # subgroup=1,
                                        model_option = "apc") {
  #
  mortality_model_fit <- mortality_models_fit[[model_option]]


  mortality_model_forecast <- forecast(mortality_model_fit, h = 1)


  # if (subgroup != 1) {
  #   tmp_exp <- data_pp$list_of_extra_exposures[[subgroup-1]]$Ext
  #
  # } else{
  #   tmp_exp <- data_pp$datahat$Ext
  # }

  tmp_exp <- data_pp$datahat$Ext

  tmp1 <- data.frame(
    y = apply(tmp_exp * mortality_model_forecast$fitted, 1, sum, na.rm = T),
    x = rownames(tmp_exp),
    group = "Group 3"
  )

  tmp_exp <- data_pp$list_of_extra_exposures[[1]]$Ext

  tmp2 <- data.frame(
    y = apply(tmp_exp * mortality_model_forecast$fitted, 1, sum, na.rm = T),
    x = rownames(tmp_exp),
    group = "Group 1"
  )

  tmp_exp <- data_pp$list_of_extra_exposures[[2]]$Ext

  tmp3 <- data.frame(
    y = apply(tmp_exp * mortality_model_forecast$fitted, 1, sum, na.rm = T),
    x = rownames(tmp_exp),
    group = "Group 2"
  )

  dt <- rbind(tmp1, tmp2, tmp3)

  p <- dt %>%
    ggplot(aes(
      x = as.numeric(x),
      y = log(y),
      color = group
    )) +
    geom_line(size = 2, alpha = .7) +
    theme_bw() +
    # scale_x_continuous(breaks = ages_breaks) +
    scale_color_manual(
      values = c(
        "Group 1" = "#4169E1",
        "Group 2" = "#2E8B57",
        "Group 3" = "#a71429"
      ),
      labels = c(expression(D[x, .]^1), expression(D[x, .]^2), expression(D[x, .]^3))
    ) +
    theme(
      text = element_text(size = 28),
      legend.position = "inside",
      # must now be "inside"
      # legend.position.inside = c(0.85, 0.85),              # numeric position
      legend.justification = c("right", "bottom"),
      legend.background = element_blank(),
      # removes legend box
      legend.key = element_blank()
    ) +
    scale_fill_discrete(guide = "none") +
    ylab("") +
    # ylab(expression(log(hat(D)[x* "."]^i ))) +
    xlab("") +
    labs(color = "")

  return(p)

}
#' @export
fill_missing_rows <- function(df) {
  # First, create a vector of ages from 0 to 110
  all_ages <- 0:110

  # Initialize an empty list to store results
  df_filled <- list()

  # For each cohort, ensure all ages are represented in the data
  for (cohort in unique(df$Cohort)) {
    # Filter data for the current cohort
    cohort_data <- df %>% filter(Cohort == cohort)

    # Get the set of ages already present for the current cohort
    present_ages <- unique(cohort_data$Age)

    # Find missing ages by subtracting present ages from all ages (0 to 110)
    missing_ages <- setdiff(all_ages, present_ages)

    # If there are missing ages, create new rows for the missing ages
    if (length(missing_ages) > 0) {
      # Create new rows for the missing ages (using data from Cohort 0 for 'qx' values)
      missing_rows <- data.frame(
        Age = missing_ages,
        Period = missing_ages + cohort,
        # Correct Period as Cohort + Age
        Cohort = cohort,
        qx = df$qx[df$Cohort == 0 &
                     df$Age %in% missing_ages]  # Match missing age rows from Cohort 0
      )

      # Combine the existing cohort rows with the missing rows (do not duplicate)
      df_filled[[length(df_filled) + 1]] <- bind_rows(
        cohort_data,
        # Keep all existing rows in cohort_data
        missing_rows  # Add only the missing rows
      )
    } else {
      # If no ages are missing, just add the current cohort's data to the list
      df_filled[[length(df_filled) + 1]] <- cohort_data
    }
  }

  # Combine all the dataframes in the list into one dataframe
  df_filled <- bind_rows(df_filled)

  return(df_filled)
}


#' @export
k_years_smoothing <- function(x, k = 5) {
  ll <- 1:length(x)

  out <- sapply(ll, average_smoother, x = x, k = k)

  return(out)

}

#' @export
average_smoother <- function(pos, x, k) {
  out <- x[max(pos - ((k - 1) / 2), 1):min(pos + ((k - 1) / 2), length(x))]
  return(mean(out))
}

#' @export
plot_smoothed_unsmoothed  <- function(data,
                                      data_pp,
                                      N_groups,
                                      mortality_models_fit,
                                      years_fit,
                                      ages_fit,
                                      bias =
                                        70,
                                      scenario = 2,
                                      model_2_plot = "lc",
                                      smoothing_years =
                                        5,
                                      forecasting_horizon) {
  out <- list()
  observed_rates <- list()

  data_pp_2 <- data_preprocessing_scenario_2(
    data = data,
    N_groups = N_groups,
    bias = bias,
    ages_fit = ages_fit,
    years_fit = c(years_fit, (last(years_fit) + 1):(last(years_fit) + forecasting_horizon))
    # scenario = scenario
  )

  for (model_option in names(mortality_models_fit)) {
    l <- list()

    mortality_model_fit <- mortality_models_fit[[model_option]]

    year.predict <- max(years_fit) + forecasting_horizon

    cv.arima.kt <- auto.arima(as.numeric(mortality_model_fit$kt), ic = "bic")

    if (model_option != "lc") {
      cv.arima.gc <- auto.arima(as.numeric(mortality_model_fit$gc), ic = "bic")
      gc.order <- unname(arimaorder(cv.arima.gc))
    } else{
      gc.order <- c(1, 1, 0)
    }

    #

    mortality_model_forecast <- forecast(
      mortality_model_fit,
      kt.method = "iarima",
      gc.order = gc.order,
      kt.order = unname(arimaorder(cv.arima.kt)),
      h = forecasting_horizon
    )

    muxt_hat <- mortality_model_forecast$fitted

    if (forecasting_horizon > 1) {
      muxt_hat_predicted <- mortality_model_forecast$rates[, as.character(last(years_fit) + forecasting_horizon)]


    } else{
      muxt_hat_predicted <- mortality_model_forecast$rates


    }

    tmp_occ <- data_pp$datahat$Dxt[as.character(ages_fit), as.character(years_fit)]
    tmp_exp <- data_pp$datahat$Ext[as.character(ages_fit), as.character(years_fit)]

    C1 <- apply(tmp_occ, 1, sum, na.rm = T) / apply(tmp_exp * muxt_hat, 1, sum, na.rm =
                                                      T)
    #
    l[['C1_unsmoothed']] <- C1
    fit1 <- rpart(C1 ~ ages_fit)
    C1 <- predict(fit1, data.frame(ages_fit))


    # C1 <- k_years_smoothing(C1, k= smoothing_years)

    Fxt_1 <- tmp_occ / tmp_exp

    tmpvar1 <- ((Fxt_1 - muxt_hat)^2 - (muxt_hat / tmp_exp))
    # tmpvar1[tmpvar1 == 1] <- NA
    varthetax_1 <- pmax(0,
                        apply(tmpvar1, 1, sum, na.rm = T) / apply(muxt_hat^2, 1, sum, na.rm = T),
                        na.rm = T)

    # varthetax_1 <- k_years_smoothing(varthetax_1, k= smoothing_years)
    l[['varthetax_1_unsmoothed']] <- varthetax_1
    fit1 <- rpart(varthetax_1 ~ ages_fit)
    varthetax_1 <- predict(fit1, data.frame(ages_fit))
    l[['varthetax_1_smoothed']] <- varthetax_1

    Z_1 <- 1 / (1 + varthetax_1 * apply(tmp_exp * muxt_hat, 1, sum, na.rm =
                                          T))

    subgroups_data <- data_pp$list_of_extra_exposures


    l[['C1_smoothed']] <- C1



    # term1 <- 2*(muxt_hat_predicted^2)*(varthetax_1)+muxt_hat_predicted^2
    #
    # term2 <- varthetax_1*(apply((tmp_exp^2) *( muxt_hat)^2, 1, sum, na.rm =
    #                   T)/ apply(tmp_exp * muxt_hat, 1, sum, na.rm =
    #                     T)^2)+(1/ apply(tmp_exp * muxt_hat, 1, sum, na.rm =
    #                                     T))
    #
    # l[['msep_1']] <- term1+term2


    # observed_rates[['muxt_actual_1']] <- data_pp_2$datahat$Dxt[, (length(years_fit) +
    #                                                                 1):(length(years_fit) + forecasting_horizon)] / data_pp_2$datahat$Ext[, (length(years_fit) +
    #                                                                                                                                            1):(length(years_fit) + forecasting_horizon)]
    observed_rates[['muxt_actual_1']] <- data_pp_2$datahat$Dxt[, (length(years_fit) + forecasting_horizon)] / data_pp_2$datahat$Ext[, (length(years_fit) + forecasting_horizon)]

    for (i in 2:N_groups) {
      tmp_occ <- subgroups_data[[i - 1]]$Dxt[as.character(ages_fit), as.character(years_fit)]
      tmp_exp <- subgroups_data[[i - 1]]$Ext[as.character(ages_fit), as.character(years_fit)]

      assign(
        paste0("C", i),
        apply(tmp_occ, 1, sum, na.rm = T) / apply(tmp_exp *
                                                    muxt_hat, 1, sum, na.rm = T)
      )


      l[[paste0('C', i, '_unsmoothed')]] <- get(paste0("C", i))
      fit1 <- rpart(get(paste0("C", i)) ~ ages_fit)
      assign(paste0("C", i), predict(fit1, data.frame(ages_fit)))

      assign(paste0("Fxt_", i), tmp_occ / tmp_exp)

      assign(paste0("tmpvar", i), ((get(
        paste0("Fxt_", i)
      ) - muxt_hat)^2 - (muxt_hat / tmp_exp)))

      tmp <- get(paste0("tmpvar", i))
      # tmp[tmp == 1] = NA
      assign(paste0("tmpvar", i), tmp)
      assign(paste0("varthetax_", i), pmax(
        0,
        apply(get(paste0("tmpvar", i)), 1, mean, na.rm =
                T) / apply(muxt_hat^2, 1, mean, na.rm =
                             T),
        na.rm = T
      ))

      l[[(paste0("varthetax_", i, "_unsmoothed"))]] <- get(paste0("varthetax_", i))
      fit1 <- rpart(get(paste0("varthetax_", i)) ~ ages_fit)
      assign(paste0("varthetax_", i), predict(fit1, data.frame(ages_fit)))
      l[[(paste0("varthetax_", i, "_smoothed"))]] <- get(paste0("varthetax_", i))


      l[[paste0('C', i, '_smoothed')]] <- get(paste0("C", i))

      assign(paste0("Z_", i), 1 / (1 + get(paste0(
        "varthetax_", i
      )) *
        apply(tmp_exp * muxt_hat, 1, sum, na.rm = T)))




      # l[[paste0("C", i)]] <- get(paste0("C", i))
      # l[[paste0("Fxt_", i)]] <- get(paste0("Fxt_", i))
      # l[[paste0("Z_", i)]] <- get(paste0("Z_", i))
      # l[[paste0("varthetax_", i)]] <- get(paste0("varthetax_", i))

      # l[[paste0("muhat", i)]] <- get(paste0("Z_", i)) * muxt_hat_predicted + (1 - get(paste0("Z_", i))) * get(paste0("C", i)) * muxt_hat_predicted
      # l[[paste0("muhat", i,"_MLE")]] <- get(paste0("C", i)) * muxt_hat_predicted


      # term1 <- 2*(muxt_hat_predicted^2)*(get(paste0("varthetax_", i)))+muxt_hat_predicted^2
      #
      # term2 <- get(paste0("varthetax_", i))*(apply((tmp_exp^2) *( muxt_hat)^2, 1, sum, na.rm =
      #                               T)/ apply(tmp_exp * muxt_hat, 1, sum, na.rm =
      #                                           T)^2)+(1/ apply(tmp_exp * muxt_hat, 1, sum, na.rm =
      #                                                             T))
      #
      # l[[paste0('msep_',i)]] <- term1+term2

      subgroups_forecast <- data_pp_2$list_of_extra_exposures[[i - 1]]

      # observed_rates[[paste0('muxt_actual_', i)]] <- subgroups_forecast$Dxt[, (length(years_fit) +
      #                                                                            1):(length(years_fit) + forecasting_horizon)] / subgroups_forecast$Ext[, (length(years_fit) +
      #                                                                                                                                                        1):(length(years_fit) + forecasting_horizon)]
      #
      observed_rates[[paste0('muxt_actual_', i)]] <- subgroups_forecast$Dxt[, (length(years_fit) + forecasting_horizon)] / subgroups_forecast$Ext[, (length(years_fit) + forecasting_horizon)]

    }

    # observed_rates[['full_pp_data']] <- data_pp_2
    # out[['actual_data']] <- observed_rates
    out[[model_option]] <- l

  }



  selected_out <- out[[model_2_plot]]



  dt_2_plot <- data.frame(
    population = rep(c(1, 2, 3), each = length(ages_fit)),
    theta_unsmoothed = c(
      selected_out$C1_unsmoothed,
      selected_out$C2_unsmoothed,
      selected_out$C3_unsmoothed
    ),
    theta_smoothed = c(
      selected_out$C1_smoothed,
      selected_out$C2_smoothed,
      selected_out$C3_smoothed
    ),
    variance_unsmoothed = c(
      selected_out$varthetax_1_unsmoothed,
      selected_out$varthetax_2_unsmoothed,
      selected_out$varthetax_3_unsmoothed
    ),
    variance_smoothed = c(
      selected_out$varthetax_1_smoothed,
      selected_out$varthetax_2_smoothed,
      selected_out$varthetax_3_smoothed
    )

  ) %>%
    mutate(
      group_variance = cumsum(c(TRUE, diff(
        variance_smoothed
      ) != 0)),
      group_theta = cumsum(c(TRUE, diff(theta_smoothed) != 0)),
      age_code = rep(as.character(ages_fit), 3)
    )

  dt_1 <- dt_2_plot %>%
    select(theta_unsmoothed, age_code, population)
  dt_2 <- dt_2_plot %>%
    select(theta_smoothed, group_theta, age_code, population)


  p_theta <- ggplot(dt_1,
                    aes(
                      x = age_code,
                      y = theta_unsmoothed,
                      color = as.character(population)
                    )) +
    geom_point(size = 2, alpha = .3) +
    geom_line(
      data = dt_2,
      aes(
        x = age_code,
        y = theta_smoothed,
        group = group_theta,
        color = as.character(population)
      ),
      size = 2,
      alpha = .7
    ) +
    ylim(0.2, 2.3) +
    theme_bw() +
    scale_x_discrete(breaks = c(20, 40, 60, 80, 100)) +
    theme(text = element_text(size = 28), legend.position = "top") +
    ylab("") +
    xlab("") +
    labs(color = "") +
    scale_color_manual(
      values = c(
        "1" = "#a71429",
        "2" = "#4169E1",
        "3" = "#006400"
      ),
      labels = c(
        expression(hat(theta)[x]^0),
        expression(hat(theta)[x]^1),
        expression(hat(theta)[x]^2)
      )
    )


  dt_1 <- dt_2_plot %>%
    select(variance_unsmoothed, age_code, population)
  dt_2 <- dt_2_plot %>%
    select(variance_smoothed, group_variance, age_code, population)

  p_var <- ggplot(dt_1,
                  aes(
                    x = age_code,
                    y = variance_unsmoothed,
                    color = as.character(population)
                  )) +
    geom_point(size = 2, alpha = .3) +
    geom_line(
      data = dt_2,
      aes(
        x = age_code,
        y = variance_smoothed,
        group = group_variance,
        color = as.character(population)
      ),
      size = 2,
      alpha = .7
    ) +
    ylim(-.0000001, 2.3) +
    theme_bw() +
    scale_x_discrete(breaks = c(20, 40, 60, 80, 100)) +
    theme(text = element_text(size = 28), legend.position = "top") +
    # ylim(0,.7)+
    ylab("") +
    xlab("") +
    labs(color = "") +
    scale_color_manual(
      values = c(
        "1" = "#a71429",
        "2" = "#4169E1",
        "3" = "#006400"
      ),
      labels = c(expression(hat(Var)(Theta[x]^0)), expression(hat(Var)(Theta[x]^1)), expression(hat(Var)(Theta[x]^2)))
    )

  return(list(variance_plot = p_var, theta_plot = p_theta))
}



# BIC ----

# Scratch in-sample BIC ---------------------------------------------------

# extract fitted in-sample rates from a StMoMo fit without using BIC()/logLik()
.get_stmomo_fitted_rates <- function(model_fit) {
  out <- NULL

  out <- tryCatch(
    stats::fitted(model_fit, type = "rates"),
    error = function(e)
      NULL
  )
  if (!is.null(out))
    return(as.matrix(out))

  out <- tryCatch(
    stats::fitted(model_fit),
    error = function(e)
      NULL
  )
  if (!is.null(out))
    return(as.matrix(out))

  # fallback used elsewhere in your code path for fitted historical rates
  out <- tryCatch({
    fc <- forecast(model_fit, h = 1)
    fc$fitted
  }, error = function(e)
    NULL)
  if (!is.null(out))
    return(as.matrix(out))

  stop("Could not extract fitted in-sample rates from the StMoMo fit.")
}

# count identifiable parameters without using BIC()/logLik()
.count_stmomo_parameters <- function(model_fit) {
  cf <- tryCatch(
    stats::coef(model_fit),
    error = function(e)
      NULL
  )

  if (!is.null(cf)) {
    k <- sum(is.finite(unlist(
      cf, recursive = TRUE, use.names = FALSE
    )))
    if (k > 0)
      return(as.integer(k))
  }

  for (nm in c("npar", "n_par", "df")) {
    if (!is.null(model_fit[[nm]]) &&
        length(model_fit[[nm]]) == 1L &&
        is.finite(model_fit[[nm]]) &&
        model_fit[[nm]] > 0) {
      return(as.integer(model_fit[[nm]]))
    }
  }

  # conservative fallback
  pieces <- c("ax", "bx", "kt", "b0x", "gc")
  k <- sum(vapply(pieces, function(nm) {
    x <- model_fit[[nm]]
    if (is.null(x))
      return(0L)
    sum(is.finite(as.numeric(x)))
  }, integer(1)))

  if (k <= 0) {
    stop("Could not determine the number of fitted parameters.")
  }

  as.integer(k)
}

# Poisson log-likelihood allowing non-integer Dxt via lgamma(Dxt + 1)
.poisson_loglik_lexis <- function(Dxt, Ext, muxt_hat, eps = 1e-12) {
  Dxt <- as.matrix(Dxt)
  Ext <- as.matrix(Ext)
  muxt_hat <- as.matrix(muxt_hat)

  valid <- is.finite(Dxt) &
    is.finite(Ext) &
    is.finite(muxt_hat) &
    (Ext > 0)

  if (!any(valid)) {
    stop("No valid cells available for scratch likelihood.")
  }

  d <- Dxt[valid]
  e <- Ext[valid]
  m <- pmax(muxt_hat[valid], eps)

  sum(d * log(e * m) - e * m - lgamma(d + 1))
}

compute_insample_bic_scratch <- function(model_fit, data_stmomo, model_label) {
  Dxt <- as.matrix(data_stmomo$Dxt)
  Ext <- as.matrix(data_stmomo$Ext)

  muxt_hat <- .get_stmomo_fitted_rates(model_fit)

  # enforce same age-year layout as the fitted data
  muxt_hat <- muxt_hat[rownames(Dxt), colnames(Dxt), drop = FALSE]

  valid <- is.finite(Dxt) &
    is.finite(Ext) &
    is.finite(muxt_hat) &
    (Ext > 0)

  n_obs <- sum(valid)
  k <- .count_stmomo_parameters(model_fit)
  loglik <- .poisson_loglik_lexis(Dxt = Dxt,
                                  Ext = Ext,
                                  muxt_hat = muxt_hat)
  bic <- -2 * loglik + log(n_obs) * k

  data.table::data.table(
    model = model_label,
    bic = bic,
    loglik = loglik,
    n_obs = n_obs,
    k = k
  )
}

#' @export
run_insample_bic_cell <- function(seed_ix,
                                  years_ix,
                                  ages_fit = 15:85,
                                  years_fit_basic = 118:140,
                                  N_groups = 3,
                                  exposure_sup = 94500,
                                  exposure_sub = c(5000, 500),

                                  global_mortality_model_candidates = c("lc", "apc", "rh"),
                                  hmd_user = "",
                                  hmd_pass = "") {
  dt_0 <- data_generator_hmd_lt(
    seed_input = seed_ix,
    hmd_username = hmd_user,
    hmd_password = hmd_pass,
    exposure_sup = exposure_sup,
    exposure_sub = exposure_sub
  )

  years_fit <- c(
    years_fit_basic,
    seq.int(max(years_fit_basic) + 1L, max(years_fit_basic) + years_ix)
  )

  data_pp <- data_preprocessing(
    data = dt_0,
    N_groups = N_groups,
    ages_fit = ages_fit,
    years_fit = years_fit
  )

  mortality_models_fit <- fit_mortality_models(
    data_pp = data_pp,
    years_fit = years_fit - min(years_fit),
    ages_fit = ages_fit,
    global_mortality_model_candidates = global_mortality_model_candidates
  )

  out <- .summarise_insample_bic(
    mortality_models_fit = mortality_models_fit,
    data_pp = data_pp
  )

  out[, `:=`(
    seed_ix = as.integer(seed_ix),
    years_ix = as.integer(years_ix),
    n_years_fit = length(years_fit)
  )]

  out[]
}



# -------------------------------------------------------------------
# Figure 5 helpers
# -------------------------------------------------------------------
#' @export
.parse_age_breakets_start <- function(x) {
  as.integer(sub("-.*$", "", x))
}
#' @export
.age_breakets_levels <- function(x) {
  unique(x[order(.parse_age_breakets_start(x))])
}
#' @export
.figure5_model_levels <- function(group_id = NULL) {
  if (is.null(group_id)) {
    c("superpopulation", "credibility", "full_mle", "separate")
  } else {
    c("total", "credibility", "full_mle", "separate")
  }
}
#' @export
summarise_age_bracket_results <- function(results_dt,
                                          predictor = "lc",
                                          group_id = NULL) {
  results_dt <- data.table::as.data.table(results_dt)

  error_col <- if (is.null(group_id))
    "error"
  else
    paste0("error_", group_id)
  dev_col <- if (is.null(group_id))
    "oos_deviance"
  else
    paste0("oos_deviance_", group_id)
  excluded_model <- if (is.null(group_id))
    "total"
  else
    "superpopulation"
  model_levels <- .figure5_model_levels(group_id)

  dt <- results_dt[predictor == ..predictor &
                     model != excluded_model, .(error = mean(get(error_col), na.rm = TRUE),
                                                oos_deviance = mean(get(dev_col), na.rm = TRUE)), by = .(seed, model, predictor, age_breakets)][, .(
                                                  error_resp = mean(error, na.rm = TRUE),
                                                  error_sd = stats::sd(error, na.rm = TRUE),
                                                  oos_deviance_resp = mean(oos_deviance, na.rm = TRUE),
                                                  oos_deviance_sd = stats::sd(oos_deviance, na.rm = TRUE)
                                                ), by = .(model, predictor, age_breakets)]

  dt[, model := factor(model, levels = model_levels)]
  dt[, age_breakets := factor(age_breakets, levels = .age_breakets_levels(age_breakets))]
  data.table::setorder(dt, age_breakets, model)
  dt[]
}

#' @export
plot_age_bracket_metric <- function(summary_dt,
                                    metric = c("mare", "deviance"),
                                    legend_position = c(0.98, 0.98),
                                    legend_justification = c(1, 1),
                                    base_size = 28) {
  metric <- match.arg(metric)
  summary_dt <- data.table::as.data.table(summary_dt)

  y_col <- if (metric == "mare")
    "error_resp"
  else
    "oos_deviance_resp"
  ysd_col <- if (metric == "mare")
    "error_sd"
  else
    "oos_deviance_sd"
  y_limits <- if (metric == "mare")
    c(0, 1)
  else
    c(0, 30)

  n_models <- length(stats::na.omit(unique(summary_dt$model)))
  model_labels <- LETTERS[seq_len(n_models)]

  ggplot2::ggplot(summary_dt,
                  ggplot2::aes(
                    x = age_breakets,
                    y = .data[[y_col]],
                    color = model,
                    group = model
                  )) +
    ggplot2::geom_line(linewidth = 1.5, na.rm = TRUE) +
    ggplot2::geom_point(na.rm = TRUE) +
    ggplot2::geom_errorbar(
      ggplot2::aes(ymin = .data[[y_col]] - .data[[ysd_col]], ymax = .data[[y_col]] + .data[[ysd_col]]),
      width = 0.2,
      alpha = 0.3,
      na.rm = TRUE
    ) +
    ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::coord_cartesian(ylim = y_limits) +
    ggplot2::labs(x = "", y = "", color = "") +
    ggplot2::scale_color_discrete(labels = model_labels) +
    ggplot2::theme(
      legend.position = legend_position,
      legend.justification = legend_justification,
      legend.key.height = grid::unit(1.5, "lines"),
      legend.background = ggplot2::element_blank(),
      legend.key = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(
        angle = 90,
        hjust = 1,
        vjust = 0.5
      )
    )
}

#' @export
write_age_bracket_plot <- function(results_dt,
                                   predictor = "lc",
                                   metric = c("mare", "deviance"),
                                   group_id = NULL,
                                   out_dir = "figures",
                                   filename = NULL,
                                   width = 8,
                                   height = 5,
                                   legend_position = c(0.98, 0.98),
                                   legend_justification = c(1, 1)) {
  metric <- match.arg(metric)

  summary_dt <- summarise_age_bracket_results(results_dt = results_dt,
                                              predictor = predictor,
                                              group_id = group_id)

  p <- plot_age_bracket_metric(
    summary_dt = summary_dt,
    metric = metric,
    legend_position = legend_position,
    legend_justification = legend_justification
  )

  if (is.null(filename)) {
    filename <- if (is.null(group_id)) {
      sprintf("%s_%s_age_classes.pdf",
              predictor,
              if (metric == "mare")
                "mare"
              else
                "deviance")
    } else {
      sprintf(
        "%s_%s_age_classes_group_%s.pdf",
        predictor,
        if (metric == "mare")
          "mare"
        else
          "deviance",
        group_id
      )
    }
  }

  fs::dir_create(out_dir)
  out <- file.path(out_dir, filename)

  ggplot2::ggsave(
    filename = out,
    plot = p,
    width = width,
    height = height
  )

  out
}



# -------------------------------------------------------------------
# Figure 6 helpers
# -------------------------------------------------------------------

#' @export
scoring_plot <- function(data,
                         metric,
                         superpopulation=FALSE,
                         subpopulation=FALSE,
                         predictor_choice) {


  if(superpopulation){

    dt_out <- data[model != "total" &
                     predictor == predictor_choice, ]

    dt_out <- dt_out[, .(error = mean(error),
                         oos_deviance = mean(oos_deviance)), by = .(seed, model, predictor, age_breakets)]


  }else{

    dt_out <- data[model != "superpopulation" &
                     predictor == predictor_choice, ]

  }

  if (subpopulation == 1) {
    dt_out <- dt_out[, .(error = mean(error_1),
                         oos_deviance = mean(oos_deviance_1)), by = .(seed, model, predictor, age_breakets)]
  }


  if (subpopulation == 2) {
    dt_out <- dt_out[, .(error = mean(error_2),
                         oos_deviance = mean(oos_deviance_2)), by = .(seed, model, predictor, age_breakets)]
  }


  if (subpopulation == 3) {
    dt_out <- dt_out[, .(error = mean(error_3),
                         oos_deviance = mean(oos_deviance_3)), by = .(seed, model, predictor, age_breakets)]
  }




  dt_out <- dt_out[, .(
    error_resp = (mean(error)),
    error_sd = (sd(error)),
    oos_deviance_resp = (mean(oos_deviance)),
    oos_deviance_sd = (sd(oos_deviance))
  ), by = .(model, predictor, age_breakets)]



  if (metric == "mse") {
    p <- ggplot(dt_out,
                aes(
                  x = age_breakets,
                  y = error_resp,
                  color = model,
                  group = model
                )) +
      geom_line(size = 2) +    # Add lines
      geom_point() +   # Add points
      geom_errorbar(
        aes(ymin = error_resp - error_sd , ymax = error_resp + error_sd),
        width = 0.2,
        size = 2,
        alpha = .3
      ) +
      theme_minimal(base_size = 28)
  } else{
    p <- ggplot(dt_out,
                aes(
                  x = age_breakets,
                  y = oos_deviance_resp,
                  color = model,
                  group = model
                )) +
      geom_line(size = 2) +    # Add lines
      geom_point() +   # Add points
      geom_errorbar(
        aes(
          ymin = oos_deviance_resp - oos_deviance_sd ,
          ymax = oos_deviance_resp + oos_deviance_sd
        ),
        width = 0.2,
        size = 2,
        alpha = .3
      ) +
      theme_minimal(base_size = 28)
  }

  return(p)

}

# -------------------------------------------------------------------
# MSEP study helpers
# -------------------------------------------------------------------

#' @export
.msep_age_index <- function(x, reference_age) {
  rn <- rownames(x)

  if (!is.null(rn)) {
    idx <- match(as.character(reference_age), rn)
    if (!is.na(idx)) {
      return(idx)
    }
  }

  idx <- match(reference_age, seq_len(nrow(x)))
  if (!is.na(idx)) {
    return(idx)
  }

  stop("Could not locate reference_age in the supplied matrix.")
}

#' @export
.msep_forecast_stmomo <- function(model_fit,
                                  model_option = "lc",
                                  forecasting_horizon = 5) {
  cv.arima.kt <- forecast::auto.arima(as.numeric(model_fit$kt), ic = "bic")

  if (model_option != "lc" && !is.null(model_fit$gc)) {
    cv.arima.gc <- forecast::auto.arima(as.numeric(model_fit$gc), ic = "bic")
    gc.order <- unname(forecast::arimaorder(cv.arima.gc))
  } else {
    gc.order <- c(1, 1, 0)
  }

  out <- tryCatch(
    forecast::forecast(
      model_fit,
      kt.method = "iarima",
      gc.order = gc.order,
      kt.order = unname(forecast::arimaorder(cv.arima.kt)),
      h = forecasting_horizon
    ),
    error = function(e) {
      forecast::forecast(model_fit, h = forecasting_horizon)
    }
  )

  out
}

#' @export
.msep_is_implausible_path <- function(muhat,
                                      global_muhat,
                                      max_ratio = 50,
                                      max_log_jump = 3.5,
                                      eps = 1e-12) {
  muhat <- as.numeric(muhat)
  global_muhat <- as.numeric(global_muhat)

  if (length(muhat) != length(global_muhat)) {
    return(TRUE)
  }

  if (any(!is.finite(muhat)) || any(muhat <= 0)) {
    return(TRUE)
  }

  log_ratio <- abs(log(pmax(muhat, eps) / pmax(global_muhat, eps)))
  if (any(!is.finite(log_ratio)) || any(log_ratio > log(max_ratio))) {
    return(TRUE)
  }

  log_muhat <- log(pmax(muhat, eps))
  if (length(log_muhat) > 1L &&
      any(abs(diff(log_muhat)) > max_log_jump, na.rm = TRUE)) {
    return(TRUE)
  }

  FALSE
}

#' @export
.msep_fit_separate_group <- function(subpop,
                                     ages_fit,
                                     years_fit,
                                     model_option = "lc",
                                     global_fit = NULL,
                                     global_future_path = NULL,
                                     forecasting_horizon = 5,
                                     max_ratio_implausible = 50,
                                     max_log_jump_implausible = 3.5) {
  mortality_model_list <- list(
    lc = StMoMo::lc(link = "log"),
    apc = StMoMo::apc(link = "log"),
    rh = StMoMo::rh(link = "log")
  )

  popspecdata <- structure(
    list(
      Dxt = subpop$Dxt[as.character(ages_fit), as.character(years_fit), drop = FALSE],
      Ext = subpop$Ext[as.character(ages_fit), as.character(years_fit), drop = FALSE],
      ages = ages_fit,
      years = years_fit,
      type = "central",
      series = "male",
      label = "subpopulation"
    ),
    class = "StMoMoData"
  )

  fit_try <- tryCatch(
    suppressWarnings(
      StMoMo::fit(
        mortality_model_list[[model_option]],
        data = popspecdata,
        years.fit = years_fit,
        ages.fit = ages_fit,
        verbose = FALSE
      )
    ),
    error = function(e) NULL
  )

  if (is.null(fit_try)) {
    return(list(
      fit = global_fit,
      used_global_fallback = TRUE
    ))
  }

  fc_try <- tryCatch(
    .msep_forecast_stmomo(
      model_fit = fit_try,
      model_option = model_option,
      forecasting_horizon = forecasting_horizon
    ),
    error = function(e) NULL
  )

  if (is.null(fc_try)) {
    return(list(
      fit = global_fit,
      used_global_fallback = TRUE
    ))
  }

  future_idx <- tail(seq_len(ncol(fc_try$rates)), forecasting_horizon)
  age_idx <- .msep_age_index(fc_try$rates, ages_fit[1L])
  sep_future_example <- as.numeric(fc_try$rates[age_idx, future_idx])

  if (!is.null(global_future_path) &&
      .msep_is_implausible_path(
        muhat = sep_future_example,
        global_muhat = global_future_path,
        max_ratio = max_ratio_implausible,
        max_log_jump = max_log_jump_implausible
      )) {
    return(list(
      fit = global_fit,
      used_global_fallback = TRUE
    ))
  }

  list(
    fit = fit_try,
    used_global_fallback = FALSE
  )
}

#' @export
.msep_plot_long <- function(dt_all_rates) {
  out <- dt_all_rates |>
    dplyr::mutate(calendar_time = 0:(nrow(dt_all_rates) - 1L)) |>
    dplyr::select(
      -true_rates
    ) |>
    tidyr::pivot_longer(
      cols = c(
        credibility_rates,
        separate_rates,
        conf1_ifreal,
        conf2_ifreal,
        conf1_cred,
        conf2_cred,
        conf1_sep,
        conf2_sep
      ),
      names_to = "quantity",
      values_to = "values"
    )

  out[["ltype"]] <- "solid"
  out$ltype[grepl("conf", out$quantity)] <- "dotted"

  out
}

#' @export
plot_msep_group <- function(group_result, ylim = NULL) {
  dt_all_rates <- group_result$plot_data
  dt_long <- .msep_plot_long(dt_all_rates)

  temporary_data <- dt_all_rates |>
    dplyr::select(true_rates) |>
    dplyr::mutate(calendar_time = 0:(nrow(dt_all_rates) - 1L))

  p <- ggplot2::ggplot(dt_long, ggplot2::aes(x = calendar_time, y = values, color = quantity)) +
    ggplot2::geom_line(
      ggplot2::aes(linetype = ltype),
      linewidth = 1
    ) +
    ggplot2::geom_point(
      data = temporary_data |>
        dplyr::rename(values = true_rates) |>
        dplyr::mutate(quantity = "true_rates"),
      ggplot2::aes(x = calendar_time, y = values),
      color = "gray"
    ) +
    ggplot2::scale_linetype_manual(values = c("solid" = "solid", "dotted" = "dotted")) +
    ggplot2::scale_color_manual(
      values = c(
        "conf1_ifreal" = "black",
        "conf2_ifreal" = "black",
        "credibility_rates" = "#4169E1",
        "separate_rates" = "#a71429",
        "conf1_cred" = "#4169E1",
        "conf2_cred" = "#4169E1",
        "conf1_sep" = "#a71429",
        "conf2_sep" = "#a71429"
      )
    ) +
    ggplot2::theme_minimal() +
    ggplot2::labs(title = "", x = "", y = "", color = "") +
    ggplot2::theme(legend.position = "none")

  if (!is.null(ylim)) {
    stopifnot(is.numeric(ylim), length(ylim) == 2, all(is.finite(ylim)))
    p <- p + ggplot2::coord_cartesian(ylim = ylim)
  }

  p
}

#' @export
write_msep_plot <- function(group_result, out, ylim = NULL) {
  fs::dir_create(dirname(out))
  p <- plot_msep_group(group_result, ylim = ylim)

  ggplot2::ggsave(
    filename = out,
    plot = p,
    width = 8,
    height = 5
  )

  out
}

#' @export
run_msep_study <- function(data,
                           data_pp,
                           mortality_models_fit,
                           credibility_model,
                           ages_fit,
                           years_fit,
                           N_groups = 3,
                           predictor = "lc",
                           reference_age = 65,
                           forecasting_horizon = 5,
                           bias = min(years_fit),
                           n_sim_sigma = 1000,
                           n_sim_obs = 1000,
                           separate_bootstrap_n = 500,
                           sigma_seed = 1,
                           poisson_seed = 1,
                           max_ratio_implausible = 50,
                           max_log_jump_implausible = 3.5,
                           if_real_alpha = 0.05,
                           competitor_alpha = 0.05) {

  model_fit <- mortality_models_fit[[predictor]]
  cred <- credibility_model[[predictor]]

  future_years <- c(
    years_fit,
    (max(years_fit) + 1):(max(years_fit) + forecasting_horizon)
  )

  data_pp_future <- data_preprocessing(
    data = data,
    N_groups = N_groups,
    bias = bias,
    ages_fit = ages_fit,
    years_fit = future_years
  )

  global_fc <- .msep_forecast_stmomo(
    model_fit = model_fit,
    model_option = predictor,
    forecasting_horizon = forecasting_horizon
  )

  age_idx_global <- .msep_age_index(global_fc$fitted, reference_age)

  fitted_idx <- seq_len(ncol(global_fc$fitted))
  future_idx <- tail(seq_len(ncol(global_fc$rates)), forecasting_horizon)

  global_fitted_age <- as.numeric(global_fc$fitted[age_idx_global, fitted_idx])
  global_future_age <- as.numeric(global_fc$rates[age_idx_global, future_idx])

  set.seed(sigma_seed)
  sigma_sim <- simulate(
    object = model_fit,
    nsim = n_sim_sigma,
    h = forecasting_horizon
  )

  sigma_age <- sigma_sim$rates[age_idx_global, , ]
  sigma2_forecast <- apply(sigma_age, 1, stats::var)

  results <- vector("list", N_groups)
  names(results) <- as.character(seq_len(N_groups))

  for (i in seq_len(N_groups)) {

    subpop_hist <- data_pp$subpopulations[[i]]
    subpop_future <- data_pp_future$subpopulations[[i]]

    age_idx_hist <- .msep_age_index(subpop_hist$Ext, reference_age)
    age_idx_future <- .msep_age_index(subpop_future$Ext, reference_age)

    hist_data_idx <- seq_len(length(years_fit))
    future_data_idx <- (length(years_fit) + 1L):(length(years_fit) + forecasting_horizon)

    tmp_exp_hist <- as.numeric(subpop_hist$Ext[age_idx_hist, hist_data_idx])
    tmp_exp_future <- as.numeric(subpop_future$Ext[age_idx_future, future_data_idx])

    true_rates <- as.numeric(
      subpop_future$Dxt[age_idx_future, seq_len(length(future_years))] /
        subpop_future$Ext[age_idx_future, seq_len(length(future_years))]
    )

    vartheta_x <- cred[[paste0("varthetax_", i)]][match(reference_age, ages_fit)]
    one_minus_Z <- cred[[paste0("Z_", i)]][match(reference_age, ages_fit)]
    theta_hat <- cred[[paste0("C", i)]][match(reference_age, ages_fit)]

    Z_manuscript <- 1 - one_minus_Z

    denom_hist <- sum(tmp_exp_hist * global_fitted_age)

    ## New Corollary 3 plug-in variance:
    ## Var(theta_hat | F_t') = Var(Theta) + 1 / sum_v E_{x,v}^i mu_{x,v}
    var_theta_hat <- vartheta_x + (1 / denom_hist)

    term1 <- sigma2_forecast * (vartheta_x + 1) +
      (global_future_age^2) * vartheta_x

    msep_cred <- term1 +
      (Z_manuscript^2) * (global_future_age^2) * var_theta_hat

    predicted_rates <- one_minus_Z * global_future_age +
      Z_manuscript * theta_hat * global_future_age

    insample_rates <- one_minus_Z * global_fitted_age +
      Z_manuscript * theta_hat * global_fitted_age

    credibility_full_path <- c(insample_rates, predicted_rates)
    full_exposure_path <- c(tmp_exp_hist, tmp_exp_future)

    # Black bands: one-SD delta-method band on the log-rate scale
    safe_rates <- pmax(credibility_full_path, .Machine$double.eps)
    safe_exposure <- pmax(full_exposure_path, .Machine$double.eps)

    poisson_sd_if_real <- sqrt(safe_rates / safe_exposure)
    poisson_se_log_if_real <- 1 / sqrt(safe_exposure * safe_rates)

    conf1_ifreal <- credibility_full_path * exp(poisson_se_log_if_real)
    conf2_ifreal <- credibility_full_path * exp(-poisson_se_log_if_real)

    sep_fit_info <- .msep_fit_separate_group(
      subpop = subpop_future,
      ages_fit = ages_fit,
      years_fit = years_fit,
      model_option = predictor,
      global_fit = model_fit,
      global_future_path = global_future_age,
      forecasting_horizon = forecasting_horizon,
      max_ratio_implausible = max_ratio_implausible,
      max_log_jump_implausible = max_log_jump_implausible
    )

    sep_fit <- sep_fit_info$fit
    sep_fc <- .msep_forecast_stmomo(
      model_fit = sep_fit,
      model_option = predictor,
      forecasting_horizon = forecasting_horizon
    )

    age_idx_sep <- .msep_age_index(sep_fc$fitted, reference_age)
    sep_future_age <- as.numeric(sep_fc$rates[age_idx_sep, future_idx])

    sep_boot <- suppressWarnings(
      bootstrap(sep_fit, nBoot = separate_bootstrap_n, type = "residual")
    )

    sep_boot_sim <- simulate(sep_boot, h = forecasting_horizon)
    sep_future_draws <- sep_boot_sim$rates[age_idx_sep, , ]

    if (is.null(dim(sep_future_draws))) {
      sep_future_draws <- matrix(
        sep_future_draws,
        nrow = forecasting_horizon,
        ncol = 1
      )
    }

    sep_var_age <- apply(sep_future_draws, 1, stats::var, na.rm = TRUE)

    conf1_sep_future <- sep_future_age + sqrt(sep_var_age)
    conf2_sep_future <- sep_future_age - sqrt(sep_var_age)

    conf1_cred_future <- predicted_rates + sqrt(msep_cred)
    conf2_cred_future <- predicted_rates - sqrt(msep_cred)

    n_hist <- length(true_rates) - forecasting_horizon
    pad_hist <- rep(NA_real_, n_hist)

    dt_all_rates <- data.frame(
      true_rates = true_rates,
      conf1_ifreal = conf1_ifreal,
      conf2_ifreal = conf2_ifreal,
      credibility_rates = c(pad_hist, predicted_rates),
      separate_rates = c(pad_hist, sep_future_age),
      conf1_cred = c(pad_hist, conf1_cred_future),
      conf2_cred = c(pad_hist, conf2_cred_future),
      conf1_sep = c(pad_hist, conf1_sep_future),
      conf2_sep = c(pad_hist, conf2_sep_future)
    )

    results[[as.character(i)]] <- list(
      reference_group = i,
      reference_age = reference_age,
      forecasting_horizon = forecasting_horizon,
      used_global_fallback_for_separate = sep_fit_info$used_global_fallback,
      theta_hat = theta_hat,
      vartheta_x = vartheta_x,
      one_minus_Z = one_minus_Z,
      Z_manuscript = Z_manuscript,
      sigma2_forecast = sigma2_forecast,
      var_theta_hat = var_theta_hat,
      msep_credibility = msep_cred,
      msep_separate = sep_var_age,
      credibility_future_rates = predicted_rates,
      separate_future_rates = sep_future_age,
      credibility_full_path = credibility_full_path,
      poisson_sd_full_path = poisson_sd_if_real,
      poisson_se_log_full_path = poisson_se_log_if_real,
      true_rates = true_rates,
      plot_data = dt_all_rates
    )
  }

  list(
    predictor = predictor,
    reference_age = reference_age,
    forecasting_horizon = forecasting_horizon,
    groups = results
  )
}




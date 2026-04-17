# _targets.R
library(targets)
library(tarchetypes)
library(credMortality)   # after you package-ify


tar_option_set(packages = c(
  "dplyr",
  "tidyr",
  "ggplot2",
  "data.table",
  "StMoMo",
  "forecast",
  "fs"
))


list(
  tar_target(ages_fit, 15:85),
  tar_target(years_fit, 118:140),
  tar_target(N_groups, 3),

  tar_target(
    hmd_user,
    {
      x <- Sys.getenv("HMD_USER", unset = "")
      if (x == "") {
        stop("Environment variable HMD_USER is not set.")
      }
      x
    }
  ),

  tar_target(
    hmd_pass,
    {
      x <- Sys.getenv("HMD_PASS", unset = "")
      if (x == "") {
        stop("Environment variable HMD_PASS is not set.")
      }
      x
    }
  ),

  tar_target(
    dt_0,
    data_generator_hmd_lt(
      seed_input = 1,
      hmd_username = hmd_user,
      hmd_password = hmd_pass,
      exposure_sup = 94500,
      exposure_sub = c(5000, 500)
    )
  ),

  tar_target(fig6_pdf, {
    fs::dir_create("figures")
    out <- "figures/mortality_rates_simulation.pdf"

    p <- dt_0 |>
      dplyr::filter(Age == 55) |>
      tidyr::pivot_longer(
        cols = c("p", "p1", "p2"),
        names_to = "group",
        values_to = "probabilities"
      ) |>
      dplyr::mutate(group = dplyr::recode(group, p = "Group 0", p1 = "Group 1", p2 =
                                            "Group 2")) |>
      ggplot2::ggplot(aes(x = Period, y = probabilities)) +
      ggplot2::geom_line(aes(colour = group), linewidth = 1.5) +
      ggplot2::theme_bw() +
      ggplot2::theme(text = element_text(size = 28), legend.position = "none") +
      ggplot2::scale_color_manual(values = c(
        "Group 0" = "#a71429",
        "Group 1" = "#4169E1",
        "Group 2" = "#2E8B57"
      )) +
      ggplot2::ylab(expression(q[xt]^i)) +
      ggplot2::xlab("") +
      ggplot2::labs(color = "")

    ggplot2::ggsave(
      filename = out,
      plot = p,
      width = 8,
      height = 5
    )
    out
  }, format = "file"),

  tar_target(
    data_pp,
    data_preprocessing(
      data = dt_0,
      ages_fit = ages_fit,
      years_fit = years_fit
    )
  ),

  tar_target(
    mortality_models_fit,
    fit_mortality_models(data_pp, years_fit - min(years_fit), ages_fit)
  ),

  tar_target(
    total_model,
    fit_and_predict_total_model(
      data = dt_0,
      N_groups = N_groups,
      mortality_models_fit = mortality_models_fit,
      years_fit = years_fit - min(years_fit),
      ages_fit = ages_fit,
      forecasting_horizon = 1,
      bias = min(years_fit)
    )
  ),

  tar_target(
    credibility_model,
    fit_and_predict_credibility_models(
      data = dt_0,
      N_groups = N_groups,
      mortality_models_fit = mortality_models_fit,
      years_fit = years_fit - min(years_fit),
      ages_fit = ages_fit,
      forecasting_horizon = 1,
      bias = min(years_fit)
    )
  ),

  tar_target(fig2a_pdf, {
    out <- "figures/lc_allthetas_zs.pdf"

    dt <- print_age_effects(
      seed_input = 1,
      hmd_username = hmd_user,
      hmd_password = hmd_pass,
      exposure_sup = 94500,
      exposure_sub = c(5000, 500)
    )
    setDT(dt)

    dt[['ages']] <- 1:nrow(dt)

    dt <- dt[ages_fit, ]

    p <- all_thetas_plotter(credibility_model, N_groups, ages_fit, predictor = "lc") +
      scale_color_manual(
        values = c(
          "C2" = "#4169E1",
          "C3" = "#006400",
          "C1" = "#a71429"
        ),
        labels = c(
          expression(hat(theta)[x]^1),
          expression(hat(theta)[x]^2),
          expression(hat(theta)[x]^3)
        )
      ) + geom_line(
        data = dt,
        aes(x = ages, y = age_eff_2),
        linewidth = 1.5,
        alpha = .4,
        color = "#006400"
      ) +
      geom_hline(
        yintercept = 1.,
        linewidth = 1.5,
        alpha = .4,
        color = "#a71429"
      ) +
      geom_line(
        data = dt,
        aes(x = ages, y = age_eff_1),
        linewidth = 1.5,
        alpha = .4,
        color = "#4169E1"
      )

    ggplot2::ggsave(
      filename = out,
      plot = p,
      width = 8,
      height = 5
    )
    out
  }, format = "file"),

  tar_target(fig2b_pdf, {
    out <- "figures/lc_variances.pdf"

    p <- variance_plotter(credibility_model,
                          ages_fit = ages_fit,
                          predictor = "lc") + scale_fill_discrete(guide =
                                                                    "none")

    ggplot2::ggsave(
      filename = out,
      plot = p,
      width = 8,
      height = 5
    )
    out
  }, format = "file"),

  tar_target(fig2c_pdf, {
    out <- "figures/lc_simulation_zs.pdf"
    p <- weights_plotter(
      credibility_model,
      N_groups,
      ages_fit,
      predictor = "lc",
      ages_breaks = c(25, 50, 75, 85)
    )

    ggplot2::ggsave(
      filename = out,
      plot = p,
      width = 8,
      height = 5
    )
    out
  }, format = "file"),

  tar_target(fig2d_pdf, {
    out <- "figures/lc_log_mortality_rates.pdf"

    p <- plot_credibility_mle_global(credibility_model, N_groups, ages_fit, predictor = "lc")

    ggplot2::ggsave(
      filename = out,
      plot = p,
      width = 8,
      height = 5
    )
    out
  }, format = "file"),
  #
  tar_target(
    separate_model,
    fit_and_predict_separate_models(
      data = dt_0,
      N_groups = N_groups,
      mortality_models_fit = mortality_models_fit,
      years_fit = years_fit - min(years_fit),
      ages_fit = ages_fit,
      forecasting_horizon = 1,
      bias = min(years_fit)
    )
  ),
  tar_target(
    ss_effect,
    sample_size_effect_on_predictions(
      data = dt_0,
      years_fit_basic =
        min(years_fit):(min(years_fit) + 4),
      number_of_rolls =
        20,
      ages_fit,
      bias = min(years_fit),
      global_mortality_model =
        "lc"
    )
  ),

  tar_target(additional_analysis_3_15_pdf, {
    out <- "figures/an315.pdf"

    p <- sample_size_effect_plotter(
      ss_effect,
      color_group =
        "#006400",
      ages_to_plot =
        15,
      label_to_plot =
        "3"
    )

    ggplot2::ggsave(
      filename = out,
      plot = p,
      width = 8,
      height = 5
    )
    out
  }, format = "file"),

  tar_target(additional_analysis_3_55_pdf, {
    out <- "figures/an355.pdf"

    p <- sample_size_effect_plotter(
      ss_effect,
      color_group = "#006400",
      ages_to_plot = 55,
      label_to_plot = "3"
    )

    ggplot2::ggsave(
      filename = out,
      plot = p,
      width = 8,
      height = 5
    )

    out
  }, format = "file"),

  tar_target(additional_analysis_3_80_pdf, {
    out <- "figures/an380.pdf"

    p <- sample_size_effect_plotter(
      ss_effect,
      color_group = "#006400",
      ages_to_plot = 80,
      label_to_plot = "3"
    )

    ggplot2::ggsave(
      filename = out,
      plot = p,
      width = 8,
      height = 5
    )

    out
  }, format = "file"),

  tar_target(additional_analysis_1_15_pdf, {
    out <- "figures/an115.pdf"

    p <- sample_size_effect_plotter(
      ss_effect,
      color_group = "#a71429",
      ages_to_plot = 15,
      label_to_plot = "1"
    )

    ggplot2::ggsave(
      filename = out,
      plot = p,
      width = 8,
      height = 5
    )

    out
  }, format = "file"),

  tar_target(additional_analysis_1_55_pdf, {
    out <- "figures/an155.pdf"

    p <- sample_size_effect_plotter(
      ss_effect,
      color_group = "#a71429",
      ages_to_plot = 55,
      label_to_plot = "1"
    )

    ggplot2::ggsave(
      filename = out,
      plot = p,
      width = 8,
      height = 5
    )

    out
  }, format = "file"),

  tar_target(additional_analysis_1_80_pdf, {
    out <- "figures/an180.pdf"

    p <- sample_size_effect_plotter(
      ss_effect,
      color_group = "#a71429",
      ages_to_plot = 80,
      label_to_plot = "1"
    )

    ggplot2::ggsave(
      filename = out,
      plot = p,
      width = 8,
      height = 5
    )

    out
  }, format = "file"),

  tar_target(additional_analysis_2_15_pdf, {
    out <- "figures/an215.pdf"

    p <- sample_size_effect_plotter(
      ss_effect,
      color_group = "#4169E1",
      ages_to_plot = 15,
      label_to_plot = "2"
    )

    ggplot2::ggsave(
      filename = out,
      plot = p,
      width = 8,
      height = 5
    )

    out
  }, format = "file"),

  tar_target(additional_analysis_2_55_pdf, {
    out <- "figures/an255.pdf"

    p <- sample_size_effect_plotter(
      ss_effect,
      color_group = "#4169E1",
      ages_to_plot = 55,
      label_to_plot = "2"
    )

    ggplot2::ggsave(
      filename = out,
      plot = p,
      width = 8,
      height = 5
    )

    out
  }, format = "file"),

  tar_target(additional_analysis_2_80_pdf, {
    out <- "figures/an280.pdf"

    p <- sample_size_effect_plotter(
      ss_effect,
      color_group = "#4169E1",
      ages_to_plot = 80,
      label_to_plot = "2"
    )

    ggplot2::ggsave(
      filename = out,
      plot = p,
      width = 8,
      height = 5
    )

    out
  }, format = "file"),



  ######### Section 4.4. simulation study
  tar_target(sim_years_fit_basic, 118:140),
  tar_target(sim_forecasting_horizon, 1),
  tar_target(sim_ageb_width, 6),
  tar_target(sim_seed_values, {
    set.seed(42)
    seeds <- sample.int(100000)
    seeds[1:30]}),
  tar_target(sim_years_ix_values, 1:6),
  tar_target(
    sim_design,
    tidyr::crossing(seed_ix = sim_seed_values, years_ix = sim_years_ix_values)
  ),
  tar_target(
    sim_age_bracket_results_branch,
    run_simulation_study_cell(
      seed_ix = sim_design$seed_ix,
      years_ix = sim_design$years_ix,
      ages_fit = ages_fit,
      years_fit_basic = sim_years_fit_basic,
      N_groups = N_groups,
      forecasting_horizon = sim_forecasting_horizon,
      ageb_width = sim_ageb_width,
      exposure_sup = 94500,
      exposure_sub = c(5000, 500),
      hmd_user = hmd_user,
      hmd_pass = hmd_pass,
      global_mortality_model = c("lc","apc","rh")
    ),
    pattern = map(sim_design),
    iteration = "list"
  ),

  tar_target(
    sim_age_bracket_results,
    data.table::rbindlist(sim_age_bracket_results_branch)
  ),
  tar_target(sim_age_bracket_results_file, {
    fs::dir_create("output")
    out <- file.path("output", "results_age_breakets_smooth.csv")
    data.table::fwrite(sim_age_bracket_results, out)
    out
  }, format = "file"),

  tar_target(
    sim_bic_table,
    {
      bic_windows <- unique(
        sim_age_bracket_results[
          ,
          .(predictor, seed_ix, years_ix, bic)
        ]
      )

      bic_by_seed <- bic_windows[
        ,
        .(bic_seed_average = mean(bic, na.rm = TRUE)),
        by = .(predictor, seed_ix)
      ]

      bic_by_seed[
        ,
        .(
          bic_average = mean(bic_seed_average, na.rm = TRUE),
          bic_sd = stats::sd(bic_seed_average, na.rm = TRUE),
          n_seeds = data.table::uniqueN(seed_ix),
          bic_div_1000 = mean(bic_seed_average, na.rm = TRUE) / 1000,
          bic_table = sprintf("%.1f", mean(bic_seed_average, na.rm = TRUE) / 1000)
        ),
        by = .(predictor)
      ][
        order(match(predictor, c("lc", "apc", "rh")))
      ]
    }
  ),

  tar_target(lc_mse_group_1_pdf, {
    out <- "figures/lc_mse_age_classes_group_1.pdf"

    p <- scoring_plot(
      data = sim_age_bracket_results,
      metric = "mse",
      subpopulation = 1,
      predictor_choice = "lc"
    )

    p <- p +
      scale_y_continuous(
        labels = function(x) x * 1e5
      ) +
      #coord_cartesian(ylim = c(0, 20e-05))+
      scale_color_discrete(labels = c("A", "B", "C", "D")) +
      labs(x = "", y = "", color = "") +
      theme(
        legend.position = "inside",
        legend.position.inside = c(0.85, 0.85),
        legend.justification = c("left", "top"),
        legend.key.height = unit(1.5, "lines"),
        legend.background = element_blank(),
        legend.key = element_blank(),
        axis.text.x = element_text(angle = 90)
      )

    ggplot2::ggsave(
      filename = out,
      plot = p,
      width = 10,
      height = 10
    )

    out
  }, format = "file"),
  tar_target(lc_mse_group_1_app, {
    out <- "figures/lc_mse_age_classes_group_1_app.pdf"

    p <- scoring_plot(
      data = sim_age_bracket_results,
      metric = "mse",
      subpopulation = 1,
      predictor_choice = "lc"
    )

    p <- p +
      scale_y_continuous(
        labels = function(x) x * 1e5
      ) +
      coord_cartesian(ylim = c(0, 20e-05))+
      scale_color_discrete(labels = c("A", "B", "C", "D")) +
      labs(x = "", y = "", color = "") +
      theme(
        legend.position = "inside",
        legend.position.inside = c(0.85, 0.85),
        legend.justification = c("left", "top"),
        legend.key.height = unit(1.5, "lines"),
        legend.background = element_blank(),
        legend.key = element_blank(),
        axis.text.x = element_text(angle = 90)
      )

    ggplot2::ggsave(
      filename = out,
      plot = p,
      width = 10,
      height = 10
    )

    out
  }, format = "file"),
  tar_target(lc_mse_group_2_pdf, {
    out <- "figures/lc_mse_age_classes_group_2.pdf"

    p <- scoring_plot(
      data = sim_age_bracket_results,
      metric = "mse",
      subpopulation = 2,
      predictor_choice = "lc"
    )

    p <- p +
      scale_y_continuous(
        labels = function(x) x * 1e5
      ) +
      coord_cartesian(ylim = c(0, 0.7e-05))+
      scale_color_discrete(labels = c("A", "B", "C", "D")) +
      labs(x = "", y = "", color = "") +
      theme(
        legend.position = "inside",
        legend.position.inside = c(0.85, 0.85),
        legend.justification = c("left", "top"),
        legend.key.height = unit(1.5, "lines"),
        legend.background = element_blank(),
        legend.key = element_blank(),
        axis.text.x = element_text(angle = 90)
      )

    ggplot2::ggsave(
      filename = out,
      plot = p,
      width = 10,
      height = 10
    )

    out
  }, format = "file"),
  tar_target(lc_mse_group_2_app, {
    out <- "figures/lc_mse_age_classes_group_2_app.pdf"

    p <- scoring_plot(
      data = sim_age_bracket_results,
      metric = "mse",
      subpopulation = 2,
      predictor_choice = "lc"
    )

    p <- p +
      scale_y_continuous(
        labels = function(x) x * 1e5
      ) +
      coord_cartesian(ylim = c(0, 20e-05))+
      scale_color_discrete(labels = c("A", "B", "C", "D")) +
      labs(x = "", y = "", color = "") +
      theme(
        legend.position = "inside",
        legend.position.inside = c(0.85, 0.85),
        legend.justification = c("left", "top"),
        legend.key.height = unit(1.5, "lines"),
        legend.background = element_blank(),
        legend.key = element_blank(),
        axis.text.x = element_text(angle = 90)
      )

    ggplot2::ggsave(
      filename = out,
      plot = p,
      width = 10,
      height = 10
    )

    out
  }, format = "file"),
  tar_target(lc_mse_group_3_pdf, {
    out <- "figures/lc_mse_age_classes_group_3.pdf"

    p <- scoring_plot(
      data = sim_age_bracket_results,
      metric = "mse",
      subpopulation = 3,
      predictor_choice = "lc"
    )

    p <- p +
      scale_y_continuous(
        labels = function(x) x * 1e5
      ) +
      # coord_cartesian(ylim = c(0, 1.2e-05))+
      scale_color_discrete(labels = c("A", "B", "C", "D")) +
      labs(x = "", y = "", color = "") +
      theme(
        legend.position = "inside",
        legend.position.inside = c(0.85, 0.85),
        legend.justification = c("left", "top"),
        legend.key.height = unit(1.5, "lines"),
        legend.background = element_blank(),
        legend.key = element_blank(),
        axis.text.x = element_text(angle = 90)
      )

    ggplot2::ggsave(
      filename = out,
      plot = p,
      width = 10,
      height = 10
    )

    out
  }, format = "file"),
  tar_target(lc_mse_group_3_app, {
    out <- "figures/lc_mse_age_classes_group_3_app.pdf"

    p <- scoring_plot(
      data = sim_age_bracket_results,
      metric = "mse",
      subpopulation = 3,
      predictor_choice = "lc"
    )

    p <- p +
      scale_y_continuous(
        labels = function(x) x * 1e5
      ) +
      coord_cartesian(ylim = c(0, 20e-05))+
      scale_color_discrete(labels = c("A", "B", "C", "D")) +
      labs(x = "", y = "", color = "") +
      theme(
        legend.position = "inside",
        legend.position.inside = c(0.85, 0.85),
        legend.justification = c("left", "top"),
        legend.key.height = unit(1.5, "lines"),
        legend.background = element_blank(),
        legend.key = element_blank(),
        axis.text.x = element_text(angle = 90)
      )

    ggplot2::ggsave(
      filename = out,
      plot = p,
      width = 10,
      height = 10
    )

    out
  }, format = "file"),

  tar_target(lc_deviance_superpop_pdf, {
    out <- "figures/lc_deviance_age_classes.pdf"

    p <- scoring_plot(
      data = sim_age_bracket_results,
      metric = "deviance",
      superpopulation = TRUE,
      predictor_choice = "lc"
    )

    p <- p +

      scale_color_discrete(labels = c("A", "B", "C", "D")) +
      labs(x = "", y = "", color = "") +
      theme(
        plot.background = element_rect(color = "red", fill = NA, linewidth = 1),
        legend.position = "inside",
        legend.position.inside = c(0.85, 0.85),
        legend.justification = c("left", "top"),
        legend.key.height = unit(1.5, "lines"),
        legend.background = element_blank(),
        legend.key = element_blank(),
        axis.text.x = element_text(angle = 90)
      )

    ggplot2::ggsave(
      filename = out,
      plot = p,
      width = 10,
      height = 10
    )
    out
  }, format = "file"),

  tar_target(lc_deviance_superpop_app, {
    out <- "figures/lc_deviance_age_classes_app.pdf"

    p <- scoring_plot(
      data = sim_age_bracket_results,
      metric = "deviance",
      superpopulation = TRUE,
      predictor_choice = "lc"
    )

    p <- p +
      coord_cartesian(ylim = c(0, 7))+
      scale_color_discrete(labels = c("A", "B", "C", "D")) +
      labs(x = "", y = "", color = "") +
      theme(
        plot.background = element_rect(color = "red", fill = NA, linewidth = 1),
        legend.position = "inside",
        legend.position.inside = c(0.85, 0.85),
        legend.justification = c("left", "top"),
        legend.key.height = unit(1.5, "lines"),
        legend.background = element_blank(),
        legend.key = element_blank(),
        axis.text.x = element_text(angle = 90)
      )

    ggplot2::ggsave(
      filename = out,
      plot = p,
      width = 10,
      height = 10
    )
    out
  }, format = "file"),

  tar_target(lc_deviance_group_1_pdf, {
    out <- "figures/lc_deviance_age_classes_group_1.pdf"

    p <- scoring_plot(
      data = sim_age_bracket_results,
      metric = "deviance",
      subpopulation = 1,
      predictor_choice = "lc"
    )

    p <- p +
      scale_color_discrete(labels = c("A", "B", "C", "D")) +
      labs(x = "", y = "", color = "") +
      theme(
        legend.position = "inside",
        legend.position.inside = c(0.85, 0.85),
        legend.justification = c("left", "top"),
        legend.key.height = unit(1.5, "lines"),
        legend.background = element_blank(),
        legend.key = element_blank(),
        axis.text.x = element_text(angle = 90)
      )

    ggplot2::ggsave(
      filename = out,
      plot = p,
      width = 10,
      height = 10
    )

    out
  }, format = "file"),


  tar_target(lc_deviance_group_1_app, {
    out <- "figures/lc_deviance_age_classes_group_1_app.pdf"

    p <- scoring_plot(
      data = sim_age_bracket_results,
      metric = "deviance",
      subpopulation = 1,
      predictor_choice = "lc"
    )

    p <- p +
      coord_cartesian(ylim = c(0, 7))+
      scale_color_discrete(labels = c("A", "B", "C", "D")) +
      labs(x = "", y = "", color = "") +
      theme(
        legend.position = "inside",
        legend.position.inside = c(0.85, 0.85),
        legend.justification = c("left", "top"),
        legend.key.height = unit(1.5, "lines"),
        legend.background = element_blank(),
        legend.key = element_blank(),
        axis.text.x = element_text(angle = 90)
      )

    ggplot2::ggsave(
      filename = out,
      plot = p,
      width = 10,
      height = 10
    )

    out
  }, format = "file"),


  tar_target(lc_deviance_group_2_pdf, {
    out <- "figures/lc_deviance_age_classes_group_2.pdf"

    p <- scoring_plot(
      data = sim_age_bracket_results,
      metric = "deviance",
      subpopulation = 2,
      predictor_choice = "lc"
    )

    p <- p +
      coord_cartesian(ylim = c(0.25, 1))+
      scale_color_discrete(labels = c("A", "B", "C", "D")) +
      labs(x = "", y = "", color = "") +
      theme(
        legend.position = "inside",
        legend.position.inside = c(0.85, 0.85),
        legend.justification = c("left", "top"),
        legend.key.height = unit(1.5, "lines"),
        legend.background = element_blank(),
        legend.key = element_blank(),
        axis.text.x = element_text(angle = 90)
      )

    ggplot2::ggsave(
      filename = out,
      plot = p,
      width = 10,
      height = 10
    )

    out
  }, format = "file"),

  tar_target(lc_deviance_group_2_app, {
    out <- "figures/lc_deviance_age_classes_group_2_app.pdf"

    p <- scoring_plot(
      data = sim_age_bracket_results,
      metric = "deviance",
      subpopulation = 2,
      predictor_choice = "lc"
    )

    p <- p +
      coord_cartesian(ylim = c(0, 7))+
      scale_color_discrete(labels = c("A", "B", "C", "D")) +
      labs(x = "", y = "", color = "") +
      theme(
        legend.position = "inside",
        legend.position.inside = c(0.85, 0.85),
        legend.justification = c("left", "top"),
        legend.key.height = unit(1.5, "lines"),
        legend.background = element_blank(),
        legend.key = element_blank(),
        axis.text.x = element_text(angle = 90)
      )

    ggplot2::ggsave(
      filename = out,
      plot = p,
      width = 10,
      height = 10
    )

    out
  }, format = "file"),

  tar_target(lc_deviance_group_3_pdf, {
    out <- "figures/lc_deviance_age_classes_group_3.pdf"

    p <- scoring_plot(
      data = sim_age_bracket_results,
      metric = "deviance",
      subpopulation = 3,
      predictor_choice = "lc"
    )

    p <- p +

      scale_color_discrete(labels = c("A", "B", "C", "D")) +
      labs(x = "", y = "", color = "") +
      theme(
        legend.position = "inside",
        legend.position.inside = c(0.85, 0.85),
        legend.justification = c("left", "top"),
        legend.key.height = unit(1.5, "lines"),
        legend.background = element_blank(),
        legend.key = element_blank(),
        axis.text.x = element_text(angle = 90)
      )

    ggplot2::ggsave(
      filename = out,
      plot = p,
      width = 10,
      height = 10
    )

    out
  }, format = "file"),

  tar_target(lc_deviance_group_3_app, {
    out <- "figures/lc_deviance_age_classes_group_3_app.pdf"

    p <- scoring_plot(
      data = sim_age_bracket_results,
      metric = "deviance",
      subpopulation = 3,
      predictor_choice = "lc"
    )

    p <- p +
      coord_cartesian(ylim = c(0, 7))+
      scale_color_discrete(labels = c("A", "B", "C", "D")) +
      labs(x = "", y = "", color = "") +
      theme(
        legend.position = "inside",
        legend.position.inside = c(0.85, 0.85),
        legend.justification = c("left", "top"),
        legend.key.height = unit(1.5, "lines"),
        legend.background = element_blank(),
        legend.key = element_blank(),
        axis.text.x = element_text(angle = 90)
      )

    ggplot2::ggsave(
      filename = out,
      plot = p,
      width = 10,
      height = 10
    )

    out
  }, format = "file"),

  ######### Section 4.4. simulation study end

  tar_target(msep_reference_age, 55),
  tar_target(msep_forecasting_horizon, 5),

  tar_target(
    lc_msep_results,
    run_msep_study(
      data = dt_0,
      data_pp = data_pp,
      mortality_models_fit = mortality_models_fit,
      credibility_model = credibility_model,
      ages_fit = ages_fit,
      years_fit = years_fit - min(years_fit),
      N_groups = N_groups,
      predictor = "lc",
      reference_age = msep_reference_age,
      forecasting_horizon = msep_forecasting_horizon,
      bias = min(years_fit),
      n_sim_sigma = 1000,
      separate_bootstrap_n = 500,
      sigma_seed = 1
    )
  ),

  tar_target(lc_msep_results_file, {
    fs::dir_create("output")
    out <- file.path("output", "lc_msep_results.rds")
    saveRDS(lc_msep_results, out)
    out
  }, format = "file"),
  tar_target(msep_ylim, NULL),
  tar_target(lc_msep_g1_pdf, {
    out <- sprintf("figures/lc_msep_%s_g1.pdf", msep_reference_age)
    write_msep_plot(lc_msep_results$groups[["1"]], out, ylim = msep_ylim)
  }, format = "file"),

  tar_target(lc_msep_g2_pdf, {
    out <- sprintf("figures/lc_msep_%s_g2.pdf", msep_reference_age)
    write_msep_plot(lc_msep_results$groups[["2"]], out, ylim = msep_ylim)
  }, format = "file"),
  tar_target(lc_msep_g3_pdf, {
    out <- sprintf("figures/lc_msep_%s_g3.pdf", msep_reference_age)
    write_msep_plot(lc_msep_results$groups[["3"]], out, ylim = msep_ylim)
  }, format = "file"),
  tar_target(msep_ylim_appendix, c(0,0.008)),

  tar_target(lc_msep_g1_app, {
    out <- sprintf("figures/lc_msep_%s_g1_appendix.pdf", msep_reference_age)
    write_msep_plot(lc_msep_results$groups[["1"]], out, ylim = msep_ylim_appendix)
  }, format = "file"),

  tar_target(lc_msep_g2_app, {
    out <- sprintf("figures/lc_msep_%s_g2_appendix.pdf", msep_reference_age)
    write_msep_plot(lc_msep_results$groups[["2"]], out, ylim = msep_ylim_appendix)
  }, format = "file"),
  tar_target(lc_msep_g3_app, {
    out <- sprintf("figures/lc_msep_%s_g3_appendix.pdf", msep_reference_age)
    write_msep_plot(lc_msep_results$groups[["3"]], out, ylim = msep_ylim_appendix)
  }, format = "file"),
  tarchetypes::tar_render(
    report_html,
    "report.Rmd"
  )



)

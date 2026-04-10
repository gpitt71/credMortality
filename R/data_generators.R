# Group and age effect -----
#' @import data.table
#' @import StMoMo
NULL
#' @export
data_generator_hmd_lt <- function(seed_input = 1,
                                  hmd_username = NULL,
                                  hmd_password = NULL,
                                  exposure_sup = 100000,
                                  exposure_sub = c(5000, 500)) {

  # --- fetch HMD data (unchanged) ---
  ita_data <- demography::hmd.mx(country = "ITA",
                                 hmd_username,
                                 hmd_password)


  E1 <- ita_data$pop[["female"]]
  D1 <- E1 * ita_data$rate[["female"]]

  periods_cols <- as.numeric(colnames(E1))
  periods_cols <- periods_cols - min(periods_cols)

  colnames(E1) <- colnames(D1) <- as.character(periods_cols)

  datahat <- structure(
    list(
      Dxt = D1,
      Ext = E1,
      ages = as.numeric(rownames(E1)),
      years = as.numeric(colnames(E1)),
      type = "central",
      series = "female",
      label = "total"
    ),
    class = "StMoMoData"
  )

  initial_data <- central2initial(datahat)

  tmp_occ <- initial_data$Dxt
  tmp_exp <- initial_data$Ext

  # --- data.table reshaping (pivot_longer + left_join replacement) ---
  tmp_occ <- as.data.frame(tmp_occ)
  tmp_exp <- as.data.frame(tmp_exp)
  tmp_occ[["Age"]] <- as.numeric(rownames(tmp_occ))
  tmp_exp[["Age"]] <- as.numeric(rownames(tmp_exp))

  data.table::setDT(tmp_occ)
  data.table::setDT(tmp_exp)

  # melt (like tidyr::pivot_longer) over the period columns
  dt1 <- data.table::melt(
    tmp_occ,
    id.vars = "Age",
    measure.vars = as.character(periods_cols),
    variable.name = "Period",
    value.name = "Dxt"
  )

  dt2 <- data.table::melt(
    tmp_exp,
    id.vars = "Age",
    measure.vars = as.character(periods_cols),
    variable.name = "Period",
    value.name = "Ext"
  )

  # left_join by Age, Period
  dt <- merge(dt1, dt2, by = c("Age", "Period"), all.x = TRUE, sort = FALSE)

  # --- mutate steps in data.table ---
  dt[, qx := Dxt / Ext]
  dt[, Period := as.numeric(as.character(Period))]
  dt[, Cohort := Period - Age]

  # mutate_all(~replace(., is.na(.) | is.infinite(.), 1))
  for (j in names(dt)) {
    v <- dt[[j]]
    if (is.numeric(v)) {
      v[is.na(v) | is.infinite(v)] <- 1
      data.table::set(dt, j = j, value = v)
    }
  }


  dt[qx / (1 - qx) <=0, intercept := 1e-08]
  dt[qx / (1 - qx) >0, intercept := log(qx / (1 - qx))]

  # mutate_all(~replace(., is.infinite(.), 0))
  for (j in names(dt)) {
    v <- dt[[j]]
    if (is.numeric(v)) {
      v[is.infinite(v)] <- 1e-08
      data.table::set(dt, j = j, value = v)
    }
  }


  # --- theta generation  ---
{  set.seed(seed = seed_input)
  age_eff_1 <- 1 - runif(length(0:110), .2, .3)
  age_eff_2 <- 1 + runif(length(0:110), .2, .3)}

  dt_theta <- data.table::data.table(
    Age = 0:110,
    Theta1x = age_eff_1,
    Theta2x = age_eff_2
  )

  # left_join(dt, dt_theta, by="Age") + mutate(...)
  dt <- merge(dt, dt_theta, by = "Age", all.x = TRUE, sort = FALSE)

  dt[, `:=`(
    predictor1 = Theta1x * exp(intercept),
    predictor2 = Theta2x * exp(intercept))][,
                                           `:=`(p1 = predictor1 / (1 + predictor1),
    p2 = predictor2 / (1 + predictor2),
    p = qx
  )]

  # select(Age,Period,Cohort,p,p1,p2)
  dt <- dt[, .(Age, Period, Cohort, p, p1, p2)]

  # arrange(Cohort, Age, Period)
  data.table::setorder(dt, Cohort, Age, Period)

  # --- cohort-wise simulations (unchanged) ---
  dt[, c("Lx", "Dx_tmp") := binomial_simulator(
    starting_exposure = exposure_sup,
    probabilities = p
  ), by = .(Cohort)]
  dt[, Dx_tmp := as.numeric(Dx_tmp)]

  dt[, c("Lx1", "Dx1_tmp") := binomial_simulator(
    starting_exposure = exposure_sub[1],
    probabilities = p1
  ), by = .(Cohort)]
  dt[, Dx1_tmp := as.numeric(Dx1_tmp)]

  dt[, c("Lx2", "Dx2_tmp") := binomial_simulator(
    starting_exposure = exposure_sub[2],
    probabilities = p2
  ), by = .(Cohort)]
  dt[, Dx2_tmp := as.numeric(Dx2_tmp)]

  # arrange again
  data.table::setorder(dt, Cohort, Age, Period)

  dt[, c("Ext", "Dx") := average_exposure_and_deaths(raw_exposure = Lx), by = .(Cohort)]
  dt[, c("Ext1", "Dx1") := average_exposure_and_deaths(raw_exposure = Lx1), by = .(Cohort)]
  dt[, c("Ext2", "Dx2") := average_exposure_and_deaths(raw_exposure = Lx2), by = .(Cohort)]

  dt
}


# Print group and age effect ----
#' @export
print_age_effects <- function(seed_input=1,
                                  hmd_username=NULL,
                                  hmd_password=NULL,
                                  exposure_sup=100000,
                                  exposure_sub=c(5000,
                                                 500)){


  # --- fetch HMD data (unchanged) ---
  ita_data <- demography::hmd.mx(country = "ITA",
                                 hmd_username,
                                 hmd_password)


  E1 <- ita_data$pop[["female"]]
  D1 <- E1 * ita_data$rate[["female"]]

  periods_cols <- as.numeric(colnames(E1))
  periods_cols <- periods_cols - min(periods_cols)

  colnames(E1) <- colnames(D1) <- as.character(periods_cols)

  datahat <- structure(
    list(
      Dxt = D1,
      Ext = E1,
      ages = as.numeric(rownames(E1)),
      years = as.numeric(colnames(E1)),
      type = "central",
      series = "female",
      label = "total"
    ),
    class = "StMoMoData"
  )

  initial_data <- central2initial(datahat)

  tmp_occ <- initial_data$Dxt
  tmp_exp <- initial_data$Ext

  # --- data.table reshaping (pivot_longer + left_join replacement) ---
  tmp_occ <- as.data.frame(tmp_occ)
  tmp_exp <- as.data.frame(tmp_exp)
  tmp_occ[["Age"]] <- as.numeric(rownames(tmp_occ))
  tmp_exp[["Age"]] <- as.numeric(rownames(tmp_exp))

  data.table::setDT(tmp_occ)
  data.table::setDT(tmp_exp)

  # melt (like tidyr::pivot_longer) over the period columns
  dt1 <- data.table::melt(
    tmp_occ,
    id.vars = "Age",
    measure.vars = as.character(periods_cols),
    variable.name = "Period",
    value.name = "Dxt"
  )

  dt2 <- data.table::melt(
    tmp_exp,
    id.vars = "Age",
    measure.vars = as.character(periods_cols),
    variable.name = "Period",
    value.name = "Ext"
  )

  # left_join by Age, Period
  dt <- merge(dt1, dt2, by = c("Age", "Period"), all.x = TRUE, sort = FALSE)

  # --- mutate steps in data.table ---
  dt[, qx := Dxt / Ext]
  dt[, Period := as.numeric(as.character(Period))]
  dt[, Cohort := Period - Age]

  # mutate_all(~replace(., is.na(.) | is.infinite(.), 1))
  for (j in names(dt)) {
    v <- dt[[j]]
    if (is.numeric(v)) {
      v[is.na(v) | is.infinite(v)] <- 1
      data.table::set(dt, j = j, value = v)
    }
  }


  dt[qx / (1 - qx) <=0, intercept := 1e-08]
  dt[qx / (1 - qx) >0, intercept := log(qx / (1 - qx))]

  # mutate_all(~replace(., is.infinite(.), 0))
  for (j in names(dt)) {
    v <- dt[[j]]
    if (is.numeric(v)) {
      v[is.infinite(v)] <- 1e-08
      data.table::set(dt, j = j, value = v)
    }
  }


  {
    set.seed(seed = seed_input)
    age_eff_1 <- 1 - runif(length(0:110), .2, .3)
    age_eff_2 <- 1 + runif(length(0:110), .2, .3)
    }


  out<-list(age_eff_1=age_eff_1,
       age_eff_2=age_eff_2)

  return(out)

}








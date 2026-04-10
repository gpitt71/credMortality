#' data_preprocessing
#'
#' Purpose: data_preprocessing for the mortality forecasting reproducibility workflow.
#' @param data Input used by `data_preprocessing`.
#' @param N_groups Input used by `data_preprocessing`.
#' @param bias Input used by `data_preprocessing`.
#' @param years_fit Input used by `data_preprocessing`.
#' @param ages_fit Input used by `data_preprocessing`.
#' @return An R object produced by this step of the workflow.
#' @export
data_preprocessing <- function(data,
                               N_groups = 3,
                               bias = 0,
                               years_fit = 110:150,
                               ages_fit = 50:65) {

  setDT(data)
  data_copy <- copy(data)

  # Dx / Ext
  tmp_occ <- wide_mat(data_copy, "Dx")
  tmp_occ <- tmp_occ[as.character(ages_fit), years_fit + bias, drop = FALSE]

  tmp_exp <- wide_mat(data_copy, "Ext")
  tmp_exp <- tmp_exp[as.character(ages_fit), years_fit + bias, drop = FALSE]

  colnames(tmp_occ) <- colnames(tmp_exp) <- years_fit - min(years_fit)

  # Dx1 / Ext1
  tmp_occ1 <- wide_mat(data_copy, "Dx1")
  tmp_occ1 <- tmp_occ1[as.character(ages_fit), years_fit + bias, drop = FALSE]

  tmp_exp1 <- wide_mat(data_copy, "Ext1")
  tmp_exp1 <- tmp_exp1[as.character(ages_fit), years_fit + bias, drop = FALSE]

  colnames(tmp_occ1) <- colnames(tmp_exp1) <- years_fit - min(years_fit)

  # Dx2 / Ext2
  tmp_occ2 <- wide_mat(data_copy, "Dx2")
  tmp_occ2 <- tmp_occ2[as.character(ages_fit), years_fit + bias, drop = FALSE]

  tmp_exp2 <- wide_mat(data_copy, "Ext2")
  tmp_exp2 <- tmp_exp2[as.character(ages_fit), years_fit + bias, drop = FALSE]

  colnames(tmp_occ2) <- colnames(tmp_exp2) <- years_fit - min(years_fit)
  datahat_tot <- structure(
    list(
      Dxt = tmp_occ + tmp_occ1 + tmp_occ2,
      Ext = tmp_exp + tmp_exp1 + tmp_exp2,
      ages = ages_fit,
      years = years_fit - min(years_fit),
      type = 'central',
      series = 'female',
      label = 'total'
    ),
    class = "StMoMoData"
  )

  l <- list()

  list_of_extra_exposures <- list()


  list_of_extra_exposures[[1]] <- list(Dxt = tmp_occ, Ext = tmp_exp)

  list_of_extra_exposures[[2]] <- list(Dxt = tmp_occ1, Ext = tmp_exp1)

  list_of_extra_exposures[[3]] <- list(Dxt = tmp_occ2, Ext = tmp_exp2)


  # l[['datahat']] <- datahat #population 0
  # l[['datahat_tot']] <- datahat_tot #superpop
  # l[['list_of_extra_exposures']] <- list_of_extra_exposures #subpopulations

  l[['superpopulation']] <- datahat_tot
  l[['subpopulations']] <- list_of_extra_exposures

  return(l)

}

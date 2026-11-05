#' Italy female mortality data from the Human Mortality Database
#'
#' Female age-specific deaths and central exposures for Italy from the Human
#' Mortality Database (HMD), stored as a [StMoMo::StMoMoData] object. This is
#' the frozen input used to generate the synthetic populations in the paper.
#'
#' @format A `StMoMoData` list with the following components:
#' \describe{
#'   \item{Dxt}{A 111 by 151 matrix of deaths, with ages 0--110 in rows and
#'   calendar years 1872--2022 in columns.}
#'   \item{Ext}{A 111 by 151 matrix of central exposures with the same
#'   dimensions and dimnames as `Dxt`.}
#'   \item{ages}{A numeric vector of ages 0--110.}
#'   \item{years}{An integer vector of calendar years 1872--2022.}
#'   \item{type}{The exposure type, `"central"`.}
#'   \item{series}{The population series, `"female"`.}
#'   \item{label}{The population label, `"Italy"`.}
#' }
#' @details The underlying HMD Italy snapshot was downloaded on 10 February
#' 2026. `ITAFemaleData` was constructed from `STATS/Mx_1x1.txt` and
#' `STATS/Exposures_1x1.txt` in that snapshot. Their respective SHA-256
#' checksums are
#' `07ddeb64683310cc61d9686ce25b8a389a0d21621e25298e307380e6f1b5fbaa`
#' and
#' `6c8469db3fcf3e0be9a5cb73c0270857909c1a1e68a5db5816fe3cd8f0e819f9`.
#' Deaths are the product of the female central exposure and mortality-rate
#' matrices, matching the input construction used for the reported results.
#' The preparation recipe is in `data-raw/ITAFemaleData.R`.
#' @source Human Mortality Database, University of California, Berkeley (USA),
#' and Max Planck Institute for Demographic Research (Germany),
#' \url{https://www.mortality.org/}.
#' @references Human Mortality Database (2026). Human Mortality Database.
#' Available at \url{https://www.mortality.org/}.
#' @seealso [StMoMo::StMoMoData], [StMoMo::EWMaleData]
#' @keywords datasets
#' @name ITAFemaleData
"ITAFemaleData"

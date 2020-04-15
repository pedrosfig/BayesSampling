#' Full Person-level Population Database
#'
#' This data set corresponds to some socioeconomic variables from 150266 people of a city in a particular year.
#'
#' @docType data
#'
#' @usage data(BigCity)
#'
#' @format A data.frame with 150266 rows and 12 variables:
#' \describe{
#' \item{HHID}{The identifier of the household. It corresponds to an alphanumeric sequence (four letters and five digits).}
#' \item{PersonID}{The identifier of the person within the household. NOTE it is not a unique identifier of a person for the whole population. It corresponds to an alphanumeric sequence (five letters and two digits).}
#' \item{Stratum}{Households are located in geographic strata. There are 119 strata across the city.}
#' \item{PSU}{Households are clustered in cartographic segments defined as primary sampling units (PSU). There are 1664 PSU and they are nested within strata.}
#' \item{Zone}{Segments clustered within strata can be located within urban or rural areas along the city.}
#' \item{Sex}{Sex of the person.}
#' \item{Income}{Per capita monthly income.}
#' \item{Expenditure}{Per capita monthly expenditure.}
#' \item{Employment}{A person's employment status.}
#' \item{Poverty}{This variable indicates whether the person is poor or not. It depends on income.}
#' }
#'
#' @references Package ‘TeachingSampling’; see \code{\link[TeachingSampling]{BigCity}}
#'
#' @source \url{https://CRAN.R-project.org/package=TeachingSampling}
#'
"BigCity"

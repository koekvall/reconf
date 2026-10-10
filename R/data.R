#' FEV1 data from the Six Cities Study
#'
#' Repeated measurements of forced expiratory volume in one second (FEV1) on
#' a random sample of 300 girls from Topeka, Kansas, in the Six Cities Study
#' of Air Pollution and Health, with one to twelve measurements per girl.
#'
#' @format A data frame with 1994 rows and 6 variables:
#' \describe{
#'   \item{id}{subject identifier.}
#'   \item{ht}{height (metres).}
#'   \item{age}{age (years).}
#'   \item{baseht}{height at the first measurement (metres).}
#'   \item{baseage}{age at the first measurement (years).}
#'   \item{logfev1}{natural logarithm of FEV1 (litres).}
#' }
#' @source Fitzmaurice, G. M., Laird, N. M., and Ware, J. H. (2011).
#'   \emph{Applied Longitudinal Analysis}, 2nd edition. Wiley.
#'   \url{https://content.sph.harvard.edu/fitzmaur/ala2e/}. Data courtesy of
#'   D. W. Dockery.
#' @references Dockery, D. W., Berkey, C. S., Ware, J. H., Speizer, F. E., and
#'   Ferris, B. G. (1983). Distribution of FVC and FEV1 in children 6 to 11
#'   years old. \emph{American Review of Respiratory Disease} 128, 405--412.
"fev1"

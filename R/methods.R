#' @importFrom generics tidy
#' @export
generics::tidy

#' Print method for reconf confidence intervals
#'
#' @param x A \code{reconf_ci} object from \code{\link{ci_lmer}} or
#'   \code{\link{ci_all_lmer}}.
#' @param digits Number of significant digits. Default 4.
#' @param ... Unused.
#' @return \code{x}, invisibly.
#' @method print reconf_ci
#' @export
print.reconf_ci <- function(x, digits = 4, ...) {
  cat(sprintf("%s%% score-based confidence intervals (%s)\n\n",
              format(100 * attr(x, "level")), attr(x, "method")))
  m <- x
  attributes(m) <- attributes(x)[c("dim", "dimnames")]
  print(signif(m, digits))
  invisible(x)
}

#' Tidy a reconf confidence interval object
#'
#' Returns the intervals as a data frame with the column names used
#' throughout the broom ecosystem, for use in pipelines and plotting.
#'
#' @param x A \code{reconf_ci} object from \code{\link{ci_lmer}} or
#'   \code{\link{ci_all_lmer}}.
#' @param ... Unused.
#' @return A data frame with columns \code{term}, \code{estimate},
#'   \code{conf.low}, and \code{conf.high}.
#' @export
tidy.reconf_ci <- function(x, ...) {
  data.frame(term = rownames(x),
             estimate = x[, "estimate"],
             conf.low = x[, "lower"],
             conf.high = x[, "upper"],
             row.names = NULL)
}

#' Print method for reconf score tests
#'
#' @param x A \code{reconf_test} object from \code{\link{score_test_all_lmer}}.
#' @param digits Number of significant digits for the statistics. Default 4.
#' @param ... Unused.
#' @return \code{x}, invisibly.
#' @method print reconf_test
#' @export
print.reconf_test <- function(x, digits = 4, ...) {
  cat(sprintf("Score tests of zero covariance parameters (%s)\n\n",
              attr(x, "method")))
  m <- x
  attributes(m) <- attributes(x)[c("dim", "dimnames")]
  out <- data.frame(statistic = signif(m[, "statistic"], digits),
                    df = m[, "df"],
                    p.value = format.pval(m[, "p.value"], digits = digits),
                    row.names = rownames(m))
  print(out)
  invisible(x)
}

#' Tidy a reconf score test object
#'
#' @param x A \code{reconf_test} object from \code{\link{score_test_all_lmer}}.
#' @param ... Unused.
#' @return A data frame with columns \code{term}, \code{statistic},
#'   \code{df}, and \code{p.value}.
#' @export
tidy.reconf_test <- function(x, ...) {
  data.frame(term = rownames(x),
             statistic = x[, "statistic"],
             df = x[, "df"],
             p.value = x[, "p.value"],
             row.names = NULL)
}

#' @importFrom generics tidy
#' @export
generics::tidy

#' Print method for confidence intervals
#'
#' @param x A \code{vc_ci} object from \code{\link{vc_ci}}.
#' @param digits Number of significant digits. Default 4.
#' @param ... Unused.
#' @return \code{x}, invisibly.
#' @method print vc_ci
#' @export
print.vc_ci <- function(x, digits = 4, ...) {
  cat(sprintf(paste("%s%% confidence intervals for variance-covariance",
                    "parameters (%s, %s)\n\n"),
              format(100 * attr(x, "level")), attr(x, "statistic"),
              attr(x, "method")))
  m <- x
  attributes(m) <- attributes(x)[c("dim", "dimnames")]
  # Fixed (non-scientific) notation as in lme4 summaries; a common format
  # across the table keeps the decimal points aligned.
  print(format(m, digits = digits, scientific = FALSE),
        quote = FALSE, right = TRUE)
  invisible(x)
}

#' Tidy confidence intervals
#'
#' Returns the intervals as a data frame with the column names used
#' throughout the broom ecosystem, for use in pipelines and plotting.
#'
#' @param x A \code{vc_ci} object from \code{\link{vc_ci}}.
#' @param ... Unused.
#' @return A data frame with columns \code{term}, \code{estimate},
#'   \code{conf.low}, and \code{conf.high}.
#' @export
tidy.vc_ci <- function(x, ...) {
  data.frame(term = rownames(x),
             estimate = x[, "estimate"],
             conf.low = x[, "lower"],
             conf.high = x[, "upper"],
             row.names = NULL)
}

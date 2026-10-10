#' Model for inference on covariance parameters
#'
#' Builds the object that \code{\link{vc_ci}}, \code{\link{vc_test}}, and
#' \code{\link{vc_loglik}} work with: the likelihood of a linear mixed model
#' as a function of its covariance parameters, with quantities precomputed
#' once, and the parameter names and estimates. Those functions also accept
#' a fitted model and build this object themselves; building it once avoids
#' repeating the precomputation across calls.
#'
#' @param x a fitted linear mixed model of class \code{lmerMod} from
#'   \code{lme4::lmer}; for \code{print}, a \code{vc_model}.
#' @param method computational method for the likelihood; see
#'   \code{\link{vc_ci}}.
#' @param digits number of significant digits printed.
#' @param ... unused.
#'
#' @return An object of class \code{"vc_model"}. Printing it shows the
#'   covariance parameters, in the order that indices in \code{vc_ci} and
#'   the other functions refer to, and their estimates.
#'
#' @details For an lme4 fit the likelihood is the restricted likelihood if
#'   the model was fitted by REML, and the estimates are those of lme4. The
#'   parameters are ordered as in
#'   \code{as.data.frame(VarCorr(x), order = "lower.tri")}, with the error
#'   variance last. Prior weights \eqn{w_i} and offsets are supported; as in
#'   \code{lme4}, observation \eqn{i} then has error variance
#'   \eqn{\psi_r / w_i}, with \eqn{\psi_r} the error variance parameter.
#'
#' @examples
#' library(lme4)
#' fit <- lmer(Reaction ~ Days + (Days | Subject), data = sleepstudy)
#' m <- vc_model(fit)
#' m
#' vc_ci(m, parm = "var_Days|Subject")
#' vc_test(m, parm = "var_Days|Subject")
#' @export
vc_model <- function(x, ...) UseMethod("vc_model")

#' @export
vc_model.default <- function(x, ...) {
  stop("vc_model() has no method for class '", class(x)[1], "'",
       call. = FALSE)
}

#' @rdname vc_model
#' @method print vc_model
#' @export
print.vc_model <- function(x, digits = 4, ...) {
  cat("Covariance parameters of a linear mixed model (",
      if (x$REML) "REML" else "ML", ")\n\n", sep = "")
  print(data.frame(estimate = x$psi_hat, row.names = x$names),
        digits = digits)
  invisible(x)
}

# The vc_model of object: object itself if it is one, otherwise
# vc_model(object, method = method). method_given says whether the caller
# set method, which a vc_model already fixes.
.as_vc_model <- function(object, method, method_given) {
  if (inherits(object, "vc_model")) {
    if (method_given) {
      stop("method is set when the model is built; see vc_model()",
           call. = FALSE)
    }
    return(object)
  }
  vc_model(object, method = method)
}

# The vc_model object, whatever the model's source: the likelihood as a
# function of psi (see make_loglik), the number r of covariance parameters,
# their names, estimates psi_hat, and structure (see .psi_structure), and
# whether the likelihood is restricted. Y, X, Z, Hlist, REML, and method are
# as for make_loglik; names default to psi1, ..., psir.
.make_setup <- function(Y, X, Z, Hlist, REML, psi_hat, names = NULL,
                        method = "auto") {
  if (is.null(X)) X <- matrix(0, length(Y), 0)
  REML    <- REML && ncol(X) > 0
  r       <- length(Hlist) + 1L
  precomp <- get_precomp(Y = Y, X = X, Z = Z, Hlist = Hlist, method = method)
  structure(
    list(ll = make_loglik(Y = Y, X = X, Z = Z, Hlist = Hlist, REML = REML,
                          precomp = precomp),
         REML = REML, precomp = precomp, psi_hat = psi_hat, r = r,
         names = if (is.null(names)) paste0("psi", seq_len(r)) else names,
         structure = .psi_structure(Hlist)),
    class = "vc_model")
}

#' Log-likelihood as a function of the covariance parameters
#'
#' Returns the log-likelihood of a linear mixed model, or its restricted
#' log-likelihood for a restricted likelihood fit, as a function of the
#' covariance parameters.
#'
#' @inheritParams vc_ci
#'
#' @return A function of \code{psi}, the covariance parameters in the order
#'   of \code{\link{vc_model}}, with logical arguments \code{get_val},
#'   \code{get_score}, and \code{get_inf} selecting what to compute and
#'   \code{expected} selecting the expected or the observed information, all
#'   \code{TRUE} by default. It returns a list with components \code{value},
#'   \code{score}, and \code{information}, the last two named by the
#'   parameters; components not requested are \code{NULL}.
#'
#' @details The likelihood is maximized over the fixed effects. Outside the
#'   parameter set, the set of covariance parameters for which the
#'   covariance matrix of the response is positive definite, the value is
#'   \code{-Inf} and the score and information are zero. With prior weights
#'   \eqn{w_i} in an lme4 fit, the value differs from \code{logLik()} by
#'   \eqn{\sum_i \log(w_i)/2}.
#'
#' @examples
#' library(lme4)
#' fit <- lmer(Reaction ~ Days + (Days | Subject), data = sleepstudy)
#' ll <- vc_loglik(fit)
#' ll(c(600, 10, 35, 650))
#' @export
vc_loglik <- function(object,
                      method = c("auto", "q_side", "n_side", "spectral")) {
  model <- .as_vc_model(object, match.arg(method), !missing(method))
  ll  <- model$ll
  r   <- model$r
  nms <- model$names
  function(psi, get_val = TRUE, get_score = TRUE, get_inf = TRUE,
           expected = TRUE) {
    .check_psi(psi, r = r)
    .check_flags(list(get_val = get_val, get_score = get_score,
                      get_inf = get_inf, expected = expected))
    out <- ll(psi, get_val = get_val, get_score = get_score,
              get_inf = get_inf, expected = expected)
    list(value = if (get_val) out$value,
         score = if (get_score) stats::setNames(as.numeric(out$score), nms),
         information = if (get_inf) {
           matrix(out$inf_mat, r, r, dimnames = list(nms, nms))
         })
  }
}

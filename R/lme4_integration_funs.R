#' Covariance structure of an lme4 fit
#'
#' The structure matrices \eqn{H_j} with \eqn{\Psi = \sum_j \psi_j H_j} as
#' in \code{make_loglik}, the estimates, and the parameter names, built from
#' the covariance structure each random-effect term declares and never
#' inferred from estimates. \eqn{H_j} has ones at the positions of
#' \eqn{\psi_j} in \eqn{\Psi} and zeros elsewhere.
#'
#' @param lmerfit An \code{lmerMod} object from \code{lme4::lmer}.
#'
#' @return A list with \code{Hlist}, \eqn{r - 1} sparse symmetric matrices,
#'   one per random-effect covariance parameter; \code{psi_hat}, the
#'   \eqn{r} estimates with the error variance last; and \code{names}. The
#'   parameters of a term are numbered in the column-major order of the
#'   lower triangle of its covariance matrix, which for an unstructured term
#'   is the order of \code{as.data.frame(VarCorr(lmerfit), order =
#'   "lower.tri")}.
#'
#' @details Terms are supported if their covariance matrix is linear in its
#'   variances and covariances; see \code{.term_pattern}. Other terms are an
#'   error. The estimates are read from \eqn{\sigma^2 \Lambda \Lambda'}, the
#'   fitted covariance matrix of the random effects, which the structure
#'   must reproduce up to rounding at every level of every term; otherwise
#'   lme4 lays out the random effects differently than assumed here, and
#'   that is an error too.
#'
#' @noRd
.lmer_structure <- function(lmerfit) {
  cnms <- lme4::getME(lmerfit, "cnms")
  Gp   <- lme4::getME(lmerfit, "Gp")
  fl   <- lme4::getME(lmerfit, "flist")
  nc   <- lengths(cnms)
  nlev <- vapply(attr(fl, "assign"), function(i) nlevels(fl[[i]]), integer(1))
  grp  <- names(lme4::VarCorr(lmerfit))
  # Term t has columns Gp[t] + 1, ..., Gp[t + 1] of Z, a block of nc[t]
  # consecutive columns per level
  if (!(length(Gp) == length(nc) + 1L && length(nlev) == length(nc) &&
        length(grp) == length(nc) && all(diff(Gp) == nc * nlev))) {
    .lmer_layout_error()
  }
  q    <- Gp[length(Gp)]
  decl <- .lmer_declared(lmerfit, nc)
  Ps   <- lapply(seq_along(nc), function(t) {
    P <- .term_pattern(nc[t], decl$cls[t], decl$hom[t])
    if (is.null(P)) {
      stop("the random-effect term for ", grp[t], " has covariance ",
           "structure ", decl$label[t], " with ", nc[t], " effects, which ",
           "is not linear in its variances and covariances. Supported: ",
           "unstructured, diag(), cs(hom = TRUE), and cs() or ar1() with ",
           "two effects", call. = FALSE)
    }
    P
  })

  # Upper triangle of the parameter map of Psi as rows (i, j, k): entry
  # (i, j) is parameter k. Each term's pattern is repeated in the blocks of
  # its levels, and the terms' parameters are numbered consecutively.
  # Entries not listed are zero.
  shift <- cumsum(c(0L, vapply(Ps, max, integer(1))))
  map <- do.call(rbind, lapply(seq_along(Ps), function(t) {
    ab   <- which(Ps[[t]] > 0 & upper.tri(Ps[[t]], diag = TRUE),
                  arr.ind = TRUE)
    base <- Gp[t] + nc[t] * (seq_len(nlev[t]) - 1L)
    cbind(as.vector(outer(ab[, 1], base, `+`)),
          as.vector(outer(ab[, 2], base, `+`)),
          rep(Ps[[t]][ab] + shift[t], nlev[t]))
  }))
  at <- unname(split(seq_len(nrow(map)),
                     factor(map[, 3], levels = seq_len(shift[length(shift)]))))
  Hlist <- lapply(at, function(s) {
    Matrix::forceSymmetric(Matrix::sparseMatrix(map[s, 1], map[s, 2], x = 1,
                                                dims = c(q, q)), uplo = "U")
  })

  # Estimates at each parameter's first position in sigma^2 Lambda Lambda'.
  # lme4 can return a fit with sigma NaN when its optimizer stops early.
  s2 <- stats::sigma(lmerfit)^2
  if (!is.finite(s2)) {
    stop("lme4 reports a residual standard deviation of ",
         stats::sigma(lmerfit), ", so the fit may be incomplete; refitting ",
         "with another optimizer, such as lmerControl(optimizer = ",
         "\"bobyqa\"), may help", call. = FALSE)
  }
  Psi_fit <- s2 * Matrix::tcrossprod(lme4::getME(lmerfit, "Lambda"))
  psi_hat <- Psi_fit[map[vapply(at, `[`, integer(1), 1L), 1:2, drop = FALSE]]
  .check_reproduces(Psi_fit, psi_hat, Hlist)
  list(Hlist = Hlist, psi_hat = c(psi_hat, s2),
       names = c(unlist(Map(.term_names, Ps, cnms, grp), use.names = FALSE),
                 "var_Residual"))
}

# Stop unless sum_j psi_hat[j] Hlist[[j]] reproduces Psi_fit, the fitted
# covariance matrix of the random effects, up to rounding. Rounding error in
# entry (a, b) of sigma^2 Lambda Lambda' is of order eps * sqrt(d_a d_b),
# with d the diagonal, so each entry is compared on that scale: a term with
# small variances is checked as tightly as one with large variances. The
# floor eps * max(d) keeps a variance that underflows to zero from tripping
# the check.
.check_reproduces <- function(Psi_fit, psi_hat, Hlist) {
  d <- Matrix::diag(Psi_fit)
  if (!all(is.finite(d))) {
    stop("the fitted covariance matrix of the random effects is not finite",
         call. = FALSE)
  }
  E   <- methods::as(Psi_fit - Reduce(`+`, Map(`*`, psi_hat, Hlist)),
                     "TsparseMatrix")
  eps <- .Machine$double.eps
  tol <- sqrt(eps) * (sqrt(d[E@i + 1L] * d[E@j + 1L]) + eps * max(d))
  if (any(abs(E@x) > tol)) .lmer_layout_error()
  invisible(TRUE)
}

# Stop when an lme4 fit does not have the random-effects layout that
# .lmer_structure assumes
.lmer_layout_error <- function() {
  stop("the random effects of this lme4 fit are laid out differently than ",
       "reconf expects; please report this, with the lme4 version, at ",
       "https://github.com/koekvall/reconf/issues", call. = FALSE)
}

# Declared covariance structure of each random-effect term of an lme4 fit
# with nc effects per term: cls, the structure ("us", "diag", "cs", "ar1",
# or NA for one reconf does not know); hom, TRUE for a variance common to
# the term's effects (NA where the class has no hom slot); and label, for
# messages. lme4 before 2.0 has only unstructured terms and no getReCovs.
.lmer_declared <- function(lmerfit, nc) {
  if (!("getReCovs" %in% getNamespaceExports("lme4"))) {
    n <- length(nc)
    return(list(cls = rep("us", n), hom = rep(NA, n), label = rep("us", n)))
  }
  covs <- getExportedValue("lme4", "getReCovs")(lmerfit)
  slot_or_na <- function(cv, s) {
    if (methods::.hasSlot(cv, s)) methods::slot(cv, s) else NA
  }
  # The declared number of effects, where there is one, must match the fit's
  nc_ok <- length(covs) == length(nc) &&
    all(vapply(seq_along(nc), function(t) {
      n <- slot_or_na(covs[[t]], "nc")
      identical(n, NA) || identical(as.integer(n), nc[[t]])
    }, logical(1)))
  if (!nc_ok) .lmer_layout_error()
  hom <- vapply(covs, function(cv) {
    h <- slot_or_na(cv, "hom")
    if (is.logical(h) && length(h) == 1L) h else NA
  }, logical(1), USE.NAMES = FALSE)
  cl  <- vapply(covs, function(cv) class(cv)[1], character(1),
                USE.NAMES = FALSE)
  cls <- unname(c(Covariance.us = "us", Covariance.diag = "diag",
                  Covariance.cs = "cs", Covariance.ar1 = "ar1")[cl])
  label <- ifelse(is.na(cls), paste0("'", cl, "'"),
                  paste0(cls, "(", ifelse(is.na(hom), "",
                                          paste0("hom = ", hom)), ")"))
  list(cls = cls, hom = hom, label = label)
}

# Pattern of the covariance matrix of one term's nc random effects under the
# declared structure cls and hom (see .lmer_declared): P[a, b] = k if entry
# (a, b) is the term's k-th covariance parameter, numbered in the
# column-major order of the lower triangle, and P[a, b] = 0 if the entry is
# zero by construction. NULL if the matrix is not linear in its variances
# and covariances or the structure is unknown.
.term_pattern <- function(nc, cls, hom) {
  if (is.na(cls)) return(NULL)
  if (nc == 1L) return(matrix(1L))  # one variance, whatever the structure
  if (cls != "us" && is.na(hom)) return(NULL)
  # With two effects, ar1 is cs, and cs with separate variances allows every
  # positive semidefinite matrix, as us does
  if (nc == 2L && cls == "ar1") cls <- "cs"
  if (nc == 2L && cls == "cs" && !hom) cls <- "us"
  P  <- matrix(0L, nc, nc)
  lt <- lower.tri(P, diag = TRUE)
  if (cls == "us") {
    P[lt] <- seq_len(sum(lt))
    P[upper.tri(P)] <- t(P)[upper.tri(P)]
  } else if (cls == "diag") {
    diag(P) <- if (hom) 1L else seq_len(nc)
  } else if (cls == "cs" && hom) {
    P[] <- 2L
    diag(P) <- 1L
  } else {
    return(NULL)
  }
  P
}

# Names of the parameters of one term with pattern P, effect names cnm, and
# group name grp. They follow lme4's names for standard deviations and
# correlations (sd_x|g, cor_y.x|g) on the variance scale: var_x|g and
# cov_y.x|g, and var_*|g and cov_*|g for a variance or covariance shared by
# the term's effects, as lme4's sd_*|g.
.term_names <- function(P, cnm, grp) {
  vapply(seq_len(max(P)), function(k) {
    ab <- which(P == k & lower.tri(P, diag = TRUE), arr.ind = TRUE)
    is_var <- ab[1, 1] == ab[1, 2]
    if (nrow(ab) > 1) {
      paste0(if (is_var) "var_*|" else "cov_*|", grp)
    } else if (is_var) {
      paste0("var_", cnm[ab[1, 1]], "|", grp)
    } else {
      paste0("cov_", cnm[ab[1, 1]], ".", cnm[ab[1, 2]], "|", grp)
    }
  }, character(1))
}

#' Structure matrices of an lme4 fit
#'
#' @param lmerfit An \code{lmerMod} object from \code{lme4::lmer}.
#'
#' @return The \code{Hlist} of \code{.lmer_structure(lmerfit)}.
#'
#' @noRd
get_Hlist_lmer <- function(lmerfit)
{
  .lmer_structure(lmerfit)$Hlist
}

# Extract Y, X, Z from an lmer fit with the offset and prior weights applied.
# Premultiplying by W^(1/2) turns the weighted model, Var(E) = psi_r W^{-1},
# into the unit-variance model the likelihood code implements, with identical
# parameters (beta, Psi, psi_r); the offset enters the mean only. All
# psi-dependent quantities (scores, information, statistics, intervals) are
# invariant under the transformation; only the log-likelihood value
# changes, by the constant 0.5 * sum(log(w)). The diagonal scaling leaves
# Z's sparsity pattern unchanged.
.lmer_matrices <- function(lmerfit) {
  Y <- lme4::getME(lmerfit, "y")
  X <- lme4::getME(lmerfit, "X")
  Z <- lme4::getME(lmerfit, "Z")
  off <- lme4::getME(lmerfit, "offset")
  if (length(off) > 0 && any(off != 0)) Y <- Y - off
  w <- stats::weights(lmerfit)
  if (length(w) > 0 && any(w != 1)) {
    if (any(w <= 0)) stop("prior weights must be positive")
    sw <- sqrt(w)
    Y <- sw * Y
    X <- sw * X
    Z <- Matrix::Diagonal(x = sw) %*% Z
  }
  list(Y = Y, X = X, Z = Z)
}

# The vc_model of an lme4 fit, with lme4's estimates and parameter names;
# the only place the inference functions depend on lme4.
#' @rdname vc_model
#' @export
vc_model.lmerMod <- function(x, method = c("auto", "q_side", "n_side",
                                           "spectral"), ...) {
  m <- .lmer_matrices(x)  # offset and prior weights applied
  s <- .lmer_structure(x)
  .make_setup(Y = m$Y, X = m$X, Z = m$Z, Hlist = s$Hlist,
              REML = lme4::getME(x, "REML") != 0,
              psi_hat = s$psi_hat, names = s$names,
              method = match.arg(method))
}

#' Estimated covariance parameters of an lme4 fit
#'
#' @param lmerfit An \code{lmerMod} object from \code{lme4::lmer}.
#'
#' @return The \code{psi_hat} of \code{.lmer_structure(lmerfit)}: the
#'   \eqn{r} covariance parameters, the error variance last.
#'
#' @noRd
get_psi_hat_lmer <- function(lmerfit)
{
  .lmer_structure(lmerfit)$psi_hat
}

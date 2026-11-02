#' Covariance structure of an lme4 fit
#'
#' Structure matrices, estimates, and names for an lme4 fit, built from the
#' covariance structure each term declares, not from the estimates.
#'
#' @param lmerfit an \code{lmerMod} object.
#'
#' @return A list with components \code{Hlist}, the \eqn{r - 1} matrices
#'   \eqn{H_j} of \code{make_loglik}; \code{psi_hat}, the estimates, error
#'   variance last; and \code{names}. A term's parameters are in
#'   column-major lower-triangle order, that of
#'   \code{as.data.frame(VarCorr(lmerfit), order = "lower.tri")} for an
#'   unstructured term.
#'
#' @noRd
.lmer_structure <- function(lmerfit) {
  cnms <- lme4::getME(lmerfit, "cnms")
  Gp   <- lme4::getME(lmerfit, "Gp")
  fl   <- lme4::getME(lmerfit, "flist")
  nc   <- lengths(cnms)
  # Levels per term; attr "assign" maps terms to grouping factors
  nlev <- vapply(attr(fl, "assign"), function(i) nlevels(fl[[i]]), integer(1))
  # Group names made unique as in VarCorr, e.g. Subject.1
  grp  <- names(lme4::VarCorr(lmerfit))
  # Term t is columns Gp[t] + 1, ..., Gp[t + 1] of Z, nc[t] per level
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

  # Rows (i, j, k): entry (i, j) of Psi, upper triangle, is parameter k;
  # entries not listed are zero
  # Term t's parameters are numbered after those of terms 1, ..., t - 1
  shift <- cumsum(c(0L, vapply(Ps, max, integer(1))))
  map <- do.call(rbind, lapply(seq_along(Ps), function(t) {
    # Positions (a, b) of the term's parameters in its upper triangle
    ab   <- which(Ps[[t]] > 0 & upper.tri(Ps[[t]], diag = TRUE),
                  arr.ind = TRUE)
    # Offset of each level's block in Psi
    base <- Gp[t] + nc[t] * (seq_len(nlev[t]) - 1L)
    cbind(as.vector(outer(ab[, 1], base, `+`)),
          as.vector(outer(ab[, 2], base, `+`)),
          rep(Ps[[t]][ab] + shift[t], nlev[t]))
  }))
  # Rows of map for each parameter
  at <- unname(split(seq_len(nrow(map)),
                     factor(map[, 3], levels = seq_len(shift[length(shift)]))))
  # H_k has ones where Psi has parameter k
  Hlist <- lapply(at, function(s) {
    Matrix::forceSymmetric(Matrix::sparseMatrix(map[s, 1], map[s, 2], x = 1,
                                                dims = c(q, q)), uplo = "U")
  })

  # NaN in some lme4 2.0-6 fits whose optimizer stopped early
  s2 <- stats::sigma(lmerfit)^2
  if (!is.finite(s2)) {
    stop("lme4 reports a residual standard deviation of ",
         stats::sigma(lmerfit), ", so the fit may be incomplete; refitting ",
         "with another optimizer, such as lmerControl(optimizer = ",
         "\"bobyqa\"), may help", call. = FALSE)
  }
  # Fitted covariance matrix of the random effects
  Psi_fit <- s2 * Matrix::tcrossprod(lme4::getME(lmerfit, "Lambda"))
  # Estimates at each parameter's first position
  first   <- vapply(at, function(s) s[1], integer(1))
  psi_hat <- Psi_fit[map[first, 1:2, drop = FALSE]]
  .check_reproduces(Psi_fit, psi_hat, Hlist)
  list(Hlist = Hlist, psi_hat = c(psi_hat, s2),
       names = c(unlist(Map(.term_names, Ps, cnms, grp), use.names = FALSE),
                 "var_Residual"))
}

# Stop unless sum_j psi_hat[j] Hlist[[j]] reproduces Psi_fit up to rounding,
# entry (a, b) on the scale sqrt(d_a d_b) of its rounding error
.check_reproduces <- function(Psi_fit, psi_hat, Hlist) {
  d <- Matrix::diag(Psi_fit)
  if (!all(is.finite(d))) {
    stop("the fitted covariance matrix of the random effects is not finite",
         call. = FALSE)
  }
  # Psi_fit - sum_j psi_hat[j] H_j as triplets (i, j, x), indices 0-based
  E   <- methods::as(Psi_fit - Reduce(`+`, Map(`*`, psi_hat, Hlist)),
                     "TsparseMatrix")
  eps <- .Machine$double.eps
  # Floor eps * max(d) for variances that underflow to zero
  tol <- sqrt(eps) * (sqrt(d[E@i + 1L] * d[E@j + 1L]) + eps * max(d))
  if (any(abs(E@x) > tol)) .lmer_layout_error()
  invisible(TRUE)
}

# Error for a fit not laid out as .lmer_structure assumes
.lmer_layout_error <- function() {
  stop("the random effects of this lme4 fit are laid out differently than ",
       "reconf expects; please report this, with the lme4 version, at ",
       "https://github.com/koekvall/reconf/issues", call. = FALSE)
}

# Declared structure of each term: cls ("us", "diag", "cs", "ar1", or NA if
# unknown), hom (NA without a hom slot), and label for messages. Before
# lme4 2.0 every term is unstructured.
.lmer_declared <- function(lmerfit, nc) {
  if (!("getReCovs" %in% getNamespaceExports("lme4"))) {
    n <- length(nc)
    return(list(cls = rep("us", n), hom = rep(NA, n), label = rep("us", n)))
  }
  covs <- getExportedValue("lme4", "getReCovs")(lmerfit)
  slot_or_na <- function(cv, s) {
    if (methods::.hasSlot(cv, s)) methods::slot(cv, s) else NA
  }
  # Declared nc, where present, must match the fit
  nc_ok <- length(covs) == length(nc) &&
    all(vapply(seq_along(nc), function(t) {
      n <- slot_or_na(covs[[t]], "nc")
      identical(n, NA) || identical(as.integer(n), nc[[t]])
    }, logical(1)))
  if (!nc_ok) .lmer_layout_error()
  # NA unless hom is a single TRUE or FALSE
  hom <- vapply(covs, function(cv) {
    h <- slot_or_na(cv, "hom")
    if (is.logical(h) && length(h) == 1L) h else NA
  }, logical(1), USE.NAMES = FALSE)
  cl  <- vapply(covs, function(cv) class(cv)[1], character(1),
                USE.NAMES = FALSE)
  # Short name of each known class, NA for others
  cls <- unname(c(Covariance.us = "us", Covariance.diag = "diag",
                  Covariance.cs = "cs", Covariance.ar1 = "ar1")[cl])
  label <- ifelse(is.na(cls), paste0("'", cl, "'"),
                  paste0(cls, "(", ifelse(is.na(hom), "",
                                          paste0("hom = ", hom)), ")"))
  list(cls = cls, hom = hom, label = label)
}

# Pattern of one term's nc x nc covariance matrix: P[a, b] = k if entry
# (a, b) is parameter k, in column-major lower-triangle order, and 0 if zero
# by construction. NULL if the structure is unknown or not linear in the
# variances and covariances.
.term_pattern <- function(nc, cls, hom) {
  if (is.na(cls)) return(NULL)
  if (nc == 1L) return(matrix(1L))  # one variance
  if (cls != "us" && is.na(hom)) return(NULL)  # structured classes need hom
  # With two effects, ar1 is cs, and cs with separate variances is us
  if (nc == 2L && cls == "ar1") cls <- "cs"
  if (nc == 2L && cls == "cs" && !hom) cls <- "us"
  P  <- matrix(0L, nc, nc)
  lt <- lower.tri(P, diag = TRUE)
  if (cls == "us") {
    # Number the lower triangle, then mirror it
    P[lt] <- seq_len(sum(lt))
    P[upper.tri(P)] <- t(P)[upper.tri(P)]
  } else if (cls == "diag") {
    diag(P) <- if (hom) 1L else seq_len(nc)
  } else if (cls == "cs" && hom) {
    # One variance and one covariance
    P[] <- 2L
    diag(P) <- 1L
  } else {
    return(NULL)
  }
  P
}

# Parameter names of one term: lme4's sd_x|g and cor_y.x|g on the variance
# scale, var_x|g and cov_y.x|g, with * for a shared parameter
.term_names <- function(P, cnm, grp) {
  vapply(seq_len(max(P)), function(k) {
    # Positions (a, b) of parameter k in the lower triangle
    ab <- which(P == k & lower.tri(P, diag = TRUE), arr.ind = TRUE)
    is_var <- ab[1, 1] == ab[1, 2]
    if (nrow(ab) > 1) {  # shared by several entries
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

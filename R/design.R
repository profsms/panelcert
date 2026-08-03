# Shared infrastructure (spec section 2.1): design summary, union-find d_K,
# two-way demeaning. Mirrors the Julia reference engine numerically.

.integer_codes <- function(unit, time) {
  if (length(unit) != length(time)) stop("unit and time must have equal length")
  if (length(unit) == 0L) stop("empty panel")
  uid <- match(unit, unique(unit))
  tid <- match(time, unique(time))
  list(uid = uid, tid = tid, N = max(uid), T = max(tid))
}

#' Fixed-effect dimension of a two-way design
#'
#' Union-find on the bipartite unit-time graph:
#' `d_K = N + T - #connected components`.
#'
#' @param uid,tid integer codes 1..N / 1..T (one per observation)
#' @param N,T number of units / periods
#' @return list with `d_K` and `ncomponents`
#' @export
fe_dimension <- function(uid, tid, N, T) {
  parent <- seq_len(N + T)
  find <- function(a) {
    while (parent[a] != a) {
      parent[a] <<- parent[parent[a]]
      a <- parent[a]
    }
    a
  }
  for (k in seq_along(uid)) {
    ra <- find(uid[k]); rb <- find(N + tid[k])
    if (ra != rb) parent[ra] <- rb
  }
  ncomp <- length(unique(vapply(seq_len(N + T), find, integer(1))))
  list(d_K = N + T - ncomp, ncomponents = ncomp)
}

.twoway_demean_codes <- function(x, uid, tid, N, T, tol = 1e-10, maxit = 10000L) {
  w <- as.numeric(x)
  ucnt <- pmax(tabulate(uid, N), 1L)
  tcnt <- pmax(tabulate(tid, T), 1L)
  converged <- FALSE
  for (it in seq_len(maxit)) {
    um <- rowsum(w, uid)[, 1L] / ucnt
    w <- w - um[uid]
    tm <- rowsum(w, tid)[, 1L] / tcnt
    delta <- max(abs(tm))
    w <- w - tm[tid]
    if (delta < tol) { converged <- TRUE; break }
  }
  if (!converged) warning("two-way demeaning did not converge within maxit")
  w
}

.control_matrix <- function(controls, n) {
  if (is.null(controls)) return(matrix(numeric(0), nrow = n, ncol = 0L))
  Z <- if (is.null(dim(controls))) matrix(as.numeric(controls), ncol = 1L) else
    as.matrix(controls)
  storage.mode(Z) <- "double"
  if (nrow(Z) != n) stop("controls has ", nrow(Z), " rows; expected n = ", n)
  if (anyNA(Z) || any(!is.finite(Z)))
    stop("controls must contain only finite values")
  Z
}

.control_basis <- function(Z, tol = 1e-10) {
  n <- nrow(Z); k <- ncol(Z)
  Q <- matrix(0, n, k); r <- 0L
  if (!k) return(list(Q = Q, rank = r))
  for (j in seq_len(k)) {
    v <- Z[, j]
    colscale <- max(sqrt(sum(v^2)), 1)
    for (pass in 1:2) if (r) for (h in seq_len(r))
      v <- v - sum(Q[, h] * v) * Q[, h]
    nv <- sqrt(sum(v^2))
    if (nv > tol * colscale) {
      r <- r + 1L
      Q[, r] <- v / nv
    }
  }
  list(Q = if (r) Q[, seq_len(r), drop = FALSE] else Q[, FALSE, drop = FALSE],
       rank = r)
}

.partial_within_codes <- function(x, uid, tid, N, T, controls = NULL,
                                  tol = 1e-10, maxit = 10000L) {
  n <- length(uid)
  if (length(x) != n) stop("x must have length n = ", n)
  xf <- .twoway_demean_codes(x, uid, tid, N, T, tol = tol, maxit = maxit)
  Z <- .control_matrix(controls, n)
  if (!ncol(Z))
    return(list(xt = xf, Q = matrix(numeric(0), n, 0L), rank = 0L))
  Zf <- vapply(seq_len(ncol(Z)), function(j)
    .twoway_demean_codes(Z[, j], uid, tid, N, T, tol = tol, maxit = maxit),
    numeric(n))
  if (is.null(dim(Zf))) Zf <- matrix(Zf, ncol = 1L)
  basis <- .control_basis(Zf, tol = tol)
  xt <- xf
  if (basis$rank) for (j in seq_len(basis$rank))
    xt <- xt - sum(basis$Q[, j] * xt) * basis$Q[, j]
  list(xt = xt, Q = basis$Q, rank = basis$rank)
}

.partial_outcome_codes <- function(y, uid, tid, N, T, Q,
                                   tol = 1e-10, maxit = 10000L) {
  yf <- .twoway_demean_codes(y, uid, tid, N, T, tol = tol, maxit = maxit)
  if (ncol(Q)) for (j in seq_len(ncol(Q)))
    yf <- yf - sum(Q[, j] * yf) * Q[, j]
  yf
}

#' Two-way within transformation by alternating projections (unbalanced-safe)
#'
#' @param x numeric vector to demean
#' @param unit,time raw identifier vectors (any type)
#' @param tol,maxit convergence controls
#' @return the two-way-demeaned numeric vector `M x`
#' @export
twoway_demean <- function(x, unit, time, tol = 1e-10, maxit = 10000L) {
  cc <- .integer_codes(unit, time)
  .twoway_demean_codes(x, cc$uid, cc$tid, cc$N, cc$T, tol = tol, maxit = maxit)
}

#' Multiway within transformation by alternating projections
#'
#' Scalable residualization for applications with more than two categorical
#' fixed-effect dimensions, including Paper A's worker--firm--year design.
#'
#' @param x numeric vector to residualize.
#' @param fe_levels list of fixed-effect identifier vectors.
#' @param tol,maxit convergence controls.
#' @return the residualized numeric vector.
#' @export
multiway_demean <- function(x, fe_levels, tol = 1e-10, maxit = 10000L) {
  x <- as.numeric(x)
  n <- length(x)
  if (!is.list(fe_levels) || !length(fe_levels))
    stop("fe_levels must be a non-empty list of id vectors")
  if (!(tol > 0)) stop("tol must be positive")
  if (maxit < 1L) stop("maxit must be positive")
  codes <- lapply(seq_along(fe_levels), function(j) {
    ids <- fe_levels[[j]]
    if (length(ids) != n)
      stop(sprintf("fixed-effect dimension %d has length %d; expected %d",
                   j, length(ids), n))
    match(ids, unique(ids))
  })
  counts <- lapply(codes, function(z) pmax(tabulate(z, max(z)), 1L))
  w <- x
  converged <- FALSE
  for (it in seq_len(maxit)) {
    max_update <- 0
    for (j in seq_along(codes)) {
      z <- codes[[j]]
      m <- rowsum(w, z)[, 1L] / counts[[j]]
      max_update <- max(max_update, max(abs(m)))
      w <- w - m[z]
    }
    if (max_update < tol) { converged <- TRUE; break }
  }
  if (!converged) warning("multiway demeaning did not converge within maxit")
  w
}

.treatment_concentration <- function(xt) {
  Vn <- sum(xt^2)
  if (!(Vn > 0))
    return(list(Vn = as.numeric(Vn), lambda_n = NULL, n_eff = NULL))
  shares <- xt^2 / Vn
  list(Vn = as.numeric(Vn), lambda_n = max(shares),
       n_eff = 1 / sum(shares^2))
}

.design_summary_codes <- function(uid, tid, N, T, x = NULL, xt = NULL) {
  if (!is.null(x) && !is.null(xt)) stop("supply raw x or residualized xt, not both")
  n <- length(uid)
  fd <- fe_dimension(uid, tid, N, T)
  within <- if (!is.null(xt)) {
    if (length(xt) != n) stop("xt must have length n = ", n)
    as.numeric(xt)
  } else if (!is.null(x)) {
    if (length(x) != n) stop("x must have length n = ", n)
    .twoway_demean_codes(x, uid, tid, N, T)
  } else NULL
  if (is.null(within)) {
    return(structure(list(n = n, N = N, T = T, d_K = fd$d_K,
                          rho = fd$d_K / n, ncomponents = fd$ncomponents,
                          tau_star2 = NULL, lambda_n = NULL, n_eff = NULL),
                     class = "DesignSummary"))
  }
  conc <- .treatment_concentration(within)
  structure(list(n = n, N = N, T = T, d_K = fd$d_K,
                 rho = fd$d_K / n, ncomponents = fd$ncomponents,
                 tau_star2 = conc$Vn, lambda_n = conc$lambda_n,
                 n_eff = conc$n_eff), class = "DesignSummary")
}

.multiway_fe_dimension <- function(fe_levels) {
  n <- length(fe_levels[[1L]])
  codes <- lapply(fe_levels, function(ids) match(ids, unique(ids)))
  levels <- vapply(codes, max, integer(1))
  offsets <- c(0L, cumsum(levels))[seq_along(levels)]
  rows <- rep(seq_len(n), times = length(codes))
  cols <- unlist(Map(function(z, off) z + off, codes, offsets), use.names = FALSE)
  D <- Matrix::sparseMatrix(i = rows, j = cols, x = 1,
                            dims = c(n, sum(levels)))
  d_K <- as.integer(Matrix::rankMatrix(D, method = "qr"))

  parent <- seq_len(sum(levels))
  find <- function(a) {
    while (parent[a] != a) {
      parent[a] <<- parent[parent[a]]
      a <- parent[a]
    }
    a
  }
  union <- function(a, b) {
    ra <- find(a); rb <- find(b)
    if (ra != rb) parent[rb] <<- ra
  }
  for (i in seq_len(n)) {
    anchor <- offsets[1L] + codes[[1L]][i]
    if (length(codes) >= 2L)
      for (j in 2:length(codes)) union(anchor, offsets[j] + codes[[j]][i])
  }
  ncomp <- length(unique(vapply(seq_along(parent), find, integer(1))))
  list(d_K = d_K, ncomponents = ncomp, levels = levels)
}

#' Design-summary primitive (spec section 2.1)
#'
#' Pre-outcome design description reused by all four inference/diagnostic modules and
#' useful standalone: n, N, T, the fixed-effect dimension `d_K` via union-find,
#' `rho = d_K/n`, and (if a regressor is supplied) its within residual
#' variation `tau_star2 = x' M x`.
#'
#' @param unit,time raw identifier vectors
#' @param x optional regressor
#' @param controls optional numeric nuisance-covariate vector or matrix. These
#'   are partialled after the fixed effects; `d_K` and `rho` still describe the
#'   fixed-effect space alone.
#' @param fe_levels optional non-empty list of fixed-effect identifier vectors;
#'   when supplied, the sparse dummy-matrix rank is used for `d_K`
#' @return object of class \code{DesignSummary}
#' @examples
#' design_summary(rep(1:20, each = 10), rep(1:10, times = 20))
#' @export
design_summary <- function(unit = NULL, time = NULL, x = NULL,
                           fe_levels = NULL, controls = NULL) {
  if (!is.null(fe_levels)) {
    if (!is.null(controls))
      stop("controls are currently supported by the two-way design_summary method only")
    if (!is.list(fe_levels) || !length(fe_levels))
      stop("fe_levels must be a non-empty list of id vectors")
    n <- length(fe_levels[[1L]])
    if (!n) stop("empty panel")
    for (j in seq_along(fe_levels))
      if (length(fe_levels[[j]]) != n)
        stop(sprintf("fixed-effect dimension %d has length %d; expected %d",
                     j, length(fe_levels[[j]]), n))
    fd <- .multiway_fe_dimension(fe_levels)
    within <- if (is.null(x)) NULL else {
      if (length(x) != n) stop("x must have length n = ", n)
      multiway_demean(x, fe_levels)
    }
    conc <- if (is.null(within))
      list(Vn = NULL, lambda_n = NULL, n_eff = NULL)
    else .treatment_concentration(within)
    return(structure(list(n = n, N = fd$levels[1L],
                          T = if (length(fd$levels) >= 2L) fd$levels[2L] else 0L,
                          d_K = fd$d_K, rho = fd$d_K / n,
                          ncomponents = fd$ncomponents,
                          tau_star2 = conc$Vn, lambda_n = conc$lambda_n,
                          n_eff = conc$n_eff), class = "DesignSummary"))
  }
  if (is.null(unit) || is.null(time))
    stop("unit and time are required when fe_levels is not supplied")
  cc <- .integer_codes(unit, time)
  if (is.null(x)) {
    if (!is.null(controls)) stop("controls require x in design_summary")
    return(.design_summary_codes(cc$uid, cc$tid, cc$N, cc$T))
  }
  partial <- .partial_within_codes(x, cc$uid, cc$tid, cc$N, cc$T,
                                   controls = controls)
  .design_summary_codes(cc$uid, cc$tid, cc$N, cc$T, xt = partial$xt)
}

#' Print a design summary
#'
#' @param x a \code{DesignSummary}
#' @param ... unused
#' @return \code{x}, invisibly
#' @export
print.DesignSummary <- function(x, ...) {
  cat("Panel Design Summary\n")
  cat(sprintf("  n = %d obs | N = %d units | T = %d periods\n", x$n, x$N, x$T))
  cat(sprintf("  d_K = %d (connected components: %d) | rho = d_K/n = %.4f",
              x$d_K, x$ncomponents, x$rho))
  if (!is.null(x$tau_star2))
    cat(sprintf("\n  within variation tau*^2 = %.6g", x$tau_star2))
  if (!is.null(x$lambda_n))
    cat(sprintf(" | lambda_n = %.6g | N_eff = %.3f", x$lambda_n, x$n_eff))
  if (x$ncomponents > 1)
    cat(sprintf("\n  NOTE: design is disconnected (%d components); within comparisons exist only inside each component.",
                x$ncomponents))
  cat("\n")
  invisible(x)
}

# Moore-Penrose pseudo-inverse via SVD (base R; no MASS dependency).
.pinv <- function(A) {
  s <- svd(A)
  tol <- max(dim(A)) * .Machine$double.eps * max(s$d)
  pos <- s$d > tol
  s$v[, pos, drop = FALSE] %*% (t(s$u[, pos, drop = FALSE]) / s$d[pos])
}

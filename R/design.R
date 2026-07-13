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

#' Design-summary primitive (spec section 2.1)
#'
#' Pre-outcome design description reused by all three diagnostic modules and
#' useful standalone: n, N, T, the fixed-effect dimension `d_K` via union-find,
#' `rho = d_K/n`, and (if a regressor is supplied) its within residual
#' variation `tau_star2 = x' M x`.
#'
#' @param unit,time raw identifier vectors
#' @param x optional regressor
#' @return object of class \code{DesignSummary}
#' @examples
#' design_summary(rep(1:20, each = 10), rep(1:10, times = 20))
#' @export
design_summary <- function(unit, time, x = NULL) {
  cc <- .integer_codes(unit, time)
  n <- length(cc$uid)
  fd <- fe_dimension(cc$uid, cc$tid, cc$N, cc$T)
  tau_star2 <- NULL
  if (!is.null(x)) {
    if (length(x) != n) stop("x must have length n = ", n)
    xt <- .twoway_demean_codes(x, cc$uid, cc$tid, cc$N, cc$T)
    tau_star2 <- sum(xt^2)
  }
  structure(list(n = n, N = cc$N, T = cc$T, d_K = fd$d_K, rho = fd$d_K / n,
                 ncomponents = fd$ncomponents, tau_star2 = tau_star2),
            class = "DesignSummary")
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

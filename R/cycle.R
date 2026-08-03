# =============================================================================
# Paper A / JBES - exact inference under concentrated identifying variation
# SOURCE: "Exact Inference in Fixed-Effect Regressions with
# Concentrated Identifying Variation".
#
# The diffuse companion's Gaussian limits require lambda_n = max_i Xt_i^2/V_n
# to vanish. When it does not, the limiting null law of the t-statistic is not fixed across
# symmetric error laws of equal variance, so no fixed distribution-free
# critical value is uniformly valid over that class, and studentization does not
# repair it. Validity must come from adapting to the unknown error law: that is
# what the contrast construction here does, exactly, in finite samples.
#
# Mirrors julia/src/cycle.jl. Both rank candidate contrasts on the RAW treatment
# (v_c'x == v_c'xtilde exactly, since M_D v_c = v_c) and enumerate unit/period
# pairs in original-label order, so the packing is deterministic and identical
# across languages.
# =============================================================================

.binary_floor_info <- function(x, alpha) {
  xr <- as.numeric(x)
  if (anyNA(xr) || any(!is.finite(xr)))
    stop("x must contain only finite values")
  is_binary <- all(xr %in% c(0, 1))
  if (!is_binary)
    return(list(is_binary = FALSE, n_treated = NA_integer_,
                binary_floor_ok = TRUE, floor = NA_real_,
                reason = "not a binary treatment; the binary granularity floor does not apply."))

  n1 <- sum(xr == 1)
  n0 <- length(xr) - n1
  cap <- min(n1, n0)
  floor <- 2^(1 - cap)
  ok <- floor <= alpha
  relation <- if (ok) "<=" else ">"
  ending <- if (ok)
    "The structural granularity floor permits a two-sided level-alpha test."
  else "Not repairable by repacking."
  reason <- sprintf(
    "binary treatment, n1 = %d treated and n0 = %d untreated observations; at most min(n1,n0) = %d contrasts can carry loading; floor 2^(1-%d) = %.6g %s alpha = %.6g. %s",
    n1, n0, cap, cap, floor, relation, alpha, ending)
  list(is_binary = TRUE, n_treated = as.integer(n1),
       binary_floor_ok = ok, floor = floor, reason = reason)
}

.ap_demean <- function(x, a, b, tol = 1e-10, maxiter = 5000) {
  .twoway_demean_codes(x, a, b, max(a), max(b),
                       tol = tol, maxit = maxiter)
}

# greedy edge-disjoint packing: digons first, then DFS-extracted cycles
.greedy_packing <- function(a, b, xr) {
  m <- length(a)
  rows <- list(); signs <- list()
  key <- paste(a, b, sep = "\r")
  used <- rep(FALSE, m)
  for (g in split(seq_len(m), key)) {
    if (length(g) < 2) next
    g <- sort(g)
    k <- 1L
    while (k + 1L <= length(g)) {
      rows[[length(rows) + 1L]] <- c(g[k], g[k + 1L])
      signs[[length(signs) + 1L]] <- c(1, -1)
      used[g[k]] <- TRUE; used[g[k + 1L]] <- TRUE
      k <- k + 2L
    }
  }

  na <- max(a); nb <- max(b)
  adj <- vector("list", na + nb)
  for (i in seq_along(adj)) adj[[i]] <- integer(0)
  ends <- matrix(0L, nrow = m, ncol = 2)
  order_v <- integer(0); seen <- rep(FALSE, na + nb)
  for (e in seq_len(m)) {
    if (used[e]) next
    va <- a[e]; vb <- na + b[e]
    adj[[va]] <- c(adj[[va]], e); adj[[vb]] <- c(adj[[vb]], e)
    ends[e, ] <- c(va, vb)
    for (v in c(va, vb)) if (!seen[v]) { seen[v] <- TRUE; order_v <- c(order_v, v) }
  }
  other <- function(e, v) if (ends[e, 1] == v) ends[e, 2] else ends[e, 1]

  peel <- function(seed) {
    st <- seed
    while (length(st)) {
      v <- st[length(st)]; st <- st[-length(st)]
      if (length(adj[[v]]) == 1L) {
        e <- adj[[v]][1]; u <- other(e, v)
        adj[[v]] <<- setdiff(adj[[v]], e)
        adj[[u]] <<- setdiff(adj[[u]], e)
        if (length(adj[[u]]) <= 1L) st <- c(st, u)
      }
    }
  }
  peel(order_v)

  ptr <- 1L
  while (ptr <= length(order_v)) {
    start <- order_v[ptr]
    if (length(adj[[start]]) < 2L) { ptr <- ptr + 1L; next }
    parent_edge <- integer(na + nb); parent_edge[] <- NA_integer_
    parent_edge[start] <- 0L
    stack <- start; cycle <- integer(0)
    while (length(stack) && !length(cycle)) {
      v <- stack[length(stack)]
      advanced <- FALSE
      for (e in sort(adj[[v]])) {
        if (!is.na(parent_edge[v]) && parent_edge[v] != 0L && e == parent_edge[v]) next
        u <- other(e, v)
        if (!is.na(parent_edge[u])) {
          path <- e; w <- v
          while (w != u) { pe <- parent_edge[w]; path <- c(path, pe); w <- other(pe, w) }
          cycle <- path
          break
        }
        parent_edge[u] <- e
        stack <- c(stack, u)
        advanced <- TRUE
        break
      }
      if (!advanced && !length(cycle)) stack <- stack[-length(stack)]
    }
    if (!length(cycle)) { ptr <- ptr + 1L; next }
    L <- length(cycle)
    rows[[length(rows) + 1L]] <- cycle
    signs[[length(signs) + 1L]] <- ifelse(seq_len(L) %% 2 == 1, 1, -1)
    touched <- integer(0)
    for (e in cycle) {
      va <- ends[e, 1]; vb <- ends[e, 2]
      adj[[va]] <- setdiff(adj[[va]], e); adj[[vb]] <- setdiff(adj[[vb]], e)
      touched <- c(touched, va, vb)
    }
    peel(touched)
  }
  list(rows = rows, signs = signs)
}

.packing_capture <- function(pk, x) {
  if (!length(pk$rows)) return(0)
  sum(vapply(seq_along(pk$rows), function(i) {
    r <- pk$rows[[i]]
    (sum(pk$signs[[i]] * x[r]) / sqrt(length(r)))^2
  }, numeric(1)))
}

.greedy_multistart <- function(a, b, x) {
  m <- length(a)
  idx <- seq_len(m)
  prime <- 2147483647
  orders <- list(
    idx,
    order(a, x, b, method = "radix"),
    order((idx * 130363) %% prime, method = "radix"),
    order((idx * 13262647) %% prime, method = "radix")
  )
  best <- list(rows = list(), signs = list())
  best_capture <- -Inf
  for (ord in orders) {
    candidate <- .greedy_packing(a[ord], b[ord], x[ord])
    candidate$rows <- lapply(candidate$rows, function(r) ord[r])
    capture <- .packing_capture(candidate, x)
    if (capture > best_capture + 1e-14 * max(capture, best_capture, 1)) {
      best <- candidate
      best_capture <- capture
    }
  }
  best
}

# structured packing: unit-pair x period-pair four-cycles, greedy by squared
# loading, then a DFS pass on the leftovers
.structured_packing <- function(a, b, xr, units, periods,
                                reverse_ties = FALSE) {
  m <- length(a)
  cell <- new.env(hash = TRUE, parent = emptyenv())
  for (e in seq_len(m)) assign(paste(a[e], b[e], sep = "\r"), e, envir = cell)
  getcell <- function(u, p) {
    k <- paste(u, p, sep = "\r")
    if (exists(k, envir = cell, inherits = FALSE)) get(k, envir = cell) else 0L
  }
  vals <- numeric(0); quads <- list()
  for (i in seq_along(units)) {
    if (i == length(units)) break
    for (j in (i + 1L):length(units)) {
      f1 <- units[i]; f2 <- units[j]
      for (k in seq_along(periods)) {
        if (k == length(periods)) break
        for (l in (k + 1L):length(periods)) {
          c1 <- periods[k]; c2 <- periods[l]
          q1 <- getcell(f1, c1); if (!q1) next
          q2 <- getcell(f2, c1); if (!q2) next
          q3 <- getcell(f2, c2); if (!q3) next
          q4 <- getcell(f1, c2); if (!q4) next
          v <- (xr[q1] - xr[q2] + xr[q3] - xr[q4]) / 2
          vals <- c(vals, v * v)
          quads[[length(quads) + 1L]] <- c(q1, q2, q3, q4)
        }
      }
    }
  }
  rows <- list(); signs <- list(); used <- rep(FALSE, m)
  if (length(vals)) {
    tie <- if (reverse_ties) -seq_along(vals) else seq_along(vals)
    ord <- order(-vals, tie, method = "radix")
    for (idx in ord) {
      q <- quads[[idx]]
      if (any(used[q])) next
      rows[[length(rows) + 1L]] <- q
      signs[[length(signs) + 1L]] <- c(1, -1, 1, -1)
      used[q] <- TRUE
    }
  }
  left <- which(!used)
  if (length(left)) {
    sub <- .greedy_packing(a[left], b[left], xr[left])
    for (i in seq_along(sub$rows)) {
      rows[[length(rows) + 1L]] <- left[sub$rows[[i]]]
      signs[[length(signs) + 1L]] <- sub$signs[[i]]
    }
  }
  list(rows = rows, signs = signs)
}

# -----------------------------------------------------------------------------
# Sparse (mobility-network) structured packing
#
# The dense routine enumerates every unit-pair x period-pair four-cycle, which is
# O(N^2 T^2) and hopeless on a matched employer-employee network (35,807 workers
# against 5,301 firms is ~10^9 candidate pairs). Mobility networks have a
# different structure to exploit: most workers contribute exactly two edges, so
#   - a worker with two edges at the SAME firm is a digon, capture d^2/2;
#   - a worker with two edges at DIFFERENT firms is a "mover" carrying the
#     within-match difference d = x[e1] - x[e2] for that firm pair; two movers
#     over the same firm pair close a four-cycle with capture (d_j - d_i)^2/4.
# Within a firm pair the optimal edge-disjoint pairing is the EXTREME (nested)
# one -- pair the largest d with the smallest and recurse inward (Lemma 5.3).
# Everything left over goes to the DFS packing. O(m log m).
#
# Mirrors _sparse_packing in julia/src/cycle.jl: workers are visited in
# ascending code order and firm pairs in lexicographic order, so the two
# languages select the same supports.
# -----------------------------------------------------------------------------
.sparse_packing <- function(a, b, xr) {
  m <- length(a)
  rows <- list(); signs <- list()
  by_worker <- split(seq_len(m), a)          # names sort ascending below
  worker_keys <- names(by_worker)[order(as.integer(names(by_worker)))]

  f1v <- integer(0); f2v <- integer(0)
  dv <- numeric(0); e1v <- integer(0); e2v <- integer(0)
  others <- integer(0)

  for (w in worker_keys) {
    es <- sort(by_worker[[w]])
    if (length(es) != 2L) {
      others <- c(others, es)
      next
    }
    e1 <- es[1L]; e2 <- es[2L]
    f1 <- b[e1]; f2 <- b[e2]
    if (f1 == f2) {                          # stayer -> digon
      rows[[length(rows) + 1L]] <- c(e1, e2)
      signs[[length(signs) + 1L]] <- c(1, -1)
    } else {
      if (f1 > f2) { tmp <- f1; f1 <- f2; f2 <- tmp
                     tmp <- e1; e1 <- e2; e2 <- tmp }
      f1v <- c(f1v, f1); f2v <- c(f2v, f2)
      dv <- c(dv, xr[e1] - xr[e2]); e1v <- c(e1v, e1); e2v <- c(e2v, e2)
    }
  }

  leftover <- integer(0)
  if (length(f1v)) {
    # firm pairs in lexicographic order; movers within a pair by ascending d
    ord <- order(f1v, f2v, dv, method = "radix")
    f1v <- f1v[ord]; f2v <- f2v[ord]; dv <- dv[ord]
    e1v <- e1v[ord]; e2v <- e2v[ord]
    grp <- cumsum(c(TRUE, (f1v[-1] != f1v[-length(f1v)]) |
                          (f2v[-1] != f2v[-length(f2v)])))
    for (g in unique(grp)) {
      idx <- which(grp == g)                 # already sorted by d
      i <- 1L; j <- length(idx)
      while (i < j) {
        hi <- idx[j]; lo <- idx[i]
        # four-cycle from the extreme pair: +e1(hi) -e2(hi) +e2(lo) -e1(lo)
        rows[[length(rows) + 1L]] <- c(e1v[hi], e2v[hi], e2v[lo], e1v[lo])
        signs[[length(signs) + 1L]] <- c(1, -1, 1, -1)
        i <- i + 1L; j <- j - 1L
      }
      if (i == j) leftover <- c(leftover, e1v[idx[i]], e2v[idx[i]])
    }
  }

  left <- sort(unique(c(leftover, others)))
  if (length(left)) {
    sub <- .greedy_packing(a[left], b[left], xr[left])
    for (i in seq_along(sub$rows)) {
      rows[[length(rows) + 1L]] <- left[sub$rows[[i]]]
      signs[[length(signs) + 1L]] <- sub$signs[[i]]
    }
  }
  list(rows = rows, signs = signs)
}

.controlled_pairing <- function(pk, x, z) {
  C <- length(pk$rows)
  load <- nuisance <- numeric(C)
  for (i in seq_len(C)) {
    q <- pk$signs[[i]] / sqrt(length(pk$rows[[i]]))
    load[i] <- sum(q * x[pk$rows[[i]]])
    nuisance[i] <- sum(q * z[pk$rows[[i]]])
  }
  tol <- 64 * .Machine$double.eps * max(sqrt(sum(z^2)), 1)
  zero <- which(abs(nuisance) <= tol)
  nonzero <- which(abs(nuisance) > tol)
  rows <- pk$rows[zero]; signs <- pk$signs[zero]
  used <- rep(FALSE, C); used[zero] <- TRUE
  if (length(nonzero) >= 2L) {
    ij <- t(utils::combn(nonzero, 2L))
    den <- nuisance[ij[, 1L]]^2 + nuisance[ij[, 2L]]^2
    num <- nuisance[ij[, 2L]] * load[ij[, 1L]] -
      nuisance[ij[, 1L]] * load[ij[, 2L]]
    val <- num^2 / den
    ord <- order(-val, ij[, 1L], ij[, 2L], method = "radix")
    for (k in ord) {
      i <- ij[k, 1L]; j <- ij[k, 2L]
      if (used[i] || used[j]) next
      scale <- sqrt(nuisance[i]^2 + nuisance[j]^2)
      a1 <- nuisance[j] / scale; a2 <- -nuisance[i] / scale
      q1 <- pk$signs[[i]] / sqrt(length(pk$rows[[i]]))
      q2 <- pk$signs[[j]] / sqrt(length(pk$rows[[j]]))
      r <- c(pk$rows[[i]], pk$rows[[j]])
      q <- c(a1 * q1, a2 * q2)
      rows[[length(rows) + 1L]] <- r
      signs[[length(signs) + 1L]] <- q * sqrt(length(r))
      used[c(i, j)] <- TRUE
    }
  }
  list(rows = rows, signs = signs)
}

.controlled_rectangles <- function(a, b, x, z) {
  key <- paste(a, b, sep = "\r")
  if (anyDuplicated(key)) return(list(rows = list(), signs = list()))
  cell <- stats::setNames(seq_along(a), key)
  N <- max(a); T <- max(b)
  if (T < 3L) return(list(rows = list(), signs = list()))
  q1 <- c(.5, -.5, 0, -.5, .5, 0)
  q2 <- c(1, 1, -2, -1, -1, 2) / sqrt(12)
  vals <- numeric(0); candidates <- list(); serial <- 0L
  for (u1 in seq_len(N - 1L)) for (u2 in (u1 + 1L):N)
    for (t1 in seq_len(T - 2L)) for (t2 in (t1 + 1L):(T - 1L))
      for (t3 in (t2 + 1L):T) {
        keys <- paste(c(u1, u1, u1, u2, u2, u2),
                      c(t1, t2, t3, t1, t2, t3), sep = "\r")
        if (!all(keys %in% names(cell))) next
        r <- unname(cell[keys])
        l1 <- sum(q1 * x[r]); l2 <- sum(q2 * x[r])
        g1 <- sum(q1 * z[r]); g2 <- sum(q2 * z[r]); gss <- g1^2 + g2^2
        if (gss > 64 * .Machine$double.eps * max(sum(z[r]^2), 1)) {
          proj <- (g1 * l1 + g2 * l2) / gss
          aa1 <- l1 - proj * g1; aa2 <- l2 - proj * g2
        } else { aa1 <- l1; aa2 <- l2 }
        anorm <- sqrt(aa1^2 + aa2^2)
        if (anorm <= 64 * .Machine$double.eps * max(sqrt(sum(x[r]^2)), 1)) next
        q <- (aa1 * q1 + aa2 * q2) / anorm
        serial <- serial + 1L
        vals <- c(vals, anorm^2)
        candidates[[serial]] <- list(rows = r, signs = q * sqrt(6), serial = serial)
      }
  ord <- order(-vals, seq_along(vals), method = "radix")
  used <- rep(FALSE, length(a)); rows <- signs <- list()
  for (k in ord) {
    r <- candidates[[k]]$rows
    if (any(used[r])) next
    rows[[length(rows) + 1L]] <- r
    signs[[length(signs) + 1L]] <- candidates[[k]]$signs
    used[r] <- TRUE
  }
  list(rows = rows, signs = signs)
}

#' Nuisance-annihilating contrast system for a two-way design
#'
#' In a two-way design the annihilating contrasts are exactly the cycle space of
#' the observation multigraph. Each contrast eliminates both sets of fixed
#' effects identically -- exactly, not asymptotically -- which is what makes the
#' sign-flip test of [signflip_test()] exact in finite samples.
#'
#' @param x treatment vector.
#' @param unit,time fixed-effect identifiers.
#' @param method `"structured"` (all unit-pair four-cycles, greedy by squared
#'   loading, then a DFS pass on leftovers; for dense panels), `"sparse"` (stayer
#'   digons then firm-pair four-cycles with the optimal extreme pairing; for
#'   MOBILITY NETWORKS, where `"structured"` is `O(N^2 T^2)` and infeasible), or
#'   `"greedy"` (digons then DFS-extracted cycles, retaining the best of four
#'   fixed traversal orders; a deterministic lower bound on achievable capture).
#' @param controls optional numeric nuisance-covariate vector or matrix. Exact
#'   packing currently supports at most one linearly independent control. With
#'   one control, structured packing uses complete 2-by-3 rectangles and local
#'   cycle-space projection; other methods combine base cycles in pairs.
#' @return A list with `rows`, `signs`, `loadings`, `V_n`, `kappa`, `max_share`.
#' @export
cycle_contrasts <- function(x, unit, time,
                            method = c("structured", "sparse", "greedy"),
                            controls = NULL) {
  method <- match.arg(method)
  n <- length(x)
  if (length(unit) != n || length(time) != n)
    stop("x, unit and time must have equal length")
  if (anyNA(unit) || anyNA(time)) stop("unit and time identifiers cannot be missing")
  xr0 <- as.numeric(x)
  if (any(!is.finite(xr0))) stop("x must contain only finite values")

  # Canonical ordering makes the selected packing invariant to input row order.
  # Supports are mapped back to the caller's observation indices below.
  ord <- order(unit, time, xr0, method = "radix")
  unitc <- unit[ord]; timec <- time[ord]; xr <- xr0[ord]
  Zr <- .control_matrix(controls, n)[ord, , drop = FALSE]
  a <- match(unitc, unique(unitc))
  b <- match(timec, unique(timec))
  partial <- .partial_within_codes(xr, a, b, max(a), max(b), controls = Zr)
  if (partial$rank > 1L)
    stop("cycle packing currently supports at most one linearly independent continuous control")
  xt <- partial$xt
  V_n <- sum(xt^2)
  if (V_n <= 1e-12 * max(sum(xr^2), 1))
    stop("regressor has no within variation")

  if (method == "structured") {
    ulab <- unitc[!duplicated(a)][order(a[!duplicated(a)])]
    tlab <- timec[!duplicated(b)][order(b[!duplicated(b)])]
    pk <- .structured_packing(a, b, xr, order(ulab), order(tlab))
    pk2 <- .structured_packing(a, b, xr, order(ulab), order(tlab),
                               reverse_ties = TRUE)
    if (.packing_capture(pk2, xr) >
        .packing_capture(pk, xr) + 1e-14 * max(V_n, 1)) pk <- pk2
  } else if (method == "sparse") {
    pk <- .sparse_packing(a, b, xr)
  } else {
    pk <- .greedy_multistart(a, b, xr)
  }
  if (method != "greedy") {
    greedy <- .greedy_multistart(a, b, xr)
    if (.packing_capture(greedy, xr) >
        .packing_capture(pk, xr) + 1e-14 * max(V_n, 1)) pk <- greedy
  }
  if (partial$rank == 1L) {
    paired <- .controlled_pairing(pk, xr, partial$Q[, 1L])
    if (method == "structured") {
      rect <- .controlled_rectangles(a, b, xr, partial$Q[, 1L])
      if (.packing_capture(rect, xr) >
          .packing_capture(paired, xr) + 1e-14 * max(V_n, 1)) pk <- rect else pk <- paired
    } else pk <- paired
  }
  C <- length(pk$rows)
  loadings <- vapply(seq_len(C), function(c)
    sum(pk$signs[[c]] * xr[pk$rows[[c]]]) / sqrt(length(pk$rows[[c]])),
    numeric(1))
  ssq <- sum(loadings^2)
  mapped_rows <- lapply(pk$rows, function(r) ord[r])
  structure(list(rows = mapped_rows, signs = pk$signs, loadings = loadings,
                 V_n = V_n, kappa = ssq / V_n,
                 max_share = if (C == 0) 0 else max(loadings^2) / max(ssq, 1e-300)),
            class = "CycleSystem")
}

#' Validated design-only annihilating contrast system
#'
#' Construct caller-supplied contrasts for any number of categorical fixed
#' effects. Supports must be disjoint and every contrast must annihilate every
#' supplied fixed effect. This covers Paper A's worker--firm--year application,
#' which is not a two-way cycle packing.
#'
#' @param x treatment vector.
#' @param fe_levels list of fixed-effect identifier vectors.
#' @param rows list of observation-index vectors, one per support.
#' @param weights list of contrast-weight vectors matching `rows`.
#' @param blocks optional dependence-block ids; incompatible supports are rejected.
#' @param tol,maxit numerical controls.
#' @return a `CycleSystem` object usable by [signflip_test()] and
#'   [signflip_interval()].
#' @export
contrast_system <- function(x, fe_levels, rows, weights, blocks = NULL,
                            tol = 1e-10, maxit = 10000L) {
  x <- as.numeric(x); n <- length(x)
  if (!is.list(fe_levels) || !length(fe_levels))
    stop("fe_levels must be a non-empty list of id vectors")
  if (!is.list(rows) || !is.list(weights) || length(rows) != length(weights) ||
      !length(rows))
    stop("rows and weights must be non-empty lists of equal length")
  for (j in seq_along(fe_levels))
    if (length(fe_levels[[j]]) != n)
      stop(sprintf("fixed-effect dimension %d has length %d; expected %d",
                   j, length(fe_levels[[j]]), n))
  codes <- lapply(fe_levels, function(ids) match(ids, unique(ids)))
  used <- rep(FALSE, n)
  out_rows <- vector("list", length(rows)); out_signs <- vector("list", length(rows))
  for (i in seq_along(rows)) {
    r <- as.integer(rows[[i]]); w <- as.numeric(weights[[i]])
    if (!length(r)) stop("support ", i, " is empty")
    if (length(r) != length(w)) stop("support ", i, " has mismatched rows/weights")
    if (anyDuplicated(r)) stop("support ", i, " repeats an observation")
    if (any(r < 1L | r > n)) stop("support ", i, " has an index outside 1:n")
    if (any(used[r])) stop("contrast supports are not disjoint")
    used[r] <- TRUE
    nw <- sqrt(sum(w^2)); if (!(nw > 0)) stop("support ", i, " has zero-norm weights")
    q <- w / nw
    for (j in seq_along(codes)) {
      bal <- rowsum(q, codes[[j]][r])[, 1L]
      if (max(abs(bal)) > tol)
        stop(sprintf("support %d does not annihilate fixed-effect dimension %d", i, j))
    }
    out_rows[[i]] <- r
    # CycleSystem consumers divide this field by sqrt(support length).
    out_signs[[i]] <- q * sqrt(length(r))
  }
  xt <- multiway_demean(x, fe_levels, tol = tol, maxit = maxit)
  V_n <- sum(xt^2)
  if (V_n <= 1e-12 * max(sum(x^2), 1))
    stop("regressor has no variation after removing the supplied fixed effects")
  loadings <- vapply(seq_along(out_rows), function(i)
    sum(out_signs[[i]] * x[out_rows[[i]]]) / sqrt(length(out_rows[[i]])), numeric(1))
  ssq <- sum(loadings^2)
  cs <- structure(list(rows = out_rows, signs = out_signs, loadings = loadings,
                       V_n = V_n, kappa = ssq / V_n,
                       max_share = if (!length(loadings)) 0 else
                         max(loadings^2) / max(ssq, 1e-300)),
                  class = "CycleSystem")
  if (!is.null(blocks)) {
    comp <- support_compatibility(cs, blocks)
    if (!comp$compatible)
      stop(sprintf("contrast supports are not unions of complete dependence blocks; %d block(s) are partial or split",
                   comp$incompatible_blocks))
  }
  cs
}

.active_supports <- function(cs) {
  tol <- 64 * .Machine$double.eps * max(sqrt(cs$V_n), 1)
  which(abs(cs$loadings) > tol)
}

#' Check dependence-block compatibility of contrast supports
#'
#' Paper A Assumption ass:sym(iii): each treatment-loaded support must contain
#' all or none of every dependence block, and a block may not be split across
#' supports. Zero-loading supports are irrelevant to the test statistic.
#'
#' @param cs a `CycleSystem`.
#' @param blocks one block id per observation.
#' @return list with `compatible`, `incompatible_blocks`, `nblocks`, and
#'   `effective_C`.
#' @export
support_compatibility <- function(cs, blocks) {
  n <- length(blocks)
  assignment <- integer(n)
  active <- .active_supports(cs)
  for (i in active) {
    r <- cs$rows[[i]]
    if (any(r < 1L | r > n)) stop("blocks is shorter than a support index")
    if (any(assignment[r] != 0L)) stop("contrast supports are not disjoint")
    assignment[r] <- i
  }
  b <- match(blocks, unique(blocks))
  bad <- vapply(split(assignment, b), function(z) length(unique(z)) > 1L, logical(1))
  list(compatible = !any(bad), incompatible_blocks = sum(bad),
       nblocks = length(bad), effective_C = length(active))
}

#' Cycle capture ratio
#'
#' `kappa_C = sum_c (v_c'x)^2 / V_n`. In Paper A's diffuse comparison it is
#' the Pitman efficiency only when both power conditions (P1) and (P2) hold.
#' In the concentrated regime it is a capture diagnostic, and
#' `1/sqrt(kappa)` is only a capture-implied signal/SE ratio. It is computable
#' from the design before any outcome is examined.
#'
#' @inheritParams cycle_contrasts
#' @return The capture ratio.
#' @export
cycle_capture <- function(x, unit, time,
                          method = c("structured", "sparse", "greedy"),
                          controls = NULL)
  cycle_contrasts(x, unit, time, match.arg(method), controls = controls)$kappa

#' Exact sign-flip randomization test
#'
#' With `blocks = NULL`, exactness is stated for independent observation-level
#' symmetric errors. For clustered errors, pass one block id per observation:
#' every treatment-loaded support is then checked to be a union of complete
#' blocks. Conditional on that check, the joint score law is invariant under
#' independent sign flips; scores need not be exchangeable or identically
#' distributed. The identity sign pattern is adjoined.
#'
#' @inheritParams cycle_contrasts
#' @param y outcome vector.
#' @param beta0 null value.
#' @param nflips number of random sign patterns.
#' @param blocks optional design-defined dependence-block ids.
#' @return A list with `p`, `C`, `effective_C`, `kappa`, `beta_tilde`, and
#'   `full_enumeration_floor = 2^(1-effective_C)`. `min_pvalue` is retained as
#'   an alias. Global sign reversal duplicates the default absolute statistic.
#' @export
signflip_test <- function(y, x, unit = NULL, time = NULL, beta0 = 0,
                          nflips = 99999,
                          method = c("structured", "sparse", "greedy"),
                          controls = NULL,
                          blocks = NULL) {
  if (inherits(x, "CycleSystem")) {
    cs <- x
  } else {
    if (is.null(unit) || is.null(time)) stop("unit and time are required")
    if (length(y) != length(x)) stop("y, x, unit and time must have equal length")
    cs <- cycle_contrasts(x, unit, time, match.arg(method), controls = controls)
  }
  if (nflips < 1L) stop("nflips must be positive")
  z <- .contrast_terms(y, cs, blocks)
  p1 <- z$p1; p2 <- z$p2; A0 <- z$A0; B0 <- z$B0
  T0 <- abs(A0 - beta0 * B0)
  orb <- .draw_orbit(p1, p2, nflips)
  Ts <- abs(orb$A - beta0 * orb$B)
  floor <- 2^(1 - z$effective_C)
  list(p = (1 + sum(Ts >= T0 - 1e-12)) / (1 + nflips),
       C = length(cs$rows), effective_C = z$effective_C, kappa = cs$kappa,
       beta_tilde = A0 / B0, full_enumeration_floor = floor,
       min_pvalue = floor, statistic = T0)
}

.contrast_terms <- function(y, cs, blocks = NULL) {
  C <- length(cs$rows)
  if (!C) stop("no admissible contrasts")
  need <- max(unlist(cs$rows, use.names = FALSE))
  if (length(y) < need) stop("y is shorter than the largest support index")
  if (!is.null(blocks)) {
    if (length(blocks) != length(y)) stop("blocks must have the same length as y")
    comp <- support_compatibility(cs, blocks)
    if (!comp$compatible)
      stop(sprintf("contrast supports are not unions of complete dependence blocks; %d block(s) are partial or split",
                   comp$incompatible_blocks))
  }
  active <- .active_supports(cs)
  if (!length(active)) stop("no treatment-loaded contrasts: every loading is zero")
  yy <- as.numeric(y)
  Vy <- vapply(active, function(i)
    sum(cs$signs[[i]] * yy[cs$rows[[i]]]) / sqrt(length(cs$rows[[i]])), numeric(1))
  b <- cs$loadings[active]
  p1 <- b * Vy; p2 <- b^2; B0 <- sum(p2)
  if (!(B0 > 0)) stop("captured treatment variation is zero")
  list(p1 = p1, p2 = p2, A0 = sum(p1), B0 = B0,
       effective_C = length(active))
}

# Draw in bounded-memory chunks. The previous nflips-by-C matrix exceeded 2 GB
# for Paper A's 2,883-support worker--firm application.
.draw_orbit <- function(p1, p2, nflips, max_cells = 2e6) {
  C <- length(p1)
  chunk <- max(1L, floor(max_cells / C))
  A <- numeric(nflips); B <- numeric(nflips)
  first <- 1L
  while (first <= nflips) {
    m <- min(chunk, nflips - first + 1L)
    S <- matrix(sample(c(-1, 1), m * C, replace = TRUE), nrow = m)
    idx <- first:(first + m - 1L)
    A[idx] <- as.numeric(S %*% p1)
    B[idx] <- as.numeric(S %*% p2)
    first <- first + m
  }
  list(A = A, B = B)
}

#' Exact confidence set by test inversion
#'
#' The contrast scores are affine in `beta0`, so the whole randomization orbit is
#' affine too and the p-value curve over an entire grid costs no more than a
#' single test.
#'
#' @inheritParams signflip_test
#' @param alpha level.
#' @param ngrid grid resolution.
#' @param span half-width of the grid in randomization-scale units.
#' @return A list with `lo`, `hi`, `contiguous`, `grid_truncated`,
#'   `beta_tilde`, `C`, `effective_C`, `kappa`, `grid`, and `pvalue`. For the
#'   unstudentized statistic the acceptance set is an interval; `contiguous` is
#'   a numerical audit. Infinite endpoints are returned when the grid boundary
#'   is reached, rather than presenting a search edge as a confidence limit.
#' @export
signflip_interval <- function(y, x, unit = NULL, time = NULL, alpha = 0.05,
                              ngrid = 24001, nflips = 99999, span = 12,
                              method = c("structured", "sparse", "greedy"),
                              controls = NULL,
                              blocks = NULL) {
  if (!(alpha > 0 && alpha < 1)) stop("alpha must lie in (0, 1)")
  if (ngrid < 3L) stop("ngrid must be at least 3")
  if (nflips < 1L) stop("nflips must be positive")
  if (!(span > 0)) stop("span must be positive")
  if (inherits(x, "CycleSystem")) {
    cs <- x
  } else {
    if (is.null(unit) || is.null(time)) stop("unit and time are required")
    if (length(y) != length(x)) stop("y, x, unit and time must have equal length")
    cs <- cycle_contrasts(x, unit, time, match.arg(method), controls = controls)
  }
  z <- .contrast_terms(y, cs, blocks)
  p1 <- z$p1; p2 <- z$p2; A0 <- z$A0; B0 <- z$B0
  beta_tilde <- A0 / B0
  orb <- .draw_orbit(p1, p2, nflips); As <- orb$A; Bs <- orb$B
  scale <- sqrt(sum(p1^2)) / B0
  if (scale <= .Machine$double.eps * max(abs(beta_tilde), 1)) {
    floor <- 2^(1 - z$effective_C)
    if (floor > alpha)
      return(list(lo = -Inf, hi = Inf, contiguous = TRUE, grid_truncated = TRUE,
                  beta_tilde = beta_tilde, C = length(cs$rows),
                  effective_C = z$effective_C, kappa = cs$kappa,
                  grid = c(beta_tilde - 1, beta_tilde, beta_tilde + 1),
                  pvalue = rep(1, 3)))
    return(list(lo = beta_tilde, hi = beta_tilde, contiguous = TRUE,
                grid_truncated = FALSE, beta_tilde = beta_tilde,
                C = length(cs$rows), effective_C = z$effective_C,
                kappa = cs$kappa, grid = beta_tilde, pvalue = 1))
  }
  grid <- seq(beta_tilde - span * scale, beta_tilde + span * scale,
              length.out = ngrid)
  # For |Bs| < B0 the comparison holds between two equality roots. Sweep
  # those intervals over the sorted grid instead of materializing B x G.
  diff <- integer(ngrid + 1L)
  orbit_tol <- 64 * .Machine$double.eps * max(abs(A0), B0, 1)
  for (j in seq_len(nflips)) {
    if (abs(abs(Bs[j]) - B0) <= orbit_tol) {
      orientation <- if (Bs[j] >= 0) 1 else -1
      if (abs(As[j] - orientation * A0) <= orbit_tol) {
        diff[1L] <- diff[1L] + 1L
        diff[ngrid + 1L] <- diff[ngrid + 1L] - 1L
        next
      }
    }
    d1 <- B0 - Bs[j]; d2 <- B0 + Bs[j]
    if (abs(d1) <= orbit_tol || abs(d2) <= orbit_tol) {
      hit <- abs(As[j] - grid * Bs[j]) >= abs(A0 - grid * B0) - 1e-12
      idx <- which(hit)
      if (length(idx)) for (g in idx) {
        diff[g] <- diff[g] + 1L; diff[g + 1L] <- diff[g + 1L] - 1L
      }
      next
    }
    roots <- c((A0 - As[j]) / d1, (A0 + As[j]) / d2)
    lo <- min(roots); hi <- max(roots)
    left <- findInterval(lo, grid) + 1L
    if (left > 1L && grid[left - 1L] >= lo) left <- left - 1L
    right <- findInterval(hi, grid)
    if (left <= right && left <= ngrid && right >= 1L) {
      left <- max(left, 1L); right <- min(right, ngrid)
      diff[left] <- diff[left] + 1L
      diff[right + 1L] <- diff[right + 1L] - 1L
    }
  }
  pv <- (1 + cumsum(diff[seq_len(ngrid)])) / (1 + nflips)
  acc <- which(pv > alpha)
  if (!length(acc))
    return(list(lo = NA_real_, hi = NA_real_, contiguous = TRUE,
                grid_truncated = FALSE, beta_tilde = beta_tilde,
                C = length(cs$rows), effective_C = z$effective_C, kappa = cs$kappa,
                grid = grid, pvalue = pv))
  first <- min(acc); last <- max(acc)
  left <- first == 1L; right <- last == length(grid)
  lo <- if (left) -Inf else grid[first]
  hi <- if (right) Inf else grid[last]
  list(lo = lo, hi = hi,
       contiguous = all(pv[first:last] > alpha),
       grid_truncated = left || right, beta_tilde = beta_tilde,
       C = length(cs$rows), effective_C = z$effective_C,
       kappa = cs$kappa, grid = grid, pvalue = pv)
}

#' Paper A adequacy report
#'
#' Reports the concentration of the identifying variation, the capture the
#' contrast system achieves, the price of exactness, and the exact confidence
#' set. The verdict answers: is the Gaussian approximation underlying
#' conventional inference trustworthy on this design?
#'
#' With `blocks = NULL`, concentration is measured across observations. With
#' dependence blocks, treatment mass is aggregated as `sum(xt_i^2)` within
#' block and realized score concentration uses squared block scores
#' `(sum(xt_i * u_i))^2`, as in Paper A's worker--firm application.
#'
#' `POINT_PASS` is a descriptive pass below the finite `lambda_max` heuristic,
#' not a theorem-level certificate; `FLAGGED` marks a concentration warning;
#' and `INCONCLUSIVE` means the two-sided full-enumeration floor
#' `2^(1-effective_C)` exceeds `alpha`.
#'
#' @inheritParams signflip_interval
#' @param delta size tolerance (reported, not used in the verdict).
#' @param interval whether to compute the exact confidence set.
#' @param lambda_max concentration warning convention; Paper A's condition is
#'   the sequence statement `lambda_n -> 0`, not a finite cutoff.
#' @return An `AdequacyReport`.
#' @export
cycle_report <- function(y, x, unit, time, alpha = 0.05, delta = 0.05,
                         method = c("structured", "sparse", "greedy"), nflips = 99999,
                         interval = TRUE, lambda_max = 0.10, controls = NULL,
                         blocks = NULL) {
  if (!(alpha > 0 && alpha < 1)) stop("alpha must lie in (0, 1)")
  if (!(delta > 0 && delta < 1 - alpha)) stop("delta must lie in (0, 1-alpha)")
  if (!(lambda_max > 0)) stop("lambda_max must be positive")
  if (nflips < 1L) stop("nflips must be positive")
  method <- match.arg(method)
  n <- length(y)
  if (length(x) != n || length(unit) != n || length(time) != n)
    stop("y, x, unit and time must have equal length")
  a <- as.integer(factor(unit, levels = unique(unit)))
  b <- as.integer(factor(time, levels = unique(time)))
  fd <- fe_dimension(a, b, max(a), max(b))
  cs <- cycle_contrasts(x, unit, time, method, controls = controls)
  C <- length(cs$rows)
  Ceff <- length(.active_supports(cs))
  if (!is.null(blocks)) {
    if (length(blocks) != n) stop("blocks must have the same length as the data")
    comp <- support_compatibility(cs, blocks)
    if (!comp$compatible)
      stop(sprintf("contrast supports are not unions of complete dependence blocks; %d block(s) are partial or split",
                   comp$incompatible_blocks))
  }
  partial <- .partial_within_codes(x, a, b, max(a), max(b), controls = controls)
  xt <- partial$xt
  block_id <- if (is.null(blocks)) NULL else match(blocks, unique(blocks))
  treatment_mass <- if (is.null(block_id)) xt^2 else
    as.numeric(rowsum(xt^2, block_id, reorder = FALSE))
  treatment_share <- treatment_mass / cs$V_n
  lambda_n <- max(treatment_share)
  H_n <- sum(treatment_share^2); n_eff <- 1 / H_n
  yt <- .partial_outcome_codes(y, a, b, max(a), max(b), partial$Q)
  beta_ols <- sum(xt * yt) / cs$V_n
  u <- yt - beta_ols * xt
  score_contribution <- if (is.null(block_id)) xt * u else
    as.numeric(rowsum(xt * u, block_id, reorder = FALSE))
  score_mass <- sum(score_contribution^2)
  if (score_mass > 0) {
    score_share <- score_contribution^2 / score_mass
    score_lambda_n <- max(score_share)
    score_H_n <- sum(score_share^2); score_n_eff <- 1 / score_H_n
  } else {
    score_lambda_n <- score_H_n <- score_n_eff <- NaN
  }
  cyc_dim <- n - (max(a) + max(b)) + fd$ncomponents - partial$rank

  itv <- if (interval && Ceff > 0)
    signflip_interval(y, cs, alpha = alpha, nflips = nflips,
                      blocks = blocks) else NULL
  min_p <- if (Ceff > 0) 2^(1 - Ceff) else 1
  binary <- .binary_floor_info(x, alpha)

  notes <- character(0)
  if (partial$rank > 0L) notes <- c(notes, sprintf(
    "every support contrast annihilates both fixed effects and %d linearly independent continuous nuisance control(s); no global outcome residualization is used in the exact test.",
    partial$rank))
  notes <- c(notes, sprintf(
    "capture kappa_C = %.4f over C = %d supports (%d treatment-loaded); the capture-implied signal/SE ratio is 1/sqrt(kappa) = %.3fx. Kappa equals Pitman efficiency only when both (P1) and diffuse-design condition (P2) hold.",
    cs$kappa, C, Ceff, 1 / sqrt(cs$kappa)))
  notes <- c(notes, sprintf(
    "cycle-space dimension = %d; the system uses %d supports (%.1f%%). Weighted edge-disjoint cycle packing is NP-hard, so the packing is a heuristic and kappa_C is a LOWER BOUND. Where a discrete treatment makes many four-cycle loadings exactly tied, this implementation ranks on the raw treatment, uses canonical label ordering, and compares fixed traversal orders. The selected supports are therefore deterministic, row-order invariant, and identical across languages.",
    cyc_dim, C, 100 * C / max(cyc_dim, 1)))
  concentrated <- lambda_n > lambda_max
  if (concentrated) {
    notes <- c(notes, sprintf(
      "CONCENTRATION WARNING: lambda_n = %.3f (N_eff = %.1f) exceeds the heuristic cutoff %.3f. Along persistently concentrated sequences the limiting t law is not fixed across symmetric error laws; at the fully concentrated boundary Paper A rules out a fixed distribution-free critical value, not every adaptive or bootstrap procedure.",
      lambda_n, n_eff, lambda_max))
  } else {
    notes <- c(notes, sprintf(
      "DESCRIPTIVE DESIGN PASS: lambda_n = %.4f (N_eff = %.1f) is below the heuristic cutoff %.3f. This is compatible with the diffuse regime but is not a finite-sample certificate; the theorem requires lambda_n -> 0 along a sequence.",
      lambda_n, n_eff, lambda_max))
  }
  if (is.finite(score_lambda_n))
    notes <- c(notes, sprintf(
      "realized %s score diagnostic: lambda_score = %.4f and N_eff,score = %.1f. This one-outcome residual diagnostic is a warning statistic, not by itself a consistent population concentration estimate.",
      if (is.null(block_id)) "observation" else "dependence-block",
      score_lambda_n, score_n_eff))
  if (cs$max_share > 0.5)
    notes <- c(notes, sprintf(
      "(P1) WARNING: one support carries %.1f%% of captured variation. Exactness is unaffected, but the Gaussian power formula and kappa efficiency interpretation are not licensed when such dominance persists.",
      100 * cs$max_share))
  notes <- c(notes, sprintf(
    "default two-sided full-enumeration floor is 2^(1-C_eff) = %.2g: global sign reversal duplicates every absolute-statistic orbit value, so level %.3g requires C_eff >= %d.",
    min_p, alpha, 1 + ceiling(log2(1 / alpha))))
  if (is.null(blocks))
    notes <- c(notes, "exactness is evaluated for singleton observation blocks. For clustered errors, pass blocks=... so compatibility is checked rather than assumed.")
  else
    notes <- c(notes, sprintf("block compatibility verified for %d design-defined dependence blocks.",
                              length(unique(blocks))))
  if (!is.null(itv) && !itv$contiguous)
    notes <- c(notes, "NUMERICAL WARNING: inversion was non-contiguous even though the unstudentized statistic has an interval acceptance set.")
  if (!is.null(itv) && itv$grid_truncated)
    notes <- c(notes, "the accepted set reached the grid boundary; the infinite endpoint is conservative and is not a finite search edge masquerading as a confidence limit.")

  verdict <- if (min_p > alpha) "INCONCLUSIVE" else
    if (concentrated) "FLAGGED" else "POINT_PASS"
  reason <- if (verdict != "INCONCLUSIVE") "" else
    if (binary$is_binary && !binary$binary_floor_ok) binary$reason else
      sprintf("selected packing has effective_C = %d; floor 2^(1-%d) = %.6g > alpha = %.6g. A stronger compatible packing may repair this.",
              Ceff, Ceff, min_p, alpha)
  if (nzchar(reason)) notes <- c(notes, reason)

  design <- .design_summary_codes(a, b, max(a), max(b), xt = xt)
  statistic <- list(kappa = cs$kappa, C = C, effective_C = Ceff,
                    cycle_dim = cyc_dim,
                    concentration_level = if (is.null(block_id)) "observation" else "block",
                    concentration_blocks = if (is.null(block_id)) n else max(block_id),
                    lambda_n = lambda_n, H_n = H_n,
                    n_eff = n_eff, score_lambda_n = score_lambda_n,
                    score_H_n = score_H_n, score_n_eff = score_n_eff,
                    lambda_score = score_lambda_n, H_score = score_H_n,
                    n_eff_score = score_n_eff,
                    beta_ols = beta_ols,
                    control_rank = partial$rank,
                    max_share = cs$max_share, se_price = 1 / sqrt(cs$kappa),
                    full_enumeration_floor = min_p, min_pvalue = min_p,
                    method = method,
                    n_treated = binary$n_treated,
                    binary_floor_ok = binary$binary_floor_ok, reason = reason,
                    beta_tilde = if (is.null(itv)) NULL else itv$beta_tilde,
                    ci_lo = if (is.null(itv)) NULL else itv$lo,
                    ci_hi = if (is.null(itv)) NULL else itv$hi,
                    ci_level = 1 - alpha,
                    ci_contiguous = if (is.null(itv)) NULL else itv$contiguous,
                    ci_grid_truncated = if (is.null(itv)) NULL else itv$grid_truncated)
  .new_AdequacyReport("cycle_inference", design, statistic, NULL, NULL, NULL,
                      NULL, verdict, alpha, delta, notes)
}

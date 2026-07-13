# Reference data now ships INSIDE the package: the four user-facing datasets
# are lazy-loaded objects (data/), and the parity-only V-Dem gate-1 vintage is
# a test fixture (tests/testthat/testdata/). The harness is therefore
# self-contained and runs under R CMD check.

# Load a bundled dataset by name (robust under both installed and load_all).
get_dataset <- function(name) {
  e <- new.env()
  utils::data(list = name, package = "panelcert", envir = e)
  e[[name]]
}

# The gate-1 V-Dem vintage (spec 7.1 design primitive) is a parity fixture,
# not a user-facing dataset -- read it from the test-fixture directory.
read_gate1 <- function() {
  utils::read.csv(testthat::test_path("testdata", "vdem_gate1.csv"),
                  stringsAsFactors = FALSE)
}

# Shim so existing read_panel("castle_panel.csv") call sites keep working:
# map the old CSV name to the bundled dataset object.
read_panel <- function(csv) {
  get_dataset(sub("_panel\\.csv$", "", csv))
}

# Dense two-way FE dummy matrix and projection, for verification on small panels.
dense_dummies <- function(uid, tid, N, T) {
  n <- length(uid)
  D <- matrix(0, n, N + T)
  D[cbind(seq_len(n), uid)] <- 1
  D[cbind(seq_len(n), N + tid)] <- 1
  D
}

pinv_ref <- function(A) {
  s <- svd(A)
  tol <- max(dim(A)) * .Machine$double.eps * max(s$d)
  pos <- s$d > tol
  s$v[, pos, drop = FALSE] %*% (t(s$u[, pos, drop = FALSE]) / s$d[pos])
}

PD <- asNamespace("panelcert")

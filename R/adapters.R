# Model-ingestion adapters (spec sections 0 and 6): the package CONSUMES an
# already-estimated model \u2014 fixest::feols, plm::plm, or stats::lm \u2014 and never
# asks the user to re-specify it. Each adapter extracts (y, x, unit, time)
# from the fitted object (respecting any observations the estimator dropped)
# and dispatches to the default method; a test asserts that the extracted
# regression reproduces the fitted coefficient exactly.

# ---- fixest ------------------------------------------------------------------

.frame_from_fixest <- function(model) {
  if (!requireNamespace("fixest", quietly = TRUE))
    stop("the fixest package is required to ingest fixest models")
  fe_vars <- model$fixef_vars
  if (length(fe_vars) != 2)
    stop("expected exactly two fixed effects (unit and time); model has ",
         length(fe_vars), " \u2014 the diagnostics are defined for two-way FE designs")
  y <- as.numeric(stats::model.matrix(model, type = "lhs"))
  X <- stats::model.matrix(model, type = "rhs")
  xcols <- setdiff(colnames(X), "(Intercept)")
  if (length(xcols) != 1)
    stop("expected exactly one (scalar) regressor of interest; model has: ",
         paste(xcols, collapse = ", "))
  FE <- stats::model.matrix(model, type = "fixef")
  list(y = y, x = as.numeric(X[, xcols]),
       unit = FE[[fe_vars[1]]], time = FE[[fe_vars[2]]],
       beta = unname(stats::coef(model)[xcols]))
}

#' @rdname eiv_adequacy
#' @export
eiv_adequacy.fixest <- function(object, ...) {
  fr <- .frame_from_fixest(object)
  .check_adapter_beta(fr, "fixest")
  eiv_adequacy.default(fr$y, fr$x, fr$unit, fr$time, ...)
}

#' @rdname leverage_report
#' @export
leverage_report.fixest <- function(object, ...) {
  fr <- .frame_from_fixest(object)
  .check_adapter_beta(fr, "fixest")
  leverage_report.default(fr$y, fr$x, fr$unit, fr$time, ...)
}

#' @rdname twfe_adequacy
#' @export
twfe_adequacy.fixest <- function(object, first_treat, ...) {
  fr <- .frame_from_fixest(object)
  if (length(first_treat) != length(fr$y))
    stop("first_treat must match the model's estimation sample (n = ",
         length(fr$y), ")")
  twfe_adequacy.default(fr$y, fr$unit, fr$time, first_treat, ...)
}

# ---- plm ---------------------------------------------------------------------

.frame_from_plm <- function(model) {
  if (!requireNamespace("plm", quietly = TRUE))
    stop("the plm package is required to ingest plm models")
  if (!(model$args$model == "within" && model$args$effect == "twoways"))
    warning("plm model is not a two-way within model (model = '",
            model$args$model, "', effect = '", model$args$effect,
            "'); the diagnostic recomputes the two-way FE design from the ",
            "panel index, which may not match the fitted specification")
  idx <- plm::index(model)
  mf <- stats::model.frame(model)
  y <- as.numeric(mf[[1]])
  xvars <- attr(stats::terms(model), "term.labels")
  if (length(xvars) != 1)
    stop("expected exactly one (scalar) regressor of interest; model has: ",
         paste(xvars, collapse = ", "))
  list(y = y, x = as.numeric(mf[[xvars]]),
       unit = as.character(idx[[1]]), time = as.character(idx[[2]]),
       beta = unname(stats::coef(model)[xvars]))
}

#' @rdname eiv_adequacy
#' @export
eiv_adequacy.plm <- function(object, ...) {
  fr <- .frame_from_plm(object)
  .check_adapter_beta(fr, "plm")
  eiv_adequacy.default(fr$y, fr$x, fr$unit, fr$time, ...)
}

#' @rdname leverage_report
#' @export
leverage_report.plm <- function(object, ...) {
  fr <- .frame_from_plm(object)
  .check_adapter_beta(fr, "plm")
  leverage_report.default(fr$y, fr$x, fr$unit, fr$time, ...)
}

#' @rdname twfe_adequacy
#' @export
twfe_adequacy.plm <- function(object, first_treat, ...) {
  fr <- .frame_from_plm(object)
  if (length(first_treat) != length(fr$y))
    stop("first_treat must match the model's estimation sample (n = ",
         length(fr$y), ")")
  twfe_adequacy.default(fr$y, fr$unit, fr$time, first_treat, ...)
}

# ---- lm ----------------------------------------------------------------------

.frame_from_lm <- function(model, x, unit, time) {
  mf <- stats::model.frame(model)
  for (v in c(x, unit, time))
    if (!v %in% names(mf))
      stop("'", v, "' is not a variable of the fitted model; model.frame has: ",
           paste(names(mf), collapse = ", "))
  list(y = as.numeric(stats::model.response(mf)), x = as.numeric(mf[[x]]),
       unit = as.character(mf[[unit]]), time = as.character(mf[[time]]),
       beta = unname(stats::coef(model)[x]))
}

#' @rdname eiv_adequacy
#' @export
eiv_adequacy.lm <- function(object, x, unit, time, ...) {
  fr <- .frame_from_lm(object, x, unit, time)
  .check_adapter_beta(fr, "lm")
  eiv_adequacy.default(fr$y, fr$x, fr$unit, fr$time, ...)
}

#' @rdname leverage_report
#' @export
leverage_report.lm <- function(object, x, unit, time, ...) {
  fr <- .frame_from_lm(object, x, unit, time)
  .check_adapter_beta(fr, "lm")
  leverage_report.default(fr$y, fr$x, fr$unit, fr$time, ...)
}

# The extracted regression must reproduce the fitted coefficient exactly \u2014
# otherwise the adapter is diagnosing a different model than the user ran.
.check_adapter_beta <- function(fr, label) {
  if (is.na(fr$beta)) return(invisible(NULL))
  cc <- .integer_codes(fr$unit, fr$time)
  xt <- .twoway_demean_codes(fr$x, cc$uid, cc$tid, cc$N, cc$T)
  yt <- .twoway_demean_codes(fr$y, cc$uid, cc$tid, cc$N, cc$T)
  beta <- sum(xt * yt) / sum(xt^2)
  if (abs(beta - fr$beta) > 1e-6 * max(abs(fr$beta), 1))
    warning("the two-way FE reproduction (beta = ", format(beta),
            ") does not match the fitted ", label, " coefficient (",
            format(fr$beta), "); the fitted model may include additional ",
            "controls or a different FE structure \u2014 the diagnostic applies ",
            "to the two-way design recomputed from the ids")
  invisible(NULL)
}

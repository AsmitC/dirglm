#' Creates a summary of a fitted \code{dirglm} object
#' @param object A fitted \code{dirglm} object
#' @param prob value between 0 and 1 indicating the desired probability
#' level for the credible intervals. Defaults to 0.95.
#' @param robust Logical indicating whether to show robust statistics for
#' \code{Estimate} and \code{Est.Error}. Defaults to \code{FALSE}. See
#' details for more information.
#'
#' @return A an object of class \code{summary.dirglm}. A list
#' containing with two elements, \code{beta} and \code{f0},
#' each containing a data frame with columns \code{Estimate},
#' \code{Est.Error}, \code{l-95% Ci}, \code{u-95% Ci}
#'
#' @details If \code{robust} is set to \code{TRUE}, the median and
#' median absolute deviation (MAD) of the posterior samples
#' are shown for \code{Estimate} and \code{Est.Error}.
#' Otherwise, the mean and the standard deviations are used.
#'
#' @method summary dirglm
#' @export
summary.dirglm <- function(object, prob=0.95, robust=FALSE) {

  fit <- object
  if (!inherits(fit, "dirglm")) stop("fit must be an object of class 'dirglm'")
  if (prob <= 0 || prob >= 1)   stop("prob must be a value between 0 and 1")

  # Get samples
  beta <- as.matrix(fit$samples$beta)

  if (is.null(colnames(beta))) {
    colnames(beta) <- paste0("beta_", seq_len(ncol(beta)) - 1L)
    colnames(beta)[1] <- "Intercept"
  }
  f0   <- as.matrix(fit$samples$f0)
  p <- ncol(beta); s <- ncol(f0)

  pct <- prob * 100
  pct_str <- if (pct == floor(pct)) as.integer(pct) else pct
  l_name <- paste0("l-", pct_str, "% CrI")
  u_name <- paste0("u-", pct_str, "% CrI")

  # Helper
  make_tbl <- function(mat) {
    qs <- c((1 - prob)/2, 1 - (1 - prob)/2)
    stats <- t(vapply(seq_len(ncol(mat)), function(j) {
      col <- mat[, j]
      if (robust) {
        est     <- median(col)
        est.err <- mad(col)
      } else {
        est     <- mean(col)
        est.err <- sd(col)
      }
      c(
        Estimate    = est,
        `Est.Error` = est.err,
        quantile(col, qs[1], names = FALSE),
        quantile(col, qs[2], names = FALSE)
      )
    }, numeric(4)))
    colnames(stats) <- c("Estimate", "Est.Error", l_name, u_name)
    df <- as.data.frame(stats, check.names = FALSE)
    rownames(df) <- colnames(mat)
    df
  }

  # Make output tables
  beta_tbl <- make_tbl(beta)
  f0_tbl   <- make_tbl(f0)

  # Arrange metadata
  lf.name <- if (!is.null(fit$link$name)) fit$link$name else "custom"
  form    <- paste(deparse(fit$formula), collapse = "")
  nobs    <- nrow(fit$data)
  burnin  <- fit$burnin
  thin    <- fit$thin
  save    <- fit$save
  iter    <- burnin + thin * save

  # Output
  meta <- list(
    Link       = lf.name,
    Formula    = form,
    nobs       = nobs,
    burnin     = burnin,
    thin       = thin,
    save       = save,
    iter       = iter,
    beta_acceptance = fit$beta_acceptance,
    f0_acceptance   = fit$f0_acceptance
  )

  out <- list(meta = meta, beta = beta_tbl, f0 = f0_tbl)
  class(out) <- "summary.dirglm"
  out
}

#' Print a summary for a \code{summary.dirglm} object
#'
#' @param x \code{summary.dirglm} object
#' @param digits The number of digits to print. Defaults to 3.
#' @param ... Additional arguments passed to \code{print}.
#'
#' @seealso \code{\link{summary.dirglm}}
#'
#' @method print summary.dirglm
#' @export
print.summary.dirglm <- function(x, digits = 3, ...) {
  if (!inherits(x, "summary.dirglm")) {
    stop("`x` must be a summary.dirglm object")
  }

  cat("Summary of DirGLM Fit\n\n")

  # Metadata
  m <- x$meta
  cat(sprintf("Formula: %s\n",                         m$Formula))
  cat(sprintf("Link:    %s\n",                         m$Link))
  cat(sprintf("Draws:   %d (burnin = %d, thin = %d, save = %d)\n",
              m$iter, m$burnin, m$thin, m$save))
  cat(sprintf("Number of observations: %d\n", m$nobs))
  cat(sprintf("Acceptance Ratio: beta = %.*f, f0 = %.*f\n\n",
              digits, m$beta_acceptance, digits, m$f0_acceptance))

  # Helper
  print_table <- function(df) {
    mat <- as.matrix(df)
    fmt <- paste0("%.", digits, "f")
    for (j in seq_len(ncol(mat))) {
      mat[, j] <- sprintf(fmt, as.numeric(mat[, j]))
    }
    colnames(mat) <- colnames(df)
    rownames(mat) <- rownames(df)
    print(mat, quote = FALSE, right = TRUE, ...)
  }

  # Coefficients
  cat("Posterior summary for beta:\n")
  print_table(x$beta)
  cat("\nPosterior summary for f0:\n")
  print_table(x$f0)

  invisible(x)
}

#' Creates a summary of a fitted \code{dpglm} object
#'
#' @param object A fitted \code{dpglm} object.
#' @param prob Value between 0 and 1 indicating the desired probability
#' level for the credible intervals. Defaults to 0.95.
#' @param robust Logical indicating whether to show robust statistics for
#' \code{Estimate} and \code{Est.Error}. Defaults to \code{FALSE}.
#'
#' @return An object of class \code{summary.dpglm}. A list containing metadata,
#' posterior summaries for \code{beta}, and CRM diagnostics.
#'
#' @details
#' If \code{robust} is \code{TRUE}, posterior medians and MADs are shown for
#' \code{Estimate} and \code{Est.Error}. Otherwise, posterior means and standard
#' deviations are used. Since the DPGLM reference distribution is represented by
#' a CRM whose atom locations and jump sizes may vary across iterations, the CRM
#' summary reports scalar diagnostics rather than a fixed-support table.
#'
#' @method summary dpglm
#' @export
summary.dpglm <- function(object, prob = 0.95, robust = FALSE) {

  fit <- object
  if (!inherits(fit, "dpglm")) stop("object must be an object of class 'dpglm'")
  if (prob <= 0 || prob >= 1) stop("prob must be a value between 0 and 1")

  beta <- as.matrix(fit$samples$beta)
  if (is.null(colnames(beta))) {
    colnames(beta) <- paste0("beta_", seq_len(ncol(beta)) - 1L)
    colnames(beta)[1] <- "Intercept"
  }

  pct <- prob * 100
  pct_str <- if (pct == floor(pct)) as.integer(pct) else pct
  l_name <- paste0("l-", pct_str, "% CrI")
  u_name <- paste0("u-", pct_str, "% CrI")
  qs <- c((1 - prob) / 2, 1 - (1 - prob) / 2)

  make_tbl <- function(mat) {
    stats <- t(vapply(seq_len(ncol(mat)), function(j) {
      col <- mat[, j]
      if (robust) {
        est <- stats::median(col, na.rm = TRUE)
        est.err <- stats::mad(col, na.rm = TRUE)
      } else {
        est <- mean(col, na.rm = TRUE)
        est.err <- stats::sd(col, na.rm = TRUE)
      }

      c(
        Estimate = est,
        `Est.Error` = est.err,
        stats::quantile(col, qs[1], names = FALSE, na.rm = TRUE),
        stats::quantile(col, qs[2], names = FALSE, na.rm = TRUE)
      )
    }, numeric(4)))

    colnames(stats) <- c("Estimate", "Est.Error", l_name, u_name)
    df <- as.data.frame(stats, check.names = FALSE)
    rownames(df) <- colnames(mat)
    df
  }

  beta_tbl <- make_tbl(beta)

  crm <- fit$samples$crm
  z_list <- crm$z.tld
  J_list <- crm$J.tld

  crm_functionals <- t(vapply(seq_along(z_list), function(i) {
    z <- z_list[[i]]
    J <- J_list[[i]]
    J <- J / sum(J)

    ord <- order(z)
    z <- z[ord]
    J <- J[ord]
    F <- cumsum(J)

    qfun <- function(p) z[which(F >= p)[1]]

    c(
      min = min(z),
      q25 = qfun(0.25),
      median = qfun(0.50),
      q75 = qfun(0.75),
      max = max(z)
    )
  }, numeric(5)))

  crm_tbl <- make_tbl(as.matrix(crm_functionals))

  lf.name <- if (!is.null(fit$link$name)) fit$link$name else "custom"
  form <- paste(deparse(fit$formula), collapse = "")
  nobs <- nrow(fit$data)
  burnin <- fit$burnin
  thin <- fit$thin
  save <- fit$save
  iter <- burnin + thin * save

  meta <- list(
    Link = lf.name,
    Formula = form,
    nobs = nobs,
    burnin = burnin,
    thin = thin,
    save = save,
    iter = iter,
    beta_acceptance = fit$beta_acceptance,
    crm_acceptance = fit$crm_acceptance
  )

  out <- list(
    meta = meta,
    beta = beta_tbl,
    crm = crm_tbl
  )

  class(out) <- "summary.dpglm"
  out
}

#' Print a summary for a \code{summary.dpglm} object
#'
#' @param x \code{summary.dpglm} object.
#' @param digits The number of digits to print. Defaults to 3.
#' @param ... Additional arguments passed to \code{print}.
#'
#' @method print summary.dpglm
#' @export
print.summary.dpglm <- function(x, digits = 3, ...) {
  if (!inherits(x, "summary.dpglm")) {
    stop("`x` must be a summary.dpglm object")
  }

  cat("Summary of DPGLM Fit\n\n")

  m <- x$meta
  cat(sprintf("Formula: %s\n", m$Formula))
  cat(sprintf("Link:    %s\n", m$Link))
  cat(sprintf("Draws:   %d (burnin = %d, thin = %d, save = %d)\n",
              m$iter, m$burnin, m$thin, m$save))
  cat(sprintf("Number of observations: %d\n", m$nobs))
  cat(sprintf("Acceptance Ratio: beta = %.*f, CRM = %.*f\n\n",
              digits, m$beta_acceptance, digits, m$crm_acceptance))

  print_table <- function(df) {
    mat <- as.matrix(df)
    fmt <- paste0("%.", digits, "f")
    for (j in seq_len(ncol(mat))) {
      mat[, j] <- sprintf(fmt, as.numeric(mat[, j]))
    }
    colnames(mat) <- colnames(df)
    rownames(mat) <- rownames(df)
    print(mat, quote = FALSE, right = TRUE, ...)
  }

  cat("Posterior summary for beta:\n")
  print_table(x$beta)

  cat("\nPosterior summary for CRM:\n")
  print_table(x$crm)

  invisible(x)
}

#' Creates a summary of a fitted \code{cdpglm} object
#'
#' @param object A fitted \code{cdpglm} object.
#' @param prob Value between 0 and 1 indicating the desired probability
#' level for the credible intervals. Defaults to 0.95.
#' @param robust Logical indicating whether to show robust statistics for
#' \code{Estimate} and \code{Est.Error}. Defaults to \code{FALSE}.
#' @param ... Additional arguments, currently unused.
#'
#' @return An object of class \code{summary.cdpglm}.
#'
#' @method summary cdpglm
#' @export
summary.cdpglm <- function(object, prob = 0.95, robust = FALSE, ...) {

  fit <- object
  if (!inherits(fit, "cdpglm")) stop("object must be an object of class 'cdpglm'")
  if (prob <= 0 || prob >= 1) stop("prob must be a value between 0 and 1")

  pct <- prob * 100
  pct_str <- if (pct == floor(pct)) as.integer(pct) else pct
  l_name <- paste0("l-", pct_str, "% CrI")
  u_name <- paste0("u-", pct_str, "% CrI")
  qs <- c((1 - prob) / 2, 1 - (1 - prob) / 2)

  make_tbl <- function(mat) {
    mat <- as.matrix(mat)

    stats <- t(vapply(seq_len(ncol(mat)), function(j) {
      col <- mat[, j]

      if (robust) {
        est <- stats::median(col, na.rm = TRUE)
        est.err <- stats::mad(col, na.rm = TRUE)
      } else {
        est <- mean(col, na.rm = TRUE)
        est.err <- stats::sd(col, na.rm = TRUE)
      }

      c(
        Estimate = est,
        `Est.Error` = est.err,
        stats::quantile(col, qs[1], names = FALSE, na.rm = TRUE),
        stats::quantile(col, qs[2], names = FALSE, na.rm = TRUE)
      )
    }, numeric(4)))

    colnames(stats) <- c("Estimate", "Est.Error", l_name, u_name)
    df <- as.data.frame(stats, check.names = FALSE)
    rownames(df) <- colnames(mat)
    df
  }

  beta <- as.matrix(fit$samples$beta)
  if (is.null(colnames(beta))) {
    colnames(beta) <- paste0("beta_", seq_len(ncol(beta)) - 1L)
  }

  beta_tbl <- make_tbl(beta)

  rho <- as.matrix(fit$samples$rho)
  colnames(rho) <- "rho"
  rho_tbl <- make_tbl(rho)

  crm <- fit$samples$crm
  z_list <- crm$z.tld
  J_list <- crm$J.tld

  reference_functionals <- t(vapply(seq_along(z_list), function(i) {
    z <- z_list[[i]]
    J <- J_list[[i]]
    J <- J / sum(J)

    ord <- order(z)
    z <- z[ord]
    J <- J[ord]
    F <- cumsum(J)

    qfun <- function(p) z[which(F >= p)[1]]

    c(
      min = min(z),
      q25 = qfun(0.25),
      median = qfun(0.50),
      q75 = qfun(0.75),
      max = max(z)
    )
  }, numeric(5)))

  reference_tbl <- make_tbl(as.matrix(reference_functionals))

  lf.name <- if (!is.null(fit$link$name)) fit$link$name else "custom"
  form <- paste(deparse(fit$formula), collapse = "")
  nobs <- nrow(fit$data)
  burnin <- fit$burnin
  thin <- fit$thin
  save <- fit$save
  iter <- burnin + thin * save

  meta <- list(
    Link = lf.name,
    Formula = form,
    nobs = nobs,
    burnin = burnin,
    thin = thin,
    save = save,
    iter = iter,
    beta_acceptance = fit$beta_acceptance,
    crm_acceptance = fit$crm_acceptance,
    rho_acceptance = fit$rho_acceptance,
    corr = fit$corr
  )

  out <- list(
    meta = meta,
    beta = beta_tbl,
    rho = rho_tbl,
    reference = reference_tbl
  )

  class(out) <- "summary.cdpglm"
  out
}

#' Print a summary for a \code{summary.cdpglm} object
#'
#' @param x \code{summary.cdpglm} object.
#' @param digits The number of digits to print. Defaults to 3.
#' @param ... Additional arguments passed to \code{print}.
#'
#' @method print summary.cdpglm
#' @export
print.summary.cdpglm <- function(x, digits = 3, ...) {
  if (!inherits(x, "summary.cdpglm")) {
    stop("`x` must be a summary.cdpglm object")
  }

  cat("Summary of CDPGLM Fit\n\n")

  m <- x$meta
  cat(sprintf("Formula: %s\n", m$Formula))
  cat(sprintf("Link:    %s\n", m$Link))
  cat(sprintf("Copula correlation: %s\n", m$corr))
  cat(sprintf("Draws:   %d (burnin = %d, thin = %d, save = %d)\n",
              m$iter, m$burnin, m$thin, m$save))
  cat(sprintf("Number of observations: %d\n", m$nobs))
  cat(sprintf("Acceptance Ratio: beta = %.*f, CRM = %.*f, rho = %.*f\n\n",
              digits, m$beta_acceptance,
              digits, m$crm_acceptance,
              digits, m$rho_acceptance))

  print_table <- function(df) {
    mat <- as.matrix(df)
    fmt <- paste0("%.", digits, "f")
    for (j in seq_len(ncol(mat))) {
      mat[, j] <- sprintf(fmt, as.numeric(mat[, j]))
    }
    colnames(mat) <- colnames(df)
    rownames(mat) <- rownames(df)
    print(mat, quote = FALSE, right = TRUE, ...)
  }

  cat("Posterior summary for beta:\n")
  print_table(x$beta)

  cat("\nPosterior summary for rho:\n")
  print_table(x$rho)

  cat("\nPosterior summary for reference distribution functionals:\n")
  print_table(x$reference)

  invisible(x)
}

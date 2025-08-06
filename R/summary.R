#' Creates a summary of a fitted \code{dirglm} object
#' @param fit A fitted \code{dirglm} object
#' @param prob value between 0 and 1 indicating the desired probability
#' level for the credible intervals. Defaults to 0.95.
#'
#' @return A an object of class \code{summary.dirglm}. A list
#' containing with two elements, \code{beta} and \code{f0},
#' each containing a data frame with columns \code{Estimate},
#' \code{Est.Error}, \code{l-95% Ci}, \code{u-95% Ci}
#'
#' @method summary dirglm
#' @export
summary.dirglm <- function(fit, prob = 0.95) {
  if (!inherits(fit, "dirglm")) stop("fit must be an object of class 'dirglm'")
  if (prob <= 0 || prob >= 1)   stop("prob must be a value between 0 and 1")

  # Get samples
  beta <- fit$samples$beta
  f0   <- fit$samples$f0
  p <- ncol(beta); s <- ncol(f0)
  colnames(beta) <- paste0("beta_", seq_len(p))
  colnames(f0)   <- paste0("f0_",   seq_len(s))

  # Build output tables
  pct <- prob * 100
  pct_str <- if (pct == floor(pct)) as.integer(pct) else pct
  l_name <- paste0("l-", pct_str, "% CI")
  u_name <- paste0("u-", pct_str, "% CI")

  # Helper
  make_tbl <- function(mat) {
    qs <- c((1 - prob)/2, 1 - (1 - prob)/2)
    stats <- t(vapply(seq_len(ncol(mat)), function(j) {
      col <- mat[, j]
      c(
        Estimate  = mean(col),
        `Est.Error` = sd(col),
        quantile(col, qs[1], names = FALSE),
        quantile(col, qs[2], names = FALSE)
      )
    }, numeric(4)))
    # now assign the dynamic names
    colnames(stats) <- c("Estimate", "Est.Error", l_name, u_name)
    df <- as.data.frame(stats, check.names = FALSE)
    rownames(df) <- colnames(mat)
    df
  }

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

  meta <- list(
    Link       = lf.name,
    Formula    = form,
    nobs       = nobs,
    burnin     = burnin,
    thin       = thin,
    save       = save,
    iter       = iter,
    p_acc_beta = fit$p_acc_beta,
    p_acc_f0   = fit$p_acc_f0
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
              digits, m$p_acc_beta, digits, m$p_acc_f0))

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

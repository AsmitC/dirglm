#' Plot DirGLM Posterior Distributions
#'
#' @param fit A fitted \code{dirglm} object.
#' @param what Specifies which posterior to plot. One of \code{"beta"}
#' or \code{"f0"}. Defaults to \code{"beta"}.
#' @param pars Specifies the names or integer indices of
#' the \code{beta} coefficients to plot. See details for more information.
#' @param top_k The number of \code{beta} coefficients to plot if \code{pars} is \code{NULL}
#' and the number of coefficients is large. Defaults to 6.
#' @param order Specifies the sorting mechanism for \code{beta} coefficients.
#' One of \code{"magnitude"}, \code{"variance"}, or \code{"none"}. Defaults to \code{"magnitude"}.
#' @param robust Logical indicating whether to plot robust statistics for
#' \code{Estimate} and \code{Est.Error}. Defaults to \code{FALSE}.
#' @param prob Value between 0 and 1 indicating the desired probability
#' level for the credible intervals. Defaults to 0.95.
#' @param dens_args A list of additional arguments passed to the
#' \code{stats::density} function for plotting the posterior density of \code{beta}.
#' @param common_xlim Logical indicating whether to share the x-axis limits
#' across all \code{beta} panels. Defaults to \code{FALSE}.
#' @param ncol Number of columns in the layout for \code{beta} panels.
#' @param col A string specifying the desired plot color.
#' @param band_alpha Specifies the transparency of the CrI band. Defaults to 0.25.
#' @param line_lwd Line width for the posterior density lines. Defaults to 2.
#' @param point_pch Point style for the \code{f0} plots.
#' @param point_cex Point size for the \code{f0} plots.
#' @param xlab Label for the x-axis.
#' @param ylab Label for the y-axis.
#' @param main Plot title.
#' @param xlim Limits for the x-axis.
#' @param ylim Limits for the y-axis.
#' @param ... Additional arguments passed to the \code{plot} function.
#'
#' @return An invisible list containing the selected coefficients or the
#' f0 support, depending on the value of \code{what}.
#'
#' @details
#' If \code{pars} is \code{NULL} and the number of \code{beta} coefficients is greater
#' than \code{top_k}, the top \code{top_k} coefficients ordered by \code{order}
#' are plotted. The default for \code{top_k} is 6 and the default sorting mechanism
#' is the coefficient magnitude. If the number of coefficients is less than \code{top_k},
#' all of them are plotted.
#'
#' @export
plot_dirglm <- function(
    fit,
    what = c("beta", "f0"),
    pars = NULL,
    top_k = 6,
    order = c("magnitude", "variance", "none"),
    robust = FALSE,
    prob = 0.95,
    # beta density options
    dens_args = list(adjust = 1, n = 512, na.rm = TRUE),
    common_xlim = FALSE,
    ncol = NULL,
    # aesthetics
    col = par("col"),
    band_alpha = 0.25,
    line_lwd = 2,
    point_pch = 16,
    point_cex = 0.8,
    xlab = NULL, ylab = NULL, main = NULL,
    xlim = NULL, ylim = NULL,
    ...
) {
  what  <- match.arg(what)
  order <- match.arg(order)

  band_col <- grDevices::adjustcolor(col, alpha.f = band_alpha) # CrI shade

  # Get summary statistics
  s <- summary(fit, prob = prob, robust = robust)
  low_key <- sprintf("l-%.0f%% CrI", 100 * prob)
  upp_key <- sprintf("u-%.0f%% CrI", 100 * prob)

  # Beta plotting (default)
  if (what == "beta") {
    beta_samples <- as.matrix(fit$samples$beta)
    btab <- as.data.frame(s$beta)

    center_all <- setNames(btab[,"Estimate"], rownames(btab))
    lower_all  <- setNames(btab[, low_key],   rownames(btab))
    upper_all  <- setNames(btab[, upp_key],   rownames(btab))

    # Choose indices to plot
    all_names <- rownames(btab)
    if (is.null(pars)) {
      if (order == "variance") ord  <- order(btab[, "Est.Error"], decreasing = TRUE)
      else if (order == "none") ord <- seq_along(all_names)
      else ord <- order(abs(center_all), decreasing = TRUE)
      keep_names <- all_names[ord][seq_len(min(top_k, length(ord)))]
    } else if (is.character(pars)) keep_names <- pars
    else keep_names <- all_names[pars]

    idx <- match(keep_names, colnames(beta_samples))
    k <- length(idx)
    if (is.null(ncol)) ncol <- ceiling(sqrt(k))
    nrow <- ceiling(k / ncol)

    if (common_xlim && is.null(xlim)) {
      xr <- range(beta_samples[, idx, drop = FALSE], finite = TRUE)
      xlim <- grDevices::extendrange(xr, f = 0.04)
    }

    if (is.null(xlab)) xlab <- "Value"
    if (is.null(ylab)) ylab <- "Density"
    if (is.null(main)) {
      main <- sprintf("Posterior %s (%.0f%% CrI shaded)",
                      if (robust) "Medians" else "Means", 100 * prob)
    }

    op <- par(no.readonly = TRUE); on.exit(par(op))
    par(mfrow = c(nrow, ncol), oma = c(0, 0, if (k > 1) 2.2 else 0, 0))

    # Plot each beta[j]
    for (j in seq_len(k)) {
      nm <- keep_names[j]
      xj <- beta_samples[, idx[j]]
      dj <- do.call(stats::density, c(list(x = xj), dens_args))
      cj <- center_all[nm]
      lj <- lower_all[nm]
      uj <- upper_all[nm]

      # CrI[j]
      L <- min(lj, uj); U <- max(lj, uj)
      xlim_j <- if (!is.null(xlim)) xlim else range(dj$x, finite = TRUE)
      ylim_j <- if (!is.null(ylim)) ylim else range(dj$y, finite = TRUE)

      plot(dj, type = "n",
           xlab = xlab, ylab = ylab,
           main = if (k == 1) main else "",
           xlim = xlim_j, ylim = ylim_j, ...)

      inside <- dj$x >= L & dj$x <= U
      if (any(inside)) {
        px <- c(dj$x[inside], rev(dj$x[inside]))
        py <- c(dj$y[inside], rep(0, sum(inside)))
        polygon(px, py, border = NA, col = band_col)
      } else warning("No density values inside the CrI range for '", nm, "'")

      lines(dj$x, dj$y, lwd = line_lwd, col = col)
      abline(v = cj, lty = 2, col = col)       # center reference
      mtext(nm, side = 3, line = 0.2, cex = 0.9)
    }

    if (k > 1 && nzchar(main)) mtext(main, outer = TRUE, line = 0.6, cex = 1.1)
    return(invisible(list(what="beta", selected=keep_names)))
  }

  # f0 plotting
  if (what == "f0") {
    x <- fit$spt
    ftab <- as.data.frame(s$f0)
    center <- ftab[,"Estimate"]
    lower  <- ftab[, low_key]
    upper  <- ftab[, upp_key]

    if (is.null(xlab)) xlab <- "Support"
    if (is.null(ylab)) ylab <- "Probability"
    if (is.null(main)) main <- sprintf("Posterior %s of f0 (%.0f%% CrI)",
                                       if (robust) "Median" else "Mean", 100 * prob)

    cap <- 1.2
    if (is.null(ylim)) {
      ymax <- min(1, cap * max(upper, na.rm=TRUE))
      ylim <- c(0, grDevices::extendrange(c(0, ymax), f = 0.04)[2])
    }
    if (is.null(xlim)) xlim <- range(x, finite = TRUE)

    plot(NA, xlim = xlim, ylim = ylim, xlab = xlab, ylab = ylab, main = main, ...)
    xs <- c(x, rev(x)); ys <- c(upper, rev(lower))
    polygon(xs, ys, border = NA, col = band_col)
    lines(x, center, lwd = line_lwd, col = col)
    points(x, center, pch = point_pch, cex = point_cex, col = col)

    return(invisible(list(what="f0")))
  }
}

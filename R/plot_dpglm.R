#' Plot DPGLM Posterior Distributions
#'
#' @param object A fitted \code{dpglm} object.
#' @param what Specifies which posterior to plot. One of \code{"beta"}
#' or \code{"crm"}. Defaults to \code{"beta"}.
#' @param pars Specifies the names or integer indices of
#' the \code{beta} coefficients to plot. See details for more information.
#' @param top_k The number of \code{beta} coefficients to plot if \code{pars}
#' is \code{NULL} and the number of coefficients is large. Defaults to 6.
#' @param order Specifies the sorting mechanism for \code{beta} coefficients.
#' One of \code{"magnitude"}, \code{"variance"}, or \code{"none"}.
#' @param robust Logical indicating whether to plot robust statistics for
#' \code{Estimate} and \code{Est.Error}. Defaults to \code{FALSE}.
#' @param prob Value between 0 and 1 indicating the desired probability
#' level for the credible intervals. Defaults to 0.95.
#' @param dens_args A list of additional arguments passed to
#' \code{stats::density}.
#' @param common_xlim Logical indicating whether to share x-axis limits across
#' beta panels. Defaults to \code{FALSE}.
#' @param ncol Number of columns in the layout for beta panels.
#' @param grid_length Number of grid points used when plotting the posterior CRM.
#' @param col Plot color.
#' @param band_alpha Transparency of the credible interval band.
#' @param line_lwd Line width.
#' @param xlab Label for the x-axis.
#' @param ylab Label for the y-axis.
#' @param main Plot title.
#' @param xlim Limits for the x-axis.
#' @param ylim Limits for the y-axis.
#' @param ... Additional arguments passed to plotting functions.
#'
#' @return An invisible list containing information about the plotted object.
#'
#' @details
#' If \code{pars = NULL} and the number of \code{beta} coefficients is greater
#' than \code{top_k}, the top \code{top_k} coefficients ordered by \code{order}
#' are plotted. The default for \code{top_k} is 6 and the default sorting mechanism
#' is the coefficient magnitude. If the number of coefficients is less than \code{top_k},
#' all of them are plotted.
#'
#' @export
plot_dpglm <- function(
    object,
    what = c("beta", "crm"),
    pars = NULL,
    top_k = 6,
    order = c("magnitude", "variance", "none"),
    robust = FALSE,
    prob = 0.95,
    dens_args = list(adjust = 1, n = 512, na.rm = TRUE),
    common_xlim = FALSE,
    ncol = NULL,
    grid_length = 200,
    col = graphics::par("col"),
    band_alpha = 0.25,
    line_lwd = 2,
    xlab = NULL, ylab = NULL, main = NULL,
    xlim = NULL, ylim = NULL,
    ...
) {
  fit <- object
  what <- match.arg(what)
  order <- match.arg(order)

  band_col <- grDevices::adjustcolor(col, alpha.f = band_alpha)

  s <- summary(fit, prob = prob, robust = robust)
  low_key <- sprintf("l-%.0f%% CrI", 100 * prob)
  upp_key <- sprintf("u-%.0f%% CrI", 100 * prob)

  if (what == "beta") {
    beta_samples <- as.matrix(fit$samples$beta)
    if (is.null(colnames(beta_samples))) {
      colnames(beta_samples) <- paste0("beta_", seq_len(ncol(beta_samples)) - 1L)
      colnames(beta_samples)[1] <- "Intercept"
    }

    btab <- as.data.frame(s$beta)

    center_all <- setNames(btab[, "Estimate"], rownames(btab))
    lower_all  <- setNames(btab[, low_key], rownames(btab))
    upper_all  <- setNames(btab[, upp_key], rownames(btab))

    all_names <- rownames(btab)
    if (is.null(pars)) {
      if (order == "variance") {
        ord <- order(btab[, "Est.Error"], decreasing = TRUE)
      } else if (order == "none") {
        ord <- seq_along(all_names)
      } else {
        ord <- order(abs(center_all), decreasing = TRUE)
      }
      keep_names <- all_names[ord][seq_len(min(top_k, length(ord)))]
    } else if (is.character(pars)) {
      keep_names <- pars
    } else {
      keep_names <- all_names[pars]
    }

    idx <- match(keep_names, colnames(beta_samples))
    k <- length(idx)

    if (anyNA(idx)) {
      stop("Some requested parameters in 'pars' were not found.", call. = FALSE)
    }

    if (is.null(ncol)) ncol <- ceiling(sqrt(k))
    nrow <- ceiling(k / ncol)

    if (common_xlim && is.null(xlim)) {
      xr <- range(beta_samples[, idx, drop = FALSE], finite = TRUE)
      xlim <- grDevices::extendrange(xr, f = 0.04)
    }

    if (is.null(xlab)) xlab <- "Value"
    if (is.null(ylab)) ylab <- "Density"
    if (is.null(main)) {
      main <- sprintf(
        "Posterior %s (%.0f%% CrI shaded)",
        if (robust) "Medians" else "Means",
        100 * prob
      )
    }

    op <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(op))
    graphics::par(
      mfrow = c(nrow, ncol),
      mar = c(3.5, 3.5, 2, 1),
      mgp = c(2, 0.7, 0)
    )

    for (j in seq_len(k)) {
      nm <- keep_names[j]
      xj <- beta_samples[, idx[j]]
      dj <- do.call(stats::density, c(list(x = xj), dens_args))

      cj <- center_all[nm]
      lj <- lower_all[nm]
      uj <- upper_all[nm]

      L <- min(lj, uj)
      U <- max(lj, uj)

      xlim_j <- if (!is.null(xlim)) xlim else range(dj$x, finite = TRUE)
      ylim_j <- if (!is.null(ylim)) ylim else range(dj$y, finite = TRUE)

      plot(
        dj,
        type = "n",
        xlab = xlab,
        ylab = ylab,
        main = if (k == 1) main else "",
        xlim = xlim_j,
        ylim = ylim_j,
        ...
      )

      inside <- dj$x >= L & dj$x <= U
      if (any(inside)) {
        px <- c(dj$x[inside], rev(dj$x[inside]))
        py <- c(dj$y[inside], rep(0, sum(inside)))
        polygon(px, py, border = NA, col = band_col)
      } else {
        warning("No density values inside the CrI range for '", nm, "'")
      }

      lines(dj$x, dj$y, lwd = line_lwd, col = col)
      abline(v = cj, lty = 2, col = col)
      mtext(nm, side = 3, line = 0.2, cex = 0.9)
    }

    if (k > 1 && nzchar(main)) {
      mtext(main, outer = TRUE, line = 0.6, cex = 1.1)
    }

    return(invisible(list(what = "beta", selected = keep_names)))
  }

  if (what == "crm") {
    crm <- fit$samples$crm

    z_list <- crm$z.tld
    J_list <- crm$J.tld

    if (is.null(xlim)) {
      xlim <- range(unlist(z_list), finite = TRUE)
    }

    grid <- seq(xlim[1], xlim[2], length.out = grid_length)
    breaks <- c(
      -Inf,
      (head(grid, -1) + tail(grid, -1)) / 2,
      Inf
    )

    crm_to_grid <- function(z, J, breaks) {
      J <- J / sum(J)
      bin <- cut(z, breaks = breaks, labels = FALSE, include.lowest = TRUE)

      out <- numeric(length(breaks) - 1)
      rs <- rowsum(J, bin, reorder = FALSE)
      out[as.integer(rownames(rs))] <- rs[, 1]
      out
    }

    grid_mass <- t(mapply(
      crm_to_grid,
      z = z_list,
      J = J_list,
      MoreArgs = list(breaks = breaks)
    ))

    if (robust) {
      center <- apply(grid_mass, 2, stats::median, na.rm = TRUE)
    } else {
      center <- colMeans(grid_mass, na.rm = TRUE)
    }

    alpha <- (1 - prob) / 2
    lower <- apply(grid_mass, 2, stats::quantile, probs = alpha, na.rm = TRUE)
    upper <- apply(grid_mass, 2, stats::quantile, probs = 1 - alpha, na.rm = TRUE)

    if (is.null(xlab)) xlab <- "Support"
    if (is.null(ylab)) ylab <- "Probability mass"
    if (is.null(main)) {
      main <- sprintf(
        "Posterior %s of CRM (%.0f%% CrI)",
        if (robust) "Median" else "Mean",
        100 * prob
      )
    }

    if (is.null(ylim)) {
      ylim <- range(c(lower, upper, center), finite = TRUE)
      ylim <- grDevices::extendrange(ylim, f = 0.04)
    }

    plot(
      NA,
      xlim = xlim,
      ylim = ylim,
      xlab = xlab,
      ylab = ylab,
      main = main,
      ...
    )

    polygon(
      c(grid, rev(grid)),
      c(lower, rev(upper)),
      border = NA,
      col = band_col
    )

    lines(grid, center, col = col, lwd = line_lwd)

    return(invisible(list(
      what = "crm",
      grid = grid,
      center = center,
      lower = lower,
      upper = upper
    )))
  }
}

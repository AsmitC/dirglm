#' Fit a Copula DP-GLM Model
#'
#'
#' @param formula An object of class \code{\link[stats]{formula}}.
#' @param data Optional data frame containing the variables in \code{formula}.
#' @param link Optional character string or link object specifying the mean link
#' function. If \code{NULL}, a default link is chosen based on the response
#' support: \code{"logit"} for responses in \code{(0, 1)}, \code{"log"} for
#' nonnegative responses, and \code{"identity"} otherwise.
#' @param group_index Vector specifying the group index for each observation.
#' Must have the same length as the response.
#' @param cdpglmControl Optional control arguments.
#' Passed as an object of class "dpglmControl", which is constructed by the
#' \code{dpglm.control} function.
#' @param thetaControl Optional control arguments for the theta update procedure.
#' Passed as an object of class "thetaControl", which is constructed by the
#' \code{theta.control} function.
#' See \code{theta.control} documentation for details.
#'
#' @return List of class "cdpglm" containing information including
#' posterior samples, Metropolis-Hastings acceptance rates, and model specification.
#' See details for more information.
#' 
#' @details 
#' The "cdpglm" class is a list of the following items.
#' \itemize{
#' \item \code{samples} A list containing the MCMC samples for
#' \code{beta}, \code{crm}, and \code{rho}.
#' \item \code{mb} Prior mean for \code{beta}.
#' \item \code{Sb} Prior variance-covariance matrix for \code{beta}.
#' \item \code{formula} Model formula.
#' \item \code{data} Model data frame.
#' \item \code{link} Link function. If a character string was passed to the
#' \code{link} argument, then this will be an object of class "link-glm".
#' Otherwise, it will be the list of three functions passed to the \code{link} argument.
#' \item \code{spt} Support of response variable \code{y} used in model fitting.
#' \item \code{mu0} Mean of the reference distribution obtained by model fit.
#' \item \code{burnin} The number of burn-in MCMC iterations.
#' \item \code{thin} Thinning parameter used for MCMC
#' \item \code{save} The number of saved MCMC iterations.
#' \item \code{rho_prior_shape} Prior parameters for dependence parameter \code{rho}.
#' \item \code{beta_acceptance} Proportion of accepted proposals for beta during MCMC.
#' \item \code{crm_acceptance} Proportion of accepted proposals for CRM during MCMC.
#' \item \code{rho_acceptance} Proportion of accepted proposals for rho during MCMC.
#' }
#'
#' @examples
#' \dontrun{
#' fit <- cdpglm(
#'   y ~ x1 + x2,
#'   data = dat,
#'   group_index = dat$id,
#'   cdpglmControl = cdpglm.control(
#'     spt = c(0, 1),
#'     save = 500,
#'     thin = 5
#'   )
#' )
#' }
#'
#' @export
cdpglm <- function(formula, data = NULL, link = NULL,
                   group_index,
                   cdpglmControl = cdpglm.control(),
                   thetaControl = theta.control()) {

  # Model
  mf <- stats::model.frame(formula, data)
  X  <- stats::model.matrix(attr(mf, "terms"), mf)

  # Warn if covariates appear substantially uncentered
  if (ncol(X) > 1L) {
    X_no_intercept <- X[, -1, drop = FALSE]

    col_means <- colMeans(X_no_intercept, na.rm = TRUE)
    col_sds <- apply(X_no_intercept, 2, stats::sd, na.rm = TRUE)

    scaled_means <- abs(col_means) / pmax(col_sds, .Machine$double.eps)
    bad <- which(is.finite(scaled_means) & scaled_means > 0.25)

    if (length(bad) > 0L) {
      warning(
        paste0(
          "Some covariates appear uncentered relative to their scale: ",
          paste(colnames(X_no_intercept)[bad], collapse = ", "),
          ". Consider centering/scaling covariates before fitting CDPGLM."
        ),
        call. = FALSE
      )
    }
  }
  
  if (qr(X)$rank < ncol(X)) {
    stop("X is singular. Revise formula or data.", call. = FALSE)
  }
  attributes(X)[c("assign", "contrasts")] <- NULL
  y <- stats::model.response(mf, type = "numeric")

  # Link
  if (is.null(link)) {
    if (all(y > 0 & y < 1, na.rm = TRUE)) {
      link <- "logit"
    } else if (all(y >= 0, na.rm = TRUE)) {
      link <- "log"
    } else {
      link <- "identity"
    }
  }

  # Group index
  if (missing(group_index)) {
    stop("cdpglm requires group_index.", call. = FALSE)
  }
  if (length(group_index) != length(y)) {
    stop("group_index must have the same length as the response.", call. = FALSE)
  }

  test.vectorized <- TRUE
  if (is.character(link)) {
    test.vectorized <- FALSE
    link <- stats::make.link(link)
  } else if (!is.list(link) ||
             !all(c("linkfun", "linkinv", "mu.eta") %in% names(link))) {
    stop("link must be a string or a list containing linkfun, linkinv, mu.eta.",
         call. = FALSE)
  }

  linkfun <- link$linkfun
  linkinv <- link$linkinv
  mu.eta  <- link$mu.eta

  if (test.vectorized) {
    is.vectorized <- function(f, data) {
      out <- tryCatch(
        f(data),
        error = function(e) e,
        warning = function(w) w
      )
      if (inherits(out, "error") || inherits(out, "warning")) return(FALSE)
      is.atomic(out) && length(out) == length(data)
    }

    linkfun.testdata <- rep(mean(y), 3)
    inveta.testdata  <- seq(-1, 1, length.out = 3)

    if (!is.vectorized(linkfun, linkfun.testdata) ||
        !is.vectorized(linkinv, inveta.testdata) ||
        !is.vectorized(mu.eta, inveta.testdata)) {
      stop("link must be vectorized.", call. = FALSE)
    }
  }

  mu0 <- if (!is.null(cdpglmControl$mu0)) cdpglmControl$mu0 else mean(y)

  eps <- cdpglmControl$eps
  spt <- cdpglmControl$spt

  y_lo <- min(y) - eps
  y_hi <- max(y) + eps

  if (is.null(spt)) {
    spt <- c(y_lo, y_hi)
  } else {
    if (!is.numeric(spt) || length(spt) != 2) {
      stop("For cdpglm, spt must be NULL or a numeric vector of length 2.",
           call. = FALSE)
    }

    if (is.na(spt[1])) spt[1] <- y_lo
    if (is.na(spt[2])) spt[2] <- y_hi

    if (spt[1] >= spt[2]) stop("Resolved support must satisfy spt[1] < spt[2].", call. = FALSE)
  }

  ## MCMC Initialization
  betaStart    <- cdpglmControl$betaStart
  varbetaStart <- cdpglmControl$varbetaStart
  thetaStart   <- cdpglmControl$thetaStart
  crmStart     <- cdpglmControl$crmStart

  if (is.null(betaStart) || is.null(crmStart) ||
      is.null(varbetaStart) || is.null(thetaStart)) {

    gfit <- gldrm(
      formula      = formula,
      data         = data,
      link         = link,
      mu0          = mu0,
      thetaControl = thetaControl
    )

    if (is.null(betaStart)) {
      betaStart <- gfit$beta
      if (any(is.na(betaStart))) {
        lmcoef <- stats::lm.fit(X, linkfun(rep(mu0, length(y))))$coef
        betaStart <- lmcoef
      }
    }

    if (is.null(varbetaStart)) varbetaStart <- gfit$varbeta
    if (is.null(thetaStart)) thetaStart <- gfit$theta

    if (is.null(crmStart)) {
      crmStart <- list(
        z.tld = as.numeric(gfit$spt),
        J.tld = as.numeric(gfit$f0)
      )
    }
  }

  init <- list(
    beta = betaStart,
    varbeta = varbetaStart,
    theta = thetaStart,
    crm = crmStart,
    rho = cdpglmControl$rhoStart
  )

  ## Fit
  fit <- cdpglmFit(
    formula       = formula,
    data          = data,
    X             = X,
    y             = y,
    group_index   = group_index,
    link          = link,
    mu0           = mu0,
    spt           = spt,
    init          = init,
    cdpglmControl = cdpglmControl,
    thetaControl  = thetaControl
  )

  out <- list(
    samples     = fit$samples,
    mb          = fit$mb,
    Sb          = fit$Sb,
    formula     = formula,
    data        = data.frame(mf),
    link        = link,
    spt         = fit$spt,
    mu0         = fit$mu0,
    burnin      = cdpglmControl$burnin,
    thin        = cdpglmControl$thin,
    save        = cdpglmControl$save,
    rho_prior_shape = cdpglmControl$rho_prior_shape,
    beta_acceptance  = fit$beta_acceptance,
    crm_acceptance   = fit$crm_acceptance,
    rho_acceptance   = fit$rho_acceptance
  )

  class(out) <- "cdpglm"
  out
}
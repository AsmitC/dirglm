#' Fits a DP-GLM model
#'
#' @import stats
#' @import mvtnorm
#' @import gldrm
#'
#' @param formula An object of class  \code{\link[stats]{formula}}.
#' @param data An optional data frame containing the variables in the model.
#' @param link Optional character string or link object specifying the mean link
#' function. If \code{NULL}, a default link is chosen based on the response
#' support: \code{"logit"} for responses in \code{(0, 1)}, \code{"log"} for
#' nonnegative responses, and \code{"identity"} otherwise.
#' @param dpglmControl Optional control arguments.
#' Passed as an object of class "dpglmControl", which is constructed by the
#' \code{dpglm.control} function.
#' See \code{dpglm.control} documentation for details.
#' @param thetaControl Optional control arguments for the theta update procedure.
#' Passed as an object of class "thetaControl", which is constructed by the
#' \code{theta.control} function.
#' See \code{theta.control} documentation for details.
#'
#' @return List of class "dpglm" containing information including
#' posterior samples, Metropolis-Hastings acceptance rates, and model specification.
#' See details for more information.
#' 
#' @details
#' The "dpglm" class is a list of the following items.
#' \itemize{
#' \item \code{samples} A list containing the MCMC samples for
#' \code{beta} and \code{crm}.
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
#' \item \code{beta_acceptance} Proportion of accepted proposals for beta during MCMC.
#' \item \code{crm_acceptance} Proportion of accepted proposals for CRM during MCMC.
#' }
#'
#' @examples
#' \dontrun{
#' set.seed(123)
#'
#' n <- 50
#' x <- rnorm(n)
#' y <- 1 + 0.5 * x + rnorm(n, sd = 0.25)
#'
#' dat <- data.frame(y = y, x = x)
#'
#' fit <- dpglm(
#'   y ~ x,
#'   data = dat,
#'   link = "identity",
#'   dpglmControl = dpglm.control(
#'     burnin = 10,
#'     thin = 1,
#'     save = 50,
#'     spt = range(y),
#'     seed = 123
#'   )
#'  )
#' }
#' 
#' @export
dpglm <- function(formula, data = NULL, link = NULL,
                  dpglmControl = dpglm.control(), thetaControl = theta.control()) {

  # Model initialization
  mf <- stats::model.frame(formula, data)
  X  <- stats::model.matrix(attr(mf, "terms"), mf)
  if (qr(X)$rank < ncol(X)) stop("X is singular. Revise formula or data.")
  attributes(X)[c("assign", "contrasts")] <- NULL
  y  <- stats::model.response(mf, type = "numeric")

  ## Extract link
  if (is.null(link)) {
    if (all(y > 0 & y < 1, na.rm = TRUE)) {
      link <- "logit"
    } else if (all(y >= 0, na.rm = TRUE)) {
      link <- "log"
    } else {
      link <- "identity"
    }
  }

  test.vectorized <- TRUE
  if (is.character(link)) {
    test.vectorized <- FALSE
    link <- stats::make.link(link)
  } else if (!is.list(link) ||
             !all(c("linkfun", "linkinv", "mu.eta") %in% names(link))) {
    stop("link must be a string or a list containing linkfun, linkinv, mu.eta")
  }

  # Check that all link functions are vectorized
  linkfun <- link$linkfun
  linkinv <- link$linkinv
  mu.eta  <- link$mu.eta

  if (test.vectorized) { # User has specified a custom link
    is.vectorized <- function(f, data) {
      out <- tryCatch(
        f(data),
        error   = function(e) e,
        warning = function(w) w
      )
      if (inherits(out, "error") || inherits(out, "warning")) return(FALSE)
      is.atomic(out) && length(out) == length(data)
    }

    linkfun.testdata <- rep(mean(y), 3)
    inveta.testdata  <- seq(-1, 1, length.out = 3)

    if (!is.vectorized(linkfun, linkfun.testdata) ||
        !is.vectorized(linkinv, inveta.testdata)  ||
        !is.vectorized(mu.eta, inveta.testdata)) {
      stop("link must be vectorized.")
    }
  }

  # Extract theoretical support
  spt <- dpglmControl$spt
  eps <- dpglmControl$eps
  y_lo <- min(y) - eps
  y_hi <- max(y) + eps

  if (is.null(spt)) {
  spt <- c(y_lo, y_hi)
    } else {
    if (!is.numeric(spt) || length(spt) != 2) {
        stop("spt must be NULL or a numeric vector of length 2.")
    }

    if (is.na(spt[1])) spt[1] <- y_lo
    if (is.na(spt[2])) spt[2] <- y_hi

    if (spt[1] >= spt[2]) {
        stop("spt must satisfy spt[1] < spt[2].")
        }
    }

  mu0 <- if (!is.null(dpglmControl$mu0)) dpglmControl$mu0 else mean(y)

  # MCMC Initialization
  betaStart    <- dpglmControl$betaStart
  varbetaStart <- dpglmControl$varbetaStart
  thetaStart   <- dpglmControl$thetaStart
  crmStart     <- dpglmControl$crmStart

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
        lmcoef <- stats::lm.fit(X, linkfun(mu0))$coef
        betaStart <- lmcoef
      }
    }

    if (is.null(varbetaStart)) varbetaStart <- gfit$varbeta
    if (is.null(thetaStart))   thetaStart   <- gfit$theta

    if (is.null(crmStart)) {
      z.tldStart <- gfit$spt
      J.tldStart <- gfit$f0
      crmStart <- list(z.tld = z.tldStart, J.tld = J.tldStart)
    }
  }

  init <- list(
    beta    = betaStart,
    varbeta = varbetaStart,
    theta   = thetaStart,
    crm     = crmStart
  )

  # Fit
  fit <- dpglmFit(
    formula      = formula,
    X            = X,
    y            = y,
    link         = link,
    spt          = spt,
    mu0          = mu0,
    init         = init,
    dpglmControl = dpglmControl,
    thetaControl = thetaControl
  )

  # Output
  out <- list(
    samples          = fit$samples,
    mb               = fit$mb,
    Sb               = fit$Sb,
    formula          = formula,
    data             = data.frame(mf),
    link             = link,
    spt              = fit$spt,
    mu0              = fit$mu0,
    burnin           = dpglmControl$burnin,
    thin             = dpglmControl$thin,
    save             = dpglmControl$save,
    beta_acceptance  = fit$beta_acceptance,
    crm_acceptance   = fit$crm_acceptance
  )

  class(out) <- "dpglm"
  out
}
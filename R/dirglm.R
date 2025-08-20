#' Fits a finite-support Dir-GLM model
#'
#' @import stats
#' @import mvtnorm
#' @import gldrm
#'
#' @param formula An object of class "formula".
#' @param data An optional data frame containing the variables in the model.
#' @param link Link function. Defaults to log.
#' Can be a character string to be passed to the
#' \code{make.link} function in the \code{stats} package (e.g. "identity",
#' "logit", or "log").
#' Alternatively, \code{link} can be a list containing three functions named
#' \code{linkfun}, \code{linkinv}, and \code{mu.eta}. The first is the link
#' function. The second is the inverse link function. The third is the derivative
#' of the inverse link function. All three functions must be vectorized.
#' @param mb Prior mean for beta. Defaults to a p-length vector whose entries are all 0.
#' @param Sb Prior variance-covariance matrix for beta.
#' Defaults to the p-dimensional identity matrix. See details for more information.
#' @param dir_pr_parm Dirichlet prior parameter for f0. Defaults to the observed
#' response frequency distribution. If specified, it should be a p-length vector
#' with positive entries.
#' @param dirglmControl Optional control arguments.
#' Passed as an object of class "dirglmControl", which is constructed by the
#' \code{dirglm.control} function.
#' See \code{dirglm.control} documentation for details.
#' @param thetaControl Optional control arguments for the theta update procedure.
#' Passed as an object of class "thetaControl", which is constructed by the
#' \code{theta.control} function.
#' See \code{theta.control} documentation for details.
#'
#' @return An S3 object of class "dirglm". See details.
#'
#' @details The arguments \code{linkfun}, \code{linkinv}, and \code{mu.eta}
#' mirror the "link-glm" class. Objects of this class can be created with the
#' \code{stats::make.link} function.
#'
#' This package currently only supports joint updates for beta
#' when \code{Sb} is non-diagonal.
#'
#'The "dirglm" class is a list of the following items.
#' \itemize{
#' \item \code{samples} A list containing the MCMC samples for \code{f0} and \code{beta}.
#' \item \code{mb} Prior mean for \code{beta}.
#' \item \code{Sb} Prior variance-covariance matrix for \code{beta}.
#' \item \code{dir_pr_parm} Dirichlet prior parameter.
#' \item \code{formula} Model formula.
#' \item \code{data} Model data frame.
#' \item \code{link} Link function. If a character string was passed to the
#' \code{link} argument, then this will be an object of class "link-glm".
#' Otherwise, it will be the list of three functions passed to the \code{link} argument.
#' \item \code{spt} Support, \eqn{\{s_j, \ j = 1, 2, ..., l\}}, of response variable \code{y}.
#' \item \code{mu0} Mean of the reference distribution \code{f0}.
#' \item \code{iter} The total number of MCMC iterations.
#' \item \code{p_acc_beta} Proportion of accepted proposals for beta during MCMC.
#' \item \code{p_acc_f0} Proportion of accepted proposals for f0 during MCMC.
#' }
#'
#' @export
dirglm <- function(formula, data=NULL, link="log", mb=NULL, Sb=NULL, dir_pr_parm=NULL,
                   dirglmControl=dirglm.control(), thetaControl=theta.control())

{
  # Model initialization
  mf <- stats::model.frame(formula, data)
  X  <- stats::model.matrix(attr(mf, "terms"), mf)
  if (qr(X)$rank < ncol(X)) stop("X is singular. Revise formula or data.")
  attributes(X)[c("assign", "contrasts")] <- NULL
  y  <- stats::model.response(mf, type = "numeric")

  ## 2. Extract link
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
      out <- tryCatch(f(data),
                      error   = function(e) e,
                      warning = function(w) w)
      if (inherits(out, "error") || inherits(out, "warning")) return(FALSE)
      is.atomic(out) && length(out) == length(data)
    }

    linkfun.testdata <- rep(mean(y), 3)
    inveta.testdata  <- seq(-1, 1, length.out=3)

    if (!is.vectorized(linkfun, linkfun.testdata) ||
        !is.vectorized(linkinv, inveta.testdata)  ||
        !is.vectorized(mu.eta, inveta.testdata)) stop("link must be vectorized.")
  }

  # Initialize (theoretical) support if not provided by the user
  spt <- dirglmControl$spt
  if (is.null(spt)) spt <- sort(unique(y)) ## Observed support
  if (is.unsorted(spt)) spt <- sort(spt)
  l <- length(spt)

  # Initialize mu0 if not provided by the user
  mu0 <- dirglmControl$mu0
  if (is.null(mu0)) mu0 <- mean(y)
  else if (mu0 <= min(spt) || mu0 >= max(spt)) {
    stop(paste0("mu0 must lie within the range of observed values. Choose a different ",
                "value or set mu0=NULL to use the default value, mean(y)."))
  }

  # MCMC Initialization
  betaStart    <- dirglmControl$betaStart
  f0Start      <- dirglmControl$f0Start

  # beta
  if (is.null(betaStart)) {
    gfit <- gldrm(formula      = formula,
                  data         = data,
                  link         = link,
                  mu0          = mu0,
                  thetaControl = thetaControl)
    betaStart <- gfit$beta
    if (any(is.na(betaStart))) {
      lmcoef <- stats::lm.fit(x, linkfun(mu0))$coef
      betaStart <- lmcoef
    }
  }

  # f0
  if (is.null(f0Start)) {
    f0      <- rep(1 / l, l)
    tht0    <- gldrm:::getTheta(
      spt       = spt,
      f0        = f0,
      mu        = mu0,
      sampprobs = NULL,
      ySptIndex = NULL
    )$theta
    f0star  <- (f0 * exp(tht0 * spt)) / sum(f0 * exp(tht0 * spt))
    f0Start <- f0star
  }

  init <- list(beta = betaStart, f0 = f0Start)

  # Fit
  fit <- dirglmFit(
    formula              = formula,
    X                    = X,
    y                    = y,
    link                 = link,
    mb                   = mb,
    Sb                   = Sb,
    dir_pr_parm          = dir_pr_parm,
    mu0                  = mu0,
    spt                  = spt,
    init                 = init,
    dirglmControl        = dirglmControl,
    thetaControl         = thetaControl
  )

  # Output
  out <- list(
    samples     = fit$samples,
    mb          = fit$mb,
    Sb          = fit$Sb,
    dir_pr_parm = fit$dir_pr_parm,
    formula     = formula,
    data        = data.frame(mf),
    link        = link,
    spt         = fit$spt,
    mu0         = fit$mu0,
    burnin      = dirglmControl$burnin,
    thin        = dirglmControl$thin,
    save        = dirglmControl$save,
    p_acc_beta  = fit$p_acc_beta,
    p_acc_f0    = fit$p_acc_f0
  )
  class(out) <- "dirglm"
  out
}

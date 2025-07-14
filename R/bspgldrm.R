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
#' @param sb Vector containing the diagonal entries of the prior variance-covariance
#' matrix for beta. Defaults to the p-dimensional identity matrix.
#' @param dir_pr_parm Dirichlet prior parameter for f0. Defaults to the observed
#' response frequency distribution. If specified, it should be a p-length vector
#' with positive entries.
#' @param bspgldrmControl Optional control arguments.
#' Passed as an object of class "bspgldrmControl", which is constructed by the
#' \code{bspgldrm.control} function.
#' See \code{bspgldrm.control} documentation for details.
#' @param thetaControl Optional control arguments for the theta update procedure.
#' Passed as an object of class "thetaControl", which is constructed by the
#' \code{theta.control} function.
#' See \code{theta.control} documentation for details.
#'
#' @return An S3 object of class "bspgldrm". See details.
#'
#' @details The arguments \code{linkfun}, \code{linkinv}, and \code{mu.eta}
#' mirror the "link-glm" class. Objects of this class can be created with the
#' \code{stats::make.link} function.
#'
#'The "bspgldrm" class is a list of the following items.
#' \itemize{
#' \item \code{samples} A list containing the MCMC samples for \code{f0} and \code{beta}.
#' \item \code{mb} Prior mean for \code{beta}.
#' \item \code{sb} Diagonal entries of the prior variance-covariance matrix for \code{beta}.
#' \item \code{dir_pr_parm} Dirichlet prior parameter.
#' \item \code{formula} Model formula.
#' \item \code{data} Model data frame.
#' \item \code{link} Link function. If a character string was passed to the
#' \code{link} argument, then this will be an object of class "link-glm".
#' Otherwise, it will be the list of three functions passed to the \code{link} argument.
#' \item \code{spt} Support, \eqn{\{s_j, \ j = 1, 2, ..., l\}}, of response variable \code{y}.
#' \item \code{mu0} Mean of the reference distribution \code{f0}.
#' }
#'
#' @examples
#' data(iris, package="datasets")
#'
#' # Fit a bspgldrm with log link
#' fit <- bspgldrm(Sepal.Length ~ Sepal.Width + Petal.Length + Petal.Width + Species,
#'                 data=iris)
#' fit
#'
#' # Fit a bspgldrm with custom link function
#' link <- list()
#' link$linkfun <- function(mu) log(mu)^3
#' link$linkinv <- function(eta) exp(eta^(1/3))
#' link$mu.eta <- function(eta) exp(eta^(1/3)) * 1/3 * eta^(-2/3)
#' fit2 <- bspgldrm(Sepal.Length ~ Sepal.Width + Petal.Length + Petal.Width + Species,
#'                  data=iris, link=link)
#' fit2
#'
#' @export
bspgldrm <- function(formula, data=NULL, link="log", mb=NULL, sb=NULL, dir_pr_parm=NULL,
                     bspgldrmControl=bspgldrm.control(), thetaControl=theta.control())

{
  ## 1. Model initialization
  mf <- stats::model.frame(formula, data)
  X  <- stats::model.matrix(attr(mf, "terms"), mf)
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

  ### 2.1 Check that all link functions are vectorized
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

  ## 3. Fit
  fit <- bspgldrmFit(
    formula              = formula,
    data                 = data,
    X                    = X,
    y                    = y,
    link                 = link,
    mb                   = mb,
    sb                   = sb,
    dir_pr_parm          = dir_pr_parm,
    bspgldrmControl      = bspgldrmControl,
    thetaControl         = thetaControl
  )

  ## 4. Output
  out <- list(
    samples     = fit$samples,
    mb          = fit$mb,
    sb          = fit$sb,
    dir_pr_parm = fit$dir_pr_parm,
    formula     = formula,
    data        = data.frame(mf),
    link        = link,
    spt         = fit$spt,
    mu0         = fit$mu0
  )
  class(out) <- "bspgldrm"
  out
}

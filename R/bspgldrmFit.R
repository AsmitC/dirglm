#' Control arguments for the \code{bspgldrm} algorithm.
#'
#' This function returns control arguments for the \code{bspgldrm} algorithm.
#' Each argument has a default value, which will be used unless a different
#' value is provided by the user.
#'
#' @param burnin Number of burn-in MCMC iterations. Defaults to 100.
#' @param thin Factor by which to thin MCMC iterations. Defaults to 10.
#' @param save Number of MCMC samples to return. Defaults to 1000.
#' @param rho MCMC update step size. A scalar in \eqn{(0, 1]}. Defaults to 0.1.
#' @param mu0 Mean of the reference distribution. The reference distribution is
#' not unique unless its mean is restricted to a specific value. This value can
#' be any number within the range of observed values, but values near the boundary
#' may cause numerical instability. This is an optional argument with \code{mean(y)}
#' being the default value.
#' @param spt Theoretical support of the response variable.
#' @param betaStart Initial value for the regression coefficients \code{beta}.
#' Defaults to the output obtained by fitting \code{gldrm}.
#' @param f0Start Initial value for the reference distribution \code{f0}.
#' Defaults to the output obtained by fitting \code{gldrm}.
#' @param joint.update Logical indicating whether to update \code{beta} jointly.
#' Defaults to \code{TRUE}.
#' @param seed Random seed. Defaults to NULL.
#'
#' @return Object of S3 class "bspgldrmControl"
#'
#' @export
bspgldrm.control <- function(burnin=100, thin=10, save=1000, rho=0.1, mu0=NULL, spt=NULL,
                             betaStart=NULL, f0Start=NULL, joint.update=TRUE, seed=NULL)
{
  if (burnin < 0 || floor(burnin) != burnin) stop("Number of burn-in samples must be an integer >= 0")
  if (thin   < 1 || floor(thin)   != thin)   stop("Thin must be an integer >= 1")
  if (save   < 1 || floor(save)   != save)   stop("Number of saved iterations must be an integer >= 1")
  if (!(rho <= 1 & rho > 0))                 stop("rho must lie in (0, 1]")
  if (!is.logical(joint.update) ||
      !joint.update %in% c(T, F))            stop("joint.update must be logical TRUE/FALSE")
  ctrl <- list(burnin       = burnin,
               thin         = thin,
               save         = save,
               rho.         = rho,
               mu0          = mu0,
               spt          = spt,
               betaStart    = betaStart,
               f0Start      = f0Start,
               joint.update = joint.update,
               seed         = seed)
  class(ctrl) <- "bspgldrmControl"
  ctrl
}

#' Main MCMC function
#'
#' This function is called by the main \code{bspgldrm} function.
#'
#' @keywords internal
bspgldrmFit <- function(formula, data, X, y,                # Data
                        link,                               # Link
                        mb, sb, dir_pr_parm,                # Priors
                        bspgldrmControl, thetaControl)      # Controls
{
  ## 1. Extract bspgldrmControl parameters
  burnin       <- bspgldrmControl$burnin
  thin         <- bspgldrmControl$thin
  save         <- bspgldrmControl$save
  rho          <- bspgldrmControl$rho
  mu0          <- bspgldrmControl$mu0
  spt          <- bspgldrmControl$spt
  joint.update <- bspgldrmControl$joint.update
  betaStart    <- bspgldrmControl$betaStart
  f0Start      <- bspgldrmControl$f0Start
  seed         <- bspgldrmControl$seed

  if (!is.null(seed)) set.seed(seed) # Set seed

  ## 1.1 Extract link
  linkfun <- link$linkfun
  linkinv <- link$linkinv
  mu.eta  <- link$mu.eta

  ## 2. Initialize (theoretical) support if not provided by the user
  if (is.null(spt)) spt <- sort(unique(y)) ## Observed support
  l <- length(spt)

  ## 3. Initialize mu0 if not provided by the user
  if (is.null(mu0)) mu0 <- mean(y)
  else if (mu0 <= min(spt) || mu0 >= max(spt)) {
    stop(paste0("mu0 must lie within the range of observed values. Choose a different ",
                "value or set mu0=NULL to use the default value, mean(y)."))
  }

  ## 4. MCMC Initialization

  ### 4.1 beta
  if (is.null(betaStart)) {
    gfit <- gldrm(formula      = formula,
                  data         = data,
                  link         = link,
                  mu0          = mu0,
                  thetaControl = thetaControl)
    betaStart <- gfit$beta
  }

  ### 4.2 f0
  if (is.null(f0Start)) {
    f0      <- rep(1 / l, l)
    tht0    <- gldrm:::getTheta(
      spt       = spt,
      f0        = f0,
      mu        = mu0,
      sampprobs = NULL,
      ySptIndex = NULL
    )$theta
    f0star  <- (f0 * exp(tht0 * spt)) %>% `/` (sum(.))
    f0Start <- f0star
  }

  init <- list(beta = betaStart, f0 = f0Start)
  X    <- as.matrix(X)
  n    <- length(y)
  p    <- ncol(X)
  iter <- burnin + thin * save

  beta_samples <- matrix(NA, nrow = save, ncol = p)
  f0_samples   <- matrix(NA, nrow = save, ncol = l)

  beta <- init$beta
  f0   <- init$f0
  beta_samples[1, ] <- beta
  f0_samples[1, ]   <- f0

  ### 4.3 Validate priors
  ### 4.3.1 Beta prior
  if (is.null(mb) || is.null(sb)) {
    if (is.null(mb)) mb <- rep(0, p)
    if (is.null(sb)) sb <- rep(1, p)
  } else if (length(mb) != p) stop("length(mb) must match the number of covariates.")
  else if   (length(sb) != p) stop("length(sb) must match the number of covariates.")
  else if   (!all(sb)    > 0) stop(paste0("Beta prior variance-covariance matrix must be positive definite. ",
                                   "Check that all(sb > 0)."))

  ### 4.3.2 Dirichlet prior
  if (is.null(dir_pr_parm)) {
    ind_mt      <- outer(y, spt, `==`) * 1
    alpha       <- 1
    dir_pr_parm <- alpha * colMeans(ind_mt)
    eps         <- 1e-6
    dir_pr_parm <- dir_pr_parm + eps
  } else if (!all(dir_pr_parm > 0) ||
             length(dir_pr_parm)   != l) stop("dir_pr_parm must be positive with K atoms.")

  ### 4.4 Theta
  mu      <- linkinv(X %*% beta)                   # Updated for general link
  out     <- tht_sol(spt, f0, mu, NULL)
  tht     <- out$tht
  btht    <- out$btht
  bpr2    <- out$bpr2
  f0_y    <- f0y(y, spt, f0)

  n_acc <- 0 # Count number of acceptences
  ## 5. MCMC loop
  for (r in 2:iter) {
    ### 5.1 Beta update
    Sig  <- Sigma_beta(X, mu, bpr2, rho, linkfun, mu.eta)
    if (joint.update) {
      b_out <- beta_update_joint(X, y, spt, beta, Sig, f0, tht,
                                 bpr2, btht, rho, linkfun, linkinv,
                                 mu.eta, mb, sb)
    } else {
      b_out <- beta_update_separate(X, y, spt, beta, Sig, f0, tht,
                                     bpr2, btht, rho, linkfun, linkinv,
                                     mu.eta, mb, sb)
    }
    beta <- b_out$cr_bt
    tht  <- b_out$cr_tht
    btht <- b_out$cr_btht
    bpr2 <- b_out$cr_bpr2
    acc  <- b_out$acc
    mu   <- linkinv(X %*% beta)                   # Updated for general link

    n_acc <- n_acc + as.integer(acc)              # Increment number of acceptances

    ### 5.2 f0 update
    propsl_dir_parm <- dir_parm(y, tht, btht, dir_pr_parm, ind_mt)
    out             <- f0_update(y, spt, f0, f0_y, propsl_dir_parm,
                                 mu, tht, bpr2, btht, dir_pr_parm, ind_mt)
    f0    <- out$cr_f0
    f0_y  <- out$cr_f0y
    tht   <- out$cr_tht
    btht  <- out$cr_btht
    bpr2  <- out$cr_bpr2

    # 5.3 Storage
    if (r > burnin & r %% thin == 0) {
      j <- (r - burnin) / thin
      beta_samples[j, ] <- beta
      f0_samples[j, ]   <- f0
    }
  }

  p_acc <- n_acc / iter # Proportion of acceptances
  if (p_acc < 0.01) warning(paste0("Markov chain did not mix well.",
                                   "Consider editing rho or increasing the total number of iterations."))

  ## 6. Tilt each f0 sample
  f0star_samples <- matrix(0, nrow = nrow(f0_samples), ncol = length(spt))
  for (iter in 1:nrow(f0_samples)) {
    wh     <- f0_samples[iter, ]
    theta0 <- gldrm:::getTheta(
      spt = spt,
      f0 = wh,
      mu = mu0,
      sampprobs = NULL,
      ySptIndex = NULL
    )$theta
    wh                     <- wh * exp(theta0 * spt)
    wh                     <- wh / sum(wh)
    f0star_samples[iter, ] <- wh
  }
  f0_samples     <- f0star_samples  # projected f0 samples

  ## 7. Output
  list(samples     = list(beta = beta_samples,
                          f0   = f0_samples),
       mb          = mb,
       sb          = sb,
       dir_pr_parm = dir_pr_parm,
       spt         = spt,
       mu0         = mu0,
       p_acc       = p_acc,
       iter        = iter)
}

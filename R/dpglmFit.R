#' Control arguments for the \code{dpglm} algorithm.
#'
#' This function returns control arguments for the \code{dpglm} algorithm.
#'
#' @param burnin Number of burn-in MCMC iterations. Defaults to 100.
#' @param thin Factor by which to thin MCMC iterations. Defaults to 10.
#' @param save Number of MCMC samples to return. Defaults to 1000.
#' @param mb Prior mean for beta. Defaults to NULL.
#' @param Sb Prior scale / diagonal variance entries for beta. Defaults to NULL.
#' @param M Number of random CRM atoms proposed at each CRM update. Defaults to 50.
#' @param alpha CRM mass parameter. Defaults to 1.
#' @param delta Tuning parameter for the \code{u} update. Defaults to 2.
#' @param c0 Bandwidth parameter used in the latent \code{z} update. Defaults to NULL.
#' @param gamma Shrinkage parameter for the default prior variance on \code{beta}.
#' Defaults to 1. Ignored if \code{Sb} is specified.
#' @param mu0 Mean of the reference distribution. Defaults to \code{NULL}, in which case
#' \code{mean(y)} is used internally.
#' @param spt Theoretical support for the response in the DP-GLM. This should be a
#' length-2 vector giving the lower and upper bounds. Defaults to \code{NULL}, in which
#' case the support is inferred internally from the data as \code{c(min(y)-eps, max(y)+eps)}.
#' @param eps Padding used when constructing default support bounds. Defaults to 1e-6.
#' @param betaStart Initial value for \code{beta}. Defaults to NULL.
#' @param varbetaStart Initial value for the proposal covariance of \code{beta}. Defaults to NULL.
#' @param thetaStart Initial value for \code{theta}. Defaults to NULL.
#' @param crmStart Initial value for the CRM state, given as a list with components
#' \code{z.tld} and \code{J.tld}. Defaults to NULL.
#' @param seed Random seed. Defaults to NULL.
#' @param robust Provides numerical stability. Defaults to \code{TRUE}. See details for more.
#'
#' @return Object of S3 class "dpglmControl".
#'
#' @details Setting \code{robust = TRUE} will tilt the CRM weights at each MCMC iteration
#' to have the desired mean \code{mu0}. If \code{robust = FALSE}, these weights can
#' vary more in magnitude and lead to downstream numerical instability.
#'
#' @return An object of class \code{"dpglmControl"}; a list of control
#' parameters for the DPGLM fitting function.
#'
#' @examples
#' \dontrun{
#' ctrl <- dpglm.control(burnin = 100,
#'                       thin = 2,
#'                       save = 500,
#'                       spt = c(0, 1),
#'                       seed = 123)
#' }
#'
#' @export
dpglm.control <- function(burnin = 100, thin = 10, save = 1000,
                          mb = NULL, Sb = NULL,
                          M = 50, alpha = 1, delta = 2, c0 = NULL,
                          gamma = 1, mu0 = NULL, spt = NULL, eps = 1e-6,
                          betaStart = NULL, varbetaStart = NULL,
                          thetaStart = NULL, crmStart = NULL,
                          seed = NULL, robust=TRUE)
{
  if (burnin < 0 || floor(burnin) != burnin) {
    stop("Number of burn-in samples must be an integer >= 0")
  }
  if (thin < 1 || floor(thin) != thin) {
    stop("Thin must be an integer >= 1")
  }
  if (save < 1 || floor(save) != save) {
    stop("Number of saved iterations must be an integer >= 1")
  }

  if (!is.null(spt)) {
    if (!is.numeric(spt) || length(spt) != 2) {
      stop("For dpglm, spt must be a numeric vector of length 2 giving the support bounds.")
    }
    if (spt[1] >= spt[2]) {
      stop("For dpglm, spt must satisfy spt[1] < spt[2].")
    }
  }

  ctrl <- list(
    burnin       = burnin,
    thin         = thin,
    save         = save,
    mb           = mb,
    Sb           = Sb,
    M            = M,
    alpha        = alpha,
    delta        = delta,
    c0           = c0,
    gamma        = gamma,
    mu0          = mu0,
    spt          = spt,
    eps          = eps,
    betaStart    = betaStart,
    varbetaStart = varbetaStart,
    thetaStart   = thetaStart,
    crmStart     = crmStart,
    seed         = seed,
    robust       = robust
  )

  class(ctrl) <- "dpglmControl"
  ctrl
}


#' Main MCMC Fitting Function for DPGLM
#'
#' Performs MCMC to fit the DPGLM model.
#'
#' @param formula Model formula supplied to \code{\link{dpglm}}.
#' @param data Data frame containing the variables appearing in \code{formula}.
#' @param X Numeric design matrix.
#' @param y Numeric response vector.
#' @param link Link object containing \code{linkfun}, \code{linkinv}, and \code{mu.eta}.
#' @param spt Numeric vector of length two giving the lower and upper bounds
#' of the response support.
#' @param mu0 Numeric target mean used to identify the reference distribution.
#' @param init List of initial values for the MCMC algorithm.
#' @param dpglmControl Object of class \code{"dpglmControl"} containing MCMC
#' and tuning parameters.
#' @param thetaControl Object of class \code{"thetaControl"} containing control
#' arguments for the theta update procedure.
#'
#' @return
#' A list containing posterior samples, acceptance rates, prior settings,
#' support information, and model specification information used by
#' \code{\link{dpglm}}.
#'
#' @keywords internal
dpglmFit <- function(formula, data, X, y,        # Data
                     link,                       # Link
                     spt, mu0, init,             # Model Specs
                     dpglmControl, thetaControl) # Controls
{
  # Extract dpglmControl parameters
  burnin       <- dpglmControl$burnin
  thin         <- dpglmControl$thin
  save         <- dpglmControl$save
  mb           <- dpglmControl$mb
  Sb           <- dpglmControl$Sb
  M            <- dpglmControl$M
  alpha        <- dpglmControl$alpha
  delta        <- dpglmControl$delta
  c0           <- dpglmControl$c0
  gamma        <- dpglmControl$gamma
  eps          <- dpglmControl$eps
  copula       <- dpglmControl$copula
  seed         <- dpglmControl$seed
  robust       <- dpglmControl$robust

  # Set theoretical support bounds
  if (is.unsorted(spt)) spt <- sort(spt)
  min_y <- spt[1]
  max_y <- spt[2]

  # MCMC Initialization
  n    <- length(y)
  X    <- as.matrix(X, nrow=n)
  p    <- dim(X)[2]
  iter <- burnin + thin * save

  if (is.null(copula) || !is.list(copula)) copula <- list()

  # Default h for the DPGLM model
  h <- h_ <- hstar <- 0

  if (!is.null(seed)) set.seed(seed)

  if (is.null(c0)) c0 <- c0_silverman(y) / 4

  # Extract link
  linkfun <- link$linkfun
  linkinv <- link$linkinv

  beta_samples <- matrix(NA, nrow = save, ncol = p)
  beta_names <- colnames(X)
  if (is.null(beta_names)) beta_names <- paste0("beta_", seq_len(p) - 1L)
  beta_names[1] <- "Intercept"
  colnames(beta_samples) <- beta_names
  theta_samples <- matrix(NA, nrow = save, ncol = n)
  u_samples <- z_samples <- matrix(NA, nrow = save, ncol = n)
  crm_samples   <- vector("list", save)
  lnlik_samples <- numeric(save)

  beta <- as.numeric(init$beta)
  #beta_samples[1, ] <- beta
  beta_mode <- beta
  beta_cov <- init$varbeta
  meanY_x <- linkinv(X %*% beta)
  z.tld <- z_tld <-  as.numeric(init$crm$z.tld)
  J.tld <- J_tld <- as.numeric(init$crm$J.tld)
  #crm_samples[[1]] <- list(z.tld = z.tld, J.tld = J.tld)

  ord <- order(J_tld)[1:M]

  RL <- z_tld[ord]
  RJ <- J_tld[ord]
  theta <- theta_ref <- init$theta
  #theta_samples[1, ] <- theta
  z <- z_sampler_unifK(y, c0, z.tld, J.tld, theta, min_y, max_y, eps)
  #z_samples[1, ] <- z

  btheta <- b_theta(theta, z.tld, J.tld)

  T_vec <- exp(btheta)
  u <- rgamma(n, shape = 1, rate = T_vec)
  #u_samples[1, ] <- u <- rgamma(n, shape = 1, rate = T_vec)

  resampled_z <- resample_zstar(z)
  zstar <- resampled_z$zstar
  nstar <- resampled_z$nstar
  Jstar <- rgamma(n = length(nstar), shape = nstar, rate = 1)

  count1 <- count2 <- 0
  burning <- TRUE

  theta0 <- getTheta(
    spt = z.tld,
    f0  = J.tld,
    mu  = mu0,
    sampprobs  = NULL,
    ySptIndex  = NULL,
    thetaStart = NULL
  )$theta

  if (any(theta0 > 50)) {
    idx <- which(theta0 > 50)
    theta0[idx] <- 50
    warning("Capped some values of Theta at 50. Consider increasing M.")
  }

  temp <- exp(theta0 * z.tld - max(theta0 * z.tld))
  Jtilt <- temp * J.tld
  W <- sum(Jtilt)
  Jtld_0 <- Jtilt / W

  if (robust) {
    J.tld <- Jtilt
    J.tld_ll <- Jtld_0
  } else {
    J.tld_ll <- Jtld_0
  }

  lnlik <- loglik(linkinv = linkinv, z = z, X = X, beta = beta, atoms = z.tld, jumps = J.tld_ll)
  # lnlik_samples[1] <- lnlik

  # Beta prior
  if (is.null(mb)) {
    mprime  <- spt[1] + (spt[2] - spt[1]) * 0.25
    Mprime  <- spt[1] + (spt[2] - spt[1]) * 0.75
    gmprime <- linkfun(mprime)
    gMprime <- linkfun(Mprime)
    mid_spt <- (gmprime + gMprime) / 2
    mb <- c(mid_spt, rep(0, p-1))
  }
  else if (length(mb) != p) stop("length(mb) must match the number of covariates.")

  if (is.null(Sb)) {
    mprime  <- spt[1] + (spt[2] - spt[1]) * 0.25
    Mprime  <- spt[1] + (spt[2] - spt[1]) * 0.75
    gmprime <- linkfun(mprime)
    gMprime <- linkfun(Mprime)
    sdX     <- c(apply(as.matrix(X[, -1], nrow=n), 2, sd))
    Sbvec   <- (gMprime - gmprime)^2 * c(100, (gamma / (2 * sdX))^2) # Re-scaling on the linear-predictor scale
    Sb <- sqrt(Sbvec) # Entries corresponding to sqrt(diagonal) of beta prior cov
  } else if (length(Sb) != p) stop("length(Sb) must match the number of betas.")
  else if (any(Sb <= 0)) stop("Sb must be positive definite.")

  # Start MCMC loop
  for(itr in seq_len(iter)){
    if (burning && itr > burnin) burning <- FALSE

    # Optimization to find the mode
    result <- tryCatch({
      optim(par = beta, fn = logpost_beta, linkinv = linkinv, z = z, X = X,
            atoms = z.tld, jumps  = J.tld, mu_beta = mb,
            sigma_beta = Sb, h = h,
            control = list(fnscale = -1), hessian = TRUE)
    }, error = function(e) {
      message(e)
      return(NULL)
    })  # If error occurs, return NULL

    # Check if `optim()` failed
    if (!is.null(result) && is.finite(det(result$hessian))) {
      beta_mode_ <- result$par

      # Safely compute Hessian inverse
      beta_cov_ <- tryCatch({
        -solve(result$hessian)
      }, error = function(e) {
        message(e)
        return(NULL)
      })  # If inversion fails, return NULL

      # Only proceed if covariance matrix is valid and positive definite
      if (!is.null(beta_cov_)) {
        # Check positive semi-definiteness
        if (all(eigen(beta_cov_, symmetric = TRUE, only.values = TRUE)$values >=
                -sqrt(.Machine$double.eps))) {

          # Force positive definiteness if needed
          beta_cov_ <- as.matrix(Matrix::nearPD(beta_cov_)$mat)

          # Sample new beta from proposal distribution
          beta_ <- as.numeric(rmvnorm(n = 1, mean = beta_mode_, sigma = beta_cov_))

          mean_z_ <- linkinv(X %*% beta_)

          if (min(mean_z_) >= min(z.tld) && max(mean_z_) <= max(z.tld)) {

            # Compute log proposal values
            pr_logprop_beta <- dmvnorm(x = beta_, mean = beta_mode, sigma = beta_cov, log = TRUE)
            cr_logprop_beta <- dmvnorm(x = beta, mean = beta_mode_, sigma = beta_cov_, log = TRUE)
            cr_logpost_beta <- logpost_beta(beta = beta, linkinv = linkinv, z = z, X = X, atoms = z.tld,
                                            jumps  = J.tld, mu_beta = mb,
                                            sigma_beta = Sb, h = h)
            pr_logpost_beta <- logpost_beta(beta = beta_, linkinv = linkinv, z = z, X = X, atoms = z.tld,
                                            jumps = J.tld, mu_beta = mb,
                                            sigma_beta = Sb, h = h_)

            # Metropolis-Hastings acceptance step
            log_acc_prob <- pr_logpost_beta - cr_logpost_beta + cr_logprop_beta - pr_logprop_beta
            if (log(runif(1)) < log_acc_prob) {
              beta <- beta_
              beta_mode <- beta_mode_
              beta_cov <- beta_cov_
              meanY_x <- mean_z <- mean_z_
              if (!burning) count1 <- count1 + 1
            }
          } else message("min(z_) not within z.tld bounds")
        } else message("beta covariance not positive semi-definite")
      } else message("beta covariance is null")
    } else message("either optim failed or det(hessian) is infinite")
    # If `optim()` fails, keep beta as it is and continue to next step.

    # theta update ------------------------------------
    theta <- getTheta(
      spt = z.tld,
      f0  = J.tld,
      mu  = meanY_x,
      sampprobs  = NULL,
      ySptIndex  = NULL,
      thetaStart = theta
    )$theta

    if (any(theta > 50)) {
      idx <- which(theta > 50)
      theta[idx] <- 50
      warning("Capped some values of Theta at 50. Consider increasing M.")
    }

    # u update ----------------------------------------
    u <- sampler_u(u, zstar, nstar, theta, alpha, delta, min_y, max_y, eps, h)

    # CRM update --------------------------------------
    crm_star <- crm_sampler(M, u, zstar, nstar, theta, alpha, min_y, max_y, eps, h)
    z.tld_star <- c(crm_star$RL, crm_star$zstar)
    J.tld_star <- c(crm_star$RJ, crm_star$Jstar)

    crm_2 <- crm_sampler(M, u, zstar, nstar, theta, alpha, min_y, max_y, eps, h)
    z.tld_2 <- c(crm_2$RL, crm_2$zstar)
    J.tld_2 <- c(crm_2$RJ, crm_2$Jstar)

    if(min(meanY_x) >= min(z.tld_star) && max(meanY_x) <= max(z.tld_star)){
      # MH step
      theta_star <- getTheta(
        spt = z.tld_star,
        f0  = J.tld_star,
        mu  = meanY_x,
        sampprobs  = NULL,
        ySptIndex  = NULL,
        thetaStart = theta
      )$theta

      if (any(theta_star > 50)) {
        idx <- which(theta_star > 50)
        theta_star[idx] <- 50
        warning("Capped some values of Theta at 50. Consider increasing M.")
      }

      crm_star_2 <- crm_sampler(M, u, zstar, nstar, theta_star, alpha, min_y, max_y, eps, h)
      z.tld_star_2 <- c(crm_star_2$RL, crm_star_2$zstar)
      J.tld_star_2 <- c(crm_star_2$RJ, crm_star_2$Jstar)

      b_thstar_crmstar <- b_theta(theta_star, z.tld_star, J.tld_star)
      b_th_crm         <- b_theta(theta,      z.tld,      J.tld)
      b_thstar_crm     <- b_theta(theta_star, z.tld,      J.tld)
      b_th_crmstar     <- b_theta(theta,      z.tld_star, J.tld_star)

      b_th_crmstar2    <- b_theta(theta_ref,      z.tld_star_2, J.tld_star_2)
      b_thstar_crmstar2 <- b_theta(theta_star, z.tld_star_2, J.tld_star_2)
      b_th_crm2        <- b_theta(theta_ref,      z.tld_2,      J.tld_2)
      b_thstar_crm2    <- b_theta(theta_star, z.tld_2,      J.tld_2)

      log_r <- sum(
        (theta_star - theta) * z -
          b_thstar_crmstar + b_th_crm -
          b_thstar_crm + b_th_crmstar
      ) +
        sum(
          -b_th_crmstar2 + b_thstar_crmstar2 +
            b_th_crm2 - b_thstar_crm2
        )

      if(log(runif(1)) < log_r){
        RL <- crm_star$RL
        RJ <- crm_star$RJ
        zstar <- crm_star$zstar
        Jstar <- crm_star$Jstar
        z.tld <- c(RL, zstar)
        J.tld <- c(RJ, Jstar)
        theta <- theta_star
        if (!burning) count2 <- count2 + 1
      }
    }

    # z update ------------------------------------------------------------------
    z <- z_sampler_unifK(y, c0, z.tld, J.tld, theta, min_y, max_y, eps)

    # zstar and nstar update ----------------------------------------------------
    resampled_z <- resample_zstar(z)
    zstar <- resampled_z$zstar
    nstar <- resampled_z$nstar

    theta0 <- getTheta(
      spt = z.tld,
      f0  = J.tld,
      mu  = mu0,
      sampprobs  = NULL,
      ySptIndex  = NULL,
      thetaStart = NULL
    )$theta

    if (any(theta0 > 50)) {
      idx <- which(theta0 > 50)
      theta0[idx] <- 50
      warning("Capped some values of Theta at 50. Consider increasing M.")
    }

    temp <- exp(theta0 * z.tld - max(theta0 * z.tld))
    Jtilt <- temp * J.tld
    W <- sum(Jtilt)
    Jtld_0 <- Jtilt / W
    if (robust) {
    J.tld <- Jtilt
    J.tld_ll <- Jtld_0
    } else {
      J.tld_ll <- Jtld_0
    }

    lnlik <- loglik(linkinv = linkinv, z = z, X = X, beta = beta, atoms = z.tld, jumps = J.tld_ll)

    # Storing MCMC simulations --------------------------------------------------
    if (itr > burnin && (itr - burnin) %% thin == 0) {
      j <- (itr - burnin) / thin

      beta_samples[j, ] <- beta
      theta_samples[j, ] <- theta
      z_samples[j, ] <- z
      u_samples[j, ] <- u
      crm_samples[[j]] <- list(z.tld = z.tld, J.tld = J.tld)
      lnlik_samples[j] <- lnlik
    }
  }

  crm_samples <- data.frame(
  z.tld = I(lapply(crm_samples, `[[`, "z.tld")),
  J.tld = I(lapply(crm_samples, `[[`, "J.tld"))
  )
  samples <- list(z = z_samples, beta = beta_samples,
                  crm = crm_samples, theta = theta_samples)

  list(samples         = samples,
       mb              = mb,
       Sb              = Sb,
       spt             = spt,
       mu0             = mu0,
       beta_acceptance = count1 / (iter - burnin),
       crm_acceptance  = count2 / (iter - burnin))
}

#' Control Arguments for CDPGLM
#'
#' Creates a list of control parameters used by \code{\link{cdpglm}} and \code{\link{cdpglmFit}}.
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
#' @param rhoStart Numeric. Starting value for the copula dependence parameter
#' \eqn{\rho}. Defaults to \code{0.5}.
#' @param corr Character string specifying copula correlation structure.
#' Options are \code{"ex"} for exchangeable correlation and \code{"ar1"} for AR(1) correlation.
#' @param rho_proposal_sd Numeric. Proposal standard deviation for the random-walk
#' update of \eqn{\rho} on the logit scale. Defaults to \code{0.1}.
#' @param rho_prior_shape Numeric vector of length two. Parameters for the
#' Beta prior on \eqn{\rho}. Defaults to \code{c(8, 2)}.
#' @param betaStart Initial value for \code{beta}. Defaults to NULL.
#' @param varbetaStart Initial value for the proposal covariance of \code{beta}. Defaults to NULL.
#' @param thetaStart Initial value for \code{theta}. Defaults to NULL.
#' @param crmStart Initial value for the CRM state, given as a list with components
#' \code{z.tld} and \code{J.tld}. Defaults to NULL.
#' @param seed Random seed. Defaults to NULL.
#' @param robust Provides numerical stability. Defaults to \code{TRUE}. See details for more.
#'
#' @return An object of class \code{"cdpglmControl"}; a list of control
#' parameters passed to the internal CDPGLM fitting routine.
#'
#' @examples
#' ctrl <- cdpglm.control(
#'   burnin = 100,
#'   thin = 5,
#'   save = 500,
#'   spt = c(0, 1),
#'   corr = "ex",
#'   seed = 123
#' )
#'
#' @export
cdpglm.control <- function(burnin = 100, thin = 10, save = 1000,
                           mb = NULL, Sb = NULL,
                           M = 50, alpha = 1, delta = 2, c0 = NULL,
                           gamma = 1, mu0 = NULL, spt = NULL, eps = 1e-6,
                           rhoStart = 0.5, corr = "ex",
                           rho_proposal_sd = 0.1,
                           rho_prior_shape = c(8, 2),
                           betaStart = NULL, varbetaStart = NULL,
                           thetaStart = NULL, crmStart = NULL,
                           seed = NULL, robust = TRUE) {
  if (burnin < 0 || floor(burnin) != burnin) {
    stop("Number of burn-in samples must be an integer >= 0.", call. = FALSE)
  }
  if (thin < 1 || floor(thin) != thin) {
    stop("Thin must be an integer >= 1.", call. = FALSE)
  }
  if (save < 1 || floor(save) != save) {
    stop("Number of saved iterations must be an integer >= 1.", call. = FALSE)
  }

  if (!(corr %in% c("ex", "ar1"))) stop("corr must be one of c('ex', 'ar1')")

  if (!is.numeric(rhoStart) || length(rhoStart) != 1L ||
      !is.finite(rhoStart) || rhoStart <= 0 || rhoStart >= 1) {
    stop("rhoStart must be a numeric value in (0, 1).", call. = FALSE)
  }

  if (!is.numeric(rho_proposal_sd) || length(rho_proposal_sd) != 1L ||
      !is.finite(rho_proposal_sd) || rho_proposal_sd <= 0) {
    stop("rho_proposal_sd must be a positive numeric value.", call. = FALSE)
  }

  if (!is.numeric(rho_prior_shape) || length(rho_prior_shape) != 2L ||
      any(!is.finite(rho_prior_shape)) || any(rho_prior_shape <= 0)) {
    stop("rho_prior_shape must be a numeric vector of length 2 with positive entries.", call. = FALSE)
  }

  ctrl <- list(
    burnin = burnin,
    thin = thin,
    save = save,
    mb = mb,
    Sb = Sb,
    M = M,
    alpha = alpha,
    delta = delta,
    c0 = c0,
    gamma = gamma,
    mu0 = mu0,
    spt = spt,
    eps = eps,
    rhoStart = rhoStart,
    corr = corr,
    rho_proposal_sd = rho_proposal_sd,
    rho_prior_shape = rho_prior_shape,
    betaStart = betaStart,
    varbetaStart = varbetaStart,
    thetaStart = thetaStart,
    crmStart = crmStart,
    seed = seed,
    robust = robust
  )

  class(ctrl) <- "cdpglmControl"
  ctrl
}

#' Main MCMC Fitting Function for CDPGLM
#'
#' Performs MCMC to fit the CDPGLM model.
#'
#' @param formula Model formula supplied to \code{\link{cdpglm}}.
#' @param data Data frame containing the variables appearing in \code{formula}.
#' @param X Numeric design matrix.
#' @param y Numeric response vector.
#' @param group_index Vector identifying the dependence group for each observation.
#' @param link Link object containing \code{linkfun}, \code{linkinv}, and \code{mu.eta}.
#' @param spt Numeric vector of length two giving the lower and upper bounds
#' of the response support.
#' @param mu0 Numeric target mean used to identify the reference distribution.
#' @param init List of initial values for the MCMC algorithm.
#' @param cdpglmControl Object of class \code{"cdpglmControl"} containing MCMC,
#' tuning, and copula-control parameters.
#' @param thetaControl Object of class \code{"thetaControl"} containing control
#' arguments for the theta update procedure.
#'
#' @return
#' A list containing posterior samples, acceptance rates, prior settings,
#' support information, copula settings, and model specification information
#' used by \code{\link{cdpglm}}.
#'
#' @keywords internal
cdpglmFit <- function(formula, data, X, y, group_index,
                      link,
                      spt, mu0, init,
                      cdpglmControl, thetaControl) {

  # Extract cdpglmControl parameters
  burnin          <- cdpglmControl$burnin
  thin            <- cdpglmControl$thin
  save            <- cdpglmControl$save
  mb              <- cdpglmControl$mb
  Sb              <- cdpglmControl$Sb
  M               <- cdpglmControl$M
  alpha           <- cdpglmControl$alpha
  delta           <- cdpglmControl$delta
  c0              <- cdpglmControl$c0
  gamma           <- cdpglmControl$gamma
  eps             <- cdpglmControl$eps
  seed            <- cdpglmControl$seed
  robust          <- cdpglmControl$robust
  corr            <- cdpglmControl$corr
  rho_proposal_sd <- cdpglmControl$rho_proposal_sd
  rho_prior_shape <- cdpglmControl$rho_prior_shape

  # Support
  if (is.unsorted(spt)) spt <- sort(spt)
  min_y <- spt[1]
  max_y <- spt[2]

  # MCMC Initialization
  n    <- length(y)
  X    <- as.matrix(X, nrow = n)
  p    <- ncol(X)
  iter <- burnin + thin * save

  if (!is.null(seed)) set.seed(seed)
  if (is.null(c0)) c0 <- c0_silverman(y) / 4

  linkfun <- link$linkfun
  linkinv <- link$linkinv

  beta_samples  <- matrix(NA, nrow = save, ncol = p)
  beta_names <- colnames(X)
  if (is.null(beta_names)) beta_names <- paste0("beta_", seq_len(p) - 1L)
  beta_names[1] <- "Intercept"
  colnames(beta_samples) <- beta_names
  theta_samples <- matrix(NA, nrow = save, ncol = n)
  h_samples     <- matrix(NA, nrow = save, ncol = n)
  u_samples     <- z_samples <- matrix(NA, nrow = save, ncol = n)
  rho_samples   <- numeric(save)
  crm_samples   <- vector("list", save)
  lnlik_samples <- numeric(save)

  beta <- as.numeric(init$beta)
  #beta_samples[1, ] <- beta
  beta_mode <- beta
  beta_cov  <- init$varbeta

  meanY_x <- as.numeric(linkinv(X %*% beta))

  z.tld <- z_tld <- as.numeric(init$crm$z.tld)
  J.tld <- J_tld <- as.numeric(init$crm$J.tld)

  theta <- theta_ref <- as.numeric(init$theta)
  #theta_samples[1, ] <- theta
  rho <- as.numeric(init$rho)
  #rho_samples[1] <- rho <- as.numeric(init$rho)

  h <- log_copula_contribution_by_obs(
    y = y,
    group_index = group_index,
    rho = rho,
    corr = corr,
    crm.atoms = z.tld,
    crm.jumps = J.tld,
    theta = theta,
    c0 = c0,
    min_y = min_y,
    max_y = max_y
  )
  # h_samples[1, ] <- h

  # crm_samples[[1]] <- list(z.tld = z.tld, J.tld = J.tld)

  ord <- order(J_tld)[1:M]
  RL <- z_tld[ord]
  RJ <- J_tld[ord]

  z <- z_sampler_unifK(y, c0, z.tld, J.tld, theta, min_y, max_y, eps)
  # z_samples[1, ] <- z

  btheta <- b_theta(theta, z.tld, J.tld)
  T_vec <- exp(btheta)
  u <- rgamma(n, shape = 1, rate = T_vec)
  # u_samples[1, ] <- u

  resampled_z <- resample_zstar(z)
  zstar <- resampled_z$zstar
  nstar <- resampled_z$nstar

  count_beta <- count_crm <- count_rho <- 0
  burning <- TRUE

  theta0 <- gldrm:::getTheta(
    spt = z.tld,
    f0  = J.tld,
    mu  = mu0,
    sampprobs = NULL,
    ySptIndex = NULL,
    thetaStart = NULL
  )$theta

  if (any(theta0 > 50)) {
    theta0[theta0 > 50] <- 50
    warning("Capped some values of theta0 at 50. Consider increasing M.")
  }

  temp <- exp(theta0 * z.tld - max(theta0 * z.tld))
  Jtilt <- temp * J.tld
  Jtld_0 <- Jtilt / sum(Jtilt)

  if (robust) {
    J.tld <- Jtilt
    J.tld_ll <- Jtld_0
  } else {
    J.tld_ll <- Jtld_0
  }

  lnlik <-  loglik(
    linkinv = linkinv,
    z = z,
    X = X,
    beta = beta,
    atoms = z.tld,
    jumps = J.tld_ll
  )
  # lnlik_samples[1] <- lnlik

  # Beta prior
  if (is.null(mb)) {
    mprime  <- spt[1] + (spt[2] - spt[1]) * 0.25
    Mprime  <- spt[1] + (spt[2] - spt[1]) * 0.75
    gmprime <- linkfun(mprime)
    gMprime <- linkfun(Mprime)
    mid_spt <- (gmprime + gMprime) / 2
    mb <- c(mid_spt, rep(0, p - 1))
  } else if (length(mb) != p) {
    stop("length(mb) must match the number of covariates.")
  }

  if (is.null(Sb)) {
    mprime  <- spt[1] + (spt[2] - spt[1]) * 0.25
    Mprime  <- spt[1] + (spt[2] - spt[1]) * 0.75
    gmprime <- linkfun(mprime)
    gMprime <- linkfun(Mprime)
    sdX     <- c(apply(as.matrix(X[, -1], nrow = n), 2, sd))
    Sbvec   <- (gMprime - gmprime)^2 * c(100, (gamma / (2 * sdX))^2)
    Sb <- sqrt(Sbvec)
  } else if (length(Sb) != p) {
    stop("length(Sb) must match the number of betas.")
  } else if (any(Sb <= 0)) {
    stop("Sb must be positive definite.")
  }

  # Start MCMC loop
  for (itr in seq_len(iter)) {
    if (burning && itr > burnin) burning <- FALSE

    # Optimization for beta mode
    result <- tryCatch({
      optim(
        par = beta,
        fn = logpost_beta,
        linkinv = linkinv,
        z = z,
        X = X,
        atoms = z.tld,
        jumps = J.tld,
        mu_beta = mb,
        sigma_beta = Sb,
        h = h,
        control = list(fnscale = -1),
        hessian = TRUE
      )
    }, error = function(e) {
      message(e)
      NULL
    })

    # Check if `optim()` failed
    if (!is.null(result) && is.finite(det(result$hessian))) {
      beta_mode_ <- result$par

      # Safely compute Hessian inverse
      beta_cov_ <- tryCatch({
        -solve(result$hessian)
      }, error = function(e) {
        message(e)
        NULL
      }) # If inversion fails, return NULL

      # Only proceed if covariance matrix is valid and positive definite
      if (!is.null(beta_cov_)) {
        # Check positive semi-definiteness
        if (all(eigen(beta_cov_, symmetric = TRUE, only.values = TRUE)$values >=
                -sqrt(.Machine$double.eps))) {

          # Force positive definiteness if needed
          beta_cov_ <- as.matrix(Matrix::nearPD(beta_cov_)$mat)

          # Sample new beta from proposal distribution
          beta_ <- as.numeric(rmvnorm(n = 1, mean = beta_mode_, sigma = beta_cov_))

          mean_z_ <- as.numeric(linkinv(X %*% beta_))

          if (min(mean_z_) >= min(z.tld) && max(mean_z_) <= max(z.tld)) {

            # Compute log-proposal values
            pr_logprop_beta <- dmvnorm(beta_, mean = beta_mode,  sigma = beta_cov,  log = TRUE)
            cr_logprop_beta <- dmvnorm(beta,  mean = beta_mode_, sigma = beta_cov_, log = TRUE)

            # Compute log-posterior values
            cr_logpost_beta <- logpost_beta(
              beta = beta,
              linkinv = linkinv,
              z = z,
              X = X,
              atoms = z.tld,
              jumps = J.tld,
              mu_beta = mb,
              sigma_beta = Sb,
              h = h
            )

            theta_ <- gldrm:::getTheta(
              spt = z.tld,
              f0  = J.tld,
              mu  = mean_z_,
              sampprobs = NULL,
              ySptIndex = NULL,
              thetaStart = theta
            )$theta

            h_ <- log_copula_contribution_by_obs(
              y = y,
              group_index = group_index,
              rho = rho,
              corr = corr,
              crm.atoms = z.tld,
              crm.jumps = J.tld,
              theta = theta_,
              c0 = c0,
              min_y = min_y,
              max_y = max_y
            )

            pr_logpost_beta <- logpost_beta(
              beta = beta_,
              linkinv = linkinv,
              z = z,
              X = X,
              atoms = z.tld,
              jumps = J.tld,
              mu_beta = mb,
              sigma_beta = Sb,
              h = h_
            )

            log_acc_prob <- pr_logpost_beta - cr_logpost_beta +
              cr_logprop_beta - pr_logprop_beta

            if (log(runif(1)) < log_acc_prob) {
              beta <- beta_
              beta_mode <- beta_mode_
              beta_cov <- beta_cov_
              meanY_x <- mean_z_
              theta <- theta_
              h <- h_
              if (!burning) count_beta <- count_beta + 1
            }
          }
        }
      }
    }

    # Theta update
    theta <- gldrm:::getTheta(
      spt = z.tld,
      f0  = J.tld,
      mu  = meanY_x,
      sampprobs = NULL,
      ySptIndex = NULL,
      thetaStart = theta
    )$theta

    if (any(theta > 50)) {
      theta[theta > 50] <- 50
      warning("Capped some values of theta at 50.")
    }

    # h update
    h <- log_copula_contribution_by_obs(
      y = y,
      group_index = group_index,
      rho = rho,
      corr = corr,
      crm.atoms = z.tld,
      crm.jumps = J.tld,
      theta = theta,
      c0 = c0,
      min_y = min_y,
      max_y = max_y
    )

    # u update
    u <- sampler_u(u, zstar, nstar, theta, alpha, delta, min_y, max_y, eps, h)

    # CRM update
    crm_star <- crm_sampler(M, u, zstar, nstar, theta, alpha, min_y, max_y, eps, h)
    z.tld_star <- c(crm_star$RL, crm_star$zstar)
    J.tld_star <- c(crm_star$RJ, crm_star$Jstar)

    crm_2 <- crm_sampler(M, u, zstar, nstar, theta, alpha, min_y, max_y, eps, h)
    z.tld_2 <- c(crm_2$RL, crm_2$zstar)
    J.tld_2 <- c(crm_2$RJ, crm_2$Jstar)

    if (min(meanY_x) >= min(z.tld_star) && max(meanY_x) <= max(z.tld_star)) {
      # MH step
      theta_star <- gldrm:::getTheta(
        spt = z.tld_star,
        f0  = J.tld_star,
        mu  = meanY_x,
        sampprobs = NULL,
        ySptIndex = NULL,
        thetaStart = theta
      )$theta

      if (any(theta_star > 50)) {
        theta_star[theta_star > 50] <- 50
        warning("Capped some values of theta_star at 50.")
      }

      h_star <- log_copula_contribution_by_obs(
        y = y,
        group_index = group_index,
        rho = rho,
        corr = corr,
        crm.atoms = z.tld_star,
        crm.jumps = J.tld_star,
        theta = theta_star,
        c0 = c0,
        min_y = min_y,
        max_y = max_y
      )

      crm_star_2 <- crm_sampler(M, u, zstar, nstar, theta_star, alpha, min_y, max_y, eps, h_star)
      z.tld_star_2 <- c(crm_star_2$RL, crm_star_2$zstar)
      J.tld_star_2 <- c(crm_star_2$RJ, crm_star_2$Jstar)

      b_thstar_crmstar <- b_theta(theta_star, z.tld_star, J.tld_star)
      b_th_crm         <- b_theta(theta,      z.tld,      J.tld)
      b_thstar_crm     <- b_theta(theta_star, z.tld,      J.tld)
      b_th_crmstar     <- b_theta(theta,      z.tld_star, J.tld_star)

      b_thref_crmstar2  <- b_theta(theta_ref,  z.tld_star_2, J.tld_star_2)
      b_thstar_crmstar2 <- b_theta(theta_star, z.tld_star_2, J.tld_star_2)
      b_thref_crm2      <- b_theta(theta_ref,  z.tld_2, J.tld_2)
      b_thstar_crm2     <- b_theta(theta_star, z.tld_2, J.tld_2)

      log_r <- sum(
        (theta_star - theta) * z +
          (h_star - h) -
          b_thstar_crmstar + b_th_crm -
          b_thstar_crm + b_th_crmstar
      ) +
        sum(
          -b_thref_crmstar2 + b_thstar_crmstar2 +
            b_thref_crm2 - b_thstar_crm2
        )

      if (log(runif(1)) < log_r) {
        RL <- crm_star$RL
        RJ <- crm_star$RJ
        zstar <- crm_star$zstar
        Jstar <- crm_star$Jstar
        z.tld <- c(RL, zstar)
        J.tld <- c(RJ, Jstar)
        theta <- theta_star
        h <- h_star
        if (!burning) count_crm <- count_crm + 1
      }
    }

    # z update
    z <- z_sampler_unifK(y, c0, z.tld, J.tld, theta, min_y, max_y, eps)

    # zstar and nstar update
    resampled_z <- resample_zstar(z)
    zstar <- resampled_z$zstar
    nstar <- resampled_z$nstar

    theta0 <- gldrm:::getTheta(
      spt = z.tld,
      f0  = J.tld,
      mu  = mu0,
      sampprobs = NULL,
      ySptIndex = NULL,
      thetaStart = NULL
    )$theta

    if (any(theta0 > 50)) {
      theta0[theta0 > 50] <- 50
      warning("Capped some values of theta0 at 50.")
    }

    # rho update
    rho_prop <- plogis(qlogis(rho) + rnorm(1, 0, rho_proposal_sd))
    # log posterior
    logpost_current <- sum(h) +
      dbeta(rho, rho_prior_shape[1], rho_prior_shape[2], log = TRUE)

    h_prop <- log_copula_contribution_by_obs(
      y = y,
      group_index = group_index,
      rho = rho_prop,
      corr = corr,
      crm.atoms = z.tld,
      crm.jumps = J.tld,
      theta = theta,
      c0 = c0,
      min_y = min_y,
      max_y = max_y
    )

    rho_pr1 <- rho_prior_shape[1]
    rho_pr2 <- rho_prior_shape[2]
    logpost_prop <- log_post_rho(rho_prop, h_prop, rho_pr1, rho_pr2)

    log_jacobian <- log(rho_prop * (1 - rho_prop)) -
      log(rho * (1 - rho))

    log_rho <- logpost_prop - logpost_current + log_jacobian

    if (log(runif(1)) < log_rho) {
      if (!burning) count_rho <- count_rho + 1
      rho <- rho_prop
      h <- h_prop
    }

    temp <- exp(theta0 * z.tld - max(theta0 * z.tld))
    Jtilt <- temp * J.tld
    Jtld_0 <- Jtilt / sum(Jtilt)

    if (robust) {
      J.tld <- Jtilt
      J.tld_ll <- Jtld_0
    } else {
      J.tld_ll <- Jtld_0
    }

    lnlik <- loglik(
      linkinv = linkinv,
      z = z,
      X = X,
      beta = beta,
      atoms = z.tld,
      jumps = J.tld_ll
    )

    # Storage
    if (itr > burnin && (itr - burnin) %% thin == 0) {
      j <- (itr - burnin) / thin

      beta_samples[j, ] <- beta
      theta_samples[j, ] <- theta
      z_samples[j, ] <- z
      u_samples[j, ] <- u
      h_samples[j, ] <- h
      rho_samples[j] <- rho
      crm_samples[[j]] <- list(z.tld = z.tld, J.tld = J.tld)
      lnlik_samples[j] <- lnlik
    }
  }

  crm_samples <- data.frame(
    z.tld = I(lapply(crm_samples, `[[`, "z.tld")),
    J.tld = I(lapply(crm_samples, `[[`, "J.tld"))
  )

  samples <- list(
    z = z_samples,
    beta = beta_samples,
    theta = theta_samples,
    h = h_samples,
    rho = rho_samples,
    crm = crm_samples
  )

  list(
    samples = samples,
    mb = mb,
    Sb = Sb,
    spt = spt,
    mu0 = mu0,
    beta_acceptance = count_beta / (iter - burnin),
    crm_acceptance = count_crm / (iter - burnin),
    rho_acceptance = count_rho / (iter - burnin)
  )
}

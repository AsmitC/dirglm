#' Control arguments for the cdpglm algorithm.
#'
#' @export
cdpglm.control <- function(burnin = 100, thin = 10, save = 1000,
                           mb = NULL, Sb = NULL,
                           M = 50, alpha = 1, delta = 2, c0 = NULL,
                           gamma = 1, mu0 = NULL, spt = NULL, eps = 1e-6,
                           robust = TRUE,
                           rhoStart = 0.5, corr = c("ex", "ar1"),
                           rho_proposal_sd = 0.1,
                           rho_prior_shape = c(8, 2),
                           betaStart = NULL, varbetaStart = NULL,
                           thetaStart = NULL, crmStart = NULL,
                           seed = NULL) {
  if (burnin < 0 || floor(burnin) != burnin) {
    stop("Number of burn-in samples must be an integer >= 0.", call. = FALSE)
  }
  if (thin < 1 || floor(thin) != thin) {
    stop("Thin must be an integer >= 1.", call. = FALSE)
  }
  if (save < 1 || floor(save) != save) {
    stop("Number of saved iterations must be an integer >= 1.", call. = FALSE)
  }

  corr <- match.arg(corr)

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
    robust = robust,
    rhoStart = rhoStart,
    corr = corr,
    rho_proposal_sd = rho_proposal_sd,
    rho_prior_shape = rho_prior_shape,
    betaStart = betaStart,
    varbetaStart = varbetaStart,
    thetaStart = thetaStart,
    crmStart = crmStart,
    seed = seed
  )

  class(ctrl) <- "cdpglmControl"
  ctrl
}

#' Main MCMC function for cdpglm
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

  beta_samples  <- matrix(NA, nrow = iter, ncol = p)
  theta_samples <- matrix(NA, nrow = iter, ncol = n)
  h_samples     <- matrix(NA, nrow = iter, ncol = n)
  u_samples     <- z_samples <- matrix(NA, nrow = iter, ncol = n)
  rho_samples   <- numeric(iter)
  crm_samples   <- list()
  lnlik_samples <- numeric(iter)

  beta_samples[1, ] <- beta <- as.numeric(init$beta)
  beta_mode <- beta
  beta_cov  <- init$varbeta

  meanY_x <- as.numeric(linkinv(X %*% beta))

  z.tld <- z_tld <- as.numeric(init$crm$z.tld)
  J.tld <- J_tld <- as.numeric(init$crm$J.tld)

  theta_samples[1, ] <- theta <- theta_ref <- as.numeric(init$theta)
  rho_samples[1] <- rho <- as.numeric(init$rho)

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
  h_samples[1, ] <- h

  crm_samples[[1]] <- list(z.tld = z.tld, J.tld = J.tld)

  ord <- order(J_tld)[1:M]
  RL <- z_tld[ord]
  RJ <- J_tld[ord]

  z_samples[1, ] <- z <- z_sampler_unifK(
    y, c0, z.tld, J.tld, theta, min_y, max_y, eps
  )

  btheta <- b_theta(theta, z.tld, J.tld)
  T_vec <- exp(btheta)
  u_samples[1, ] <- u <- rgamma(n, shape = 1, rate = T_vec)

  resampled_z <- resample_zstar(z)
  zstar <- resampled_z$zstar
  nstar <- resampled_z$nstar

  count_beta <- count_crm <- count_rho <- 0

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

  lnlik_samples[1] <- loglik(
    linkinv = linkinv,
    z = z,
    X = X,
    beta = beta,
    atoms = z.tld,
    jumps = J.tld_ll
  )

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
  for (itr in 2:iter) {
    if (itr %% 50 == 0) cat("\nStarting iteration:", itr)

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
              count_beta <- count_beta + 1
              beta <- beta_
              beta_mode <- beta_mode_
              beta_cov <- beta_cov_
              meanY_x <- mean_z_
              theta <- theta_
              h <- h_
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
        count_crm <- count_crm + 1
        RL <- crm_star$RL
        RJ <- crm_star$RJ
        zstar <- crm_star$zstar
        Jstar <- crm_star$Jstar
        z.tld <- c(RL, zstar)
        J.tld <- c(RJ, Jstar)
        theta <- theta_star
        h <- h_star
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
      count_rho <- count_rho + 1
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
    z_samples[itr, ] <- z
    u_samples[itr, ] <- u
    beta_samples[itr, ] <- beta
    theta_samples[itr, ] <- theta
    h_samples[itr, ] <- h
    rho_samples[itr] <- rho
    crm_samples[[itr]] <- list(z.tld = z.tld, J.tld = J.tld)
    lnlik_samples[itr] <- lnlik
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
    beta_acceptance = count_beta / iter,
    crm_acceptance = count_crm / iter,
    rho_acceptance = count_rho / iter
  )
}
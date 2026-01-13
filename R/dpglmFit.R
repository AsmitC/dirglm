#' Control arguments for the \code{dpglm} algorithm.
#'
#' This function returns control arguments for the \code{dpglm} algorithm.
#' Each argument has a default value, which will be used unless a different
#' value is provided by the user.
#'
#' @param burnin Number of burn-in MCMC iterations. Defaults to 100.
#' @param thin Factor by which to thin MCMC iterations. Defaults to 10.
#' @param save Number of MCMC samples to return. Defaults to 1000.
#' @param rho MCMC update step size. Either a single number or a vector matching
#' the length of \code{beta}.
#' @param mb Prior mean for beta. Defaults to a p-length vector whose entries are all 0.
#' @param Sb Vector containing the diagonal entries in the prior variance-covariance matrix for beta.
#' @param gamma Shrinkage parameter for the (default) prior variance on \code{beta}.
#' Defaults to 1. Will not be used if \code{Sb} is specified in \code{dirglm}.
#' @param mu0 Mean of the reference distribution. The reference distribution is
#' not unique unless its mean is restricted to a specific value. This value can
#' be any number within the range of observed values, but values near the boundary
#' may cause numerical instability. This is an optional argument with \code{mean(y)}
#' being the default value.
#' @param eps Padding for the theoretical support. Defaults to 1e-6.
#' @param spt Theoretical support of the response variable. Defaults to the
#' empirical distribution of \code{y}.
#' @param betaStart Initial value for the regression coefficients \code{beta}.
#' Defaults to the output obtained by fitting \code{gldrm}.
#' @param crmStart Initial value for the reference distribution.
#' Defaults to the output obtained by fitting \code{gldrm}.
#' @param seed Random seed. Defaults to NULL.
#'
#' @return Object of S3 class "dpglmControl"
#'
#' @export
dpglm.control <- function(burnin=100, thin=10, save=1000,
                          mb=NULL, Sb=NULL,
                          M=20, alpha=1, delta=2, c0=NULL,
                          gamma=1, mu0=NULL, spt=NULL, H=NULL, flag=c("dp", "cdp", "ods"), eps=1e-6,
                          betaStart=NULL, varbetaStart=NULL, thetaStart=NULL, crmStart=NULL,
                          seed=NULL) #, add H, flag 
{
  if (burnin < 0 || floor(burnin) != burnin) stop("Number of burn-in samples must be an integer >= 0")
  if (thin   < 1 || floor(thin)   != thin)   stop("Thin must be an integer >= 1")
  if (save   < 1 || floor(save)   != save)   stop("Number of saved iterations must be an integer >= 1")
  #if (!(rho <= 1 & rho > 0))                 stop("rho must lie in (0, 1]")
  #if (length(spt) != 2)                      stop("Support must be a vector of length 2")
  ctrl <- list(burnin       = burnin,
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
               H            = H,
               flag         = flag,
               eps          = eps,
               betaStart    = betaStart,
               varbetaStart = varbetaStart,
               thetaStart   = thetaStart,
               crmStart     = crmStart,
               seed         = seed)
  class(ctrl) <- "dpglmControl"
  ctrl
}


#' Main MCMC function
#' This function is called by the main \code{dirglm} function.
#' @keywords internal
dpglmFit <- function(formula, data, X, y,        # Data
                     link,                       # Link
                     mu0, spt, init,             # Specs
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
  H            <- dpglmControl$H
  flag         <- dpglmControl$flag
  eps          <- dpglmControl$eps
  seed         <- dpglmControl$seed

  #if (!is.null(H)) {
    #if (class(H) != "function") stop("H must be a.")
  #} else {
    # Set H based on flag
    # if flag is closest to  "dpglm" then make H = 0 for all inputs
    # Else, define H as some default
  #}

  if (is.null(flag)) flag <- "dpglm"
  flag <- match.arg(as.character(flag), choices = c("dpglm", "copula", "ods"))

  if (is.null(H) || flag == "dpglm") H <- function(...) 0
  else { # H is non-null and model is copula or ods
    if (class(H) != "function") stop("H must be a function.")
    # Now that H is definitely a function, handle copula/ods cases separately
  }

  if (!is.null(seed)) set.seed(seed)
  
  if(is.null(c0)) c0 <- c0_silverman(y) / 4
  
  # Extract link
  linkfun <- link$linkfun
  linkinv <- link$linkinv
  mu.eta  <- link$mu.eta

  # MCMC Initialization
  n    <- length(y)
  X    <- as.matrix(X, nrow=n)
  p    <- dim(X)[2]
  l    <- length(spt)
  iter <- burnin + thin * save

  # Extract support bounds
  if (is.unsorted(spt)) spt <- sort(spt)
  r <- diff(range(spt))
  min_y <- spt[1] - eps # A: make eps later?
  max_y <- spt[l] + eps # A: same here

  # Verify MCMC step size
  #if (!is.null(rho)) {
    #rho <- as.numeric(rho)
    #if (length(rho) == 1) rho <- rep(rho, p) # Scalar rho
    #else if (length(rho) != p) stop("length(rho) must match the number of betas.")
  #} else rho <- rep(1, p) # Defaults to no step size scaling

  beta_samples <- matrix(NA, nrow = iter, ncol = p)
  theta_samples <- matrix(NA, nrow = iter, ncol = n)
  u_samples <- z_samples <- matrix(NA, nrow = iter, ncol = n)
  crm_samples   <- list()
  lnlik_samples <- numeric(iter)

  # gldrm_fit <- gldrm(y ~ X[, -1], link = link) # E old
  beta_samples[1, ] <- beta <- as.numeric(init$beta)
  beta_mode <- beta
  #beta_cov <- gldrm_fit$varbeta # What to do about this
  beta_cov <- init$varbeta # this
  meanY_x <- linkinv(X %*% beta)
  z.tld <- z_tld <-  as.numeric(init$crm$z.tld)
  J.tld <- J_tld <- as.numeric(init$crm$J.tld)
  crm_samples[[1]] <- list(z.tld = z.tld, J.tld = J.tld)

  ord <- order(J_tld)[1:M]
  
  RL <- z_tld[ord]
  RJ <- J_tld[ord]
  #theta_samples[1, ] <- theta <- true_theta #gldrm_fit$theta %>% as.numeric() # What to do about this
  theta_samples[1, ] <- theta <- init$theta # this
  z_samples[1, ] <- z <- z_sampler_unifK(y, c0, z.tld, J.tld, theta, min_y, max_y, eps)

  btheta <- b_theta(theta, z.tld, J.tld)
  
  T_vec <- exp(btheta)
  u_samples[1, ] <- u <- rgamma(n, shape = 1, rate = T_vec)
  
  resampled_z <- resample_zstar(z)
  zstar <- resampled_z$zstar
  nstar <- resampled_z$nstar
  Jstar <- rgamma(n = length(nstar), shape = nstar, rate = 1)
  
  count1 <- count2 <- 0
  #scale_factors <- c(1.5, 2.5) %>% sqrt()

  theta0 <- gldrm:::getTheta(
    spt = z.tld,
    f0  = J.tld,
    mu  = mu0,
    sampprobs  = NULL,
    ySptIndex  = NULL,
    thetaStart = NULL
  )$theta
  #Jtld_0 <- exp(theta0 * z.tld) * J.tld / sum(exp(theta0 * z.tld) * J.tld)
  temp <- exp(theta0 * z.tld - max(theta0 * z.tld))
  # Adding H, something like
  # tilt <- H(...)
  # term <- theta0 * z.tld + tilt
  #temp <- exp(term - max(term))
  Jtld_0 <- (temp * J.tld) / sum(temp * J.tld)
  lnlik_samples[1] <- lnlik <- loglik(linkinv = linkinv, z = z, X = X, beta = beta, atoms = z.tld, jumps = Jtld_0)

  # Beta prior
  if (is.null(mb)) {
    mprime  <- spt[1] + (spt[2] - spt[1]) * 0.25
    Mprime  <- spt[l - 1] + (spt[l] - spt[l - 1]) * 0.75
    gmprime <- linkfun(mprime)
    gMprime <- linkfun(Mprime)
    mid_spt <- (gmprime + gMprime) / 2
    mb <- c(mid_spt, rep(0, p-1))
  }
  else if (length(mb) != p) stop("length(mb) must match the number of covariates.")

  Sbdiag <- TRUE
  if (is.null(Sb)) {
    mprime  <- spt[1] + (spt[2] - spt[1]) * 0.25
    Mprime  <- spt[l - 1] + (spt[l] - spt[l - 1]) * 0.75
    gmprime <- linkfun(mprime)
    gMprime <- linkfun(Mprime)
    sdX     <- c(apply(as.matrix(X[, -1], nrow=n), 2, sd))
    Sbvec   <- (gMprime - gmprime)^2 * c(100, (gamma / (2 * sdX))^2) # Re-scaling on the linear-predictor scale
    #Sb      <- diag(Sbvec)
    Sb <- Sbvec # Entries corresponding to the diagonal of beta prior cov
  } else if (length(Sb) != p) stop("length(Sb) must match the number of betas.")
  else if (any(Sb <= 0)) stop("Sb must be positive definite.")

  ### A: below few comments are really old (dirglm), probably entirely unnecessary
  #else if (!all(dim(Sb)   == c(p, p))) stop("dim(Sb) must match the number of covariates.")
  #else if   (!all(diag(Sb)) > 0)         stop("Sb must be positive definite.")
  #else if (!all(Sb == diag(diag(Sb))))   Sbdiag <- FALSE

  #if (!Sbdiag) {
    #joint.update <- TRUE
    #warning("Beta prior variance-covariance matrix is non-diagonal. Forcing joint update.")
    #}
  #if (!joint.update) Sb <- diag(Sb)


  for(itr in 2:iter){
    message(sprintf("Starting iter: %d", itr))
    # Optimization to find the mode
    result <- tryCatch({
      optim(par = beta, fn = logpost_beta, linkinv = linkinv, z = z, X = X,
            atoms = z.tld, jumps  = J.tld, mu_beta = mb, # Used to be mu_beta
            sigma_beta = Sb, # Used to be sigma_beta
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
      }, error = function(e) return(NULL))  # If inversion fails, return NULL
      
      # Only proceed if covariance matrix is valid and positive definite
      if (!is.null(beta_cov_)) {
        #beta_cov_ <- diag(scale_factors) %*% beta_cov_ %*% diag(scale_factors) # E old
        #beta_cov_ <- diag(rho) %*% beta_cov_ %*% diag(rho)
        # Check positive semi-definiteness
        if (all(eigen(beta_cov_, symmetric = TRUE, only.values = TRUE)$values >= 
                -sqrt(.Machine$double.eps))) {
          
          # Force positive defineness if needed
          beta_cov_ <- Matrix::nearPD(beta_cov_)$mat %>% as.matrix()
          
          # Sample new beta from proposal distribution
          beta_ <- rmvnorm(n = 1, mean = beta_mode_, sigma = beta_cov_) %>% as.numeric()
          
          #mean_z_ <- exp(X %*% beta_) / (1 + exp(X %*% beta_))
          # mean_z_ <- plogis(X %*% beta_) # E old
          mean_z_ <- linkinv(X %*% beta_)
          
          if (min(mean_z_) >= min(z.tld) && max(mean_z_) <= max(z.tld)) {
            
            # Compute log proposal values
            pr_logprop_beta <- dmvnorm(x = beta_, mean = beta_mode, sigma = beta_cov, log = TRUE)
            cr_logprop_beta <- dmvnorm(x = beta, mean = beta_mode_, sigma = beta_cov_, log = TRUE)
            
            # Compute log-posterior values
            cr_logpost_beta <- logpost_beta(beta = beta, linkinv = linkinv, z = z, X = X, atoms = z.tld, 
                                            jumps  = J.tld, mu_beta = mb, # Used to be mu_beta
                                            sigma_beta = Sb) # Used to be sigma_beta
            pr_logpost_beta <- logpost_beta(beta = beta_, linkinv = linkinv, z = z, X = X, atoms = z.tld, 
                                            jumps  = J.tld, mu_beta = mb, # Used to be mu_beta
                                            sigma_beta = Sb) # Used to be sigma_beta
            
            # Metropolis-Hastings acceptance step
            log_acc_prob <- pr_logpost_beta - cr_logpost_beta + cr_logprop_beta - pr_logprop_beta
            if (log(runif(1)) < log_acc_prob) {
              count1 <- count1 + 1
              beta <- beta_
              beta_mode <- beta_mode_
              beta_cov <- beta_cov_
              meanY_x <- mean_z <- mean_z_
            }
          }
        } 
      }
    } else message("  `optim()` failed on this iteration")
    # If `optim()` fails, keep beta as it is and continue to next step.
    # theta update ------------------------------------
    theta_tilde <- gldrm:::getTheta(
      spt = z.tld,
      f0  = J.tld,
      mu  = meanY_x,
      sampprobs  = NULL,
      ySptIndex  = NULL,
      thetaStart = theta
    )$theta
    
    theta <- theta_tilde
    message("  theta updated successfully")
    # u update ----------------------------------------
    u <- sampler_u(u, zstar, nstar, theta, alpha, delta, min_y, max_y, eps)
    message("  u updated successfully")
    
    # CRM update --------------------------------------
    crm_star <- crm_sampler(M, u, zstar, nstar, theta, alpha, min_y, max_y, eps)
    z.tld_star <- c(crm_star$RL, crm_star$zstar)
    J.tld_star <- c(crm_star$RJ, crm_star$Jstar)
    
    if(min(meanY_x) >= min(z.tld_star) && max(meanY_x) <= max(z.tld_star)){
      # MH step
      theta_tilde_star <- gldrm:::getTheta(
        spt = z.tld_star,
        f0  = J.tld_star,
        mu  = meanY_x,
        sampprobs  = NULL,
        ySptIndex  = NULL,
        thetaStart = theta
      )$theta
      
      b1 <- b_theta(theta_tilde_star, z.tld_star, J.tld_star)
      b2 <- b_theta(theta_tilde, z.tld, J.tld)
      b3 <- b_theta(theta_tilde_star, z.tld, J.tld)
      b4 <- b_theta(theta_tilde, z.tld_star, J.tld_star)
      # log_r <- log(exp(sum(2*(theta_tilde_star - theta_tilde)*z - b1 + b2 - b3 + b4))) E old
      log_r <- sum(2*(theta_tilde_star - theta_tilde)*z - b1 + b2 - b3 + b4)
      
      if(log(runif(1)) < log_r){
        count2 <- count2 + 1
        RL <- crm_star$RL
        RJ <- crm_star$RJ
        zstar <- crm_star$zstar
        Jstar <- crm_star$Jstar
        z.tld <- c(RL, zstar)
        J.tld <- c(RJ, Jstar)
        theta <- theta_tilde_star
      }
    }
    
    message("  crm updated successfully")
    # z update ------------------------------------------------------------------
    z <- z_sampler_unifK(y, c0, z.tld, J.tld, theta, min_y, max_y, eps)
    message("  z updated successfully")
    # zstar and nstar update ----------------------------------------------------
    resampled_z <- resample_zstar(z)
    zstar <- resampled_z$zstar
    nstar <- resampled_z$nstar
    
    theta0 <- gldrm:::getTheta(
      spt = z.tld,
      f0  = J.tld,
      mu  = mu0, # Used to be m0
      sampprobs  = NULL,
      ySptIndex  = NULL,
      thetaStart = NULL
    )$theta
    
    temp <- exp(theta0 * z.tld - max(theta0 * z.tld))
    Jtld_0 <- (temp * J.tld) / sum(temp * J.tld)
    lnlik <- loglik(linkinv = linkinv, z = z, X = X, beta = beta, atoms = z.tld, jumps = Jtld_0)
    message("  loglik calculated successfully")
    
    # Storing MCMC simulations --------------------------------------------------
    z_samples[itr, ] <- z
    u_samples[itr, ] <- u
    beta_samples[itr,] <- beta
    theta_samples[itr, ] <- theta
    message("max(theta) = ", max(theta))
    crm_samples[[itr]] <- list(z.tld = z.tld, J.tld = J.tld)
    lnlik_samples[itr] <- lnlik
  }
  
  crm_samples <- data.frame(
  z.tld = I(lapply(crm_samples, `[[`, "z.tld")),
  J.tld = I(lapply(crm_samples, `[[`, "J.tld")))
  samples <- list(z = z_samples, beta = beta_samples, crm = crm_samples)
  
  # A: Add priors to this later?
  # Control params?
  # iter? So people dont have to calculate?
  list(samples    = samples,
       mb         = mb,
       Sb         = Sb,
       spt        = spt,
       mu0        = mu0,
       p_acc_beta = count1 / iter,
       p_acc_crm  = count2 / iter)
}
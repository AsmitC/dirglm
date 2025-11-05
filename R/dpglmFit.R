#' Control arguments for the \code{dpglm} algorithm.
#'
#' This function returns control arguments for the \code{dpglm} algorithm.
#' Each argument has a default value, which will be used unless a different
#' value is provided by the user.
#'
#' @param burnin Number of burn-in MCMC iterations. Defaults to 100.
#' @param thin Factor by which to thin MCMC iterations. Defaults to 10.
#' @param save Number of MCMC samples to return. Defaults to 1000.
#' @param rho MCMC update step size. A scalar in \eqn{(0, 1]}. Defaults to 0.1.
#' @param mb Prior mean for beta. Defaults to a p-length vector whose entries are all 0.
#' @param Sb Prior variance-covariance matrix for beta.
#' Defaults to the p-dimensional identity matrix. See details for more information.
#' @param gamma Shrinkage parameter for the (default) prior variance on \code{beta}.
#' Defaults to 1. Will not be used if \code{Sb} is specified in \code{dirglm}.
#' @param mu0 Mean of the reference distribution. The reference distribution is
#' not unique unless its mean is restricted to a specific value. This value can
#' be any number within the range of observed values, but values near the boundary
#' may cause numerical instability. This is an optional argument with \code{mean(y)}
#' being the default value.
#' @param spt Theoretical support of the response variable. Should be a vector of length 2.
#' Defaults to \code{c(min(y), max(y))}. 
#' @param betaStart Initial value for the regression coefficients \code{beta}.
#' Defaults to the output obtained by fitting \code{gldrm}.
#' @param crmStart Initial value for the reference distribution.
#' Defaults to the output obtained by fitting \code{gldrm}.
#' @param seed Random seed. Defaults to NULL.
#'
#' @return Object of S3 class "dpglmControl"
#'
#' @export
dpglm.control <- function(burnin=100, thin=10, save=1000, rho=0.1,
                          M=20, alpha=1, delta=2, c0=NULL,
                          gamma=1, mu0=NULL, spt=NULL,
                          betaStart=NULL, varbetaStart=NULL, thetaStart=NULL, crmStart=NULL,
                          seed=NULL)
{
  if (burnin < 0 || floor(burnin) != burnin) stop("Number of burn-in samples must be an integer >= 0")
  if (thin   < 1 || floor(thin)   != thin)   stop("Thin must be an integer >= 1")
  if (save   < 1 || floor(save)   != save)   stop("Number of saved iterations must be an integer >= 1")
  if (!(rho <= 1 & rho > 0))                 stop("rho must lie in (0, 1]")
  if (length(spt) != 2)                      stop("Support must be a vector of length 2")
  if(spt[1] > spt[2])                        stop("spt[1] must be less than or equal to spt[2]")
  ctrl <- list(burnin       = burnin,
               thin         = thin,
               save         = save,
               rho          = rho,
               M            = M,
               alpha        = alpha,
               delta        = delta,
               c0           = c0,
               gamma        = gamma,
               mu0          = mu0,
               spt          = spt,
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
  rho          <- dpglmControl$rho
  M            <- dpglmControl$M
  alpha        <- dpglmControl$alpha
  delta        <- dpglmControl$delta
  c0           <- dpglmControl$c0
  gamma        <- dpglmControl$gamma
  seed         <- dpglmControl$seed

  if (!is.null(seed)) set.seed(seed)
  
  if(is.null(c0)) c0 <- c0_silverman(y) / 4
  
  # Extract link
  linkfun <- link$linkfun
  linkinv <- link$linkinv
  mu.eta  <- link$mu.eta

  # Extract spt
  min_y <- spt[1]
  max_y <- spt[2]

  # MCMC Initialization
  X    <- as.matrix(X, nrow=n)
  n    <- length(y)
  p    <- dim(X)[2]
  iter <- burnin + thin * save

  beta_samples <- matrix(NA, nrow = iter, ncol = p)
  theta_samples <- matrix(NA, nrow = iter, ncol = n)
  u_samples <- z_samples <- matrix(NA, nrow = iter, ncol = n)
  crm_samples   <- list()
  lnlik_samples <- numeric(iter)

  gldrm_fit <- gldrm(y ~ X[, -1], link = link) # E old
  beta_samples[1, ] <- beta <- as.numeric(init$beta)
  beta_mode <- beta
  #beta_cov <- gldrm_fit$varbeta # What to do about this
  beta_cov <- init$varbeta # this
  meanY_x <- linkinv(X %*% beta)
  z.tld <- z_tld <-  as.numeric(init$crmStart$z.tld)
  J.tld <- J_tld <- as.numeric(init$crmStart$J.tld)
  crm_samples[[1]] <- list(z.tld = z.tld, J.tld = J.tld)

  ord <- order(J_tld)[1:M]
  
  RL <- z_tld[ord]
  RJ <- J_tld[ord]
  #theta_samples[1, ] <- theta <- true_theta #gldrm_fit$theta %>% as.numeric() # What to do about this
  theta_samples[1, ] <- theta <- init$theta # this
  z_samples[1, ] <- z <- z_sampler_unifK(y, c0, z.tld, J.tld, theta, min_y, max_y)

  btheta <- b_theta(theta, z.tld, J.tld)
  
  T_vec <- exp(btheta)
  u_samples[1, ] <- u <- rgamma(n, shape = 1, rate = T_vec)
  
  resampled_z <- resample_zstar(z)
  zstar <- resampled_z$zstar
  nstar <- resampled_z$nstar
  Jstar <- rgamma(n = length(nstar), shape = nstar, rate = 1)
  
  count1 <- count2 <- 0
  #scale_factors <- c(1.5, 2.5) %>% sqrt() # A: Using rho instead?
  # Could allow rho to be either a scalar (specify rep(rho, p))
  # Or a vector of length(p). Need to relax (0, 1) requirement then.

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
  Jtld_0 <- (temp * J.tld) / sum(temp * J.tld)
  lnlik_samples[1] <- lnlik <- loglik(linkinv = linkinv, z = z, X = X, beta = beta, atoms = z.tld, jumps = Jtld_0)

  for(itr in 2:iter){
    # Optimization to find the mode
    result <- tryCatch({
      optim(par = beta, fn = logpost_beta, linkinv = linkinv, z = z, X = X,
            atoms = z.tld, jumps  = J.tld, mu_beta = mu_beta,
            sigma_beta = sigma_beta,
            control = list(fnscale = -1), hessian = TRUE)
    }, error = function(e) return(NULL))  # If error occurs, return NULL
    
    # Check if `optim()` failed
    if (!is.null(result) && is.finite(det(result$hessian))) {
      beta_mode_ <- result$par
      
      # Safely compute Hessian inverse
      beta_cov_ <- tryCatch({
        -solve(result$hessian)
      }, error = function(e) return(NULL))  # If inversion fails, return NULL
      
      # Only proceed if covariance matrix is valid and positive definite
      if (!is.null(beta_cov_)) {
        beta_cov_ <- diag(scale_factors) %*% beta_cov_ %*% diag(scale_factors)
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
                                            jumps  = J.tld, mu_beta = mu_beta, 
                                            sigma_beta = sigma_beta)
            pr_logpost_beta <- logpost_beta(beta = beta_, linkinv = linkinv, z = z, X = X, atoms = z.tld, 
                                            jumps  = J.tld, mu_beta = mu_beta, 
                                            sigma_beta = sigma_beta)
            
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
    }
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
    
    # u update ----------------------------------------
    u <- sampler_u(u, zstar, nstar, theta, alpha, delta, min_y, max_y)
    
    # CRM update --------------------------------------
    crm_star <- crm_sampler(M, u, zstar, nstar, theta, alpha, min_y, max_y)
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
      log_r <- log(exp(sum(2*(theta_tilde_star - theta_tilde)*z - b1 + b2 - b3 + b4)))
      
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
    

    # z update ------------------------------------------------------------------
    z <- z_sampler_unifK(y, c0, z.tld, J.tld, theta, min_y, max_y)
    
    # zstar and nstar update ----------------------------------------------------
    resampled_z <- resample_zstar(z)
    zstar <- resampled_z$zstar
    nstar <- resampled_z$nstar
    
    theta0 <- gldrm:::getTheta(
      spt = z.tld,
      f0  = J.tld,
      mu  = m0,
      sampprobs  = NULL,
      ySptIndex  = NULL,
      thetaStart = NULL
    )$theta
    
    temp <- exp(theta0 * z.tld - max(theta0 * z.tld))
    Jtld_0 <- (temp * J.tld) / sum(temp * J.tld)
    lnlik <- loglik(linkinv = linkinv, z = z, X = X, beta = beta, atoms = z.tld, jumps = Jtld_0)
    
    # Storing MCMC simulations --------------------------------------------------
    z_samples[itr, ] <- z
    u_samples[itr, ] <- u
    beta_samples[itr,] <- beta
    theta_samples[itr, ] <- theta
    crm_samples[[itr]] <- list(z.tld = z.tld, J.tld = J.tld)
    lnlik_samples[itr] <- lnlik
  }

  samples <- list(z = z_samples, beta = beta_samples, crm = crm_samples)
  
  # A: Add priors to this later?
  # Control params?
  # iter? So people dont have to calculate?
  dpglm_fit <- list(samples    = samples,
                    spt        = spt,
                    mu0        = mu0,
                    p_acc_beta = count1 / iter,
                    p_acc_crm  = count2 / iter)
  
  
  out <- list(data = dat, dpglm = dpglm_fit)
  return(out)
}
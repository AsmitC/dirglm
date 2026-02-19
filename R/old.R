fit_func <- function(dat, iter, truth, mu_beta, sigma_beta, scale_factors, m0, ebp) {
  y <- dat[, 1]
  min_y <- 0
  max_y <- 1
  true_theta <- dat[, 2]
  X <- dat[, -c(1, 2)]
  n <- length(y)
  # Tuning Parameters --------------------------------------------------------------------
  rho <- 1
  M <- 20
  alpha <- 1
  delta <- 2
  c0 <- old_c0_silverman(y) / 4 #0.25
  # mu_beta <- truth$beta
  # sigma_beta <- 2.5
  # sigma_beta <- rho * sigma_beta
  
  # Data Preparation ---------------------------------------------------------------------
  X <- as.matrix(X)
  y <- as.numeric(y)
  n <- length(y)
  p <- dim(X)[2]
  
  # Link Function -----------------------------------------------------------------------
  link <- 'logit'
  
  ## Initialization -----------------------------------------------------------------------
  beta_samples <- matrix(NA, nrow = iter, ncol = p)
  theta_samples <- matrix(NA, nrow = iter, ncol = n)       
  u_samples <- z_samples <- matrix(NA, nrow = iter, ncol = n)
  crm_samples   <- list()
  lnlik_samples <- numeric(iter)
  gldrm_fit <- gldrm(y ~ X[, -1], link = link)
  beta_samples[1, ] <- beta <- as.numeric(truth$beta) #gldrm_fit$beta %>% as.numeric()
  beta_mode <- beta
  beta_cov <- gldrm_fit$varbeta
  #meanY_x <- exp(X %*% beta) / (1 + exp(X %*% beta))
  meanY_x <- plogis(X %*% beta)
  z.tld <- z_tld <-  as.numeric(truth$f0[, 1]) #gldrm_fit$spt %>% as.numeric() 
  J.tld <- J_tld <- as.numeric(truth$f0[, 2]) #gldrm_fit$f0 %>% as.numeric()   
  crm_samples[[1]] <- list(z.tld = z.tld, J.tld = J.tld)
  
  ord <- order(J_tld)[1:M]
  
  RL <- z_tld[ord]
  RJ <- J_tld[ord]
  theta_samples[1, ] <- theta <- true_theta #gldrm_fit$theta %>% as.numeric() #
  z_samples[1, ] <- z <- old_z_sampler_unifK(y, c0, z.tld, J.tld, theta, min_y, max_y)
  #z_samples[1, ] <- z <- z_sampler_triK(y, c0, z.tld, J.tld, theta)
  #z_samples[1, ] <- z <- z_sampler_epanK(y, c0, z.tld, J.tld, theta)
  
  btheta <- old_b_theta(theta, z.tld, J.tld)
  
  T_vec <- exp(btheta)
  u_samples[1, ] <- u <- rgamma(n, shape = 1, rate = T_vec)
  
  resampled_z <- old_resample_zstar(z)
  zstar <- resampled_z$zstar
  nstar <- resampled_z$nstar
  Jstar <- rgamma(n = length(nstar), shape = nstar, rate = 1)
  
  count1 <- count2 <- 0
  #scale_factors <- c(1.5, 2.5) %>% sqrt() 
  #beta_cov <- 0.5 * gldrm_fit$varbeta
  
  theta0 <- gldrm:::getTheta(
    spt = z.tld,
    f0  = J.tld,
    mu  = m0,
    sampprobs  = NULL,
    ySptIndex  = NULL,
    thetaStart = NULL
  )$theta

  temp  <- exp(theta0 * z.tld - max(theta0 * z.tld))
  Jtilt <- temp * J.tld
  W     <- sum(Jtilt)
  Jtld_0 <- Jtilt / W

  if (ebp == 0L) {
    J.tld_ll <- Jtld_0
  } else if (ebp == 1L) {
    J.tld <- Jtilt
    J.tld_ll <- J.tld
  } else if (ebp == 2L) {
    J.tld <- Jtld_0
    J.tld_ll <- J.tld
  } else stop("ebp must be one of {0L, 1L, 2L}")

  lnlik_samples[1] <- lnlik <- old_loglik(z = z, X = X, beta = beta, atoms = z.tld, jumps = J.tld_ll)
  #itr <- 1
  for(itr in 2:iter){
    #message(sprintf("Starting iter: %d", itr))
    #itr <- itr + 1 
    # Optimization to find the mode
    result <- tryCatch({
      optim(par = beta, fn = old_logpost_beta, z = z, X = X,
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
          
          # Optional: Force positive defineness if needed
          beta_cov_ <- as.matrix(Matrix::nearPD(beta_cov_)$mat)
          
          # Sample new beta from proposal distribution
          beta_ <- as.numeric(rmvnorm(n = 1, mean = beta_mode_, sigma = beta_cov_))
          
          
          #mean_z_ <- exp(X %*% beta_) / (1 + exp(X %*% beta_))
          mean_z_ <- plogis(X %*% beta_)
          
          if (min(mean_z_) >= min(z.tld) && max(mean_z_) <= max(z.tld)) {
            
            # Compute log proposal values
            
            pr_logprop_beta <- dmvnorm(x = beta_, mean = beta_mode, sigma = beta_cov, log = TRUE)
            cr_logprop_beta <- dmvnorm(x = beta, mean = beta_mode_, sigma = beta_cov_, log = TRUE)
            
            # Compute log-posterior values
            cr_logpost_beta <- old_logpost_beta(beta = beta, z = z, X = X, atoms = z.tld, 
                                            jumps  = J.tld, mu_beta = mu_beta, 
                                            sigma_beta = sigma_beta)
            pr_logpost_beta <- old_logpost_beta(beta = beta_, z = z, X = X, atoms = z.tld, 
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
          } else message("min(z_) not within z.tld bounds")
        } else message("beta covariance not positive semi-definite")
      } else message("beta covariance is null")
    } else message("either optim failed or det(hessian) is infinite")
    # If `optim()` fails, keep beta as it is and continue to next step.    
    
    # theta update ------------------------------------
    theta <- gldrm:::getTheta(
      spt = z.tld,
      f0  = J.tld,
      mu  = meanY_x,
      sampprobs  = NULL,
      ySptIndex  = NULL,
      thetaStart = theta
    )$theta
    
    # u update ----------------------------------------
    u <- old_sampler_u(u, zstar, nstar, theta, alpha, delta, min_y, max_y)
    
    # CRM update --------------------------------------
    sd_theta <- rep(1, n) # Initial proposal standard deviation for theta_tilde
    burn_in <- min(200, floor(iter / 2))
    if(itr == burn_in + 1) sd_theta <- apply(theta_samples[2:burn_in, ], 2, sd) # sd could be replaced by sd(diff)
    crm_star <- old_crm_sampler(M, u, zstar, nstar, theta, sd_theta, alpha, min_y, max_y, itr)
    z.tld_star <- c(crm_star$RL, crm_star$zstar)
    J.tld_star <- c(crm_star$RJ, crm_star$Jstar)
    theta_tilde <- crm_star$theta_tilde
    
    if(min(meanY_x) >= min(z.tld_star) && max(meanY_x) <= max(z.tld_star)){
      # MH step
      theta_star <- gldrm:::getTheta(
        spt = z.tld_star,
        f0  = J.tld_star,
        mu  = meanY_x,
        sampprobs  = NULL,
        ySptIndex  = NULL,
        thetaStart = theta
      )$theta
      
      b1 <- old_b_theta(theta_star, z.tld_star, J.tld_star)
      b2 <- old_b_theta(theta, z.tld, J.tld)
      b3 <- old_b_theta(theta_tilde, z.tld, J.tld)
      b4 <- old_b_theta(theta_tilde, z.tld_star, J.tld_star)
      #log_r <- sum(2*(theta_star - theta)*z - b1 + b2 - b3 + b4)
      log_r <- sum((theta_star - theta)*z - b1 + b2 - b3 + b4) +
        sum(dnorm(theta_tilde, mean = theta_star, sd = sd_theta, log = TRUE) - 
        dnorm(theta_tilde, mean = theta, sd = sd_theta, log = TRUE))
      
      if(log(runif(1)) < log_r){
        count2 <- count2 + 1
        RL <- crm_star$RL
        RJ <- crm_star$RJ
        zstar <- crm_star$zstar
        Jstar <- crm_star$Jstar
        z.tld <- c(RL, zstar)
        J.tld <- c(RJ, Jstar)
        theta <- theta_star
      }
    }
    

    # z update ------------------------------------------------------------------
    z <- old_z_sampler_unifK(y, c0, z.tld, J.tld, theta, min_y, max_y)
    
    # zstar and nstar update ----------------------------------------------------
    resampled_z <- old_resample_zstar(z)
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
    
    temp  <- exp(theta0 * z.tld - max(theta0 * z.tld))
    Jtilt <- temp * J.tld
    W     <- sum(Jtilt)
    Jtld_0 <- Jtilt / W

    if (ebp == 0L) {
      # EBP 0: keep chain J.tld as-is; likelihood uses tilted+normalized
      J.tld_ll <- Jtld_0
    } else if (ebp == 1L) {
      # EBP 1: chain uses tilt-only (unnormalized); likelihood uses same
      J.tld <- Jtilt
      J.tld_ll <- J.tld
    } else if (ebp == 2L) {
      # EBP 2: chain uses tilt+renorm; likelihood uses same
      J.tld <- Jtld_0
      J.tld_ll <- J.tld
    } else stop("ebp must be one of {0L, 1L, 2L}")

    if (itr %% 250 == 0) {
      Jsum <- sum(J.tld)
      Jnorm <- J.tld / Jsum
      cat(
        "itr", itr,
        "| ebp=", ebp,
        "| sumJ=", format(Jsum, scientific=TRUE),
        "| Jmax_share=", round(max(Jnorm), 3),
        "| ess(J)=", round(1/sum(Jnorm^2), 3),
        "| W=", format(W, scientific=TRUE),
        "| mean(Jtld_0)=", round(sum(z.tld * Jtld_0), 6),
        "| mu0=", round(m0, 6),
        "\n\n"
      )
    }

    lnlik <- old_loglik(z = z, X = X, beta = beta, atoms = z.tld, jumps = J.tld_ll)
    
    # Storing MCMC simulations --------------------------------------------------
    z_samples[itr, ] <- z
    u_samples[itr, ] <- u
    beta_samples[itr,] <- beta
    theta_samples[itr, ] <- theta
    theta_samples[itr, ] <- theta
    crm_samples[[itr]] <- list(z.tld = z.tld, J.tld = J.tld)
    lnlik_samples[itr] <- lnlik
  }
  
  dpglm_fit <- list(z = z_samples, beta = beta_samples,
                    theta = theta_samples,
                    crm = crm_samples, 
                    beta_acceptance = count1 / iter,
                    crm_acceptance = count2 / iter)
  
  
  out <- list(data = dat, dpglm = dpglm_fit)
  return(out)
}

# -----------------------------------------------------------------------------#
#                               Update z                                      #
# -----------------------------------------------------------------------------#

old_z_sampler_unifK <- function(y, c0, crm.atoms, crm.jumps, tht, min_y, max_y) {
  n <- length(y)
  eps <- 1e-6
  lower <- crm.atoms - c0
  upper <- crm.atoms + c0
  # lower[lower < 0] <- 0
  # upper[upper > 1] <- 1
  lower[lower < min_y] <- min_y + eps  
  upper[upper > max_y] <- max_y - eps
  width <- upper - lower
  z <- numeric(n)
  for (i in 1:n) {
    indices <- y[i] >= lower & y[i] <= upper
    indx <- which(indices)
    log_prob <- -log(width[indx]) + (tht[i] * crm.atoms[indx]) + log(crm.jumps[indx])
    prob <- exp(log_prob - max(log_prob))
    if(sum(indices) == 1){
      z[i] <- crm.atoms[indices]
    } else {
      z[i] <- sample(crm.atoms[indx], 1, prob = prob)
    }
  }
  return(z)
}

old_c0_silverman <- function(y) {
  n <- length(y)
  sd_y <- sd(y)
  iqr_y <- IQR(y)
  return(0.9 * min(sd_y, iqr_y / 1.34) * n^(-1/5))
}

# -----------------------------------------------------------------------------#
#                               Resample z                                     #
# -----------------------------------------------------------------------------#

old_resample_zstar <- function(z){
  z_table <- table(z)
  zstar   <- as.numeric(names(z_table))
  nstar   <- as.numeric(z_table)
  ## Write code for resampling zstar to avoid the ‘sticky clusters effect’
  
  return(list(zstar = zstar, nstar = nstar))
}

old_logpost_beta <- function(beta, z, X, atoms, jumps, mu_beta, sigma_beta){
  n <- length(z)
  #mu <- exp(X %*% beta) / (1 + exp(X %*% beta))
  mu <- plogis(X %*% beta)
  theta <- gldrm:::getTheta(spt = atoms, f0 = jumps, mu = mu, 
                            ySptIndex = NULL, sampprobs = NULL)$theta
  btheta <- old_b_theta(theta, atoms, jumps)
  log_post <- sum(theta * z - btheta) - sum((beta - mu_beta)^2 / (2 * sigma_beta^2)) # others const in beta
  return(log_post)
}

old_loglik <- function(z, X, beta, atoms, jumps){
  n <- length(z)
  #mu <- exp(X %*% beta) / (1 + exp(X %*% beta))   
  mu <- plogis(X %*% beta)
  theta <- gldrm:::getTheta(spt = atoms, f0 = jumps, mu = mu, ySptIndex = NULL, 
                            sampprobs = NULL)$theta
  
  btheta <- old_b_theta(theta, atoms, jumps) 
  
  f0_z <- numeric(n)
  for (i in 1:n) {
    f0_z[i] <- sum(jumps[z[i] == atoms])
  }
  loglik <- sum(theta * z - btheta + log(f0_z))
  return(loglik)
}


old_log_post_u <- function(u, zstar, nstar, theta, alpha, min_y, max_y) {
  # Number of grid points for the continuous part
  R <- 250
  eps <- 1e-6
  
  # Construct a grid for integration over the continuous part (G_0, uniform on (0,1))
  #z_grid <- seq(eps, 1 - eps, length.out = R)
  z_grid <- seq(min_y + eps, max_y - eps, length.out = R)
  diff_z <- diff(z_grid)[1]  # uniform grid spacing
  
  # ---- Continuous part ----
  # For each grid point v in z_grid, compute in a stabilized manner:
  #    log(1 + sum_i u_i * exp(theta_i * v))
  cont_vals <- sapply(z_grid, function(v) {
    A <- log(u) + theta * v
    max_A <- max(A)
    S <- exp(max_A) * sum(exp(A - max_A))
    log1p(S)  # log(1+S) computed in a numerically stable way
  })
  
  # Riemann-sum approximation of the integral over the continuous part:
  integral_continuous <- sum(cont_vals) * diff_z
  
  # ---- Discrete part ----
  # For each unique atom in zstar, with multiplicity given by nstar,
  # compute: nstar * log(1 + sum_i u_i * exp(theta_i * zstar))
  disc_vals <- sapply(zstar, function(v) {
    A <- log(u) + theta * v
    max_A <- max(A)
    S <- exp(max_A) * sum(exp(A - max_A))
    log1p(S)
  })
  
  sum_discrete <- sum(disc_vals * nstar)
  
  # Combine the continuous and discrete contributions
  neg_log_post <- alpha * integral_continuous + sum_discrete
  
  return(-neg_log_post)
}

old_sampler_u <- function(u, zstar, nstar, theta, alpha, delta, min_y, max_y) {
  n <- length(u) # Get the length of u
  
  for (i in seq_len(n)) {
    u_star <- u
    # Sample from Gamma distribution: shape = delta, scale = u[i] / delta
    u_star[i] <- rgamma(1, shape = delta, scale = u[i] / delta)
    
    # Compute logQ_ratio: log(q(ui | ui_star)) - log(q(ui_star | ui))
    logQ_ratio <- dgamma(u[i], shape = delta, scale = u_star[i] / delta, log = TRUE) - 
      dgamma(u_star[i], shape = delta, scale = u[i] / delta, log = TRUE)
    
    # Compute logratio
    logratio <- old_log_post_u(u_star, zstar, nstar, theta, alpha, min_y, max_y) - 
      old_log_post_u(u, zstar, nstar, theta, alpha, min_y, max_y) + logQ_ratio
    
    # Metropolis-Hastings acceptance step
    if (log(runif(1)) < logratio) {
      u[i] <- u_star[i]
    } # else, u[i] remains unchanged
  }
  
  return(u)
}
# --------------------------------------------------------------------------- #
#                          For beta update                                   #
# --------------------------------------------------------------------------- #

old_b_theta <- function(theta, spt, f0) {
  # # Compute the matrix of exponentiated outer products
  # exp_matrix <- exp(outer(theta, spt, "*"))
  # 
  # # Compute the dot product row-wise and apply the log
  # result <- log(exp_matrix %*% f0)
  
  # Compute the vector for each spt:= log(f0) + log(grid spacing)
  log_f0_dz <- log(f0) #+ log(diff(spt)[1])
  
  # Create the log weights matrix: each element is theta[i]*spt[j] + log(f0[j]) + log(diff(spt)[1])
  log_weights <- outer(theta, spt, "*") +
    matrix(log_f0_dz, nrow = length(theta), ncol = length(spt), byrow = TRUE)
  
  # For each row, compute the log-sum-exp
  result <- apply(log_weights, 1, function(x) {
    m <- max(x)
    m + log(sum(exp(x - m)))
  })
  
  
  return(as.vector(result))  # Ensure the output is a vector
}

# --------------------------------------------------------------------------- #
#                               CRM update                                   #
# --------------------------------------------------------------------------- #

old_crm_sampler <- function(M, u, zstar, nstar, tht_, sd_, alpha, min_y, max_y, itr){
  N <- 3001
  R <- 3001
  eps <- 1e-6
  s <- -log(seq(exp(-eps), exp(-5e-4), length.out = N))

  #tht <- rnorm(1, mean = tht_, sd = sd_)
  tht <- rnorm(length(tht_), mean = tht_, sd = sd_)
  
  # Sorted, ascending order, needed in RL
  #z <- seq(eps, 1 - eps, length.out = R)
  z <- seq(min_y + eps, max_y - eps, length.out = R)
  
  # Assume u is a vector and z is the grid where you want to compute psi_z.
  # Here we compute psi_z for each grid point z[j] by summing over u.
  
  # Compute a matrix of log-terms:
  log_terms <- outer(log(u), z, function(lu, z_val) lu + tht * z_val)
  
  # Now compute psi_z using the log-sum-exp trick along the u dimension (rows):
  psi_z <- apply(log_terms, 2, function(log_vec) {
    max_val <- max(log_vec)
    exp(max_val + log(sum(exp(log_vec - max_val))))
  })

  if (itr %% 250 == 0) {
  cat("  psi_z range:", paste(round(range(psi_z), 3), collapse=","), "\n")
}
  
  # Assume:
  #   psi_z: vector of length R computed as before over the z grid
  #   s: vector of length N representing the s grid
  #   eps: small positive number for numerical boundaries
  dz <- (1 - 2 * eps) / (R - 1)   # Grid spacing for the uniform measure
  
  # Preallocate the result vector
  fnS <- numeric(length(s))
  
  for (j in seq_along(s)) {
    # For a fixed s[j], the integrand at each z is:
    #    f(z, s[j]) = exp( -(1+psi_z) * s[j] ) / s[j]
    # Taking logs gives:
    #    log f(z, s[j]) = -(1+psi_z) * s[j] - log(s[j])
    log_vals <- -(1 + psi_z) * s[j] - log(s[j])
    
    # Use the log-sum-exp trick for the summation over z grid:
    max_val <- max(log_vals)
    log_sum_exp <- max_val + log(sum(exp(log_vals - max_val)))
    
    # Multiply by the grid spacing to approximate the integral:
    log_integral <- log(dz) + log_sum_exp
    
    # Store the stabilized value (exponentiating back)
    fnS[j] <- exp(log_integral)
  }
  
  
  ds <- diff(s)
  h <- (fnS[-N] + fnS[-1]) / 2
  Nv <- rev(cumsum(rev(ds * h)))
  Nv <- Nv * alpha
  Nv <- c(Nv, 0)
  
  
  #Generate random jumps RJ and random locations RL
  xi <- cumsum(rexp(M, rate = 1.0))
  RJ <- numeric(M)
  iNv <- N - 1
  
  for (i in seq_len(M)) {
    while (iNv > 0 && Nv[iNv] < xi[i]) {
      iNv <- iNv - 1
    }
    RJ[i] <- s[iNv + 1]
  }
  
  
  RL <- numeric(M)
  for (m in seq_len(M)) {
    xi_rl <- runif(1)
    # Compute log probabilities: log_temp = -(1 + psi_z) * RJ[m]
    log_temp <- -(1 + psi_z) * RJ[m]
    
    # Stabilize by subtracting the maximum log value
    max_log <- max(log_temp)
    
    # Convert to probabilities in a numerically stable way
    p <- exp(log_temp - max_log)
    p <- p / sum(p)
    
    # Compute the cumulative probabilities
    cumsum_p <- cumsum(p)
    
    # Select the index where the cumulative sum exceeds xi_rl
    RL[m] <- z[min(which(cumsum_p > xi_rl))]
  }
  
  
  # Second Part: random jumps [for fixed locations]
  # Suppose u is a vector and zstar is a vector.
  # We want to compute, for each fixed zstar[j]:
  #   psi_star[j] = sum(u_i * exp(tht * zstar[j]))
  # We'll compute it in a numerically stable way:
  
  # Create a matrix of log-terms: each element is log(u_i) + tht * zstar[j]
  log_terms <- outer(log(u), zstar, function(lu, z) lu + tht * z)
  
  # For each column (corresponding to a fixed zstar), use the log-sum-exp trick
  psi_star <- sapply(seq_along(zstar), function(j) {
    max_val <- max(log_terms[, j])
    exp(max_val + log(sum(exp(log_terms[, j] - max_val))))
  })
  
  
  Jstar <- rgamma(length(zstar), shape = nstar, rate = psi_star + 1)
  
  return(list(
    RL = RL,
    RJ = RJ,
    zstar = zstar,
    Jstar = Jstar,
    theta_tilde = tht
  ))
}
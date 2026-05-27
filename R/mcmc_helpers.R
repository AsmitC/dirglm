## Required packages -----
## gldrm, mvtnorm

## Computes log(sum(exp(x))) with better precision
logSumExp <- function(x)
{
  i <- which.max(x)
  m <- x[i]
  lse <- log1p(sum(exp(x[-i]-m))) + m
  lse
}

#' Title: function for finding theta, btheta & bprime2
#'
#' @param spt (vector) support, \{s_j, j = 1(1)l\}, of response variable y
#' @param f0 (vector) \{f0(y = s_j), j=1(1)l\}: free parameter, f0(.) centering density
#' @param mu (vector) \{E(y_i | x_i), i=1(1)n\}
#' @param thtst (vector) initial value in theta solving iterative algorithm
#'
#' @return (list) containing:
#' 1. bpr2 (vector) \{second derivative of b(theta_i), i=1(1)n\}
#' 2. tht (vector) derived parameter theta = \{theta_i = function(beta, f0, x_i), i=1(1)n\}
#' 3. btht (vector)  normalizing constant b(theta) = \{b(theta_i), i=1(1)n\}
#' @keywords internal

tht_sol <- function(spt, f0, mu, thtst) {
  out <- gldrm:::getTheta(
    spt = spt,
    f0 = f0,
    mu = mu,
    sampprobs = NULL,
    ySptIndex = NULL,
    thetaStart = thtst,
    thetaControl = gldrm:::theta.control(logsumexp=TRUE)
  )
  tht <- out$theta
  bpr2 <- out$bPrime2
  lw   <- sweep(tcrossprod(out$theta, spt), 2, log(f0), "+")
  btht <- apply(lw, 1, logSumExp)
  return(list(bpr2 = bpr2, tht = tht, btht = btht))
}

#' Title: function for finding dirichlet proposal parameter
#'
#' @param y (vector) response variable values \{y_i, i=1(1)n\}
#' @param tht (vector) derived parameter, \{theta_i = function(beta, f0, x_i), i=1(1)n\}
#' @param btht (vector)  normalizing constant b(theta) = \{b(theta_i), i=1(1)n\}
#' @param dir_pr_parm (vector) dirichlet prior parameter
#' @param ind_mt indicator matrix with (i,j)th element: 1(y_i == s_j)
#'
#' @return (vector) dirichlet proposal parameter
#' @keywords internal

dir_parm <- function(y, tht, btht, dir_pr_parm, ind_mt) {
  wgt  <- (1 / exp(tht * y - btht)) / sum(1 / exp(tht * y - btht))
  wgt  <- as.numeric(wgt)
  parm <- dir_pr_parm + colSums(ind_mt * wgt)

  return(parm)
}

#' Title: function for log Dir pdf for x ~ Dir(a)
#'
#' @param x (vector) value at which to evaluate pdf
#' @param a (vector) Dirichlet parameter
#'
#' @return (vector) log Dirichlet pdf
#' @keywords internal
lDir = function(x, a){
  a0 = sum(a)
  logB = sum(lgamma(a))-lgamma(a0)
  f = sum( (a-1)*log(x) ) - logB
  return(f)
}

#' Title: function for calculating \{f0(y_i), i = 1(1)n\}
#'
#' @param y (vector) response variable
#' @param spt (vector) support, \{s_j, j = 1(1)l\}, of response variable y
#' @param f0 (vector) \{f0(s_j), j=1(1)l\}: free parameter, f0(.) centering density
#'
#' @return (vector) \{f0(y_i), i=1(1)n\}
#' @keywords internal

f0y <- function(y, spt, f0) {
  ind_mt <- outer(spt, y, `==`)
  f0_y <- colSums(f0 * ind_mt)
  return(f0_y)
}

#' Title: function for finding variance-covariance matrix of beta
#'
#' @param X design (matrix)
#' @param mu (vector) \{E(y_i | x_i), i=1(1)n\}
#' @param bpr2 (vector) \{second derivative of b(theta_i), i=1(1)n\}
#' @param rho MCMC update step size, a (scalar) in (0, 1]
#'
#' @return variance-covariance (matrix) of beta
#' @keywords internal

Sigma_beta <- function(X, mu, bpr2, rho, linkfun, mu.eta) {
  eta <- linkfun(mu)
  # gpr <- as.numeric(1 / mu)              # due to log link: g'(mu) = 1 / mu
  gpr <- as.numeric(1 / mu.eta(eta))       # Updated for general link
  gprsq_bpr2 <- gpr ^ 2 * bpr2
  Xstar <- X / sqrt(gprsq_bpr2)
  info_mt <- t(Xstar) %*% Xstar
  Sigma <- rho ^ 2 * solve(info_mt)

  return(Sigma)
}

#' Title: function for updating f0
#'
#' @param y response variable
#' @param spt y support
#' @param cr_f0 current f0 = \{f0(s_j), j=1(1)l\}
#' @param cr_f0y current f0(y) = \{f0(y_i), i=1(1)n\}
#' @param cr_dir_parm current dirichlet proposal parameter
#' @param cr_mu current mu
#' @param cr_tht current theta
#' @param cr_bpr2 current \{second derivative of b(theta_i), i=1(1)n\}
#' @param cr_btht current b(theta)
#' @param dir_pr_parm dirichlet prior parameter
#' @param ind_mt indicator matrix with (i,j)th element: 1(y_i == s_j)
#'
#' @importFrom extraDistr ddirichlet
#' @keywords internal
f0_update <- function(y,
                      spt,
                      cr_f0,
                      cr_f0y,
                      cr_dir_parm,
                      cr_mu,
                      cr_tht,
                      cr_bpr2,
                      cr_btht,
                      dir_pr_parm,
                      ind_mt) {
  n <- length(y)
  l <- length(spt)
  eps <- .Machine$double.eps

  for (j in 1:l) {
    pr_f0 <- cr_f0
    pr_f0[j] <- rbeta(1, cr_dir_parm[j], sum(cr_dir_parm[-j]))
    pr_f0[-j] <- pr_f0[-j] * (1 - pr_f0[j]) / sum(pr_f0[-j])
    c1 <- cr_mu > spt[1] + eps
    c2 <- cr_mu < spt[l] - eps

    if (all(c1 & c2)) {
      out <- tht_sol(spt, pr_f0, cr_mu, cr_tht)
      pr_tht <- out$tht
      pr_btht <- out$btht
      pr_bpr2 <- out$bpr2
      pr_dir_parm <- dir_parm(y, pr_tht, pr_btht, dir_pr_parm, ind_mt)
      pr_f0y <- f0y(y, spt, pr_f0)

      pr_llik <- sum(pr_tht * y - pr_btht + log(pr_f0y))
      cr_llik <- sum(cr_tht * y - cr_btht + log(cr_f0y))

      R_pf0 = sum((dir_pr_parm-1)*log(pr_f0/cr_f0))

      cr_qf0  <- lDir(cr_f0, pr_dir_parm)
      pr_qf0  <- lDir(pr_f0, cr_dir_parm)
      alp <- min(0, (pr_llik - cr_llik + R_pf0 + cr_qf0 - pr_qf0))

      if (log(runif(1)) < alp) {
        cr_f0   <- pr_f0
        cr_tht  <- pr_tht
        cr_btht <- pr_btht
        cr_f0y  <- pr_f0y
        cr_bpr2 <- pr_bpr2
        acc_f0  <- TRUE
      } else acc_f0 <- FALSE
    } else acc_f0 <- FALSE # Closes line 141
  } # Closes line 129
  return(
    list(
      cr_f0   = cr_f0,
      cr_f0y  = cr_f0y,
      cr_tht  = cr_tht,
      cr_btht = cr_btht,
      cr_bpr2 = cr_bpr2,
      acc_f0  = acc_f0
    )
  )
}



#' Title : function for updating beta (jointly)
#'
#' @param X design matrix
#' @param y response variable
#' @param spt support of y
#' @param cr_bt current beta
#' @param cr_Sig variance-covariance matrix for current beta
#' @param cr_f0 current f0 = \{f0(s_j), j=1(1)l\}
#' @param cr_tht current theta
#' @param cr_bpr2 current \{second derivative of b(theta_i), i=1(1)n\}
#' @param cr_btht current b(theta)
#' @param rho MCMC update step size, a (scalar) in (0, 1]
#'
#' @keywords internal
beta_update_joint <- function(X,
                              y,
                              spt,
                              cr_bt,
                              cr_Sig,
                              cr_f0,
                              cr_tht,
                              cr_bpr2,
                              cr_btht,
                              rho,
                              linkfun,
                              linkinv,
                              mu.eta,
                              mb,
                              Sb) {
  n <- dim(X)[1]
  l <- length(spt)

  pr_bt <- as.vector(mvtnorm::rmvnorm(1, mean = cr_bt, sigma = cr_Sig))
  pr_mu <- as.numeric(linkinv(X %*% pr_bt))

  if (sum(spt[1] <= pr_mu & pr_mu <= spt[l]) == n) {
    out <- tht_sol(spt, cr_f0, pr_mu, cr_tht)
    pr_bpr2 <- out$bpr2
    pr_tht <- out$tht
    pr_btht <- out$btht
    pr_Sig <- Sigma_beta(X, pr_mu, pr_bpr2, rho, linkfun, mu.eta)

    pr_llik <- sum(pr_tht * y - pr_btht)
    cr_llik <- sum(cr_tht * y - cr_btht)
    pr_pbt <- mvtnorm::dmvnorm(pr_bt,
                      mean = mb,
                      sigma = Sb,
                      log = T)
    cr_pbt <- mvtnorm::dmvnorm(cr_bt,
                      mean = mb,
                      sigma = Sb,
                      log = T)
    cr_qbt <- mvtnorm::dmvnorm(cr_bt,
                      mean = pr_bt,
                      sigma = pr_Sig,
                      log = T)
    pr_qbt <- mvtnorm::dmvnorm(pr_bt,
                      mean = cr_bt,
                      sigma = cr_Sig,
                      log = T)

    alp <- min(0, pr_llik - cr_llik + pr_pbt - cr_pbt + cr_qbt - pr_qbt)

    if (log(runif(1)) < alp) {
      cr_bt <- pr_bt
      cr_tht <- pr_tht
      cr_btht <- pr_btht
      cr_bpr2 <- pr_bpr2
      cr_Sig <- pr_Sig
      acc_beta <- TRUE
    } else {
      acc_beta <- FALSE
    }
  } else { # Closes line 226
    acc_beta <- FALSE
  }
  return(list(
    cr_bt    = cr_bt,
    cr_tht   = cr_tht,
    cr_btht  = cr_btht,
    cr_bpr2  = cr_bpr2,
    acc_beta = acc_beta
  ))
}



#' Title: function for updating beta (one at a time)
#'
#' @param X design matrix
#' @param y response variable
#' @param spt support of y
#' @param cr_bt current beta
#' @param cr_Sig variance-covariance matrix for current beta
#' @param cr_f0 current f0 = \{f0(s_j), j=1(1)l\}
#' @param cr_tht current theta
#' @param cr_bpr2 current \{second derivative of b(theta_i), i=1(1)n\}
#' @param cr_btht current b(theta)
#' @param rho MCMC update step size, a (scalar) in (0, 1]
#'
#' @keywords internal
beta_update_separate <- function(X,
                                 y,
                                 spt,
                                 cr_bt,
                                 cr_Sig,
                                 cr_f0,
                                 cr_tht,
                                 cr_bpr2,
                                 cr_btht,
                                 rho,
                                 linkfun,
                                 linkinv,
                                 mu.eta,
                                 mb,
                                 Sb) {
  n <- dim(X)[1]
  p <- dim(X)[2]
  l <- length(spt)
  cr_sd <- sqrt(diag(cr_Sig))
  for (j in 1:p) {
    pr_bt <- cr_bt
    pr_bt[j] <- cr_bt[j] + rnorm(1, mean = 0, sd = cr_sd[j])
    pr_mu <- as.numeric(linkinv(X %*% pr_bt))        # Updated for general link

    if (sum(spt[1] <= pr_mu & pr_mu <= spt[l]) == n) {
      out <- tht_sol(spt, cr_f0, pr_mu, cr_tht)
      pr_tht <- out$tht
      pr_btht <- out$btht
      pr_bpr2 <- out$bpr2
      pr_Sig <- Sigma_beta(X, pr_mu, pr_bpr2, rho, linkfun, mu.eta)
      pr_sd <- sqrt(diag(pr_Sig))

      pr_llik <- sum(pr_tht * y - pr_btht)
      cr_llik <- sum(cr_tht * y - cr_btht)
      pr_pbt <- dnorm(pr_bt[j],
                      mean = mb[j],
                      sd = sqrt(Sb[j]),
                      log = T)
      cr_pbt <- dnorm(cr_bt[j],
                      mean = mb[j],
                      sd = sqrt(Sb[j]),
                      log = T)
      cr_qbt <- dnorm(cr_bt[j],
                      mean = pr_bt[j],
                      sd = pr_sd[j],
                      log = T)
      pr_qbt <- dnorm(pr_bt[j],
                      mean = cr_bt[j],
                      sd = cr_sd[j],
                      log = T)

      alp <- min(0, pr_llik - cr_llik + pr_pbt - cr_pbt + cr_qbt - pr_qbt)

      if (log(runif(1)) < alp) {
        cr_bt    <- pr_bt
        cr_tht   <- pr_tht
        cr_btht  <- pr_btht
        cr_bpr2  <- pr_bpr2
        acc_beta <- TRUE
      } else {
        acc_beta <- FALSE
      }
    } else { # Closes line 316
      acc_beta <- FALSE
    }
  }
  return(list(
    cr_bt    = cr_bt,
    cr_tht   = cr_tht,
    cr_btht  = cr_btht,
    cr_bpr2  = cr_bpr2,
    acc_beta = acc_beta
  ))
}


#####################################
# ALL FUNCTIOMS BELOW ARE FOR DPGLM #
#####################################
# -----------------------------------------------------------------------------#
#                               Update z                                      #
# -----------------------------------------------------------------------------#

z_sampler_unifK <- function(y, c0, crm.atoms, crm.jumps, tht, min_y, max_y, eps) {
  n <- length(y)
  lower <- crm.atoms - c0
  upper <- crm.atoms + c0
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

c0_silverman <- function(y) {
  n <- length(y)
  sd_y <- sd(y)
  iqr_y <- IQR(y)
  return(0.9 * min(sd_y, iqr_y / 1.34) * n^(-1/5))
}

z_sampler_triK <- function(y, c0, crm.atoms, crm.jumps, tht) {
  n <- length(y)
  z <- numeric(n)
  
  for (i in 1:n) {
    # Compute triangular kernel weights
    tri_weight <- pmax(0, 1 - abs((y[i] - crm.atoms) / c0))
    
    # Compute log probabilities
    log_prob <- log(tri_weight) + (tht[i] * crm.atoms) + log(crm.jumps)
    
    # Normalize and sample
    prob <- exp(log_prob - max(log_prob))
    prob <- prob / sum(prob)
    
    # Sample z[i] based on triangular kernel weighting
    z[i] <- sample(crm.atoms, 1, prob = prob)
  }
  
  return(z)
}

z_sampler_epanK <- function(y, c0, crm.atoms, crm.jumps, tht) {
  n <- length(y)
  z <- numeric(n)
  
  for (i in 1:n) {
    # Compute Epanechnikov kernel weights
    epan_weight <- (3/4) * (1 - ((y[i] - crm.atoms) / c0)^2)
    epan_weight[abs(y[i] - crm.atoms) > c0] <- 0  # Ensure weights are 0 outside the range
    
    # Compute log probabilities
    log_prob <- log(epan_weight) + (tht[i] * crm.atoms) + log(crm.jumps)
    
    # Normalize and sample
    prob <- exp(log_prob - max(log_prob))
    prob <- prob / sum(prob)  
    
    # Sample z[i] based on Epanechnikov kernel weighting
    z[i] <- sample(crm.atoms, 1, prob = prob)
  }
  
  return(z)
}


# -----------------------------------------------------------------------------#
#                               Resample z                                     #
# -----------------------------------------------------------------------------#

resample_zstar <- function(z){
  z_table <- table(z)
  zstar   <- as.numeric(names(z_table))
  nstar   <- as.numeric(z_table)
  ## Write code for resampling zstar to avoid the ‘sticky clusters effect’
  
  return(list(zstar = zstar, nstar = nstar))
}

logpost_beta <- function(beta, linkinv, z, X, atoms, jumps, mu_beta, sigma_beta, h){
  n <- length(z)
  mu <- linkinv(X %*% beta)
  theta <- gldrm:::getTheta(spt = atoms, f0 = jumps, mu = mu, 
                            ySptIndex = NULL, sampprobs = NULL)$theta
  btheta <- b_theta(theta, atoms, jumps)
  log_post <- sum(theta * z + h - btheta) - sum((beta - mu_beta)^2 / (2 * sigma_beta^2)) # others const in beta
  return(log_post)
}

loglik <- function(linkinv, z, X, beta, atoms, jumps){
  n <- length(z)
  #mu <- exp(X %*% beta) / (1 + exp(X %*% beta))   
  # mu <- plogis(X %*% beta) # E old
  mu <- linkinv(X %*% beta)
  theta <- gldrm:::getTheta(spt = atoms, f0 = jumps, mu = mu, ySptIndex = NULL, 
                            sampprobs = NULL)$theta
  
  btheta <- b_theta(theta, atoms, jumps) 
  
  f0_z <- numeric(n)
  for (i in 1:n) {
    f0_z[i] <- sum(jumps[z[i] == atoms])
  }
  loglik <- sum(theta * z - btheta + log(f0_z))
  return(loglik)
}


log_post_u <- function(u, zstar, nstar, theta, alpha, min_y, max_y, eps, h) {
  # Number of grid points for the continuous part
  R <- 250
  
  # Construct a grid for integration over the continuous part (G_0, uniform on (0,1))
  #z_grid <- seq(eps, 1 - eps, length.out = R)
  z_grid <- seq(min_y + eps, max_y - eps, length.out = R)
  diff_z <- diff(z_grid)[1]  # uniform grid spacing
  
  # ---- Continuous part ----
  # For each grid point v in z_grid, compute in a stabilized manner:
  #    log(1 + sum_i u_i * exp(theta_i * v))
  cont_vals <- sapply(z_grid, function(v) {
    A <- log(u) + theta * v + h
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
    A <- log(u) + theta * v + h
    max_A <- max(A)
    S <- exp(max_A) * sum(exp(A - max_A))
    log1p(S)
  })
  
  sum_discrete <- sum(disc_vals * nstar)
  
  # Combine the continuous and discrete contributions
  neg_log_post <- alpha * integral_continuous + sum_discrete
  
  return(-neg_log_post)
}

sampler_u <- function(u, zstar, nstar, theta, alpha, delta, min_y, max_y, eps, h) {
  n <- length(u) # Get the length of u
  
  for (i in seq_len(n)) {
    u_star <- u
    # Sample from Gamma distribution: shape = delta, scale = u[i] / delta
    u_star[i] <- rgamma(1, shape = delta, scale = u[i] / delta)
    
    # Compute logQ_ratio: log(q(ui | ui_star)) - log(q(ui_star | ui))
    logQ_ratio <- dgamma(u[i], shape = delta, scale = u_star[i] / delta, log = TRUE) - 
      dgamma(u_star[i], shape = delta, scale = u[i] / delta, log = TRUE)
    
    # Compute logratio
    logratio <- log_post_u(u_star, zstar, nstar, theta, alpha, min_y, max_y, eps, h) - 
      log_post_u(u, zstar, nstar, theta, alpha, min_y, max_y, eps, h) + logQ_ratio
    
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

#' Compute b(theta) for discrete base measure
#'
#' Calculates \eqn{b(\theta) = \log \left( \int \exp(\theta z) dG_0(z) \right)} for a discrete base measure \eqn{G_0}
#' with support points \code{spt} and weights \code{f0}, using log-sum-exp trick for numerical stability.
#'
#' @param theta Numeric vector of parameter values.
#' @param spt Numeric vector of support points for the base measure.
#' @param f0 Numeric vector of weights or density values at each support point.
#'
#' @return Numeric vector of \eqn{b(\theta)} values, one for each \code{theta}.
#' @export

b_theta <- function(theta, spt, f0) { # Add arg H, type?
  log_f0_dz <- log(f0) #+ log(diff(spt)[1]) # log density values
  
  # Create the log weights matrix of dimensions length(theta) x length(spt)
  # (i, j)th element is theta[i]*spt[j] + log(f0[j])
  log_weights <- outer(theta, spt, "*") +
    matrix(log_f0_dz, nrow = length(theta), ncol = length(spt), byrow = TRUE)

  # Adding H, something like
  #log_weights <- outer(theta, spt, "*")  + H(...),
    #matrix(log_f0, nrow = length(theta), ncol = length(spt), byrow = TRUE)
  
  # Use log-sum-exp trick to compute b(theta) for each theta
  result <- apply(log_weights, 1, function(x) {
    m <- max(x)
    m + log(sum(exp(x - m)))
  })
  
  return(as.vector(result))  # Ensure the output is a vector
}

crm_sampler <- function(M, u, zstar, nstar, theta, alpha,
                        min_y, max_y, eps, h = 0) {
  N <- 3001
  R <- 3001
  s <- -log(seq(exp(-eps), exp(-5e-4), length.out = N))

  tht <- theta

  if (length(h) == 1L) {
    h <- rep(h, length(u))
  }

  z <- seq(min_y + eps, max_y - eps, length.out = R)

  # log_terms[i, j] = log(u_i) + h_i + theta_i * z_j
  log_terms <- sweep(tht %o% z, 1, log(u) + h, "+")

  psi_z <- apply(log_terms, 2, function(log_vec) {
    max_val <- max(log_vec)
    exp(max_val + log(sum(exp(log_vec - max_val))))
  })

  dz <- (max_y - min_y - 2 * eps) / (R - 1)

  fnS <- numeric(length(s))

  for (j in seq_along(s)) {
    log_vals <- -(1 + psi_z) * s[j] - log(s[j])

    max_val <- max(log_vals)
    log_sum_exp <- max_val + log(sum(exp(log_vals - max_val)))

    log_integral <- log(dz) + log_sum_exp
    fnS[j] <- exp(log_integral)
  }

  ds <- diff(s)
  h_trap <- (fnS[-N] + fnS[-1]) / 2
  Nv <- rev(cumsum(rev(ds * h_trap)))
  Nv <- Nv * alpha
  Nv <- c(Nv, 0)

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

    log_temp <- -(1 + psi_z) * RJ[m]
    max_log <- max(log_temp)

    p <- exp(log_temp - max_log)
    p <- p / sum(p)

    cumsum_p <- cumsum(p)
    RL[m] <- z[min(which(cumsum_p > xi_rl))]
  }

  # Fixed-location jumps
  log_terms <- sweep(tht %o% zstar, 1, log(u) + h, "+")

  psi_star <- sapply(seq_along(zstar), function(j) {
    max_val <- max(log_terms[, j])
    exp(max_val + log(sum(exp(log_terms[, j] - max_val))))
  })

  Jstar <- rgamma(length(zstar), shape = nstar, rate = psi_star + 1)

  list(
    RL = RL,
    RJ = RJ,
    zstar = zstar,
    Jstar = Jstar
  )
}

# ------------------------------------------------------------------
# Compute observation-specific marginal CDFs F_{x_j}(y_j) using a
# CRM-induced mixture of truncated uniform kernels.
#
# The mixture weights are covariate-dependent and vary across
# observations through theta[j].
#
# Arguments:
#   y         : numeric vector of outcomes (length n)
#   crm.atoms : vector of CRM atoms (z_l)
#   crm.jumps : vector of CRM jump sizes (J_l)
#   theta   : numeric vector of covariate effects (length n)
#   c0        : half-width of the uniform kernel
#   min_y     : global lower bound of the outcome support
#   max_y     : global upper bound of the outcome support
#
# Returns:
#   Fy        : numeric vector of marginal CDF values F_{x_j}(y_j)
# ------------------------------------------------------------------
marginal_cdf_unif_kde <- function(y, crm.atoms, crm.jumps, theta, c0,
                         min_y, max_y) {

  eps <- 1e-6
  n <- length(y)

  ## 1. Construct truncated kernel bounds (shared across observations)
  lower <- crm.atoms - c0
  upper <- crm.atoms + c0

  lower[lower < min_y] <- min_y + eps
  upper[upper > max_y] <- max_y - eps

  width <- upper - lower

  ## 2. Evaluate CDF observation by observation
  Fy <- numeric(n)

  for (j in seq_len(n)) {

    ## Observation-specific mixture weights (log-scale)
    logw <- log(crm.jumps) + theta[j] * crm.atoms
    logw <- logw - max(logw)   # log-sum-exp stabilization
    w <- exp(logw)
    w <- w / sum(w)

    ## Mixture CDF at y_j
    yj <- y[j]
    Fy[j] <- sum(
      w * ifelse(
        yj < lower, 0,
        ifelse(
          yj > upper, 1,
          (yj - lower) / width
        )
      )
    )
  }

  Fy
}



# ------------------------------------------------------------------
# Compute the log-density of a Gaussian copula with an exchangeable
# (compound symmetry) correlation structure.
#
# The copula density is evaluated at a vector of marginal probabilities
# u = (u_1, ..., u_n), where u_j = F_x(y_ij).
#
# Arguments:
#   u   : numeric vector of marginal CDF values in (0, 1)
#   rho : dependence parameter (pairwise correlation)
#
# Returns:
#   logc : scalar value of log c_rho(u)
# ------------------------------------------------------------------
log_copula_gaussian <- function(u, rho, corr) {
  n <- length(u)

  # avoid boundary
  eps <- 1e-10
  u <- pmin(pmax(u, eps), 1 - eps)

  ## Transform to latent Gaussian scale
  z <- qnorm(u)

  # Exchangeable correlation structure
  # Sigma <- matrix(rho, n, n)
  # diag(Sigma) <- 1

  # AR(1) correlation
  # Sigma <- outer(1:n, 1:n, function(i, j) rho^abs(i - j))

  if (corr == "ex") {
    Sigma <- matrix(rho, n, n)
    diag(Sigma) <- 1
  } else {
    Sigma <- outer(1:n, 1:n, function(i, j) rho^abs(i - j))
  }

  Sigma_inv <- solve(Sigma)
  logdet <- determinant(Sigma, logarithm = TRUE)$modulus

  ## Gaussian copula log-density
  # log c(u) = -1/2 log|Sigma| - 1/2 z' (Sigma^{-1} - I) z
  val <- -0.5 * logdet -
    0.5 * t(z) %*% (Sigma_inv - diag(n)) %*% z

  as.numeric(val)
}

# ------------------------------------------------------------------
# Compute observation-level Gaussian copula log-likelihood contributions.

#
# Marginal CDFs F_j = F_{x_j}(y_j) are computed internally via a
# CRM-induced mixture of truncated uniform kernels using
# marginal_cdf_unif_kde().

#
# Observations are coupled within groups using a Gaussian copula with
# exchangeable correlation parameter rho.
#
# For a group i with n_i observations, each observation j in that group
# receives an equal share of the group copula log-density:
#
#   log c_rho(F_{i1}, ..., F_{in_i}) / n_i
#
# Arguments:
#   y           : numeric vector of outcomes
#   group_index : vector indicating group membership for each y
#   rho         : copula dependence parameter
#   crm.atoms   : vector of CRM atoms (z_l)
#   crm.jumps   : vector of CRM jump sizes (J_l)
#   theta       : numeric vector of covariate effects (length = length(y))
#   c0          : half-width of the uniform kernel
#   min_y       : global lower bound of the outcome support
#   max_y       : global upper bound of the outcome support
#
# Returns:
#   logcop_obs  : numeric vector of observation-level copula
#                 log-likelihood contributions (same length as y)
# ------------------------------------------------------------------
log_copula_contribution_by_obs <- function(y, group_index, rho, corr,
                                           crm.atoms, crm.jumps, theta, c0,
                                           min_y, max_y) {
  
  ## 1. Compute marginal CDFs internally
  Fy <- marginal_cdf_unif_kde(y = y,
                     crm.atoms = crm.atoms, crm.jumps = crm.jumps, theta = theta, c0 = c0,
                     min_y = min_y, max_y = max_y)

  ## 2. Initialize observation-level output
  logcop_obs <- numeric(length(y))
  
  ## 3. Compute copula contribution per group, then distribute to observations
  groups <- unique(group_index)

  for (g in groups) {
    idx <- which(group_index == g)
    u <- Fy[idx]
    n_i <- length(u)

    logc_group <- log_copula_gaussian(u, rho, corr)
    
    # Assign equal share to each observation in the group
    logcop_obs[idx] <- logc_group / n_i
  }

  logcop_obs
}

# ------------------------------------------------------------------
# Log posterior kernel for the copula dependence parameter rho
# (simplified interface).
#
# This function evaluates the log posterior of rho under the
# simplified joint distribution
#
#   pi(rho | y, lambda, beta) ∝ ∏_i c(y_i, lambda, beta, rho) p(rho),
#
# where all marginal likelihood terms are handled elsewhere and the
# input 'logc_obs' already contains the log copula density contributions
# evaluated at the current value of rho.
#
# Arguments:
#   rho      : scalar copula dependence parameter, must lie in (0, 1)
#   logc_obs : numeric vector of log copula contributions whose sum
#              equals log ∏_i c(y_i, lambda, beta, rho)
#
# Returns:
#   Scalar log posterior value log pi(rho | ·), up to an additive constant
# ------------------------------------------------------------------
log_post_rho <- function(rho, logc_obs, shape1, shape2) {
  
  ## Enforce support of the copula parameter
  if (rho <= 0 || rho >= 1) return(-Inf)
  
  ## Copula log-likelihood contribution
  loglik <- sum(logc_obs)
  
  ## Beta prior on rho
  logprior <- dbeta(rho, shape1 = shape1, shape2 = shape2, log = TRUE)
  
  ## Log posterior kernel
  loglik + logprior
}

# calc_exceed <- function(z_atoms, J_jumps, thr_vec) {
#     Jsum <- sum(J_jumps)
#     if (!is.finite(Jsum) || Jsum <= 0) return(rep(NA_real_, length(thr_vec)))
#     w <- J_jumps / Jsum  # temporary normalization for probability calc
#     vapply(thr_vec, function(t) sum(w[z_atoms >= t]), numeric(1))
#   }
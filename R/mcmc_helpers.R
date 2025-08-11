## Required packages -----
## gldrm, mvtnorm

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
    thetaStart = thtst
  )
  tht <- out$theta
  bpr2 <- out$bPrime2
  btht <- apply(exp(outer(tht, spt, "*")), 1, function(row)
    log(sum(row * f0)))
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


#' Title: function for calculating \{f0(y_i), i = 1(1)n\}
#'
#' @param y (vector) response variable
#' @param spt (vector) support, \{s_j, j = 1(1)l\}, of response variable y
#' @param f0 (vector) \{f0(s_j), j=1(1)l\}: free parameter, f0(.) centering density
#'
#' @return (vector) \{f0(y_i), i=1(1)n\}
#' @keywords internal

f0y <- function(y, spt, f0) {
  #n <- length(y)
  #f0_y <- numeric(n)
  #for (i in 1:n) {
    #f0_y[i] <- sum(f0[y[i] == spt])
  #}
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

  for (j in 1:l) {
    pr_f0 <- cr_f0
    pr_f0[j] <- rbeta(1, cr_dir_parm[j], sum(cr_dir_parm[-j]))
    pr_f0[-j] <- pr_f0[-j] * (1 - pr_f0[j]) / sum(pr_f0[-j])

    if (sum(pr_f0 < 1e-3) == 0) {
      out <- tht_sol(spt, pr_f0, cr_mu, cr_tht)
      pr_tht <- out$tht
      pr_btht <- out$btht
      pr_bpr2 <- out$bpr2
      pr_dir_parm <- dir_parm(y, pr_tht, pr_btht, dir_pr_parm, ind_mt)
      pr_f0y <- f0y(y, spt, pr_f0)

      pr_llik <- sum(pr_tht * y - pr_btht + log(pr_f0y))
      cr_llik <- sum(cr_tht * y - cr_btht + log(cr_f0y))

      pr_pf0 <- extraDistr::ddirichlet(pr_f0, dir_pr_parm, log = T)                   # prior probability (proposal)
      cr_pf0 <- extraDistr::ddirichlet(cr_f0, dir_pr_parm, log = T)                   # prior probability (current)
      cr_qf0 <- extraDistr::ddirichlet(cr_f0, pr_dir_parm, log = T)
      pr_qf0 <- extraDistr::ddirichlet(pr_f0, cr_dir_parm, log = T)

      alp <- min(0, (pr_llik - cr_llik + pr_pf0 - cr_pf0 + cr_qf0 - pr_qf0))

      if (log(runif(1)) < alp) {
        cr_f0   <- pr_f0
        cr_tht  <- pr_tht
        cr_btht <- pr_btht
        cr_f0y  <- pr_f0y
        cr_bpr2 <- pr_bpr2
        acc_f0  <- TRUE
      } else {
        acc_f0 <- FALSE
      }
    } else {
      acc_f0 <- FALSE
    } # Closes line 134
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
                              sb) {
  n <- dim(X)[1]
  l <- length(spt)

  pr_bt <- as.vector(mvtnorm::rmvnorm(1, mean = cr_bt, sigma = cr_Sig))
  pr_mu <- as.numeric(linkinv(X %*% pr_bt))          # Updated for general link

  if (sum(spt[1] <= pr_mu & pr_mu <= spt[l]) == n) {
    out <- tht_sol(spt, cr_f0, pr_mu, cr_tht)
    pr_bpr2 <- out$bpr2
    pr_tht <- out$tht
    pr_btht <- out$btht
    pr_Sig <- Sigma_beta(X, pr_mu, pr_bpr2, rho, linkfun, mu.eta)

    pr_llik <- sum(pr_tht * y - pr_btht)
    cr_llik <- sum(cr_tht * y - cr_btht)
    # pr_pbt <- dmvnorm(pr_bt, log = T) # Prior probability (for proposal)
    pr_pbt <- dmvnorm(pr_bt,            # Updated
                      mean = mb,
                      sigma = diag(sb),
                      log = T)
    # cr_pbt <- dmvnorm(cr_bt, log = T) # Prior probability (for current)
    cr_pbt <- dmvnorm(cr_bt,            # Updated
                      mean = mb,
                      sigma = diag(sb),
                      log = T)
    cr_qbt <- dmvnorm(cr_bt,
                      mean = pr_bt,
                      sigma = pr_Sig,
                      log = T)
    pr_qbt <- dmvnorm(pr_bt,
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
                                 sb) {
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
      # pr_pbt <- dnorm(pr_bt[j], log = T)  # Prior probability (for proposal)?
      pr_pbt <- dnorm(pr_bt[j],             # Updated
                      mean = mb[j],
                      sigma = sqrt(sb[j]),
                      log = T)
      # cr_pbt <- dnorm(cr_bt[j], log = T)  # Prior probability (for current)?
      cr_pbt <- dnorm(cr_bt[j],             # Updated
                      mean = mb[j],
                      sigma = sqrt(sb[j]),
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

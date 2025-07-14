simData <- function(n, p, betaMax, link, ySim) {
  beta <- c(0, seq(-betaMax, betaMax, length.out = p))
  X    <- cbind(1, matrix(rnorm(n * p), nrow = n))
  colnames(X) <- paste0("X", seq_len(ncol(X)))
  eta  <- as.vector(X %*% beta)
  mu   <- link$linkinv(eta)
  y    <- ySim(n, mu)
  list(x = X, y = y, eta = eta, mu = mu, beta = beta)
}


test_that("bspgldrm and bspgldrmFit match", {
  set.seed(100)
  formula <- y ~ X2
  lf <- stats::make.link("log")
  ySim <- function(n, mu) rpois(n, mu)
  sim  <- simData(n=50, p=1, betaMax=0.7, link=lf, ySim=ySim)
  data <- data.frame(sim$x, y = sim$y)

  # Setup for bspgldrmFit
  mf <- stats::model.frame(formula, data)
  X  <- stats::model.matrix(attr(mf, "terms"), mf)
  attributes(X)[c("assign", "contrasts")] <- NULL
  y  <- stats::model.response(mf, type = "numeric")

  set.seed(100)
  ctrl <- bspgldrm:::bspgldrm.control(mu0=mean(sim$y))
  m1   <- bspgldrm(formula=formula, data=data,
                   bspgldrmControl=ctrl)

  set.seed(100)
  m2   <- bspgldrm:::bspgldrmFit(formula=formula, data=data, X=X, y=y, link=lf,
                                 mb=NULL, sb=NULL, dir_pr_parm=NULL,
                                 bspgldrmControl=ctrl,
                                 thetaControl=gldrm:::theta.control())

  expect_equal(m1$samples, m2$samples, tolerance=1e-3, ignore_attr=TRUE)
})

test_that("bspgldrm matches intercept-only (empirical distribution) model", {
  set.seed(100)
  ySim <- function(n, mu) rpois(n, 1)
  lf <- stats::make.link("identity")
  sim <- simData(n=100, p=0, betaMax=0, link=lf, ySim=ySim)
  data <- data.frame(sim$x, y=sim$y)

  m1 <- as.vector(table(data$y)) / length(data$y)
  m2 <- bspgldrm(y ~ X1 - 1, data=data, link="identity",
                 bspgldrmControl=bspgldrm:::bspgldrm.control(save=10000))  # link function doesn't matter with no covariates

  ## this is an intercept-only model, so all observations have fitted mean equal to mean(y)
  expect_equal(mean(data$y), mean(m2$samples$beta), # Compare to posterior mean
               tolerance=1e-2, ignore_attr=TRUE)
  ## f0 should match response frequency table
  expect_equal(m1, colMeans(m2$samples$f0), tolerance=1e-2) # Compare to posterior mean
})

test_that("bspgldrm matches logistic regression", {
  set.seed(100)
  ySim <- function(n, mu) rbinom(n, 1, mu)
  sim <- simData(n=100, p=5, betaMax=0.5, make.link("logit"), ySim)
  data <- data.frame(sim$x, y=sim$y)
  m1 <- glm(y ~ X2 + X3 + X4 + X5 + X6 - 1, data=data, family=binomial(link="logit"))
  m2 <- bspgldrm(y ~ X2 + X3 + X4 + X5 + X6 - 1, data=data, link="logit",
                 bspgldrmControl=bspgldrm.control(burnin=10000, thin=20, save=25000))

  ## BSPGLDRM should match logistic regression coefficient estimates
  ## (semiparametric model is identical to fully parametric in this case)
  expect_equal(as.vector(coef(m1)), colMeans(m2$samples$beta), tolerance=1e-1, ignore_attr=TRUE)
})

#test_that("Can handle muHat on boundary of spt", {
  #n <- 10
  #y <- rep(c(0, 1), each=n/2)
  #x <- cbind(1, y)

  #m1 <- bspgldrm(y ~ x - 1, data=NULL, link="identity")
  #expect_equal(colMeans(m1$samples$beta), c(0, 1), tolerance=1e-2, ignore_attr=TRUE)

  #lf <- stats::make.link("logit")
  #m2 <- bspgldrm(y ~ x - 1, data=NULL, link=lf)
  #eta <- x %*% colMeans(m2$samples$beta)
  #mu <- lf$linkinv(eta)
  #expect_equal(mu, y, tolerance=1e-2)
#})

test_that("Can handle singular covariate matrix", {
  n <- 10
  y <- rep(c(0, 1), each=n/2)
  x <- matrix(1, nrow=n, ncol=2)

  ctrl <- bspgldrm::bspgldrm.control(burnin=1000, thin=20, save=5000)
  m1 <- bspgldrm(y ~ x - 1, data=NULL, link="identity")
  expect_equal(colMeans(m1$samples$beta), c(.5, NA), tolerance=1e-2, ignore_attr=TRUE)

  m2 <- bspgldrm(y ~ x - 1, data=NULL, link="logit")
  expect_equal(colMeans(m2$samples$beta), c(0, NA), tolerance=1e-2, ignore_attr=TRUE)
})

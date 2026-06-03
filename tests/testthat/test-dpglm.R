test_data_dpglm <- function() {
  data.frame(
    y = c(0.42, 0.45, 0.48, 0.50, 0.52, 0.55, 0.58, 0.60),
    x = c(-1.5, -1, -0.5, 0, 0.25, 0.5, 1, 1.5)
  )
}

test_control_dpglm <- function(save = 3, burnin = 2, thin = 1) {
  z <- seq(0.01, 0.99, length.out = 25)
  J <- rep(1 / length(z), length(z))
  dat <- test_data_dpglm()
  list(
    burnin = burnin,
    thin = thin,
    save = save,
    spt = c(0, 1),
    c0 = 0.15,
    betaStart = c(0, 0),
    varbetaStart = diag(0.01, 2),
    thetaStart = rep(0, nrow(dat)),
    crmStart = list(z.tld = z, J.tld = J),
    seed = 832
  )
}

test_that("dpglm.control returns expected control object", {
  ctrl <- dpglm.control(
    burnin = 2,
    thin = 1,
    save = 3,
    spt = c(0, 1),
    seed = 832
  )

  expect_s3_class(ctrl, "dpglmControl")
  expect_equal(ctrl$burnin, 2)
  expect_equal(ctrl$thin, 1)
  expect_equal(ctrl$save, 3)
  expect_equal(ctrl$spt, c(0, 1))
  expect_equal(ctrl$seed, 832)
})

test_that("dpglm fits and returns expected object structure", {
  dat <- test_data_dpglm()
  ctrl <- do.call(dpglm.control, test_control_dpglm())

  fit <- suppressWarnings(
    dpglm(
      y ~ x,
      data = dat,
      link = "logit",
      dpglmControl = ctrl
    )
  )

  expect_s3_class(fit, "dpglm")
  expect_named(fit$samples, c("z", "beta", "crm", "theta"))

  expect_equal(dim(fit$samples$beta), c(3, 2))
  expect_equal(dim(fit$samples$theta), c(3, nrow(dat)))
  expect_equal(dim(fit$samples$z), c(3, nrow(dat)))

  expect_equal(colnames(fit$samples$beta), c("Intercept", "x"))

  expect_equal(nrow(fit$samples$crm), 3)
  expect_true(all(c("z.tld", "J.tld") %in% colnames(fit$samples$crm)))

  expect_true(is.numeric(fit$beta_acceptance))
  expect_true(is.numeric(fit$crm_acceptance))
  expect_true(fit$beta_acceptance >= 0 && fit$beta_acceptance <= 1)
  expect_true(fit$crm_acceptance >= 0 && fit$crm_acceptance <= 1)
})

test_that("dpglm respects burnin, thin, and save in stored samples", {
  dat <- test_data_dpglm()
  ctrl <- do.call(dpglm.control, test_control_dpglm(save = 4, burnin = 2, thin = 2))

  fit <- suppressWarnings(
    dpglm(y ~ x, data = dat, link = "logit", dpglmControl = ctrl)
  )

  expect_equal(nrow(fit$samples$beta), 4)
  expect_equal(nrow(fit$samples$theta), 4)
  expect_equal(nrow(fit$samples$z), 4)
  expect_equal(nrow(fit$samples$crm), 4)

  expect_false(anyNA(fit$samples$beta))
  expect_false(anyNA(fit$samples$theta))
  expect_false(anyNA(fit$samples$z))
})

test_that("dpglm supports burnin equal to zero", {
  dat <- test_data_dpglm()
  ctrl <- do.call(dpglm.control, test_control_dpglm(save = 3, burnin = 0, thin = 1))

  fit <- suppressWarnings(
    dpglm(y ~ x, data = dat, link = "logit", dpglmControl = ctrl)
  )

  expect_equal(nrow(fit$samples$beta), 3)
  expect_false(anyNA(fit$samples$beta))
})

test_that("summary.dpglm returns expected components", {
  dat <- test_data_dpglm()
  ctrl <- do.call(dpglm.control, test_control_dpglm())

  fit <- suppressWarnings(
    dpglm(y ~ x, data = dat, link = "logit", dpglmControl = ctrl)
  )

  s <- summary(fit)

  expect_s3_class(s, "summary.dpglm")
  expect_named(s, c("meta", "beta", "crm"))

  expect_equal(rownames(s$beta), c("Intercept", "x"))
  expect_equal(rownames(s$crm), c("min", "q25", "median", "q75", "max"))

  expect_true(all(c("Estimate", "Est.Error", "l-95% CrI", "u-95% CrI") %in% colnames(s$beta)))
  expect_true(all(c("Estimate", "Est.Error", "l-95% CrI", "u-95% CrI") %in% colnames(s$crm)))
})

test_that("plot_dpglm runs for beta and crm", {
  dat <- test_data_dpglm()
  ctrl <- do.call(dpglm.control, test_control_dpglm())

  fit <- suppressWarnings(
    dpglm(y ~ x, data = dat, link = "logit", dpglmControl = ctrl)
  )

  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)

  expect_silent(plot_dpglm(fit, what = "beta", pars = "Intercept"))
  expect_silent(plot_dpglm(fit, what = "crm"))
})

test_that("dpglm.control rejects invalid MCMC settings", {
  expect_error(dpglm.control(burnin = -1), "burn-in")
  expect_error(dpglm.control(burnin = 1.5), "burn-in")

  expect_error(dpglm.control(thin = 0), "Thin")
  expect_error(dpglm.control(thin = 1.5), "Thin")

  expect_error(dpglm.control(save = 0), "saved")
  expect_error(dpglm.control(save = 1.5), "saved")

  expect_error(dpglm.control(spt = 1), "spt")
  expect_error(dpglm.control(spt = c(1, 0)), "spt")
  expect_error(dpglm.control(spt = c("a", "b")), "spt")
})

test_that("dpglm rejects mb and Sb with invalid dimensions", {
  dat <- test_data_dpglm()

  expect_error(
    suppressWarnings(dpglm(
      y ~ x,
      data = dat,
      link = "logit",
      dpglmControl = do.call(dpglm.control, c(
        test_control_dpglm(),
        list(mb = c(0, 0, 0))
      ))
    )),
    "length\\(mb\\) must match"
  )

  expect_error(
    suppressWarnings(dpglm(
      y ~ x,
      data = dat,
      link = "logit",
      dpglmControl = do.call(dpglm.control, c(
        test_control_dpglm(),
        list(Sb = c(1, 1, 1))
      ))
    )),
    "length\\(Sb\\) must match"
  )

  expect_error(
    suppressWarnings(dpglm(
      y ~ x,
      data = dat,
      link = "logit",
      dpglmControl = do.call(dpglm.control, c(
        test_control_dpglm(),
        list(Sb = c(1, 0))
      ))
    )),
    "Sb must be positive"
  )
})

test_that("dpglm rejects singular design matrices", {
  dat <- data.frame(
    y = c(0.42, 0.45, 0.48, 0.50, 0.52, 0.55),
    x = rep(1, 6)
  )

  ctrl <- do.call(dpglm.control, test_control_dpglm())

  expect_error(
    suppressWarnings(
      dpglm(y ~ x, data = dat, link = "logit", dpglmControl = ctrl)
    ),
    "singular"
  )
})

test_that("dpglm rejects invalid custom link objects", {
  dat <- test_data_dpglm()
  ctrl <- do.call(dpglm.control, test_control_dpglm())

  bad_link <- list(
    linkfun = function(mu) mu,
    linkinv = function(eta) eta
  )

  expect_error(
    suppressWarnings(
      dpglm(y ~ x, data = dat, link = bad_link, dpglmControl = ctrl)
    ),
    "link"
  )
})

test_that("dpglm rejects non-vectorized custom links", {
  dat <- test_data_dpglm()
  ctrl <- do.call(dpglm.control, test_control_dpglm())

  bad_link <- list(
    linkfun = function(mu) mu[1],
    linkinv = function(eta) eta,
    mu.eta = function(eta) rep(1, length(eta))
  )

  expect_error(
    suppressWarnings(
      dpglm(y ~ x, data = dat, link = bad_link, dpglmControl = ctrl)
    ),
    "vectorized"
  )
})

test_that("dpglm and dpglmFit agree", {
  dat <- test_data_dpglm()
  ctrl <- do.call(dpglm.control, test_control_dpglm())

  fit_public <- suppressWarnings(
    dpglm(y ~ x, data = dat, link = "logit", dpglmControl = ctrl)
  )

  mf <- stats::model.frame(y ~ x, dat)
  X <- stats::model.matrix(attr(mf, "terms"), mf)
  attributes(X)[c("assign", "contrasts")] <- NULL
  y <- stats::model.response(mf, type = "numeric")

  init <- list(
    beta = ctrl$betaStart,
    varbeta = ctrl$varbetaStart,
    theta = ctrl$thetaStart,
    crm = list(
      z.tld = ctrl$crmStart$z.tld,
      J.tld = ctrl$crmStart$J.tld
    )
  )

  fit_internal <- suppressWarnings(
    dpglmFit(
      formula = y ~ x,
      data = dat,
      X = X,
      y = y,
      link = stats::make.link("logit"),
      spt = ctrl$spt,
      mu0 = mean(dat$y),
      init = init,
      dpglmControl = ctrl,
      thetaControl = theta.control()
    )
  )

  expect_equal(fit_public$samples$beta, fit_internal$samples$beta)
  expect_equal(fit_public$samples$theta, fit_internal$samples$theta)
  expect_equal(fit_public$beta_acceptance, fit_internal$beta_acceptance)
  expect_equal(fit_public$crm_acceptance, fit_internal$crm_acceptance)
})

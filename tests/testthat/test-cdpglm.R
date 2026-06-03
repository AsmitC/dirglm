test_data_cdpglm <- function() {
  data.frame(
    y = c(0.42, 0.45, 0.48, 0.50, 0.52, 0.55, 0.58, 0.60),
    x = c(-1.5, -1, -0.5, 0, 0.25, 0.5, 1, 1.5),
    id = rep(1:4, each = 2)
  )
}

test_control_cdpglm <- function(save = 3, burnin = 2, thin = 1) {
  dat <- test_data_cdpglm()

  z <- seq(0.01, 0.99, length.out = 25)
  J <- rep(1 / length(z), length(z))

  list(
    burnin = burnin,
    thin = thin,
    save = save,
    spt = c(0, 1),
    mu0 = mean(dat$y),
    c0 = 0.15,
    M = 5,
    betaStart = c(0, 0),
    varbetaStart = diag(0.01, 2),
    thetaStart = rep(0, nrow(dat)),
    crmStart = list(z.tld = z, J.tld = J),
    rhoStart = 0.3,
    rho_proposal_sd = 0.05,
    rho_prior_shape = c(8, 2),
    corr = "ex",
    seed = 832
  )
}

test_that("cdpglm.control returns expected control object", {
  ctrl <- cdpglm.control(
    burnin = 2,
    thin = 1,
    save = 3,
    spt = c(0, 1),
    rhoStart = 0.3,
    corr = "ex",
    seed = 832
  )

  expect_s3_class(ctrl, "cdpglmControl")
  expect_equal(ctrl$burnin, 2)
  expect_equal(ctrl$thin, 1)
  expect_equal(ctrl$save, 3)
  expect_equal(ctrl$spt, c(0, 1))
  expect_equal(ctrl$rhoStart, 0.3)
  expect_equal(ctrl$corr, "ex")
  expect_equal(ctrl$seed, 832)
})

test_that("cdpglm fits and returns expected object structure", {
  dat <- test_data_cdpglm()
  ctrl <- do.call(cdpglm.control, test_control_cdpglm())

  fit <- suppressWarnings(
    cdpglm(
      y ~ x,
      data = dat,
      group_index = dat$id,
      link = "logit",
      cdpglmControl = ctrl
    )
  )

  expect_s3_class(fit, "cdpglm")
  expect_named(fit$samples, c("z", "beta", "theta", "h", "rho", "crm"))

  expect_equal(dim(fit$samples$beta), c(3, 2))
  expect_equal(dim(fit$samples$theta), c(3, nrow(dat)))
  expect_equal(dim(fit$samples$h), c(3, nrow(dat)))
  expect_equal(dim(fit$samples$z), c(3, nrow(dat)))
  expect_equal(length(fit$samples$rho), 3)

  expect_equal(colnames(fit$samples$beta), c("Intercept", "x"))
  expect_equal(nrow(fit$samples$crm), 3)
  expect_true(all(c("z.tld", "J.tld") %in% colnames(fit$samples$crm)))

  expect_true(is.numeric(fit$beta_acceptance))
  expect_true(is.numeric(fit$crm_acceptance))
  expect_true(is.numeric(fit$rho_acceptance))

  expect_true(fit$beta_acceptance >= 0 && fit$beta_acceptance <= 1)
  expect_true(fit$crm_acceptance >= 0 && fit$crm_acceptance <= 1)
  expect_true(fit$rho_acceptance >= 0 && fit$rho_acceptance <= 1)
})

test_that("cdpglm respects burnin, thin, and save in stored samples", {
  dat <- test_data_cdpglm()
  ctrl <- do.call(cdpglm.control, test_control_cdpglm(save = 4, burnin = 2, thin = 2))

  fit <- suppressWarnings(
    cdpglm(
      y ~ x,
      data = dat,
      group_index = dat$id,
      link = "logit",
      cdpglmControl = ctrl
    )
  )

  expect_equal(nrow(fit$samples$beta), 4)
  expect_equal(nrow(fit$samples$theta), 4)
  expect_equal(nrow(fit$samples$h), 4)
  expect_equal(nrow(fit$samples$z), 4)
  expect_equal(nrow(fit$samples$crm), 4)
  expect_equal(length(fit$samples$rho), 4)

  expect_false(anyNA(fit$samples$beta))
  expect_false(anyNA(fit$samples$theta))
  expect_false(anyNA(fit$samples$h))
  expect_false(anyNA(fit$samples$z))
  expect_false(anyNA(fit$samples$rho))
})

test_that("summary.cdpglm returns expected components", {
  dat <- test_data_cdpglm()
  ctrl <- do.call(cdpglm.control, test_control_cdpglm())

  fit <- suppressWarnings(
    cdpglm(
      y ~ x,
      data = dat,
      group_index = dat$id,
      link = "logit",
      cdpglmControl = ctrl
    )
  )

  s <- summary(fit)

  expect_s3_class(s, "summary.cdpglm")
  expect_named(s, c("meta", "beta", "rho", "reference"))

  expect_equal(rownames(s$beta), c("Intercept", "x"))
  expect_equal(rownames(s$rho), "rho")
  expect_equal(rownames(s$reference), c("min", "q25", "median", "q75", "max"))

  expect_true(all(c("Estimate", "Est.Error", "l-95% CrI", "u-95% CrI") %in% colnames(s$beta)))
  expect_true(all(c("Estimate", "Est.Error", "l-95% CrI", "u-95% CrI") %in% colnames(s$rho)))
  expect_true(all(c("Estimate", "Est.Error", "l-95% CrI", "u-95% CrI") %in% colnames(s$reference)))
})

test_that("plot_cdpglm runs for beta, crm, and rho", {
  dat <- test_data_cdpglm()
  ctrl <- do.call(cdpglm.control, test_control_cdpglm())

  fit <- suppressWarnings(
    cdpglm(
      y ~ x,
      data = dat,
      group_index = dat$id,
      link = "logit",
      cdpglmControl = ctrl
    )
  )

  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)

  expect_silent(plot_cdpglm(fit, what = "beta", pars = "Intercept"))
  expect_silent(plot_cdpglm(fit, what = "crm"))
  expect_silent(plot_cdpglm(fit, what = "rho"))
})

test_that("cdpglm.control rejects invalid MCMC settings", {
  expect_error(cdpglm.control(burnin = -1), "burn-in")
  expect_error(cdpglm.control(burnin = 1.5), "burn-in")

  expect_error(cdpglm.control(thin = 0), "Thin")
  expect_error(cdpglm.control(thin = 1.5), "Thin")

  expect_error(cdpglm.control(save = 0), "saved")
  expect_error(cdpglm.control(save = 1.5), "saved")

  expect_error(cdpglm.control(corr = "bad"), "corr")
  expect_error(cdpglm.control(rhoStart = 0), "rho")
  expect_error(cdpglm.control(rhoStart = 1), "rho")
  expect_error(cdpglm.control(rho_proposal_sd = 0), "rho_proposal_sd")
  expect_error(cdpglm.control(rho_prior_shape = c(8, -2)), "rho_prior_shape")
  expect_error(cdpglm.control(rhoStart = NA_real_), "rhoStart")
  expect_error(cdpglm.control(rho_proposal_sd = NA_real_), "rho_proposal_sd")
  expect_error(cdpglm.control(rho_prior_shape = c(8, NA_real_)), "rho_prior_shape")
})

test_that("cdpglm rejects singular design matrices", {
  dat <- data.frame(
    y = c(0.42, 0.45, 0.48, 0.50, 0.52, 0.55, 0.58, 0.60),
    x = rep(1, 8),
    id = rep(1:4, each = 2)
  )

  ctrl <- do.call(cdpglm.control, test_control_cdpglm())

  expect_error(
    suppressWarnings(
      cdpglm(
        y ~ x,
        data = dat,
        group_index = dat$id,
        link = "logit",
        cdpglmControl = ctrl
      )
    ),
    "singular"
  )
})

test_that("cdpglm rejects invalid group_index", {
  dat <- test_data_cdpglm()
  ctrl <- do.call(cdpglm.control, test_control_cdpglm())

  expect_error(
    suppressWarnings(
      cdpglm(
        y ~ x,
        data = dat,
        group_index = dat$id[-1],
        link = "logit",
        cdpglmControl = ctrl
      )
    ),
    "group_index"
  )
})

test_that("cdpglm rejects invalid custom link objects", {
  dat <- test_data_cdpglm()
  ctrl <- do.call(cdpglm.control, test_control_cdpglm())

  bad_link <- list(
    linkfun = function(mu) mu,
    linkinv = function(eta) eta
  )

  expect_error(
    suppressWarnings(
      cdpglm(
        y ~ x,
        data = dat,
        group_index = dat$id,
        link = bad_link,
        cdpglmControl = ctrl
      )
    ),
    "link"
  )
})

test_that("cdpglm rejects non-vectorized custom links", {
  dat <- test_data_cdpglm()
  ctrl <- do.call(cdpglm.control, test_control_cdpglm())

  bad_link <- list(
    linkfun = function(mu) mu[1],
    linkinv = function(eta) eta,
    mu.eta = function(eta) rep(1, length(eta))
  )

  expect_error(
    suppressWarnings(
      cdpglm(
        y ~ x,
        data = dat,
        group_index = dat$id,
        link = bad_link,
        cdpglmControl = ctrl
      )
    ),
    "vectorized"
  )
})

test_that("cdpglm rejects mb and Sb with invalid dimensions", {
  dat <- test_data_cdpglm()

  expect_error(
    suppressWarnings(cdpglm(
      y ~ x,
      data = dat,
      group_index = dat$id,
      link = "logit",
      cdpglmControl = do.call(cdpglm.control, c(
        test_control_cdpglm(),
        list(mb = c(0, 0, 0))
      ))
    )),
    "length\\(mb\\) must match"
  )

  expect_error(
    suppressWarnings(cdpglm(
      y ~ x,
      data = dat,
      group_index = dat$id,
      link = "logit",
      cdpglmControl = do.call(cdpglm.control, c(
        test_control_cdpglm(),
        list(Sb = c(1, 1, 1))
      ))
    )),
    "length\\(Sb\\) must match"
  )

  expect_error(
    suppressWarnings(cdpglm(
      y ~ x,
      data = dat,
      group_index = dat$id,
      link = "logit",
      cdpglmControl = do.call(cdpglm.control, c(
        test_control_cdpglm(),
        list(Sb = c(1, 0))
      ))
    )),
    "Sb must be positive"
  )
})

test_that("cdpglm and cdpglmFit agree", {
  dat <- test_data_cdpglm()
  ctrl <- do.call(cdpglm.control, test_control_cdpglm())

  fit_public <- suppressWarnings(
    cdpglm(
      y ~ x,
      data = dat,
      group_index = dat$id,
      link = "logit",
      cdpglmControl = ctrl
    )
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
    ),
    rho = ctrl$rhoStart
  )

  fit_internal <- suppressWarnings(
    cdpglmFit(
      formula = y ~ x,
      data = dat,
      X = X,
      y = y,
      group_index = dat$id,
      link = stats::make.link("logit"),
      spt = ctrl$spt,
      mu0 = ctrl$mu0,
      init = init,
      cdpglmControl = ctrl,
      thetaControl = theta.control()
    )
  )

  expect_equal(fit_public$samples$beta, fit_internal$samples$beta)
  expect_equal(fit_public$samples$theta, fit_internal$samples$theta)
  expect_equal(fit_public$samples$h, fit_internal$samples$h)
  expect_equal(fit_public$samples$rho, fit_internal$samples$rho)

  expect_equal(fit_public$beta_acceptance, fit_internal$beta_acceptance)
  expect_equal(fit_public$crm_acceptance, fit_internal$crm_acceptance)
  expect_equal(fit_public$rho_acceptance, fit_internal$rho_acceptance)
})

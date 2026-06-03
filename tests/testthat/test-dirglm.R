test_data <- function() {
  data.frame(
    y = c(0, 0.2, 0.4, 0.6, 0.8, 1),
    x = c(-2, -1, 0, 1, 2, 3)
  )
}

test_that("dirglm.control returns expected control object", {
  ctrl <- dirglm.control(
    burnin = 2,
    thin = 1,
    save = 3,
    seed = 832
  )

  expect_s3_class(ctrl, "dirglmControl")
  expect_equal(ctrl$burnin, 2)
  expect_equal(ctrl$thin, 1)
  expect_equal(ctrl$save, 3)
  expect_equal(ctrl$seed, 832)
})

test_that("dirglm fits and returns expected object structure", {
  dat <- test_data()

  fit <- suppressWarnings(
    dirglm(
      y ~ x,
      data = dat,
      link = "logit",
      dirglmControl = dirglm.control(
        burnin = 2,
        thin = 1,
        save = 3,
        spt = sort(unique(dat$y)),
        seed = 832
      )
    )
  )

  expect_s3_class(fit, "dirglm")
  expect_named(fit$samples, c("beta", "f0"))
  expect_equal(dim(fit$samples$beta), c(3, 2))
  expect_equal(dim(fit$samples$f0), c(3, length(unique(dat$y))))

  expect_equal(colnames(fit$samples$beta), c("Intercept", "x"))
  expect_equal(colnames(fit$samples$f0), paste0("f0_", seq_along(unique(dat$y))))

  expect_true(is.numeric(fit$beta_acceptance))
  expect_true(is.numeric(fit$f0_acceptance))
  expect_true(fit$beta_acceptance >= 0 && fit$beta_acceptance <= 1)
  expect_true(fit$f0_acceptance >= 0 && fit$f0_acceptance <= 1)
})

test_that("dirglm respects burnin, thin, and save in stored samples", {
  dat <- test_data()

  fit <- suppressWarnings(dirglm(
    y ~ x,
    data = dat,
    link = "logit",
    dirglmControl = dirglm.control(
      burnin = 2,
      thin = 2,
      save = 4,
      spt = sort(unique(dat$y)),
      seed = 832
    )
  ))

  expect_equal(nrow(fit$samples$beta), 4)
  expect_equal(nrow(fit$samples$f0), 4)
})

test_that("summary.dirglm returns expected components", {
  dat <- test_data()

  fit <- suppressWarnings(dirglm(
    y ~ x,
    data = dat,
    link = "logit",
    dirglmControl = dirglm.control(
      burnin = 2,
      thin = 1,
      save = 3,
      spt = sort(unique(dat$y)),
      seed = 832
    )
  ))

  s <- summary(fit)

  expect_s3_class(s, "summary.dirglm")
  expect_named(s, c("meta", "beta", "f0"))

  expect_equal(rownames(s$beta), c("Intercept", "x"))
  expect_equal(rownames(s$f0), paste0("f0_", seq_along(unique(dat$y))))

  expect_true(all(c("Estimate", "Est.Error", "l-95% CrI", "u-95% CrI") %in% colnames(s$beta)))
  expect_true(all(c("Estimate", "Est.Error", "l-95% CrI", "u-95% CrI") %in% colnames(s$f0)))
})

test_that("plot_dirglm runs for beta and f0", {
  dat <- test_data()

  fit <- suppressWarnings(dirglm(
    y ~ x,
    data = dat,
    link = "logit",
    dirglmControl = dirglm.control(
      burnin = 2,
      thin = 1,
      save = 3,
      spt = sort(unique(dat$y)),
      seed = 832
    )
  ))

  grDevices::pdf(NULL)
  expect_silent(plot_dirglm(fit, what = "beta", pars = "Intercept"))
  expect_silent(plot_dirglm(fit, what = "f0"))
  grDevices::dev.off()
})

test_that("dirglm.control rejects invalid MCMC settings", {
  expect_error(dirglm.control(burnin = -1), "burn-in")
  expect_error(dirglm.control(burnin = 1.5), "burn-in")

  expect_error(dirglm.control(thin = 0), "Thin")
  expect_error(dirglm.control(thin = 1.5), "Thin")

  expect_error(dirglm.control(save = 0), "saved")
  expect_error(dirglm.control(save = 1.5), "saved")
})

test_that("dirglm.control rejects invalid tuning parameters", {
  expect_error(dirglm.control(rho = 0), "rho")
  expect_error(dirglm.control(rho = 1.5), "rho")

  expect_error(dirglm.control(joint.update = "yes"), "joint.update")
})

test_that("dirglm rejects singular design matrices", {
  # x constant, linear dependence
  dat <- data.frame(
    y = c(0, 0.2, 0.4, 0.6, 0.8, 1),
    x = rep(1, 6)
  )

  expect_error(
    suppressWarnings(dirglm(
      y ~ x,
      data = dat,
      link = "logit",
      dirglmControl = dirglm.control(
        burnin = 2,
        thin = 1,
        save = 3,
        spt = sort(unique(dat$y)),
        seed = 832
      ))
    ),
    "singular"
  )
})

test_that("dirglm rejects invalid mu0 values", {
  # mu0 not in [0, 1]
  dat <- test_data()

  expect_error(
    suppressWarnings(dirglm(
      y ~ x,
      data = dat,
      link = "logit",
      dirglmControl = dirglm.control(
        burnin = 2,
        thin = 1,
        save = 3,
        spt = sort(unique(dat$y)),
        mu0 = 1.5,
        seed = 832
      ))
    ),
    "mu0"
  )
})

test_that("dirglm rejects invalid custom link objects", {
  dat <- test_data()

  # Missing mu.eta
  bad_link <- list(
    linkfun = function(mu) mu,
    linkinv = function(eta) eta
  )

  expect_error(
    suppressWarnings(dirglm(
      y ~ x,
      data = dat,
      link = bad_link,
      dirglmControl = dirglm.control(
        burnin = 2,
        thin = 1,
        save = 3,
        spt = sort(unique(dat$y)),
        seed = 832
      ))
    ),
    "link"
  )
})

test_that("dirglm rejects non-vectorized custom links", {
  dat <- test_data()

  bad_link <- list(
    linkfun = function(mu) mu[1],
    linkinv = function(eta) eta,
    mu.eta = function(eta) rep(1, length(eta))
  )

  expect_error(
    suppressWarnings(dirglm(
      y ~ x,
      data = dat,
      link = bad_link,
      dirglmControl = dirglm.control(
        burnin = 2,
        thin = 1,
        save = 3,
        spt = sort(unique(dat$y)),
        seed = 832
      ))
    ),
    "vectorized"
  )
})

test_that("dirglm warns and forces joint update for non-diagonal Sb", {
  dat <- test_data()

  # Non diagonal prior covariance
  Sb <- matrix(c(1, 0.25, 0.25, 1), nrow = 2)

  expect_warning(
    withCallingHandlers(
      dirglm(
        y ~ x,
        data = dat,
        link = "logit",
        dirglmControl = dirglm.control(
          burnin = 2,
          thin = 1,
          save = 3,
          spt = sort(unique(dat$y)),
          mb = c(0, 0),
          Sb = Sb,
          joint.update = FALSE,
          seed = 832
        )
      ),
      warning = function(w) {
        if (grepl("Acceptance ratios are high", conditionMessage(w))) {
          invokeRestart("muffleWarning")
        }
      }
    ),
    "Forcing joint update"
  )
})

test_that("dirglm rejects mb with wrong length", {
  dat <- test_data()

  expect_error(
    dirglm(
      y ~ x,
      data = dat,
      link = "logit",
      dirglmControl = dirglm.control(
        burnin = 2,
        thin = 1,
        save = 3,
        spt = sort(unique(dat$y)),
        mb = c(0, 0, 0),
        seed = 832
      )
    ),
    "length\\(mb\\) must match"
  )
})

test_that("dirglm rejects Sb with wrong dimension", {
  dat <- test_data()

  expect_error(
    dirglm(
      y ~ x,
      data = dat,
      link = "logit",
      dirglmControl = dirglm.control(
        burnin = 2,
        thin = 1,
        save = 3,
        spt = sort(unique(dat$y)),
        Sb = diag(3),
        seed = 832
      )
    ),
    "dim\\(Sb\\) must match"
  )
})

test_that("dirglm rejects Sb with nonpositive diagonal entries", {
  dat <- test_data()

  expect_error(
    dirglm(
      y ~ x,
      data = dat,
      link = "logit",
      dirglmControl = dirglm.control(
        burnin = 2,
        thin = 1,
        save = 3,
        spt = sort(unique(dat$y)),
        Sb = diag(c(1, 0)),
        seed = 832
      )
    ),
    "Sb must be positive definite"
  )
})

test_that("dirglm rejects invalid dir_pr_parm", {
  dat <- test_data()

  expect_error(
    dirglm(
      y ~ x,
      data = dat,
      link = "logit",
      dirglmControl = dirglm.control(
        burnin = 2,
        thin = 1,
        save = 3,
        spt = sort(unique(dat$y)),
        dir_pr_parm = c(1, 1),
        seed = 832
      )
    ),
    "dir_pr_parm must be positive"
  )

  expect_error(
    dirglm(
      y ~ x,
      data = dat,
      link = "logit",
      dirglmControl = dirglm.control(
        burnin = 2,
        thin = 1,
        save = 3,
        spt = sort(unique(dat$y)),
        dir_pr_parm = rep(0, length(unique(dat$y))),
        seed = 832
      )
    ),
    "dir_pr_parm must be positive"
  )
})

test_that("dirglm and dirglmFit agree", {
  dat <- test_data()

  ctrl <- dirglm.control(
    burnin = 2,
    thin = 1,
    save = 3,
    spt = sort(unique(dat$y)),
    mu0 = mean(dat$y),
    betaStart = c(0, 0),
    f0Start = rep(1 / length(unique(dat$y)), length(unique(dat$y))),
    seed = 832
  )

  fit_public <- suppressWarnings(
    dirglm(y ~ x, data = dat, link = "logit", dirglmControl = ctrl)
  )

  mf <- stats::model.frame(y ~ x, dat)
  X <- stats::model.matrix(attr(mf, "terms"), mf)
  attributes(X)[c("assign", "contrasts")] <- NULL
  y <- stats::model.response(mf, type = "numeric")

  fit_internal <- suppressWarnings(
    dirglmFit(
      formula = y ~ x,
      data = dat,
      X = X,
      y = y,
      link = stats::make.link("logit"),
      mu0 = ctrl$mu0,
      spt = ctrl$spt,
      init = list(
        beta = ctrl$betaStart,
        f0 = ctrl$f0Start
      ),
      dirglmControl = ctrl,
      thetaControl = theta.control()
    )
  )

  expect_equal(fit_public$samples$beta, fit_internal$samples$beta)
  expect_equal(fit_public$samples$f0, fit_internal$samples$f0)
  expect_equal(fit_public$beta_acceptance, fit_internal$beta_acceptance)
  expect_equal(fit_public$f0_acceptance, fit_internal$f0_acceptance)
})

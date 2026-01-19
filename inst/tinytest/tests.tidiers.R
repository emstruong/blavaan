# Tests for tidy and glance methods

if (requireNamespace("rstan", quietly = TRUE)) {
  load(system.file("testdata", "sysdata.rda", package = "blavaan"))
  library("lavaan", quietly = TRUE)

  # Test tidy.blavaan with Stan model
  tidy_stan <- tidy.blavaan(fitstan)
  expect_true(inherits(tidy_stan, "data.frame"))
  expect_true("term" %in% names(tidy_stan))
  expect_true("estimate" %in% names(tidy_stan))
  expect_true("std.error" %in% names(tidy_stan))
  expect_true("conf.low" %in% names(tidy_stan))
  expect_true("conf.high" %in% names(tidy_stan))
  expect_true("rhat" %in% names(tidy_stan))
  expect_true("ess" %in% names(tidy_stan))
  expect_true("prior" %in% names(tidy_stan))
  expect_true(nrow(tidy_stan) > 0)

  # Test tidy without confidence intervals
  tidy_noci <- tidy.blavaan(fitstan, conf.int = FALSE)
  expect_false("conf.low" %in% names(tidy_noci))
  expect_false("conf.high" %in% names(tidy_noci))

  # Test tidy without rhat and ess
  tidy_nodiag <- tidy.blavaan(fitstan, rhat = FALSE, ess = FALSE)
  expect_false("rhat" %in% names(tidy_nodiag))
  expect_false("ess" %in% names(tidy_nodiag))

  # Test tidy with standardized estimates
  tidy_std <- tidy.blavaan(fitstan, standardized = TRUE)
  expect_true("std.all" %in% names(tidy_std))

  # Test glance.blavaan with Stan model
  glance_stan <- glance.blavaan(fitstan)
  expect_true(inherits(glance_stan, "data.frame"))
  expect_equal(nrow(glance_stan), 1)
  expect_true("npar" %in% names(glance_stan))
  expect_true("nobs" %in% names(glance_stan))
  expect_true("ngroups" %in% names(glance_stan))
  expect_true("estimator" %in% names(glance_stan))
  expect_true("nchains" %in% names(glance_stan))

  # Check that fit measures are included
  expect_true("ppp" %in% names(glance_stan) || "looic" %in% names(glance_stan))

  # Test glance with fit indices (this is computationally expensive)
  # Only test if blavFitIndices works
  bfi_works <- tryCatch({
    blavFitIndices(fitstan, fit.measures = "BRMSEA")
    TRUE
  }, error = function(e) FALSE)

  if (bfi_works) {
    glance_bfi <- glance.blavaan(fitstan, fit.indices = "BRMSEA")
    expect_true("BRMSEA" %in% names(glance_bfi))
  }

  # Test with classic Stan model
  tidy_stanc <- tidy.blavaan(fitstanc)
  expect_true(inherits(tidy_stanc, "data.frame"))
  expect_true(nrow(tidy_stanc) > 0)

  glance_stanc <- glance.blavaan(fitstanc)
  expect_true(inherits(glance_stanc, "data.frame"))
  expect_equal(nrow(glance_stanc), 1)

  # Test with JAGS model if available
  if (requireNamespace("runjags", quietly = TRUE)) {
    tidy_jags <- tidy.blavaan(fitjags)
    expect_true(inherits(tidy_jags, "data.frame"))
    expect_true(nrow(tidy_jags) > 0)

    glance_jags <- glance.blavaan(fitjags)
    expect_true(inherits(glance_jags, "data.frame"))
    expect_equal(nrow(glance_jags), 1)
  }
}

res <- simulate_ms_pva(
 leaving_prob = 0.2,
 reaching_prob = 0.5,
 growth_rate_means = c(0.043, -0.002, 0),
 growth_rate_vars = c(0.051, 0.041, 0.051),
 initial_pops = c(70, 26, 33),
 growth_rate_corrs = {matrix(
   c(	1.000,	0.995,   0.896,
   0.995,	1.000,   0.938,
   0.896,	0.938,   1.000),
   nrow = 3,
   ncol = 3,
   byrow = TRUE
 )},
 K = c(286, 60, 58),
 quasi_extinction_thresholds = c(20, 20, 20),
 with_progress_bar = FALSE
)

test_that("Verifying agreement with the MatLab code for Multi-site code",{

  expect_equal(
    res$lam0,
    1.0320,
    tolerance = .001
  )

  expect_equal(
    mean(res$stochLam),
    0.99875,
    tolerance = .001
  )

  expect_equal(
    mean(res$logLam),
    -0.0012787,
    tolerance = .005
  )

  expect_equal(
    sd(res$logLam),
    0.0070218,
    tolerance = .01
  )


})

test_that("Parameter loading from file",{

  expect_warning(
    params <- calculate_params_from_file(
      system.file("extdata", "PolarBear_Stirling2004.csv", package = "msPVA")
    )
  )

  expect_equal(
    c(params$growth_rate_corrs),
    c(1.0000000, -0.6207679, -0.6207679,  1.0000000),
    tolerance = .0001
  )

  expect_equal(
    params$initial_pops,
    c(201, 140),
    check.attributes = FALSE
  )

  attributes(params$growth_rate_means) <- NULL

  expect_equal(
    params$growth_rate_means,
    c(0.02362944, 0.01572059),
    tolerance = .001
  )

  attributes(params$growth_rate_vars) <- NULL

  expect_equal(
    params$growth_rate_vars,
    c(0.1292640, 0.1639641),
    tolerance = .0001
  )

})
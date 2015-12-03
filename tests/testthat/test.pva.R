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

test_that("Our results agree with the MatLab code for Multi-site code",{

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

test_that("Parameters are correctly loading from a file",{

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

test_that("Our results agree with the old R code for single-site simulations",{
  # Old R code was run for 10 000 iterations

  res <- simulate_ss_pva(
    lambdas = c(
      0.808510638,
      0.828947368,
      1,
      1.047619048,
      0.833333333,
      1.777777778
    ),
    initial_pop = 136,
    n_years = 100,
    n_runs = 3000
  )

  attributes(res$decline_risk) <- NULL
  expect_equal(
    res$decline_risk,
    c(0.49566),
    tolerance = .05
  )

  attributes(res$extinction_risk) <- NULL
  expect_equal(
    res$extinction_risk,
    c(0.24153),
    tolerance = .06
  )

})

test_that("K is correctly accounted for in single-site simulations",{
  res <- simulate_ss_pva(
    lambdas = c(
      0.808510638,
      0.828947368,
      1,
      1.047619048,
      0.833333333,
      1.777777778
    ),
    initial_pop = 136,
    n_years = 100,
    n_runs = 2000,
    K = 200
  )

  expect_lte(
    max(res$final_pops),
    200
  )

  res <- simulate_ss_pva(
    lambdas = c(
      0.808510638,
      0.828947368,
      1,
      1.047619048,
      0.833333333,
      1.777777778
    ),
    initial_pop = 136,
    n_years = 100,
    n_runs = 2000
  )

  expect_gte(
    max(res$final_pops),
    200
  )
})

test_that("The multiple ways of calling simulate_ss_pva are correctly verified",{

  expect_error(
    res <- simulate_ss_pva(
      initial_pop = 20
    )
  )

  expect_error(
    res <- simulate_ss_pva(
      initial_pop = 20,
      lambdas = c(1,2,3),
      log_lambdas = c(1,2,3)
    )
  )

})

test_that("Stochastic model for a single population",{

  res <- simulate_ss_pva(
    growth_rate_means = 0.015,
    growth_rate_vars = .041,
    initial_pop = 26,
    quasi_extinction_thresholds = 20,
    n_years = 50,
    n_runs = 1000
  )

  # popbio returns 0.7624 for this :
  # max(ex<-extCDF(0.015, 0.041, Nc=26, Ne=20))

  expect_equal(
    res$extinction_risk,
    0.7624,
    tolerance = .05
  )

  # x <- rnorm(n = 20, mean = 0.015, sd = sqrt(0.041))

  res <- simulate_ss_pva(
    log_lambdas = c(-0.0503626618483076, -0.0316522478682412, -0.205890697055539,
                    -0.0407897021414208, 0.151024474883104, -0.141017433696716, 0.105149579850484,
                    0.104087724782143, 0.18297223483855, -0.0114362475217544, 0.358832191971109,
                    0.144782148474211, 0.176846281421632, -0.34330593611314, 0.337176958392917,
                    -0.178642571388776, -0.298701626025117, -0.0984130919939361,
                    -0.214594528811423, 0.431114720192186),
    initial_pop = 26,
    quasi_extinction_thresholds = 20,
    n_years = 50,
    n_runs = 1000
  )

  expect_equal(
    res$extinction_risk,
    0.7624,
    tolerance = .04
  )

})
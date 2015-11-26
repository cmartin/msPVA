context("Verifying agreement with the MatLab code")

res <- simulate_pva(
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
 quasi_extinction_thresholds = c(20, 20, 20)
)

expect_equal(
  res$lam0,
  1.0320,
  tolerance = .01
)

expect_equal(
  mean(res$stochLam),
  0.99875,
  tolerance = .01
)
expect_equal(
  res$lam0,
  1.0320,
  tolerance = .01
)
expect_equal(
  mean(res$logLam),
  -0.0012787,
  tolerance = .01
)
expect_equal(
  sd(res$logLam),
  0.0070218,
  tolerance = .01
)
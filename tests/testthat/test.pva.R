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

all.equal(
  res$lam0,
  1.0320,
  tolerance = .001
)

all.equal(
  mean(res$stochLam),
  0.99875,
  tolerance = .001
)

all.equal(
  res$lam0,
  1.0320,
  tolerance = .001
)

all.equal(
  mean(res$logLam),
  -0.0012787,
  tolerance = .005
)

all.equal(
  sd(res$logLam),
  0.0070218,
  tolerance = .005
)

context("Parameter loading from file")

params <- calculate_params_from_file(
  system.file("extdata", "PolarBear_Stirling2004.csv", package = "PopulationViabilityAnalysis")
)

all.equal(
  c(params$growth_rate_corrs),
  c(1.0000000, -0.6207679, -0.6207679,  1.0000000),
  tolerance = .0001
)

all.equal(
  params$initial_pops,
  c(201, 140),
  check.attributes = FALSE
)

all.equal(
  params$growth_rate_means,
  c( 0.02362944, 0.01572059),
  tolerance = .0001,
  check.attributes = FALSE
)

all.equal(
  params$growth_rate_vars,
  c(0.1292640, 0.1639641),
  tolerance = .0001,
  check.attributes = FALSE
)

## @knitr ssExample

# From a precalculated mean log-lambdas

res <- simulate_ss_pva(
 growth_rate_means = 0.043,
 growth_rate_vars = 0.051,
 initial_pops = 70,
 K = 286,
 quasi_extinction_thresholds = 20,
 n_years = 50,
 n_runs = 100
)

print(res)
hist(res)

# From a vector of log-lambdas
res <- simulate_ss_pva(
  log_lambdas = c(-0.0503626618483076, -0.0316522478682412, -0.205890697055539,
                  -0.0407897021414208, 0.151024474883104, -0.141017433696716, 0.105149579850484,
                  0.104087724782143, 0.18297223483855),
  initial_pops = 70,
  K = 286,
  quasi_extinction_thresholds = 20,
  n_years = 50,
  n_runs = 100
)

# Or from a vector of lambdas
res <- simulate_ss_pva(
  lambdas = c(
    0.808510638,
    0.828947368,
    1,
    1.047619048,
    0.833333333,
    1.777777778
  ),
  initial_pops = 70,
  K = 286,
  quasi_extinction_thresholds = 20,
  n_years = 50,
  n_runs = 100
)

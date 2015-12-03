## @knitr ssExample

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

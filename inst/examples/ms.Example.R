## @knitr msExample

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
 n_years = 50,
 n_runs = 100,
 with_progress_bar = FALSE

)

print(res)
hist(res)
plot(res)

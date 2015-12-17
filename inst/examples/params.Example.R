## @knitr paramsExample

params <- calculate_params_from_file(
 system.file("extdata", "PolarBear_Stirling2004.csv", package = "msPVA")
)
res <- do.call(
 "simulate_ms_pva",
 c(
   params,
   list(
     K = c(300,200),
     leaving_prob = 0.1,
     reaching_prob = 0.7,
     quasi_extinction_thresholds = c(20, 20),
     n_years = 50,
     n_runs = 100,
     with_progress_bar = FALSE
   )
 )
)
print(res)

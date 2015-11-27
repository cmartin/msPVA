#' Run a stochastic simulation for a count-based multi-site population viability analysis
#'
#' See Morris & Doak 2002, p. 423-431
#' @param n_years A number. How many years should we simulate?
#' @param n_runs A number. How many simulations should we do?
#' @param leaving_prob A number. The probability of an individual leaving a pop. for another pop.
#' @param reaching_prob A number. The probability of an individual reaching safely another pop.
#' @param growth_rate_means A vector of mean log(lambdas), one for each pop.
#' @param growth_rate_vars A vector of variance of log(lambdas), one for each pop.
#' @param growth_rate_corrs A correlation matrix of log(lambdas) between pops.
#' @param initial_pops A vector of initial population sizes.
#' @param K A vector of maximum numbers in each pop. (carrying capacity)
#' @param quasi_extinction_thresholds A vector of near extinction threshold for each pops.
#' @return A list-based S3 object of class \code{PVARes} containing elements CDFExt, lam0, logLam and stochLam.
#' @export
#' @examples
#' res <- simulate_ms_pva(
#'  leaving_prob = 0.2,
#'  reaching_prob = 0.5,
#'  growth_rate_means = c(0.043, -0.002, 0),
#'  growth_rate_vars = c(0.051, 0.041, 0.051),
#'  initial_pops = c(70, 26, 33),
#'  growth_rate_corrs = {matrix(
#'    c(	1.000,	0.995,   0.896,
#'    0.995,	1.000,   0.938,
#'    0.896,	0.938,   1.000),
#'    nrow = 3,
#'    ncol = 3,
#'    byrow = TRUE
#'  )},
#'  K = c(286, 60, 58),
#'  quasi_extinction_thresholds = c(20, 20, 20),
#'  n_years = 50,
#'  n_runs = 100
#' )
#'
#' print(res)
#' hist(res)
#' plot(res)
simulate_ms_pva <- function(
  n_years = 100,
  n_runs = 1000,
  leaving_prob = 0,
  reaching_prob = 1,
  growth_rate_means,
  growth_rate_vars,
  initial_pops,
  growth_rate_corrs,
  K,
  quasi_extinction_thresholds
) {

  if (!all(
    c(
      length(growth_rate_means),
      length(growth_rate_vars),
      length(initial_pops),
      dim(growth_rate_corrs),
      length(K)
    ) == length(quasi_extinction_thresholds)
  )) stop("Please ensure all population arguments (growth_rate_means, K, initial_pops, etc.) are of the same dimensions.")

  ################# SETUP ###################

  n_pops <- length(growth_rate_means) 		# number of local populations
  np <- n_pops + 2  		    # number of vital rates
  # adding movement to the growth_rate_means vectors

  growth_rate_means <- c(growth_rate_means,leaving_prob,reaching_prob)
  growth_rate_vars <- c(growth_rate_vars,0,0)

  initial_pops <- matrix(initial_pops, ncol = 1)

  # adding movement to the vital rates correlation matrix
  growth_rate_corrs <- rbind(
    cbind(growth_rate_corrs, matrix(0,ncol = 2, nrow = n_pops)),
    matrix(0, nrow = 2, ncol = np)
  )
  diag(growth_rate_corrs) <- 1

  D <- matrix(0, nrow = np, ncol = np)
  diag(D) <- rev(eigen(growth_rate_corrs)$values)

  W <- eigen(growth_rate_corrs)$vectors
  W <- W[, np:1]

  M12 <- W %*% (sqrt(abs(D))) %*% t(W)

  vrs <- c(
    growth_rate_means[1:n_pops] + 0.5*growth_rate_vars[1:n_pops],
    growth_rate_means[(n_pops + 1):(n_pops + 2)]) # the addition of 0.5*variances

  mx <- .makemx(vrs, leaving_prob, reaching_prob)

  lam1 <- matrix(0, nrow = n_pops, ncol = n_pops)
  diag(lam1) <- rev(eigen(mx)$values)
  uu <- eigen(mx)$vectors
  uu <- uu[, n_pops:1]

  lam1s <- apply(lam1, 2, max)
  lam0 <- max(lam1s)

  PrExt = rep(0,n_years) # extinction time tracker
  logLam = rep(0, n_runs) # tracker of stoch-log-lambda values
  stochLam = rep(0, n_runs) # tracker of stochastic lambda values

  pb <- utils::txtProgressBar(
    min = 0,
    max = n_runs,
    initial = 0,
    char = "=",
    width = 75,
    style = 3
  )

  ################# SIMULATION LOOPS ###################

  for (xx in 1:n_runs) {

    utils::setTxtProgressBar(pb, xx)

    nt <- initial_pops; # start at initial population size vector
    extinct <- rep(0, n_pops) # vector of extinction recorders

    for (tt in 1:n_years) {
      m <- rnorm(np)

      rawelems <- t(M12 %*% m) #correlated str. normals
      vrs <- growth_rate_means + sqrt(growth_rate_vars)*rawelems
      mx <- .makemx(vrs, leaving_prob, reaching_prob)

      nt = mx %*% nt	 # multiply by the population vector
      nt = pmin(nt, K) # applying population cap

      # these two ifs check for extinction and count it
      if (sum(extinct) < n_pops	) {
        for (nn in 1:n_pops) {
          if (nt[nn] <= quasi_extinction_thresholds[nn]) extinct[nn] <- 1
        }
        if (sum(extinct) == n_pops) PrExt[tt] <- PrExt[tt] + 1
      } # if

    } # tt

    logLam[xx] <- (1/n_years)*log(sum(nt)/sum(initial_pops))
    stochLam[xx] <- (sum(nt)/sum(initial_pops)) ^ (1/n_years)

  } # xx

  l = list(
    CDFExt = cumsum(PrExt/n_runs),
    lam0 = lam0,
    logLam = logLam,
    stochLam = stochLam
  )
  class(l) <- append(class(l),"PVARes")
  return(l)

}

.makemx = function(vrs, leaving_prob, reaching_prob, n_pops = (length(vrs) - 2)) {
  m <- matrix(leaving_prob*reaching_prob, ncol = n_pops, nrow = n_pops)
  diag(m) <- (1 - leaving_prob)*exp(vrs[1:n_pops])
  return(m)
}

#' @export
print.PVARes <- function(x, ...) {
    cat(paste('This is the deterministic lambda value : ', x$lam0))
    cat(paste('\r\nAnd this is the mean stochastic lambda : ', mean(x$stochLam)))
    cat('\r\nBelow is mean and standard deviation of log lambda :\r\n')
    cat(paste(mean(x$logLam), sd(x$logLam)))
}

#' @export
plot.PVARes <- function(x, ...) {
  plot(
    x$CDFExt,
    main = "Extinction time CDF",
    xlab = "Years",
    ylab = "Cumulative probability of quasi-extinction",
    type = "l"
  )
}

#' @export
hist.PVARes <- function(x, ...) {
  hist(x$logLam, main = "logLams")
}

#' Extract population parameters from a data file
#'
#' This function expects a n_pop+1 columns layout with column headers as the first line.
#' First column should be years, then counts for each population. E.g. :\cr\cr
#' Year;Pop1;Pop2\cr
#' 2010;35;26\cr
#' 2011;26;41\cr
#' 2012;33;35
#'
#' This function calculates log(lambda) means and variances the standard way,
#' unless there are gaps in the time series.
#' In the latter case, regression analysis is used.
#'
#' See Morris & Doak 2002, p.64-69 for calculation details
#'
#' @param data_file The path to the file
#' @return A \code{list} containing initial_pops, growth_rate_corrs,
#' growth_rate_means and growth_rate_vars parameters. Please see \link{simulate_ms_pva} for a description of each parameter
#' @export
#' @examples
#' params <- calculate_params_from_file(
#'  system.file("extdata", "PolarBear_Stirling2004.csv", package = "msPVA")
#' )
#' do.call(
#'  "simulate_ms_pva",
#'  c(
#'    params,
#'    list(
#'      K = c(300,200),
#'      leaving_prob = 0.1,
#'      reaching_prob = 0.7,
#'      quasi_extinction_thresholds = c(20, 20),
#'      n_years = 50,
#'      n_runs = 100
#'    )
#'  )
#' )
calculate_params_from_file <- function(data_file) {
  df <- read.csv(
    data_file,
    sep = ";",
    stringsAsFactors = FALSE
  )

  # Calculate yearly growth rates (log-lambdas; Morris & Doak 2002, p.64-65)
  # This calculation assumes regular time intervals

  log_lambdas <- apply(df[, 2:ncol(df)], 2, function(pop_sizes){
    log(pop_sizes[2:length(pop_sizes)]) - log(pop_sizes[1:(length(pop_sizes) - 1)])
  })

  if (isTRUE(all.equal(df[1:(nrow(df) - 1),1], df[2:nrow(df),1] - 1))) {

    # Calculate means and variances
    growth_rate_means <- apply(log_lambdas, 2, mean)
    growth_rate_vars <- apply(log_lambdas, 2, var)

  } else {

    warning("There were gaps in the time series.
  Mean and variance were computed by linear regression.
  Correlation matrix was calculated assuming regular intervals.")

    # Morris & Doak 2002, p. 66-69
    time_intervals <- sqrt(df[2:nrow(df),1] - df[1:(nrow(df) - 1),1])

    dependent_var <- apply(df[, 2:ncol(df)], 2, function(pop_sizes){
      log(pop_sizes[2:length(pop_sizes)]) - log(pop_sizes[1:(length(pop_sizes) - 1)])
    }) / time_intervals

    growth_rate_means <- apply(dependent_var, 2, function(y){
      summary(
        lm(y ~ 0 + time_intervals)
      )$coefficients[1]
    })

    growth_rate_vars <- apply(dependent_var, 2, function(y){
      summary(
        lm(y ~ 0 + time_intervals)
      )$sigma ^ 2
    })

  }

  list(
    # Correlation between growth rates of each pops
    growth_rate_corrs = cor(log_lambdas),
    # Use last year of data as initial population sizes
    initial_pops = unlist(df[nrow(df), 2:ncol(df)]),

    growth_rate_means = growth_rate_means,
    growth_rate_vars = growth_rate_vars
  )


}

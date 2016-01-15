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
#' @param with_progress_bar A boolean value, to show a text-based progress bar while the simulation runs
#' @return A list-based S3 object of class \code{msPVARes} containing elements CDFExt, lam0, logLam and stochLam.
#' @export
#' @example /inst/examples/ms.Example.R
simulate_ms_pva <- function(
  n_years = 100,
  n_runs = 1000,
  leaving_prob = 0,
  reaching_prob = 1,
  with_progress_bar = TRUE,
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

  if (with_progress_bar) {
    pb <- utils::txtProgressBar(
      min = 0,
      max = n_runs,
      initial = 0,
      char = "=",
      width = 75,
      style = 3
    )
  }

  ################# SIMULATION LOOPS ###################

  for (xx in 1:n_runs) {

    if (with_progress_bar) {
      utils::setTxtProgressBar(pb, xx)
    }

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
  class(l) <- append(class(l),"msPVARes")
  return(l)

}

.makemx = function(vrs, leaving_prob, reaching_prob, n_pops = (length(vrs) - 2)) {
  m <- matrix(leaving_prob*reaching_prob, ncol = n_pops, nrow = n_pops)
  diag(m) <- (1 - leaving_prob)*exp(vrs[1:n_pops])
  return(m)
}

#' @export
print.msPVARes <- function(x, ...) {
    cat(paste('This is the deterministic lambda value : ', x$lam0))
    cat(paste('\r\nAnd this is the mean stochastic lambda : ', mean(x$stochLam)))
    cat('\r\nBelow is mean and standard deviation of log lambda :\r\n')
    cat(paste(mean(x$logLam), sd(x$logLam)))
}

#' @export
plot.msPVARes <- function(x, ...) {
  plot(
    x$CDFExt,
    main = "Extinction time CDF",
    xlab = "Years",
    ylab = "Cumulative probability of quasi-extinction",
    type = "l"
  )
}

#' @export
hist.msPVARes <- function(x, ...) {
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
#' @example /inst/examples/params.Example.R
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

#' Run a stochastic simulation for a count-based single-site population viability analysis
#'
#' This function uses a stochastic algorithm if log-lambdas or mean log-lambdas and variance
#' are provided.
#'
#' This function is not yet vectorized, so provide a single population per run.
#'
#' @param n_years A number. How many years should we simulate?
#' @param n_runs A number. How many simulations should we do?
#' @param initial_pops Number. Initial population sizes.
#' @param K A number. Maximum population size (carrying capacity).
#' @param quasi_extinction_thresholds A number. Near extinction threshold for the population.
#' @param ... Either log_lambdas (a vector of log(lambdas) calculated as in Morris & Doak 2002, p.64-65) or growth_rate_means and growth_rate_vars (Numbers: Mean and variance of log(lambdas))
#' @return A list-based S3 object of class \code{ssPVARes} containing elements final_pops (vector), n_years, n_runs, initial_pops, decline_risk and extinction_risk.
#' @export
#' @example /inst/examples/ss.Example.R
#' @export
simulate_ss_pva <- function(
  initial_pops,
  n_years = 100,
  n_runs = 1000,
  K = NA,
  quasi_extinction_thresholds = 0,
  ...
) {

  args = list(...)

  if (sum(
    is.null(args$log_lambdas),
    is.null(args$growth_rate_means)
  ) != 1) {
    stop("You must provide either log-lambdas or growth-rate means and variances, but only one of them")
  }

  if (is.null(args$growth_rate_means)) {
    growth_rate_means <- mean(args$log_lambdas)
    growth_rate_vars <- var(args$log_lambdas)
  } else {
    if (is.null(args$growth_rate_vars)) {
      stop("You must also provide growth rate variances (growth_rate_vars argument)")
    }
    growth_rate_means <- args$growth_rate_means
    growth_rate_vars <- args$growth_rate_vars
  }

  results = c()

  for (iteration in 1:n_runs) {
    population = initial_pops
    for (year in 1:n_years) {
      population <- population * exp(rnorm(n = 1, mean = growth_rate_means, sd = sqrt(growth_rate_vars)))

      population = floor(population)

      if (!(is.na(K))  && (population > K)) population <- min(population,K)

      # If population goes extinct, go to next iteration (it cannot get better...)
      if (population < quasi_extinction_thresholds) {
        break;
      }

    }
    results = append(results,population)
  }

  ss <- list(
    final_pops = results,
    n_years = n_years,
    n_runs = n_runs,
    initial_pops = initial_pops,
    decline_risk = sum(results < initial_pops) / n_runs,
    extinction_risk = sum(results <= quasi_extinction_thresholds) / n_runs
  )

  class(ss) <- "ssPVARes"
  return(ss)

}

#' @export
print.ssPVARes <- function(x, ...) {
  cat(paste('Over a', x$n_years, 'years span, the extinction risk of this population is', x$extinction_risk))
  cat(paste('\r\nAnd the risk of decline is', x$decline_risk))
}

#' @export
hist.ssPVARes <- function(x, ...) {
  hist(res$final_pops, xlab = "Population size", main = paste("Population size at the end of", res$n_runs, "simulations"))
}
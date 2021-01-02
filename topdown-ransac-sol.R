# a)
############## PSEUDO - CODE ##################
## do not execute ##

# apply RANSAM model to data for a certain model
ransaclm(formula, 
         data, 
         min_obs,
         error_threshold, 
         inlier_threshold,
         stop_conditions,
         iterations) {
  # check the inputs
  
  first_sample <- sample_from_data(data, min_obs) # no separate function
  first_model <- fit_model(formula, first_sample)
  
  ## start loop 
  inliers <- detect_inliers(model, 
                            all_observations, 
                            sample/inliers, 
                            error_threshold)
  inlier_proportion <- inliers/datapoints
  ## check if stop_conditions are fulfilled:
  # iterations are reached OR
  # inlier_threshold is met OR
  # BOTH is fulfilled
  ## if yes, END
  
  # otherwise:
  new_model <- fit_model(formula, inliers)
  ## end loop
}

# fit model
fit_model(formula, data) {
  # apply formula in linear model
  # what is important? 
  # - coefficients
  # - evtl. responses (instead of manual computation)
  # - evtl. residuals (instead of manual computation) 
}

# find the inliers of the data (with respect to a model)
detect_inliers(model, data, sample, error_threshold) {
  prediction <- predict_observations(model, sample)
  residuals <- compute_residuals(prediction, data, sample)
  inliers <- compare(residuals, error_threshold) # no separate function
}

predict_observations(model, data) {
  # Three possibilities: A, B, C
  
  # A: Either use coefficients to predict on data
  
  # OR 
  
  # B: Merge model responses with predictions of data not used to 
  # compute the model
  
  # OR
  
  # C: only predict on data not used in the model
}

compute_residuals(prediction, data, (model)) {
  # A, B: Either compute errors purely from prediction
  
  # OR 
  
  # C: Merge residuals from model with manually computed errors
  # from prediction of data points, which are not in the model
}


## Concerning Parallelisation:
# since the results of the iterations are based on each other (one uses the
# inliers from the previous iteration for modeling), it is not possible to 
# iterate parallel 
# since Hinweis 2 states that it is not of much use anyway, I won't consider it

# b)
############## IMPLEMENTATION ##################
# chose variant C for the implementation

# function to check the input of ransaclm to be correct
check_input_ransaclm <- function(formula,
                                 data,
                                 min_obs,
                                 error_threshold,
                                 inlier_threshold,
                                 stop_when_reached,
                                 iterations,
                                 seed) {
  # check the inputs
  checkmate::assert_class(formula, classes = "formula")

  # extract attribues of formula
  formula_attr <- attributes(terms.formula(formula, data = data))
  # look-up the number of parameters
  n_parameters <- length(formula_attr$term.labels)

  # data must at least contain one feature and one target variable
  # and there must be at least as many rows as n_parameters or 2 if it is lower
  valid_data <-
    checkmate::check_data_frame(data,
      min.rows = n_parameters + 1,
    )

  if (!is.logical(valid_data)) {
    stop("dimensions of 'data' do not suit together with 'formula'. Please check if p > n.")
  }

  valid_data <-
    checkmate::check_data_frame(data,
      min.cols = 2,
    )

  if (!is.logical(valid_data)) {
    stop("'data' needs at least two columns in order to compute a linear model")
  }

  # min_obs must be a integerish value
  # not greater than the number of observations
  # but at least n_parameters or 2 if it is lower
  valid_min_obs <-
    checkmate::check_integerish(min_obs,
      lower = n_parameters + 1,
      upper = nrow(data),
      any.missing = FALSE,
      len = 1
    )

  if (!is.logical(valid_min_obs)) {
    stop("'min_obs' must be at least 'n_parameters' + 1 and at most 'data's number of rows")
  }

  # error threshold must be a positive numerical value
  checkmate::assert_number(error_threshold,
    lower = 0
  )

  # inlier_threshold must be a positive integerish value
  # not greater than the number of observations
  checkmate::assert_integerish(inlier_threshold,
    lower = 0,
    upper = nrow(data),
    all.missing = FALSE,
    null.ok = FALSE,
    len = 1
  )

  # stop_when_reached must be  if (inlier_threshold > 0) {
  checkmate::assert_flag(stop_when_reached)


  # iterations must be a strictly positive integerish value
  # must be defined to avoid non-terminating models
  checkmate::assert_count(iterations,
    positive = TRUE
  )

  # seed must be a number (allow doubles like in the set.seed function)
  checkmate::assert_number(seed, null.ok = TRUE)

  # is needed later
  return(n_parameters)
}


# find the inliers of the data (with respect to a model)
detect_inliers <- function(model, target, data, sample_indizes, error_threshold) {
  # compute the residuals for observations used in model
  residuals_im <- residuals(model)
  # first check for them if they are within the error threshold
  inliers <- sample_indizes[abs(residuals_im) <= error_threshold]
  # then compute the residuals for observations not used in the model
  predictions_nim <- predict(model, data[-sample_indizes, ])
  residuals_nim <- predictions_nim - data[-sample_indizes, paste(target)]
  # also add the respective observations to inliers
  inliers <-
    c(inliers, seq_len(nrow(data))[-sample_indizes][abs(residuals_nim) <= error_threshold])
  inliers
}

ransaclm <- function(formula,
                     data,
                     min_obs = nrow(data),
                     error_threshold,
                     inlier_threshold = 0,
                     stop_when_reached = FALSE,
                     iterations = 100,
                     seed = NULL) {
  # do input checks
  n_parameters <-
    check_input_ransaclm(
      formula,
      data,
      min_obs,
      error_threshold,
      inlier_threshold,
      stop_when_reached,
      iterations,
      seed
    )

  # set seed
  set.seed(seed)

  # get target variable from formula
  target <- as.character(formula)[2]

  # generate first sample
  inliers <- sample(nrow(data), size = min_obs, replace = FALSE)

  # remember if inlier_threshold is fulfilled or not
  inlier_threshold_reached <- FALSE

  i <- 1
  while (i <= iterations) {
    # check if model can be computed with the current number of inliers
    enough_inliers <- checkmate::check_true(length(inliers) >= n_parameters + 1)
    if (!is.logical(enough_inliers)) {
      stop("inliers of current iteration are not sufficient for the specified linear model")
    }
    # compute a model with the chosen observations
    model <- lm(formula, data[inliers, ])
    # detect the inliers for the current model
    inliers_new <- detect_inliers(model, target, data, inliers, error_threshold)
    # if there is no change in the inliers stop looping
    if (identical(inliers_new, inliers)) {
      inliers <- inliers_new
      break
    }
    inliers <- inliers_new
    # check if the inlier_threshold is reached
    if (length(inliers) >= inlier_threshold) {
      inlier_threshold_reached <- TRUE
    }
    # if so, end the loop if stop_when_reached is set to TRUE
    if (stop_when_reached && inlier_threshold_reached) {
      break
    }
    # use the inliers for the next model
    i <- i + 1
  }

  output <- list()
  # if inlier_threshold was not reached return NULL instead of model
  if (inlier_threshold_reached) {
    output$model <- model
    # add the consensus set (chosen observations for the output model)
    data$.consensus_set <- FALSE
    data[inliers, ".consensus_set"] <- TRUE
    output$data <- data
  }
  else {
    output$model <- NULL
    output$data <- data
    warning("finding a model suiting the given criteria failed")
  }
  # return output
  output
}


###### Testing #######

# generate toy data with very clear inlier / outlier distinction to
# validate ransaclm. intercept is 0, coefs have value 1.
# inputs: as named -- number  of observations, number of coefs, fraction of inliers
# output: a data.frame with columns y, x.1, ..., x.<n_coef>, inlier
make_ransac_data <- function(n_obs, n_coef, inlier_fraction = 0.7) {
  coef <- rep(1, n_coef)
  design <- matrix(runif(n_obs * n_coef), n_obs, n_coef)
  inlier <- sample(
    c(TRUE, FALSE)[c(
      rep(1, n_obs * inlier_fraction),
      rep(2, n_obs - n_obs * inlier_fraction)
    )]
  )
  # E(y) = design %*% coef if inlier, else random draw from
  # uniform dist. on [-10 * n_coef, -2 * n_coef] and [2 * n_coef, 10 * n_coef].
  # error variance = n_coef * Var(U[0,1])
  inlier_expected_y <- design %*% coef
  outlier_expected_y <- sample(c(-1, 1), n_obs, replace = TRUE) *
    runif(n_obs, min = 2 * n_coef, max = 10 * n_coef)
  y <- ifelse(inlier, inlier_expected_y, outlier_expected_y) +
    rnorm(n_obs, sd = sqrt(1 / 12 * n_coef))
  data.frame(y, x = design, inlier)
}

#-------------------------------------------------------------------------------
# summarize & visualize ransaclm() results on data from "make_ransac_data". Only
# works if ransaclm's output list has a "data" entry with columns "inlier" as
# generated by make_ransac_data() as well as a column ".consensus_set" flagging
# the observations in the consensus set found by RANSAC.
validate_ransac <- function(ransacmodel, plot = TRUE) {
  checkmate::assert_class(ransacmodel[["model"]], "lm")
  checkmate::assert_data_frame(ransacmodel[["data"]])
  data <- ransacmodel[["data"]]
  checkmate::assert_logical(data[, "inlier"], any.missing = FALSE)
  checkmate::assert_logical(data[, ".consensus_set"], any.missing = FALSE)

  consensus_set <- data[, ".consensus_set"]
  true_inliers <- data[, "inlier"]
  cat("Inlier:\n")
  print(table(true = true_inliers, estimated = consensus_set))
  cat("\nCoefficients: (should be intercept ~0, effects ~1)\n")
  print(summary(ransacmodel[["model"]])$coefficients)
  if (plot &
    (!all(c("x", "y") %in% names(data)))) {
    warning("Can't plot this data, expecting columns 'x' and 'y'.")
    plot <- FALSE
  }
  if (plot) {
    plot_ransac(ransacmodel, data)
  }
  invisible(ransacmodel)
}

plot_ransac <- function(ransacmodel, data,
                        colors = c(
                          black = rgb(0, 0, 0, .5),
                          red   = rgb(1, 0, 0, .5),
                          blue  = rgb(0, 0, 1, .5)
                        )) {
  # scatterplot with true inliers in red, estimated inliers in blue
  # true regression line in black, RANSAC model in blue
  with(
    data,
    plot(x, y, col = colors[inlier + 1], pch = 19, bty = "n")
  )
  abline(c(0, 1), lwd = 2, col = colors[1])
  abline(ransacmodel$model, col = colors[3], lty = 2, lwd = 2)
  with(
    subset(data, .consensus_set),
    points(x, y, pch = 4, cex = 1.5, col = colors[3])
  )
  legend("top",
    lty = c(NA, NA, 1, 2, NA), lwd = c(NA, NA, 2, 2, NA),
    pch = c(19, 19, NA, NA, 4), col = colors[c(1, 2, 1, 3, 3)],
    legend = c(
      "'true' outliers", "'true' inliers", "true regression line",
      "RANSAC estimate", "RANSAC consensus set"
    ),
    cex = .7, bg = NA, inset = c(0, -.15), ncol = 3, xpd = NA
  )
}

# immer set.seed() um Ergebnisse reproduzierbar zu machen...
set.seed(1874374111)
data_simple <- make_ransac_data(
  n_obs = 100,
  n_coef = 1,
  inlier_fraction = 0.7
)

# univariate example:
ransac_simple <- ransaclm(y ~ . - inlier,
  data = data_simple,
  min_obs = nrow(data_simple),
  error_threshold = 2,
  inlier_threshold = 70,
  stop_when_reached = FALSE,
  iterations = 100,
  seed = 20171111
)
validate_ransac(ransac_simple)

### Own tests

testthat::test_that("ransaclm can solve an easy example with more than one feature", {
  set.seed(1234)
  ransac_data <- make_ransac_data(
    n_obs = 100,
    n_coef = 10,
    inlier_fraction = 0.7
  )

  ransac_object <- ransaclm(y ~ . - inlier,
    data = ransac_data,
    min_obs = nrow(ransac_data),
    error_threshold = 2,
    inlier_threshold = 20,
    stop_when_reached = FALSE,
    iterations = 100,
    seed = 20171111
  )
  testthat::expect_equal(class(ransac_object[[1]]), "lm")
  testthat::expect_vector(ransac_object[[2]][, ".consensus_set"],
    ptype = logical()
  )
})

testthat::test_that("ransaclm produces a warning if it fails", {
  set.seed(1234)
  ransac_data <- make_ransac_data(
    n_obs = 100,
    n_coef = 10,
    inlier_fraction = 0.8
  )

  testthat::expect_warning(ransaclm(y ~ . - inlier,
    data = ransac_data,
    min_obs = nrow(ransac_data),
    error_threshold = 5,
    inlier_threshold = 100,
    stop_when_reached = FALSE,
    iterations = 100,
    seed = 20171111
  ))
})

testthat::test_that("ransaclm produces an error if there are too few inliers for a new model at some point", {
  set.seed(1234)
  ransac_data <- make_ransac_data(
    n_obs = 100,
    n_coef = 10,
    inlier_fraction = 0.7
  )

  testthat::expect_error(ransaclm(y ~ . - inlier,
    data = ransac_data,
    min_obs = nrow(ransac_data),
    error_threshold = 0,
    inlier_threshold = 90,
    stop_when_reached = FALSE,
    iterations = 100,
    seed = 20171111
  ))
})

testthat::test_that("ransaclm produces an error if min_obs is not high enough", {
  set.seed(1234)
  ransac_data <- make_ransac_data(
    n_obs = 100,
    n_coef = 10,
    inlier_fraction = 0.7
  )

  testthat::expect_error(ransaclm(y ~ . - inlier,
    data = ransac_data,
    min_obs = 5,
    error_threshold = 10000000,
    inlier_threshold = 50,
    stop_when_reached = FALSE,
    iterations = 5000,
    seed = 20171111
  ))
})

testthat::test_that("ransaclm stops early if there is no change in the inliers", {
  set.seed(1234)
  ransac_data <- make_ransac_data(
    n_obs = 100,
    n_coef = 5,
    inlier_fraction = 0.7
  )

  short_loop <- ransaclm(y ~ . - inlier,
    data = ransac_data,
    min_obs = 50,
    error_threshold = 2,
    inlier_threshold = 50,
    stop_when_reached = FALSE,
    iterations = 10,
    seed = 20171111
  )
  long_loop <- ransaclm(y ~ . - inlier,
    data = ransac_data,
    min_obs = 50,
    error_threshold = 2,
    inlier_threshold = 50,
    stop_when_reached = FALSE,
    iterations = 500000,
    seed = 20171111
  )

  # compare results
  testthat::expect_equal(short_loop, long_loop)

  # compare running times
  compare_times <-
    microbenchmark::microbenchmark(
      ransaclm(y ~ . - inlier,
        data = ransac_data,
        min_obs = 50,
        error_threshold = 2,
        inlier_threshold = 50,
        stop_when_reached = FALSE,
        iterations = 10,
        seed = 20171111
      ),
      ransaclm(y ~ . - inlier,
        data = ransac_data,
        min_obs = 50,
        error_threshold = 2,
        inlier_threshold = 50,
        stop_when_reached = FALSE,
        iterations = 500000,
        seed = 20171111
      )
    )
  mean_times <- summary(compare_times)$mean
  # mean times should be approximately the same
  testthat::expect_true(mean_times[2] > 0.95 * mean_times[1] &&
    mean_times[2] < 1.05 * mean_times[1])
})


testthat::test_that("stop_when_reached saves running time for short inlier_thresholds and many iterations", {
  set.seed(1234)
  ransac_data <- make_ransac_data(
    n_obs = 500,
    n_coef = 5,
    inlier_fraction = 0.7
  )

  # compare running times
  compare_times <-
    microbenchmark::microbenchmark(
      ransaclm(y ~ . - inlier,
        data = ransac_data,
        min_obs = 50,
        error_threshold = 2,
        inlier_threshold = 5,
        stop_when_reached = FALSE,
        iterations = 10000,
        seed = 20171111
      ),
      ransaclm(y ~ . - inlier,
        data = ransac_data,
        min_obs = 50,
        error_threshold = 2,
        inlier_threshold = 5,
        stop_when_reached = TRUE,
        iterations = 10000,
        seed = 20171111
      )
    )
  mean_times <- summary(compare_times)$mean
  # function call with stop_when_reached = TRUE (index [2]) should be much faster
  print(mean_times[2] / mean_times[1])
  testthat::expect_true(mean_times[2] < mean_times[1])
})

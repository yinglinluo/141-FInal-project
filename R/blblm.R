#' @import purrr
#' @import stats
#' @import utils
#' @importFrom magrittr %>%
#' @aliases NULL
#' @details
#' Regression models with Little Bag of Bootstraps
"_PACKAGE"

## quiets concerns of R CMD check re: the .'s that appear in pipelines
# from https://github.com/jennybc/googlesheets/blob/master/R/googlesheets.R
utils::globalVariables(c("."))

#' find the linear regression with bag of little bootstraps
#' blblm is used to fit linear regression  models. It can be used to carry out regression,
#' single stratum analysis of variance and analysis of covariance
#' @param formula an object of class "formula"
#' @param  data an data frame, list or environment containing the variables in the model.
#' @param  m split data into m parts，default number is 10
#' @param  B times of bootstrap
#' @export

blblm <- function(formula, data, m = 10, B = 5000) {
  data_list <- split_data(data, m)
  estimates <- map(
    data_list,
    ~ lm_each_subsample(formula = formula, data = ., n = nrow(data), B = B))
  res <- list(estimates = estimates, formula = formula)
  class(res) <- "blblm"
  invisible(res)
}

#' find the linear regression with bag of little bootstraps useing future_map function
#' par_blblm is used to fit linear regression  models. It can be used to carry out regression,
#' single stratum analysis of variance and analysis of covariance
#' It can speed up the calculation
#' @param formula an object of class "formula" (or one that can be coerced to that class)
#' @param  data an data frame, list or environment containing the variables in the model.
#' @param  m split data into m parts，default number is 10
#' @param  B times of bootstrap
#'
#' @export

par_blblm <- function(formula, data, m = 10, B = 5000) {
  if(class(data) == "character"){
    data_list <- read_data(data)
  }
  else{
    data_list <- split_data(data, m)
  }
  estimates <- future_map(
    data_list,
    ~ lm_each_subsample(formula = formula, data = ., n = nrow(.), B = B))
  res <- list(estimates = estimates, formula = formula)
  class(res) <- "blblm"
  invisible(res)
}

#' find the logistic regression with bag of little bootstraps
#' blbglm is used to fit logistic models. It can be used to carry out regression,
#' single stratum analysis of variance and analysis of covariance
#' @param formula an object of class "formula" (or one that can be coerced to that class)
#' @param  data an data frame, list or environment containing the variables in the model.
#' @param  m split data into m parts，default number is 10
#' @param  B eachoots run B times,dedefault number is 5000
#' @export


blbglm <- function(formula, data, m = 10, B = 5000, family) {
  if(class(data) == "character"){
    data_list <- read_data(data)
  }
  else{
    data_list <- split_data(data, m)
  }
  estimates <- map(
    data_list,
    ~ glm_each_subsample(formula = formula, data = ., n = nrow(.), B = B, family))
  res <- list(estimates = estimates, formula = formula)
  class(res) <- "blbglm"
  invisible(res)
}

#' find the logistics regression with bag of little bootstraps useing future_map function
#' par_blblm is used to fit linear regression  models. It can be used to carry out regression,
#' single stratum analysis of variance and analysis of covariance
#' It can speed up the calculation
#' @param formula an object of class "formula" (or one that can be coerced to that class)
#' @param  data an data frame, list or environment containing the variables in the model.
#' @param  m split data into m parts，default number is 10
#' @param  B times of bootstrap
#'
#' @export

par_blbglm <- function(formula, data, m = 10, B = 5000, family) {
  if(class(data) == "character"){
    data_list <- read_data(data)
  }
  else{
    data_list <- split_data(data, m)
  }
  estimates <- future_map(
    data_list,
    ~ glm_each_subsample(formula = formula, data = ., n = nrow(.), B = B, family))
  res <- list(estimates = estimates, formula = formula)
  class(res) <- "blbglm"
  invisible(res)
}

#' split data into m parts of approximated equal sizes
#' @param data an data frame containing the variables in the function.
#' @param m split data into m parts，default number is 10

split_data <- function(data, m) {
  idx <- sample.int(m, nrow(data), replace = TRUE)
  data %>% split(idx)
}

#’ print the formula
#' @param x fittef model
#' @param ...extra conditions
#' @export
#' @method print blblm

print.blblm <- function(x,...) {
  cat(class(x),"model:",capture.output(x$formula))
  cat("\n")
}

#' compute the estimates
#' compute the estimates for each subsample and fit linear regression
#' model for each subsmple B times

#' @param formula an object of class "formula"
#' @param data the sub_data from original data set
#' @param n the size of original data set
#' @param B the subsample repeat B times
#' @export


lm_each_subsample <- function(formula, data, n, B) {
  replicate(B, lm_each_boot(formula, data, n), simplify = FALSE)
}


#' compute the regression estimates for a blb dataset.
#'@param formula an object of class "formula" (or one that can be coerced to that class)
#'@param data the sub_data from original data set
#'@param n the size of original data set
#'@export


lm_each_boot <- function(formula, data, n) {
  freqs <- rmultinom(1, n, rep(1, nrow(data)))
  lm1(formula, data, freqs)
}


#' estimate the regression estimates based on given the number of repetitions
#' @param formula an object of class "formula" (or one that can be coerced to that class)
#' @param data the sub_data from original data set
#' @param freqs the frequence from lm_each_boot function
#' @export



lm1 <- function(formula, data, freqs) {
  # drop the original closure of formula,
  # otherwise the formula will pick a wront variable from the global scope.
  environment(formula) <- environment()
  fit <- lm(formula, data, weights = freqs)
  list(coef = blbcoef(fit), sigma = blbsigma(fit))
}


#' compute the coefficients from fit
#' compute the coefficients from each model
#' @param fit the linear regression model for each subsample
#' @export blbcoef

blbcoef <- function(fit) {
  coef(fit)
}


#' compute sigma from fit
#'
#' @param fit the linear regression model for each subsample
#' @export blbsigma
blbsigma <- function(fit) {
  p <- fit$rank
  y <- model.extract(fit$model, "response")
  e <- fitted(fit) - y
  w <- fit$weights
  sqrt(sum(w * (e^2)) / (sum(w) - p))
}

#' compute the GLM estimates for a blb dataset
#' @param formula regression formula to fit
#' @param data data
#' @param n number of random vectors to draw
#' @param family family in glm

glm_each_boot <- function(formula, data, n, family) {
  freqs <- rmultinom(1, n, rep(1, nrow(data)))
  glm1(formula, data, freqs, family)
}

#' compute the estimates for GLM
#' @param formula regression formula to fit
#' @param data data
#' @param n data
#' @param B times bootstrap
#' @param family family in glm


glm_each_subsample <- function(formula, data, n, B, family) {
  replicate(B, glm_each_boot(formula, data, n, family), simplify = FALSE)
}




#' estimate the logistic regression estimates based on given number of repetitions
#' @param formula an object of class "formula"
#' @param data the sub_data from original data set
#' @param freqs the frequence from glm_each_boot function
#' @param family family in glm

glm1 <- function(formula, data, freqs, family) {
  # drop the original closure of formula,
  # otherwise the formula will pick wrong variables from a parent scope.
  environment(formula) <- environment()
  fit <- glm(formula, data, family = family, weights = freqs)
  list(coef = coef(fit), sigma = sigma(fit))
}

#' @export  print blblm

print.blblm <- function(x, ...) {
  cat("blblm model:", capture.output(x$formula))
  cat("\n")
}


#' sigma for blb regression model
#' @export sigma blblm

sigma.blblm <- function(object, confidence = FALSE, level = 0.95, ...) {
  est <- object$estimates
  sigma <- mean(map_dbl(est, ~ mean(map_dbl(., "sigma"))))
  if (confidence) {
    alpha <- 1 - 0.95
    limits <- est %>%
      map_mean(~ quantile(map_dbl(., "sigma"), c(alpha / 2, 1 - alpha / 2))) %>%
      set_names(NULL)
    return(c(sigma = sigma, lwr = limits[1], upr = limits[2]))
  } else {
    return(sigma)
  }
}


#' coefficients for blb regression model
#' @param object fitted model
#' @param ... extra conditions
#' @export coef.blblm

coef.blblm <- function(object, ...) {
  est <- object$estimates
  map_mean(est, ~ map_cbind(., "coef") %>% rowMeans())
}



#' @export  confint blblm
confint.blblm <- function(object, parm = NULL, level = 0.95, ...) {
  if (is.null(parm)) {
    parm <- attr(terms(fit$formula), "term.labels")
  }
  alpha <- 1 - level
  est <- object$estimates
  out <- map_rbind(parm, function(p) {
    map_mean(est, ~ map_dbl(., list("coef", p)) %>% quantile(c(alpha / 2, 1 - alpha / 2)))
  })
  if (is.vector(out)) {
    out <- as.matrix(t(out))
  }
  dimnames(out)[[1]] <- parm
  out
}


#' @export
#' @method predict.blblm

predict.blblm <- function(object, new_data, confidence = FALSE, level = 0.95, ...) {
  est <- object$estimates
  X <- model.matrix(reformulate(attr(terms(object$formula), "term.labels")), new_data)
  if (confidence) {
    map_mean(est, ~ map_cbind(., ~ X %*% .$coef) %>%
               apply(1, mean_lwr_upr, level = level) %>%
               t())
  } else {
    map_mean(est, ~ map_cbind(., ~ X %*% .$coef) %>% rowMeans())
  }
}

mean_lwr_upr <- function(x, level = 0.95) {
  alpha <- 1 - level
  c(fit = mean(x), quantile(x, c(alpha / 2, 1 - alpha / 2)) %>% set_names(c("lwr", "upr")))
}

map_mean <- function(.x, .f, ...) {
  (map(.x, .f, ...) %>% reduce(`+`)) / length(.x)
}

map_cbind <- function(.x, .f, ...) {
  map(.x, .f, ...) %>% reduce(cbind)
}

map_rbind <- function(.x, .f, ...) {
  map(.x, .f, ...) %>% reduce(rbind)
}


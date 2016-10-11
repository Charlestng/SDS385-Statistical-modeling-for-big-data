#### Statistical modeling for big data ####

### Exercise 5 ###

## Sparcity ##

# question B1 #

# define a vector theta

theta <- floor(11 * runif(1000))
cache <- rbinom(n = 1000, size = 1, prob = 0.01)
theta <- theta * cache
theta <- as.vector(theta)
theta_sample <- rnorm(1000, mean = theta)

crossprod(theta - theta_estimate,theta - theta_estimate)

theta_estimate <- soft_threshold(theta_sample,3)

plot(x = theta)
lines(x = theta_estimate, col = "red" )

# get the diabetes data
diabetesX <- read.csv(url("https://raw.githubusercontent.com/jgscott/SDS385/master/data/diabetesX.csv"))
diabetesY <- read.csv(url("https://raw.githubusercontent.com/jgscott/SDS385/master/data/diabetesY.csv"))

diabetesX <- diabetesX[1:441,]
diabetesX <- as.matrix(diabetesX)
diabetesY <- as.matrix(diabetesY)

test <- glmnet(x = diabetesX, y = diabetesY, family = "gaussian")

# extract coefficients at a single value of lambda
coefficients <- matrix(ncol = 66)
for (e in seq(from = 0, to = 10, by = 0.01)){
  coefficients <- rbind(coefficients, cbind(e,t(as.matrix(coef(test, s = e)))))
  colnames(coefficients)[1] <- c("lambda")
}

MSE <- matrix(ncol = 2)
for (i in 1:nrow(coefficients)){
  error <- 1 / length(diabetesY) * crossprod(diabetesY - 
                                               diabetesX %*% coefficients[i,3:66] - 
                                               coefficients[i,"(Intercept)"] * 
                                               matrix(1, nrow = length(diabetesY)),
                                             diabetesY - 
                                               diabetesX %*% coefficients[i,3:66] - 
                                               coefficients[i,"(Intercept)"] * 
                                               matrix(1, nrow = length(diabetesY)))
  MSE <- rbind(MSE,
               cbind(coefficients[i,"lambda"],error))
}

MOOSE <- function(target, data, split_ratio, lambda){
  n_obs <- nrow(data)
  n_features <- ncol(data)
  training_indexes <- sample(1:n_obs, floor(n_obs * split_ratio))
  training_set <- data[training_indexes,]
  training_target <- target[training_indexes,]
  test_set <- data[!training_indexes,]
  test_target <- target[!training_indexes,]
  model <- glmnet(x = training_set, y = training_target, family = "gaussian")
  coefficients <- t(as.matrix(coef(model, s = lambda)))
  test_n_obs <- length(test_target)
  error <- 1 / test_n_obs * crossprod(test_target - 
                                      test_set %*% coefficients[2:n_features] - 
                                      coefficients["(Intercept)"] * 
                                      matrix(1, nrow = test_n_obs),
                                      test_target - 
                                      test_set %*% coefficients[2:n_features] - 
                                      coefficients["(Intercept)"] * 
                                      matrix(1, nrow = test_n_obs))
  return(error)
}

Mallows_Cp <- function(target, data, lambda){
  n_obs <- nrow(data)
  n_features <- ncol(data)
  model <- glmnet(x = data, y = target, family = "gaussian")
  coefficients <- t(as.matrix(coef(model, s = lambda)))
  MSE <- 1 / n_obs * crossprod(target -
                                 data %*% coefficients[2:n_features] - 
                                 coefficients["(Intercept)"] * 
                                 matrix(1, nrow = n_obs),
                                 target - 
                                 data %*% coefficients[2:n_features] - 
                                 coefficients["(Intercept)"] * 
                                 matrix(1, nrow = n_obs))
  S_lambda <- sum(coefficients != 0)
  OLS_model <- lm(target ~ data)
  sigma_hat_2 <- 1 / (n_obs - n_features) * sum(residuals(OLS_model)^2)
  CP_stat <- MSE + 2 * S_lambda / n_obs * sigma_hat_2
  return(CP_stat)
}

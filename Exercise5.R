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
  training_target <- target[training_indexes]
  test_set <- data[-training_indexes,]
  test_target <- target[-training_indexes]
  model <- glmnet(x = training_set, y = training_target, family = "gaussian")
  coefficients <- as.matrix(coef(model, s = lambda))
  test_n_obs <- length(test_target)
  residuals_vector <- test_target - 
    test_set %*% coefficients[2:(n_features+1)] - 
    coefficients[1] * 
    matrix(1, nrow = test_n_obs)
  error <- 1 / test_n_obs * crossprod(residuals_vector,residuals_vector)
  return(error)
}

Mallows_Cp <- function(target, data, lambda){
  n_obs <- nrow(data)
  n_features <- ncol(data)
  model <- glmnet(x = data, y = target, family = "gaussian")
  coefficients <- as.matrix(coef(model, s = lambda))
  MSE <- 1 / n_obs * crossprod(target -
                                 data %*% coefficients[2:(n_features+1)] - 
                                 coefficients[1] * 
                                 matrix(1, nrow = n_obs),
                               target - 
                                 data %*% coefficients[2:(n_features+1)] - 
                                 coefficients[1] * 
                                 matrix(1, nrow = n_obs))
  S_lambda <- sum(coefficients != 0)
  OLS_model <- lm(target ~ data)
  sigma_hat_2 <- 1 / (n_obs - n_features) * sum(residuals(OLS_model)^2)
  CP_stat <- MSE + 2 * S_lambda / n_obs * sigma_hat_2
  return(CP_stat)
}

cross_validation <- function(target, data, lambda_min, lambda_max, step, repeats){
  range <- seq(from = lambda_min, to = lambda_max, by = step)
  MOOSE_matrix <- matrix(data = range, nrow = length(range))
  for (i in 1:repeats){
    temporary_matrix <- matrix(ncol = 1)
    for (lambda in range){
      temporary_matrix <- rbind(temporary_matrix, t(as.matrix(c(
                                           MOOSE(target = diabetesY, 
                                                 data = diabetesX, 
                                                 split_ratio = 0.7, 
                                                 lambda)))))
    }
    MOOSE_matrix <- cbind(MOOSE_matrix, temporary_matrix[-1,])
  }

  return(MOOSE_matrix)
}

test <- cross_validation(target = diabetesY, data = diabetesX,
                         lambda_min = 0, lambda_max = 2, step = 0.01, repeats = 20)

average_MOOSE <- matrix(ncol = 2, nrow = 201)
for(i in 1:201){
  average_MOOSE[i,1]<- test[i,1]
  average_MOOSE[i,2]<- mean(test[i,-1])
}

test <- MOOSE(target = diabetesY, data = diabetesX, split_ratio = 0.7, lambda = 1)
test <- Mallows_Cp(target = diabetesY, data = diabetesX, lambda = 1)

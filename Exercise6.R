# get the diabetes data
diabetesX <- read.csv(url("https://raw.githubusercontent.com/jgscott/SDS385/master/data/diabetesX.csv"))
diabetesY <- read.csv(url("https://raw.githubusercontent.com/jgscott/SDS385/master/data/diabetesY.csv"))

diabetesX <- diabetesX[1:441,]
diabetesX <- as.matrix(diabetesX)
diabetesY <- as.matrix(diabetesY)

grad <- function(data, target, beta){
  t(data) %*% data %*% beta - crossprod(data, target)
}

prox_gradient <- function(data, target, beta, max_iter, gamma, lambda){
  n_obs <- nrow(data)
  n_features <- ncol(data)
  
  # initialise the proxi
  proxi <- beta*0
  
  # bookeeping variable
  beta_hist <- t(beta)
  
  for (iter in 1:max_iter){
    proxi <- beta - gamma * grad(data, target, beta)
    for (feat in 1:n_features){
      beta[feat] <- sign(proxi[feat]) * max(abs(proxi[feat]) - lambda * gamma, 0)
    }
    beta_hist <- rbind (beta_hist, t(beta))
  }
  return(beta_hist)
}

beta <- diabetesX[1,]*0
test <- prox_gradient(data = diabetesX,
                      target = diabetesY,
                      beta,
                      max_iter = 10000, gamma = 0.01, lambda = 10)
result <- diabetesY - diabetesX %*% test[10000,]
crossprod(result,result)

accelerated_prox_gradient <- function(data, target, beta, max_iter, gamma, lambda, s0){
  n_obs <- nrow(data)
  n_features <- ncol(data)
  
  # initialise the proxi
  proxi <- beta*0
  
  # initiate the accelerator
  accel <- beta*0
  
  # bookeeping variable
  beta_hist <- t(beta)
  
  # initiate the momentum
  s <- c(s0)
  
  #loop over the number of iterations
  for (iter in 1:max_iter){
    #update the proxi of beta
    proxi <- accel - gamma * grad(data, target, beta)
    
    # update the momentum
    s <- c(s, (1 + sqrt(1 + 4 * s[iter]^2)) / 2)
    print((s[iter] - 1) / s[iter + 1])
    for (feat in 1:n_features){
      beta[feat] <- sign(proxi[feat]) * max(abs(proxi[feat]) - lambda * gamma, 0)
    }
    
    beta_hist <- rbind (beta_hist, t(beta))
    
    # update the accelerated beta
    accel <- beta + (s[iter] - 1) / s[iter + 1] * (beta - beta_hist[iter,])
    
  }
  return(beta_hist)
}

test <- accelerated_prox_gradient(data = diabetesX,
                      target = diabetesY,
                      beta,
                      max_iter = 1000, gamma = 0.01, lambda = 10, s0 = 100)
result <- diabetesY - diabetesX %*% test[1000,]
crossprod(result,result)

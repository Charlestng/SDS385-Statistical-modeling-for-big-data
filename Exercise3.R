#### Statistical modeling for big data ####

### Exercise 3 ###

## Better learning ##

# let's get the data for the exercise

data <- read.csv(url("https://raw.githubusercontent.com/jgscott/SDS385/master/data/wdbc.csv"))
data <- data.frame(data)

# let's define our y and X
y <- data$M
y <- data.frame(y)
# let's convert y to a logical format
y <- y=="M"
y <- data.frame(y)
y <- y*1
y <- data.matrix(y)

# let's include an intercept column in the data set
X <- data[,3:12]
data$intercept <- 1
X <- cbind(data$intercept,X)
X <- data.matrix(X)

# log likelihood function
# inputs:
#    y, the target vector
#    X, the explanatory variables
#    beta, the parameters to be estimated
#    m, the number or trials for binomial distribution
# outputs:
#    log_l, the log likelihood for the binomial distribution at point beta (scalar)
log_likelihood <- function(y, X, beta, m){
  log_l <- 0
  for (i in 1:length(y)){
    log_l <- log_l - lchoose(m[i],y[i]) - (y[i] - m[i])*X[i,]%*%beta +
      m[i]*log(1 + exp(-X[i,]%*%beta))
  }
  return(log_l)
}


# gradient function
# inputs:
#    y, the target vector
#    X, the explanatory variables
#    beta, the parameters to be estimated
#    m, the number or trials for binomial distribution
# outputs:
#    grad, the beta-gradient of the log likelihood at point beta (column matrix)

gradient <- function(y, X, beta, m){
  
  # create grad, a column vectors the size of beta
  grad <- matrix(data=beta*0, nrow=ncol(X), ncol=1)
  
  # loop over every observation to calculate the gradient
  for (i in 1:length(y)){
    grad <- grad - ((y[i] - m[i]) + m[i]*(1 - 1/(1 + exp(-X[i,]%*%beta))))*X[i,] 
  }
  return(grad)
}

# backtracking line search function
# inputs:
#    the target vector, y
#    the data set containing the dependent variables, X
#    the vector of parameters to be optimised, beta
#    the vector giving the number of trials for each yi, m
#    the list of indexes for the batch, indexes
# outputs:
#    the new value of beta after the line search in the direction of the batch

backtracking <- function(y, X, beta, m){
  
  # initialise the step size
  alpha <- 10
  
  # the constant for the sufficient descent Wolfe condition
  c <- runif(1, min = 0, max = 1)
  
  # the decreasing factor for the search of a good alpha
  rho <- runif(1, min = 0, max = 1)
  
  # the initial log_likelyhood
  log_l <- log_likelihood(y, X, beta, m)
  
  # the search direction
  direction <- gradient(y, X, beta, m) / 
    sqrt(sum(gradient(y, X, beta, m)^2))
  
  grad <-gradient(y, X, beta, m)
  
  beta_update <- beta - alpha * direction
  
  # let's define the two components of the Wolfe condition
  sufficient <- log_l - c * alpha * t(grad) %*% direction
  new_likelihood <- log_likelihood(y, X, beta_update, m)
  
  # loop until a value of alpha that satisfies the Wolfe condition is found
  while (new_likelihood > sufficient){
    # update the value of alpha
    alpha <- alpha * rho
    # update the value of beta
    beta_update <- beta - alpha * direction
    
    # update the two components of the condition
    sufficient <- log_l - c * alpha * t(grad) %*% direction
    new_likelihood <- log_likelihood(y, X, beta_update, m)
  }
  # return the new value of beta after the line search is done
  return(beta_update)
}

# create a batch gradient descent function that uses line search
# inputs:
#    y, the explained variable
#    X, the dependent variables
#    beta, the estimated parameters
#    m, the vector containing the number of trials for each yi
#    batch_size, the number of data points we use for each iteration
# outputs:
#    beta, the estimated parameters
#    a convergence plot
gradient_descent_line_search <- function(y, X, beta, m){
  
  # initialise the number of iterations
  iteration <- 0
  
  # initialise the convergence indicator
  convergence <- 100
  
  # list of the betas calculated at each iteration
  beta_list <- data.frame(t(beta))
  
  # loop until convergence or limit number of iterations reached
  while(iteration < 100){  
    
    # update beta with line search
    beta <- backtracking(y, X, beta, m)
    
    # update the list of betas
    beta_list <- rbind(beta_list, data.frame(t(beta)))
    
    # add one iteration to the counter
    iteration <- iteration + 1
  }
  return(beta_list)
}

test2 <- gradient_descent_line_search(y, X, beta, m)
log_likelihood(y, X, t(as.matrix(test2[100,])), m)

# Line search algorithm for Quasi Newton method

# bactracking for any direction where you can also adapt the log likelihood type
# inputs:
#    the target vector, y
#    the data set containing the dependent variables, X
#    the vector of parameters to be optimised, beta
#    the vector giving the number of trials for each yi, m
#    direction, the direction of the line search
#    C, a constant strictly between 0 and 1 (not too close to 1)
#    rho, a constant strictly between 0 and 1 (not too close to 1)
#    alpha, a constant greater than 0 (better have it too big than too small)

# outputs:
#    the new value of beta after the line search in the direction of the batch

backtracking_general <- function(y, X, beta, m, direction, C, rho, alpha){
  
  # the initial log_likelyhood
  log_l <- log_likelihood(y, X, beta, m)
  
  # the gradient of the log likelihood
  grad <- gradient(y, X, beta, m)
  
  # initialise beta
  beta_update <- beta - alpha * direction
  
  # let's define the two components of the Wolfe condition
  sufficient <- log_l - C * alpha * t(grad) %*% direction
  new_likelihood <- log_likelihood(y, X, beta_update, m)
  
  # loop until a value of alpha that satisfies the Wolfe condition is found
  while (new_likelihood > sufficient){
    
    # update the value of alpha
    alpha <- alpha * rho
    
    # update the value of beta
    beta_update <- beta - alpha * direction
    
    # update the two components of the condition
    sufficient <- log_l - C * alpha * t(grad) %*% direction
    new_likelihood <- log_likelihood(y, X, beta_update, m)
  }
  
  # return the new value of beta after the line search is done
  return(beta_update)
}

# Quasi Newton algorithm

Quasi_Newton_line_search <- function(y, X, beta, m, C, rho, alpha, max_iter){
  
  # initialise the number of iterations
  iteration <- 0
  
  # create the identity matrix
  Identity <- diag(x = 1, nrow = length(beta), ncol = length(beta))
  
  # Initialise the inverse of the Hessian estimate
  Hk <- Identity
  
  # Initialise convergence condition
  convergence <- FALSE
  
  # create beta list
  beta_list <- t(beta)
  
  # initialise the gradient
  grad_new <- gradient(y, X, beta, m)
  
  # loop until the maximum number of iterations is reached of convergence achieved
  while (iteration < max_iter & convergence==FALSE){
    
    # update iteration
    iteration <- iteration + 1
    
    # initial gradient of log likelihood
    grad <- grad_new
    
    # initial direction
    direction <- - Hk %*% grad_new

    norm_direction <- as.numeric (sqrt((crossprod(direction, direction))))

    direction <- direction / norm_direction
    
    # search for optimal beta for this search direction
    beta_new <- backtracking_general(y, X, beta, m, direction, C, rho, alpha)
    
    # create sk
    sk <- beta_new - beta
    
    #update beta
    beta <- beta_new
    
    # update the list of betas
    beta_list <- rbind(beta_list, t(beta))
    
    # create new gradient
    grad_new <- gradient(y, X, beta, m)
    
    # create yk
    yk <- grad_new - grad
    
    # create rhok
    rhok <- as.numeric(1 / crossprod(yk, sk))
    
    # create sk yk T
    skykT <- tcrossprod(yk, sk)
    
    # update the inverse of the Hessian
    Hk <- (Identity - rhok * skykT) %*% Hk %*% (Identity - rhok * skykT) + rhok * sk %*% t(sk)
    
    # update convergence status
  }
  return(beta_list)
}

test <- Quasi_Newton_line_search(y, X, beta, m, C=0.7, rho=0.5, alpha=10, max_iter=1000)

######################################
###                                ###
### Functions for SDS385 exercises ###
###                                ###
######################################

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
#    grad, the beta-gradient of the negative log likelihood at point beta (column matrix)

gradient <- function(y, X, beta, m){
  
  # create grad, a column vectors the size of beta
  grad <- matrix(data=beta*0, nrow=ncol(X), ncol=1)
  
  # loop over every observation to calculate the gradient
  for (i in 1:length(y)){
    grad <- grad - (y[i] - m[i]*(1/(1 + exp(-X[i,]%*%beta))))*X[i,] 
  }
  return(grad)
}

# gradient descent function
# inputs:
#    y, the target vector
#    X, the explanatory variables
#    beta, the parameters to be estimated
#    m, the number or trials for binomial distribution
#    step, the step taken in the direction of the negative gradient at each iteration
# outputs:
#    beta, the estimated parameters (column matrix)
#    iteration, the number of iterations to convergence
#    l_likelihood, the list of the log likelihoods

gradient_descent <- function(y, X, beta, m, step, max_iter){
  
  # initialise the list of betas
  beta_list <- t(beta)
  
  # initialise the number of iterations
  iteration <- 0
  
  # initialise the convergence metric
  convergence <- 100
  
  # loop until convergence criterion is met or maximum number of iterations is reached
  while (iteration < max_iter){
    
    # update the number of iterations
    iteration <- iteration + 1
    
    # update beta along the negative gradient direction
    beta <- beta - step*gradient(y, X, beta, m)
    
    # calculate the log likelihood for the new beta and add it to the list
    beta_list <- rbind(beta_list, t(beta))
  }  
  return(beta_list)
}

# let's write a function that returns the local gradient component for an index
# inputs:
#    y, the target vector
#    X, the explanatory variables
#    beta, the parameters to be estimated
#    m, the number or trials for binomial distribution
#    index, the observation for which we calculate the contribution to the gradient
# outputs:
#    grad, the contribution to the gradient from the observation "index"

gradient_local <- function(y, X, beta, m, index){
  grad <- - ((y[index] - m[index]) + m[index]*(1 - 1/(1 + exp(-X[index,]%*%beta))))*X[index,] 
  return(grad)
}

# let's write a function that return the gradient for a batch of indexes
# inputs:
#    y, the target vector
#    X, the explanatory variables
#    beta, the parameters to be estimated
#    m, the number or trials for binomial distribution
#    indexes, the list of observations for which we calculate the contribution to the gradient
# outputs:
#    grad, the contribution to the gradient from the observation "index"

gradient_batch <- function(y, X, beta, m, indexes){
  
  # initialise the value of the partial gradient to 0
  grad <- 0
  
  # loop over the observations in the list
  for (index in indexes){
    grad <- grad + gradient_local(y, X, beta, m, index)
  }
  return(grad)
}

# The stochastic gradient descent function has the following characteristics
#   The number of iteration for the algorithm is limited to 10000
#   The convergence is measured by the a criterion of growth rate for the exponential
# moving average of the local contributions to the log likelihood
#   At each iteration in the while loop a new index between 1 and n is picked at random in 
# the sample, beta is upgraded according to that local component of the gradient, the number
# of iterations is upgraded, the list of the log likelihood is updated, and we calculate the
# absolute value of the log likelihood's growth rate until convergence.
#   The function returns:
#      The final betas
#      The list of log_likelihoods
#      The number of iterations needed for convergence
stochastic_gradient_descent <- function(y, X, beta, m, step, max_iter){
  
  # set the number of observations that can contribute
  sample_size <- length(y)
  
  # initialise the number of iterations
  iteration <- 0
  
  # initialise the list of betas
  beta_list <- t(beta)
  
  while(iteration < max_iter){ 
    
    # randomly pick an observation in the data set
    index <- sample(1:sample_size,1)
    
    # update beta according to the negative gradient direction
    beta <- beta - step * (iteration + 1)^(-5/9) * gradient_local(y, X, beta, m, index)
    
    # update the list of betas
    beta_list <- rbind(beta_list, (t(beta)))
    
    # update the number of iterations
    iteration <- iteration + 1
  }
  return(beta_list)
}

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
  
  # initialise condition
  condition <- TRUE
  
  # loop until a value of alpha that satisfies the Wolfe condition is found
  while (condition){
    
    # update the value of alpha
    alpha <- alpha * rho
    
    # update the value of beta
    beta_update <- beta - alpha * direction
    
    # update the two components of the condition
    sufficient <- log_l - C * alpha * t(grad) %*% direction
    new_likelihood <- log_likelihood(y, X, beta_update, m)
    
    condition <- (new_likelihood > sufficient)
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
gradient_descent_line_search <- function(y, X, beta, m, C, rho, alpha, max_iter){
  
  # initialise the number of iterations
  iteration <- 0
  
  # list of the betas calculated at each iteration
  beta_list <- t(beta)
  
  # loop until convergence or limit number of iterations reached
  while(iteration < max_iter){  
    
    # update gradient of log likelihood to current version
    grad <- gradient(y, X, beta, m)
    
    # define search direction
    direction <- grad
    
    # update beta with line search
    beta <- backtracking_general(y, X, beta, m, direction, C, rho, alpha)
    
    # update the list of betas
    beta_list <- rbind(beta_list, t(beta))
    
    # add one iteration to the counter
    iteration <- iteration + 1
  }
  return(beta_list)
}

# Quasi Newton algorithm line search
# inputs:
#    y, the target vector
#    X, the explanatory variables
#    beta, parameters to be estimated
#    m, number of trials for each data point (binomial probs)
#    C, strictly between 0 and 1 Wolfe condition
#    rho, decreasing step factor strictly between 0 and 1
#    alpha, initial step length greater than 0
#    max_iter, maximal number of iterations
# outputs:
#    beta_list, list of estimated betas

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
    
    # gradient of log likelihood
    grad <- grad_new
    
    # initial direction
    direction <- Hk %*% grad_new
    
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
    if(crossprod(yk, sk) == 0){break}
    rhok <- as.numeric(1 / crossprod(yk, sk))
    
    # create sk yk T
    skykT <- tcrossprod(yk, sk)
    
    # update the inverse of the Hessian
    Hk <- (Identity - rhok * skykT) %*% Hk %*% (Identity - rhok * skykT) + rhok * sk %*% t(sk)
    
    # update convergence status
  }
  return(beta_list)
}

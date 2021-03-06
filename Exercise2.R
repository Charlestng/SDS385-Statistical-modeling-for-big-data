#### Statistical modeling for big data ####

### Exercise 2 ###

## Stochastic gradient descent ##

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

X <- data[,3:12]
data$intercept <- 1
X <- cbind(data$intercept,X)
X <- data.matrix(X)

# pick an initial value for beta, scale matrix X and define m
beta <- X[1,]*0+0.5
beta <- matrix(beta,nrow=length(beta))
m <- y*0+1
for (i in 2:ncol(X)){
  X[,i] <- scale(X[,i])
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

# let's write a function that returns the log likelihood
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

# let's write a function that returns the individual contribution to the log likelihood
# for an index
# local log likelihood function
# inputs:
#    y, the target vector
#    X, the explanatory variables
#    beta, the parameters to be estimated
#    m, the number or trials for binomial distribution
#    index, the observation fo
# outputs:
#    log_l, the log likelihood for the binomial distribution at point beta (scalar)
local_l_likelihood <- function(y, X, beta, m, index){
  L <- -lchoose(m[index],y[index]) - 
    (y[index] - m[index])*X[index,]%*%beta + 
    m[index]*log(1 + exp(-X[index,]%*%beta))
  return(L)
}

# let's write a function that returns the batch's contribution to the log likelihood
batch_l_likelihood <- function(y, X, beta, m, indexes){
  L <- 0
  for (index in indexes){
    L <- L + local_l_likelihood(y, X, beta, m, index)
  }
  return(L)
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
stochastic_gradient_descent <- function(y, X, beta, m, step, converged){
  
  # set the number of observations that can contribute
  sample_size <- length(y)
  
  # initialise the number of iterations
  iteration <- 0
  
  # initialise the convergence metric
  convergence <- 100
  
  # initialise the list of betas
  beta_list <- data.frame(t(beta))
  
  while(iteration < 10000){ 
    
    # randomly pick an observation in the data set
    index <- sample(1:sample_size,1)
    
    # update beta according to the negative gradient direction
    beta <- beta - step * (iteration + 1)^(-5/9) * gradient_local(y, X, beta, m, index)
    
    # update the list of betas
    beta_list <- rbind(beta_list, data.frame(t(beta)))
    
    # update the number of iterations
    iteration <- iteration + 1
  }
  return(c("beta",
         beta))
}

fin <- stochastic_gradient_descent(y, X, beta, m, step=1, converged=0.000001, decay=0.1)

# let's now write a batch version of the stochastic gradient descent
batch_gradient_descent <- function(y, X, beta, m, step, converged, decay, batch_size){
  sample_size <- length(y)
  iteration <- 0
  convergence <- 100
  beta_list <- c(beta)
  while(iteration < 10000){ 
    indexes <- sample(1:sample_size, batch_size)
    beta <- beta - step * (iteration + 1)^(-5/9) *  gradient_batch(y, X, beta, m, indexes)
    beta_list <- c(beta_list, beta)
    iteration <- iteration + 1
  }
  return(c("beta",
           beta))
}

test <- batch_gradient_descent(y, X, beta, m, step=0.1, converged=0.1, decay=0.1, batch_size=20)

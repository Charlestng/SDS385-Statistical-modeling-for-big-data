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

# let's include an intercept column in the data set
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

# let's write a function that returns the log likelihood
log_likelihood <- function(y, X, beta, m){
  log_l <- 0
  for (i in 1:length(y)){
    log_l <- log_l - lchoose(m[i],y[i]) - (y[i] - m[i])*X[i,]%*%beta +
      m[i]*log(1 + exp(-X[i,]%*%beta))
  }
  return(log_l)
}

# let's write a function that returns the local gradient component for an index
gradient_local <- function(y, X, beta, m, index){
  grad <- - ((y[index] - m[index]) + m[index]*(1 - 1/(1 + exp(-X[index,]%*%beta))))*X[index,] 
  return(grad)
}

# let's write a function that return the gradient for a batch of indexes
gradient_batch <- function(y, X, beta, m, indexes){
  grad <- 0
  for (index in indexes){
    grad <- grad + gradient_local(y, X, beta, m, index)
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

bactracking <- function(y, X, beta, m, indexes){
  
  # initialise the step size
  alpha <- 1
  
  # the constant for the sufficient descent Wolfe condition
  c <- runif(1, min = 0, max = 1)
  
  # the decreasing factor for the search of a good alpha
  rho <- runif(1, min = 0, max = 1)
  
  # the initial log_likelyhood
  log_l <- log_likelihood(y, X, beta, m)
  
  # the search direction
  direction <- gradient_batch(y, X, beta, m, indexes) / 
               abs(gradient_batch(y, X, beta, m, indexes))
  
  # loop until a value of alpha that satisfies the Wolfe condition is found
  while (
    log_l + c * aplha * gradient_batch(y, X, beta, m, indexes) <
    log_likelihood(y, X, beta + alpha * direction)
    )
    {
    # update the value of alpha
    alpha <- alpha * rho
    
    # update the value of beta
    beta <- beta + alpha * direction
  }
  # return the new value of beta after the line search is done
  return(beta)
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
batch_gradient_descent_line_search <- function(y, X, beta, m, batch_size){
  # number of data points
  sample_size <- length(y)
  
  # initialise the number of iterations
  iteration <- 0
  
  # initialise the convergence indicator
  convergence <- 100
  
  # list of the betas calculated at each iteration
  beta_list <- data.frame(beta)
  
  # loop until convergence or limit number of iterations reached
  while(iteration < 10000){ 
    
    # pick a random batch of data points
    indexes <- sample(1:sample_size, batch_size)
    
    # update beta with line search
    beta <- bactracking(y, X, beta, m, indexes)
    
    # update the list of betas
    beta_list <- cbind(beta_list, data.frame(beta))
    
    # add one iteration to the counter
    iteration <- iteration + 1
  }
  return(beta_list)
}

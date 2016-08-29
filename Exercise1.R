#### Statistical modeling for big data ####

### Exercise 1 ###

## Linear regression ##

# C #

# inversion method

# let's create a random matrix that contains our explanatory variables
# we are going to call it X, like in the exercise 1's regression problem
# y = Xb + e

inversion_method <- function(number_of_obs, number_of_variables){
  X <- rnorm(number_of_obs*number_of_variables, mean=0, sd=1)
  X <- matrix(X, nrow=number_of_obs, ncol=number_of_variables)

  # let's now generate the vector y to be explained
  y <- rnorm(number_of_obs, mean=0, sd=1)
  y <- matrix(y, nrow=number_of_obs, ncol=1)

  # we have everything we need in order to solve brutally the least square
  # problem we are going to call beta b here
  Inverse <- solve(t(X)%*%X)

  b <- Inverse %*% t(X) %*% y
  return(b)
}

ptm <- proc.time()
inversion_method(5000,1000)
proc.time() - ptm
# results
# user    system  elapsed 
# 18.28   0.12    18.49 

ptm <- proc.time()
inversion_method(10000,2000)
proc.time() - ptm
# results
# user      system  elapsed 
# 168.77    0.68    173.58 

# QR method

QR_method <- function(number_of_obs, number_of_variables){
  # let's create X and y like previously
  X <- rnorm(number_of_obs*number_of_variables, mean=0, sd=1)
  X <- matrix(X, nrow=number_of_obs, ncol=number_of_variables)
  
  y <- rnorm(number_of_obs, mean=0, sd=1)
  y <- matrix(y, nrow=number_of_obs, ncol=1)
  
  # let's now solve the least square problem using the QR decomposition
  b <- qr.solve(X,y)
  
  return(b)
}

ptm <- proc.time()
QR_method(5000,1000)
proc.time() - ptm
# results
# user    system  elapsed 
# 8.42    0.05    8.49

ptm <- proc.time()
QR_method(10000,2000)
proc.time() - ptm
# results
#user   system  elapsed 
#75.22  0.45    78.12

# D #

# let's install the "Matrix" package, that provides R prefered algebra routines
# for sparse matrices

library("Matrix", lib.loc="C:/Program Files/R/R-3.0.2/library")

N = 2000
P = 500
X = matrix(rnorm(N*P), nrow=N)
mask = matrix(rbinom(N*P,1,0.05), nrow=N)
X = mask*X

# inversion method

inversion_sparse_method <- function(number_of_obs, number_of_variables){
  X <- rnorm(number_of_obs*number_of_variables, mean=0, sd=1)
  X <- matrix(X, nrow=number_of_obs, ncol=number_of_variables)
  mask = matrix(rbinom(number_of_obs*number_of_variables,1,0.05), nrow=number_of_obs)
  X = mask*X
  
  # let's convert our matrix X to the sparse class to take advantage of its sparsity
  X <- Matrix(X, sparse=TRUE)
  
  # let's now generate the vector y to be explained
  y <- rnorm(number_of_obs, mean=0, sd=1)
  y <- matrix(y, nrow=number_of_obs, ncol=1)
  
  # we have everything we need in order to solve brutally the least square
  # problem we are going to call beta b here
  Inverse <- solve(t(X)%*%X, sparse=TRUE)
  
  b <- Inverse %*% t(X) %*% y
  return(b)
}

ptm <- proc.time()
inversion_sparse_method(5000,1000)
proc.time() - ptm
# results
# user  system elapsed 
# 8.82    0.11    9.00

ptm <- proc.time()
inversion_sparse_method(10000,2000)
proc.time() - ptm
# results
# user  system elapsed 
# 78.13    1.31   85.76

# QR method

QR_sparse_method <- function(number_of_obs, number_of_variables){
  # let's create X and y like previously
  X <- rnorm(number_of_obs*number_of_variables, mean=0, sd=1)
  X <- matrix(X, nrow=number_of_obs, ncol=number_of_variables)
  mask = matrix(rbinom(number_of_obs*number_of_variables,1,0.05), nrow=number_of_obs)
  X = mask*X
  
  # let's convert our matrix X to the sparse class to take advantage of its sparsity
  X <- Matrix(X, sparse=TRUE)
  
  y <- rnorm(number_of_obs, mean=0, sd=1)
  y <- matrix(y, nrow=number_of_obs, ncol=1)
  
  # let's now solve the least square problem using the QR decomposition
  b <- qr.solve(X,y)
  
  return(b)
}

ptm <- proc.time()
QR_sparse_method(5000,1000)
proc.time() - ptm
# results
#    user  system elapsed 
# 11.72    0.53   12.95 

ptm <- proc.time()
QR_sparse_method(10000,2000)
proc.time() - ptm
# results
# user  system elapsed 
# 87.89    0.90   92.39

## Generalised Linear Models ##

# B #

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

gradient <- function(y, X, beta, m){
  grad <- matrix(data=beta*0, nrow=ncol(X), ncol=1)
  for (i in 1:length(y)){
    grad <- grad - ((y[i] - m[i]) + m[i]*(1 - 1/(1 + exp(-X[i,]%*%beta))))*X[i,] 
  }
  return(grad)
}

log_likelyhood <- function(y, X, beta, m){
  log_l <- 0
  for (i in 1:length(y)){
    log_l <- log_l - lchoose(m[i],y[i]) - (y[i] - m[i])*X[i,]%*%beta +
      m[i]*log(1 + exp(-X[i,]%*%beta))
  }
  return(log_l)
}

gradient_descent <- function(y, X, beta, m, step){
  l_likelyhood <- c(log_likelyhood(y, X, beta, m))
  step_n <- 0
  convergence <- 100
  while (step_n < 100 & convergence > 0.001){
    step_n <- step_n + 1
    beta <- beta - step*gradient(y, X, beta, m)
    l_likelyhood <- c(l_likelyhood,log_likelyhood(y, X, beta, m))
    convergence <- abs((l_likelyhood[step_n+1]-l_likelyhood[step_n])/l_likelyhood[step_n])
  }
  return(c("beta",beta,"step",step_n,"log likelyhood",l_likelyhood))
}

# let's test our functions now on the given data
beta <- X[1,]*0+0.5
beta <- matrix(beta,nrow=length(beta))
m <- y*0+1
for (i in 2:ncol(X)){
  X[,i] <- scale(X[,i])
}
test <- gradient_descent(y, X, beta, m, step=0.01)

# D #

# let's create a function that computes the Hessian matrix for any beta

Hessian <- function(X, beta, m){
  H <- matrix(0,nrow=length(beta),ncol=length(beta))
  for (row in 1:length(beta)){
    for (col in 1:length(beta)){
      for (obs in 1:nrow(X)){
        H[row,col] <- m[obs]*X[obs,row]*X[obs,col]*exp(-X[i,]%*%beta)/
          (1 + exp(-X[i,]%*%beta))^2
      }
    }
  }
  return(H)
}

Newton_method <- function(y, X, beta, m){
  l_likelyhood <- c(log_likelyhood(y, X, beta, m))
  step_n <- 0
  convergence <- 100
  while (step_n < 100 & convergence > 0.001){
    step_n <- step_n + 1
    beta <- beta - solve(Hessian(X, beta, m))%*%gradient(y, X, beta, m)
    l_likelyhood <- c(l_likelyhood,log_likelyhood(y, X, beta, m))
    convergence <- abs((l_likelyhood[step_n+1]-l_likelyhood[step_n])/l_likelyhood[step_n])
  }
  return(c("beta",beta,"step",step_n,"log likelyhood",l_likelyhood))
}

#let's remove some columns in X that are too correlated
X <- X[,c(1,3,5,6,9,10,11)]
beta <- matrix(as.numeric(test[c(2,4,6,7,10,11,12)],nrow=length(beta)))
m <- y*0+1
test <- Newton_method(y, X, beta, m)



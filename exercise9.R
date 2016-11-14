#### Statistical modeling for big data ####

### Exercise 9 ###

## Matrix Factorization ##

# function for penalized rank K=1 single value decomposition

# penalized single value decomposition for K = 1
# input:
#    data, the matrix to be decomposed
#    lambda1, penalizer for u (must be between 1 and sqrt(n_obs))
#    lambda2, penalizer for v (must be between 1 and sqrt(n_feat))
#    max_iter, maximum number of iterations for the algorithm
#    convergence, convergence criteria
K1_PSVD <- function(data, lambda1, lambda2, max_iter, convergence = 1e-3){
  # number of observations
  n_obs <- nrow(data)
  
  # number of features 
  n_feat <- ncol(data)
  
  # test if lambda1 and lambda2 are appropriate
  if(lambda1<1 | lambda1>sqrt(n_obs)){
    stop("lambda1 must be between 1 and sqrt(nrow(data))")
  }
  if(lambda1<1 | lambda1>sqrt(n_feat)){
    stop("lambda1 must be between 1 and sqrt(ncol(data))")
  }
  
  # Initialise v
  v <- matrix(data = 1/sqrt(n_feat), nrow = n_feat)
  
  # Initialise objective function
  objective <- c(0)
  
  # Loop over the number of iterations
  for (iter in 1:max_iter){
    
    # Pre cache Xv
    Xv <- data %*% v
    
    # Compute u
    u <- Xv / sqrt(crossprod(Xv, Xv))
    
    # Pre cache L1 of u
    uL1 <- sum(abs(u))
    
    # Shrink to L1 = lambda1
    if ( uL1 > lambda1){
      
      # the penalized L1 norm behaves like sqrt
      delta1 <- max(Xv) * (1 - (lambda1/uL1)^2)
      SXv <- soft_threshold(Xv, delta1)
      u <- SXv / sqrt(crossprod(SXv, SXv))
    }
    
    # Same process for v
    Xu <- t(data) %*% u
    v <- Xu / sqrt(crossprod(Xu, Xu))
    vL1 <- sum(abs(v))
    if ( vL1 > lambda2){
      delta2 <- max(Xu) * (1 - (lambda2/vL1)^2)
      SXu <- soft_threshold(Xu, delta1)
      v <- SXu / sqrt(crossprod(SXu, SXu))
    }
    
    # compute d
    d <- t(u) %*% data %*% v
    
    objective <- c(objective, d)
    change <- abs((objective[iter+1] - objective[iter])/objective[iter])
    if(change < convergence){
      return(c(u,v,d))
      break
    }
  }
  return(c(u,v,d))
}


# the function S(X, delta)L1 / S(X, delta)L2 behaves sort of like a mirrored offset sqrt
vec <- matrix(1:100,nrow = 100)
C <- c()
for (i in seq(from = 0, to = 100, by = 0.1)){
  u <- soft_threshold(vec, i)
  v <- sum(abs(u))
  w <- sqrt(crossprod(u,u))
  x <- v/w
  C <- c(C,x)
}
plot(C)
comparaison <- sqrt(100 - seq(from = 0, to = 100, by = 0.1)) * C[1]/10
lines(comparaison)


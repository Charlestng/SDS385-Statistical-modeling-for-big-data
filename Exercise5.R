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

library(Matrix)
library(Rcpp)
library(RcppEigen)
sourceCpp('C:/Users/Charles/Documents/McCombs/Statistical modeling for big data/C++/MyAdagradLogit.cpp')

processSVM("C:/Users/Charles/Documents/McCombs/Statistical modeling for big data/exercise 4/url_svmlight1/")
processSVM("C:/Users/Charles/Documents/McCombs/Statistical modeling for big data/exercise 4/url_svmlight2/")

# Read in serialized objects
system.time(X <- readRDS('url_X.rds'))
system.time(y <- readRDS('url_y.rds'))

n = length(y)
p = ncol(X)

# column-oriented storage of each observation
tX = t(X)


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
X <- t(X)

# pick an initial value for beta, scale matrix X and define m
beta <- X[1,]*0+0.5
beta <- matrix(beta,nrow=length(beta))
m <- y*0+1
for (i in 2:ncol(X)){
  X[,i] <- scale(X[,i])
}

init_beta = rep(0.0, p)
system.time(sgd1 <- sparsesgd_logit(X, y, m, npass=1, beta0 = beta, lambda=.1))


system.time(sgd1 <- sparsesgd_logit(tX[1:10, 1:1000], y[1:1000], rep(1,n), npass=1, beta0 = init_beta[1:10], lambda=.1))



system.time(sgd1 <- sparsesgd_logit(tX, y, rep(1,n), npass=1, beta0 = init_beta, lambda=1e-8, discount = 0.001))
names(sgd1)

rafalib::splot(1:p, sort(sgd1$beta))
rafalib::splot(1:length(sgd1$nll_tracker), sgd1$nll_tracker)

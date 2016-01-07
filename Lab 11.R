library(truncnorm)
library(MASS)

data <- read.csv('Lab 11 - data.csv')
y <- data$y
n <- length(y)
x <- as.matrix(data)
x = x[,1:3]
x = cbind(rep(1, n), x)
x = t(x)

#Priors
b <- matrix(c(0, 0, 0, 0), 4)
B <- diag(4)
#Gibbs
y.star <- list(rep(1, n))
beta <- list(c(1, 1, 1, 1))
sum.xx <- Reduce('+', lapply(seq(1, n), function(i) {
  return(x[,i] %*% t(x[,i]))
}))
for(k in 2:10010) {
  y.star[[k]] = sapply(seq(1, n), function(i) {
    return(ifelse(y[i] == 1, 
                  rtruncnorm(1, a = 0, b = Inf, mean = t(x[,i]) %*% matrix(beta[[k-1]], 4), sd = 1), 
                  rtruncnorm(1, a = -Inf, b = 0, mean = t(x[,i]) %*% matrix(beta[[k-1]], 4), sd = 1)))
  })
  sum.xystar = Reduce('+', lapply(seq(1, n), function(i) {
    return(x[,i]*y.star[[k]][i])
  }))
  beta[[k]] = mvrnorm(1, solve(sum.xx + solve(B)) %*% (sum.xystar + (solve(B) %*% b)), solve(sum.xx + solve(B)))
}
#Trace plots
beta.0 <- sapply(beta, function(x) {
  return(x[1])
})
beta.1 <- sapply(beta, function(x) {
  return(x[2])
})
beta.2 <- sapply(beta, function(x) {
  return(x[3])
})
beta.3 <- sapply(beta, function(x) {
  return(x[4])
})
plot(seq(1, 10010), beta.0, type = 'l', main = 'Trace plot for beta.0', xlab = 'iteration')
plot(seq(1, 10010), beta.1, type = 'l', main = 'Trace plot for beta.1', xlab = 'iteration')
plot(seq(1, 10010), beta.2, type = 'l', main = 'Trace plot for beta.2', xlab = 'iteration')
plot(seq(1, 10010), beta.3, type = 'l', main = 'Trace plot for beta.3', xlab = 'iteration')
y.star.i.trace <- function(i) {
  y.star.i = sapply(y.star, function(x) {
    return(x[i])
  })
  plot(seq(1, 10010), y.star.i, type = 'l', main = paste0('Trace plot for y.star.', i), xlab = 'iteration', ylab = paste0('y.star.', i))
}

#Posterior predictive for TA
predictive <- sapply(beta, function(x) {
  #return(matrix(c(1, 26, 0, 1), 1) %*% matrix(x, 4))
  return(1 - pnorm(0, matrix(c(1, 26, 0, 1), 1) %*% matrix(x, 4), 1))
})
hist(predictive[11:10010], breaks = 20, main = 'Posterior predictive probability of accident for 26-year-old TA', xlab = 'Probability')
mean(predictive[11:10010])

#Mercedes question
predictive1 <- sapply(beta, function(x) {
  return(1 - pnorm(0, matrix(c(1, 17, 1, 0), 1) %*% matrix(x, 4), 1))
})
predictive2 <- sapply(beta, function(x) {
  return(1 - pnorm(0, matrix(c(1, 18, 1, 1), 1) %*% matrix(x, 4), 1))
})
mean(predictive1[11:10010])
mean(predictive2[11:10010])

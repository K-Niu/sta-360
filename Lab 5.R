library(data.table)
library(ggplot2)

beta = 0.05
N = 50
lambda = 25
y = 20
for(i in 2:100010) {
  N[i] = rpois(1, lambda*(1 - beta[i-1])) + y
  beta[i] = rbeta(1, y + 1, N[i] - y + 1)
}
data <- data.table(N = N, beta = beta)
ggplot(data[1:10,], aes(x = beta, y = N)) + geom_point() + geom_path() + labs(title = "2D trace plot for first 10 draws of Gibbs sampler")

#Throw out burn-ins
N = N[11:100010]
beta = beta[11:100010]

#90% credible interval for beta
beta = sort(beta)
c(beta[5000], beta[95000])

#P(N = 20 | y)
mean(N == 20)

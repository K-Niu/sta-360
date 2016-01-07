library(data.table)

#Commented code is for more hypotheses (narrower bins)

#Function that simulates fishing in a pond for a collective sum_ti hours
#theta is the true rate of bites, a and b are the parameters for a Gamma prior for theta
#Prior is truncated for each hypothesis
#Returns posterior probabilites under each hypothesis
simulate <- function(theta, a, b, sum_ti) {
  data = rpois(sum_ti, theta)
  sum_yi = sum(data)
  #cdf of Gamma(a, b)
  F1 <- function(c) {
    pgamma(c, a, rate = b)
  }
  #cdf of Gamma(a + sum_yi, b + sum_ti)
  F2 <- function(c) {
    pgamma(c, a + sum_yi, rate = b + sum_ti)
  }
  #Calculates marginal likehood of a hypothesis given the upper and lower bounds
  #If upper bound is infinity, omit upper bound parameter
  #Outputs marginal likehood without constant factor
  ML <- function(lower, upper) {
    if(missing(upper)) {
      return((1 - F2(lower))/(1 - F1(lower)))
    } else {
      return((F2(upper) - F2(lower))/(F1(upper) - F1(lower)))
    }
  }
  post_prob_H0 = 1/(1 + (ML(1, 5) + ML(5))/ML(0, 1))
  post_prob_H1 = 1/(1 + (ML(0, 1) + ML(5))/ML(1, 5))
  post_prob_H2 = 1/(1 + (ML(0, 1) + ML(1, 5))/ML(5)) 
  #post_prob_H0 = 1/(1 + (ML(1, 2) + ML(2, 3) + ML(3, 4) + ML(4, 5) + ML(5))/ML(0, 1))
  #post_prob_H1 = 1/(1 + (ML(0, 1) + ML(2, 3) + ML(3, 4) + ML(4, 5) + ML(5))/ML(1, 2))
  #post_prob_H2 = 1/(1 + (ML(0, 1) + ML(1, 2) + ML(3, 4) + ML(4, 5) + ML(5))/ML(2, 3))
  #post_prob_H3 = 1/(1 + (ML(0, 1) + ML(1, 2) + ML(2, 3) + ML(4, 5) + ML(5))/ML(3, 4))
  #post_prob_H4 = 1/(1 + (ML(0, 1) + ML(1, 2) + ML(2, 3) + ML(3, 4) + ML(5))/ML(4, 5))
  #post_prob_H5 = 1/(1 + (ML(0, 1) + ML(1, 2) + ML(2, 3) + ML(3, 4) + ML(4, 5))/ML(5))
  return(c(post_prob_H0, post_prob_H1, post_prob_H2))
  #return(c(post_prob_H0, post_prob_H1, post_prob_H2, post_prob_H3, post_prob_H4, post_prob_H5))
}

simulate_specified <- function(theta) {
  simulate(theta, 2, 2/3, 30) #Change last three arugments to desired prior and sum_ti
}

#Testing how approach does
theta <- rep(0.5, 1000) #Change first argument to desired theta
results <- lapply(theta, simulate_specified)
results_table <- data.table(matrix(unlist(results), nrow = 1000, ncol = 3, byrow = TRUE))
#results_table <- data.table(matrix(unlist(results), nrow = 1000, ncol = 6, byrow = TRUE))
correct <- results_table[V1 > V2 & V1 > V3,] #Change depending on hypothesis that is true
#correct <- results_table[V1 > V2 & V1 > V3 & V1 > V4 & V1 > V5 & V1 > V6,] #Change depending on hypothesis that is true

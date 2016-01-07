library(data.table)
library(ggplot2)

#Function that calculates the Bayes factor in favor of H0
#Under the null, theta ~ beta(a, b)
#Under the alternative, theta1 ~ beta(c, d), theta2 ~ beta(e, f)
#The control group is represent by A, the exposed group is represented by B
BF <- function(a, b, c, d, e, f, tumors_A, sample_size_A, tumors_B, sample_size_B) {
  (beta(a + tumors_A + tumors_B, b + sample_size_A + sample_size_B - tumors_A - tumors_B)*beta(c, d)*beta(e, f))/(beta(a, b)*beta(c + tumors_A, d + sample_size_A - tumors_A)*beta(e + tumors_B, f + sample_size_B - tumors_B))
}

#Function that calculates the posterior probability P(M = 1 | data)
#Takes in BF calculated with the BF function
post_prob <- function(bayesfactor) {
  1/(1 + bayesfactor)
}

#Trying different priors
prior_1.BF <- BF(1, 299, 1, 299, 0.5, 0.5, 0, 30, 2, 30)
prior_1.prob <- post_prob(prior_1.BF)
prior_2.BF <- BF(0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0, 30, 2, 30)
prior_2.prob <- post_prob(prior_2.BF)
prior_3.BF <- BF(1, 2499, 1, 2499, 2, 28, 0, 30, 2, 30)
prior_3.prob <- post_prob(prior_3.BF)

#Fisher's exact test
data <- data.table(Tumor = c(0, 2), NoTumor = c(30, 28))
fisher.test(data)

#Functions to simulate a trial under the null being true or the alternative being true
simulation_null <- function(theta, a, b, c, d, e, f) {
  control = sum(rbinom(30, 1, theta))
  exposed = sum(rbinom(30, 1, theta))
  prob = post_prob(BF(a, b, c, d, e, f, control, 30, exposed, 30))
  data_trial = data.table(Tumor = c(control, exposed), NoTumor = c(30 - control, 30 - exposed))
  pvalue = fisher.test(data_trial)$p.value
  return(c(prob, pvalue))
}
simulation_alt <- function(theta1, theta2, a, b, c, d, e, f) {
  control = sum(rbinom(30, 1, theta1))
  exposed = sum(rbinom(30, 1, theta2))
  prob = post_prob(BF(a, b, c, d, e, f, control, 30, exposed, 30))
  data_trial = data.table(Tumor = c(control, exposed), NoTumor = c(30 - control, 30 - exposed))
  pvalue = fisher.test(data_trial)$p.value
  return(c(prob, pvalue))
}

results = data.table(prob = numeric(), pvalue = numeric(), hypothesis = character(), prior = character())
#Simulation for prior 1
for(i in 1:100) {
  trial.null = simulation_null(1/2500, 1, 299, 1, 299, 0.5, 0.5)
  trial.alt = simulation_alt(1/2500, 2/30, 1, 299, 1, 299, 0.5, 0.5)
  results = rbindlist(list(results, data.table(prob = trial.null[1], pvalue = trial.null[2], hypothesis = "null", prior = "1"), data.table(prob = trial.alt[1], pvalue = trial.alt[2], hypothesis = "alternative", prior = "1")))
}
#Simulation for prior 2
for(i in 1:100) {
  trial.null = simulation_null(1/2500, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5)
  trial.alt = simulation_alt(1/2500, 2/30, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5)
  results = rbindlist(list(results, data.table(prob = trial.null[1], pvalue = trial.null[2], hypothesis = "null", prior = "2"), data.table(prob = trial.alt[1], pvalue = trial.alt[2], hypothesis = "alternative", prior = "2")))
}
#Simulation for prior 3
for(i in 1:100) {
  trial.null = simulation_null(1/2500, 1, 2499, 1, 2499, 2, 28)
  trial.alt = simulation_alt(1/2500, 2/30, 1, 2499, 1, 2499, 2, 28)
  results = rbindlist(list(results, data.table(prob = trial.null[1], pvalue = trial.null[2], hypothesis = "null", prior = "3"), data.table(prob = trial.alt[1], pvalue = trial.alt[2], hypothesis = "alternative", prior = "3")))
}

#Visualize results
ggplot(results[hypothesis == "null",], aes(x = prob, group = prior)) + geom_density(aes(fill = prior), alpha = 0.4) + labs(title = "Posterior probabilities for different groups of priors (null true)")
ggplot(results[hypothesis == "null",], aes(x = pvalue, group = prior)) + geom_density(aes(fill = prior), alpha = 0.4, adjust = 0.000001) + labs(title = "p-values for different groups of priors (null true)")
ggplot(results[hypothesis == "alternative",], aes(x = prob, group = prior)) + geom_density(aes(fill = prior), alpha = 0.4) + labs(title = "Posterior probabilities for different groups of priors (alternative true)")
ggplot(results[hypothesis == "alternative",], aes(x = pvalue, group = prior)) + geom_density(aes(fill = prior), alpha = 0.4) + labs(title = "p-values for different groups of priors (alternative true)")

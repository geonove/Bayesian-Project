# This code is combining the following two papers
# Bayesian Methods for Hidden Markov Models: Recursive Computing in the 21st Century - Author(s): Steven L. Scott
# Gibbs Sampling for Bayesian Finite Mixtures of Normal Distributions - Cristina Mollica, Luca Tardella
graphics.off()
rm(list=ls())
cat("\014")
set.seed(42)

# PARAMETERS
C <- 3 # number of possible hidden states
alpha_0 <- runif(n = C, min = 1, max = 4)
mu_0 <- 0
tau_0 <- 0.02 # precision
a_0 <- 5
b_0 <- 5

LENGTH <- 500 # length of the chain
niter <- 1000 # number of Gibbs Sampling iterations
burnin <- 200

# Defining the unknown mixture
mu_real <- rnorm(C, mu_0, sqrt(1 / tau_0))
tau_real <- rgamma(C, shape = a_0, rate = b_0)
cat("\nmu_real\n")
print(mu_real)
cat("\ntau_real\n")
print(tau_real)

# Transition probability matrix 
# (for finite mixture models all rows are equal)
Q_real <- c()
for (i in 1:C) {
  Q_real <- rbind(Q_real, gtools::rdirichlet(1, alpha_0))
}
cat("\nQ_real\n")
print(Q_real)


# Generating an Hidden Markov Chain
HMC <- function(Q) {
  h <- numeric(LENGTH)
  h[1] <- 1
  for (i in 2:LENGTH) {
    h[i] <- sample(1:C, prob = Q[h[i-1], ], size = 1)
  }
  return(h)
}
h_real <- HMC(Q_real)


# Sampling from the unknown distribution
rmix <- function(h_real, mu_real, tau_real) {
  d <- rnorm(length(h_real), mu_real[h_real], sqrt(1 / tau_real[h_real]))
  return(d)
}
d <- rmix(h_real, mu_real, tau_real)


# Plotting the Hidden Markov Model and samples
library(rgl)
plot_HMM_samples <- function(x_seq, d, h_real) {
  norms <- data.frame(x_seq)
  norms <- c()
  for (i in 1:C) {
    norms <- rbind(norms, dnorm(x_seq, mu_real[i], sqrt(1 / tau_real[i])))
  }

  for (i in 1:LENGTH) {
    plot3d(x_seq, norms[h_real[i], ], rep(i, length(x_seq)), type = "l", lwd = 4, col = h_real[i], zlim = c(1, LENGTH))
    plot3d(d[i], 0, i, size = 4, col = h_real[i], zlim = c(1, LENGTH))
  }
  #grid3d(c("x", "y+", "z"))
}
x_seq <- seq(min(d) - 2*max(sqrt(1/tau_real)), max(d)+ 2*max(sqrt(1/tau_real)), by = 0.001)
#plot_HMM_samples(x_seq, d, h_real)


# Full conditionals
sample_Q <- function(alpha_0, h) {
  Q <- matrix(0, nrow = C, ncol = C)
  NN <- matrix(0, nrow = C, ncol = C)
  for (i in 2:LENGTH) {
    NN[h[i - 1], h[i]] <- NN[h[i - 1], h[i]] + 1
  }
  for (c in 1:C) {
    Q[c, ] <- gtools::rdirichlet(1, alpha_0 + NN[c, ])
  }
  return(Q)
}

sample_tau <- function(mu, z, x, a_0, b_0, N) {
  tau <- numeric(C)
  for (c in 1:C) {
    summation <- (z == c) %*% (x - mu[c])^2
    tau[c] <- rgamma(1, shape = a_0 + N[c] / 2, rate = b_0 + summation / 2)
  }
  return(tau)
}

sample_mu <- function(tau, z, x, tau_0, mu_0, N) {
  mu <- numeric(C)
  for (c in 1:C) {
    # to avoid division by 0
    if (N[c] == 0) {
      N[c] <- 1
    }
    delta <- N[c] * tau[c] / (tau_0 + N[c] * tau[c])
    x_bar <- ((z == c) %*% x) / N[c]
    mu[c] <- rnorm(1, delta * x_bar + (1 - delta) * mu_0, sqrt(1 / (tau_0 + N[c] * tau[c])))
  }
  return(mu)
}


# Forward-Backward 
sample_h <- function(d, Q, mu, tau) {
  h <- numeric(LENGTH)
  # Forward recursion
  P <- array(0, dim = c(C, C, LENGTH))
  pi <- matrix(0, nrow = LENGTH, ncol = C)
  pi[1, 1] <- 1
  for (t in 2:LENGTH) {
    P[, , t] <- exp(log(pi[t - 1, ]) + sweep(log(Q), 2, log(dnorm(d[t], mu, sqrt(1 / tau))), "+"))
    summation <- sum(P[, , t])

    # to reconcile proportionality
    # to avoid division by 0
    if (summation == 0) {
      P[, , t] <- matrix(1 / C^2, nrow = C, ncol = C)
    } else {
      P[, , t] <- P[, , t] / summation
    }

    pi[t, ] <-  colSums(P[, , t])
  }

  # Backward recursion
  h[LENGTH] <- sample(1:C, prob = pi[LENGTH, ], size = 1)
  for (i in (LENGTH - 1):1) {
    # to avoid division by 0
    if (sum(P[, h[i + 1], i + 1]) == 0) {
      prob <- rep(1 / C, C)
    } else {
      prob <- P[, h[i + 1], i + 1]
    }
    h[i] <- sample(1:C, prob = prob, size = 1)
  }
  return(h)
}


# Gibbs Sampler
gibbs <- function(d, niter, alpha_0, mu_0, tau_0, a_0, b_0) {
  cat("\n\nGibbs Sampler\n")
  w <- gtools::rdirichlet(1, alpha_0)
  Q <- matrix(nrow = C, ncol = C)
  for (c in 1:C) {
    Q[c, ] <- w
  }
  tau <- rgamma(C, shape = a_0, rate = b_0)
  mu <- rnorm(C, mu_0, sqrt(1 / tau_0))
  h <- HMC(Q)

  mu_GS <- matrix(, nrow = niter, ncol = C)
  mu_GS[1, ] <- mu

  cat(0, "/", niter, "\n")
  for (i in 1:niter) {
    if (i %% 1 == 0) {
      cat("\r", i, "/", niter)
    }
    
    N <- tabulate(h)
    Q <- sample_Q(alpha_0, h)
    tau <- sample_tau(mu, h, d, a_0, b_0, N)
    mu <- sample_mu(tau, h, d, tau_0, mu_0, N)
    h <- sample_h(d, Q, mu, tau)
    mu_GS[i, ] <- mu
  }
  cat("\nmu\n")
  print(mu)
  cat("\nQ\n")
  print(Q)
  return(mu_GS)
}

mu_GS <- gibbs(d, niter, alpha_0, mu_0, tau_0, a_0, b_0)

x11()
matplot(mu_GS, main="Markov Chain for mu", type = 'l', xlim = c(0, niter), lty = 1, lwd = 2)

mu_hat <- apply(mu_GS[burnin:niter, ], 2, mean)
cat("mu_real: \n")
print(mu_real)
cat("mu_Gs: \n")
print(mu_hat)
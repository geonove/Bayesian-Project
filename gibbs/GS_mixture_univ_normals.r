# This code is implementing the paper:
# Gibbs Sampling for Bayesian Finite Mixtures of Normal Distributions - Cristina Mollica, Luca Tardella
# https://web.uniroma1.it/memotef/sites/default/files/file%20lezioni/Lezione3_CMollica.pdf

set.seed(42)

# PARAMETERS
C <- 4 # number of classes
alpha_0 <- runif(n = C, min = 1, max = 4)
mu_0 <- 0
tau_0 <- 0.02 # precision
a_0 <- 5
b_0 <- 5

N_SAMPLES <- 2000 # number of samples
niter <- 300 # number of Gibbs Sampling iterations


# Defining the unknown mixture
w_real <- gtools::rdirichlet(1, alpha_0)
mu_real <- rnorm(C, mu_0, sqrt(1 / tau_0))
tau_real <- rgamma(C, shape = a_0, rate = b_0)
cat("w_real :", w_real, "\n")
cat("mu_real :", mu_real, "\n")
cat("tau_real :", tau_real, "\n")


# Sampling from the unknown distribution
rmix <- function(w_real, mu_real, tau_real) {
  z <- sample(1:C, prob = w_real, size = N_SAMPLES, replace = TRUE)
  x <- rnorm(N_SAMPLES, mu_real[z], sqrt(1 / tau_real[z]))
  return(x)
}
x <- rmix(w_real, mu_real, tau_real)


# Plotting
plot_mixture <- function(x_seq, w, mu, tau, title) {
  norms <- data.frame(x_seq)
  col_names <- c("x_seq")
  for (i in 1:C) {
    norms <- cbind(norms, dnorm(x_seq, mu[i], sqrt(1 / tau[i])) * w[i])
    col_names <- c(col_names, paste("norm", i, sep = ""))
  }

  norms_matrix <- data.matrix(norms[2:length(norms)])
  mixture <- rowSums(norms_matrix)
  norms <- cbind(norms, mixture)

  col_names <- c(col_names, "mixture")
  colnames(norms) <- col_names
  library("tidyr")
  long_norms <- pivot_longer(norms, cols = all_of(col_names[2:length(col_names)]))

  library(ggplot2)
  line_sizes <- rep(0.6, C + 1)
  line_sizes[1] <- 2
  line_colors <- rep(1, C + 1)
  line_colors[1] <- 2
  ggp <- ggplot(long_norms) + 
        geom_line(aes(x_seq, value, col = name, size = name)) + 
        scale_size_manual(values = line_sizes) +
        scale_color_manual(values = line_colors) +
        ggtitle(title) +
        theme(plot.title = element_text(size = 20, face = "bold"))

  # x_seq <- data.frame(x_seq)
  # ggp <- ggp +
  #   geom_histogram(data = x_seq, aes(x_seq = x_seq, y = after_stat(density)), alpha = 0.3,
  #                 bins = 300, position = "identity", lwd = 0.2) +
  #   ggtitle("Unknown mixture of gaussians + samples")

  # ggp <- x_seq[[1]]

  x11(type = "cairo")
  plot(ggp)
}
x_seq <- seq(min(x), max(x), by = 0.001)
title <- "Unknown mixture of gaussians"
plot_mixture(x_seq, w_real, mu_real, tau_real, title)


# Full conditionals
sample_w <- function(alpha_0, N) {
  alpha_new <- alpha_0 + N
  w <- gtools::rdirichlet(1, alpha_new)
  return(w)
}

sample_tau <- function(mu, z, x, a_0, b_0, N) {
  tau <- numeric(C)
  for (c in 1:C) {
    summation <- z[, c] %*% (x - mu[c])^2
    a_new <- a_0 + N[c] / 2
    b_new <- b_0 + summation / 2
    tau[c] <- rgamma(1, shape = a_new, rate = b_new)
  }
  return(tau)
}

sample_mu <- function(tau, z, x, tau_0, mu_0, N) {
  mu <- numeric(C)
  for (c in 1:C) {
    if (N[c] == 0) {
      N[c] <- 1
    }
    delta <- N[c] * tau[c] / (tau_0 + N[c] * tau[c])
    x_bar <- (z[, c] %*% x) / N[c]
    mu[c] <- rnorm(1, delta * x_bar + (1 - delta) * mu_0, sqrt(1 / (tau_0 + N[c] * tau[c])))
  }
  return(mu)
}

sample_z <- function(mu, tau, w, x) {
  z <- matrix(, nrow = length(x), ncol = length(w))
  for (i in 1:N_SAMPLES) {
    prob <- c()
    summation <- 0
    for (c in 1:C) {
      summation <- summation + w[c] * dnorm(x[i], mu[c], sqrt(1 / tau[c]))
    }
    
    # avoid division by 0
    if (summation != 0) {
      for (c in 1:C) {
        prob <- c(prob, w[c] * dnorm(x[i], mu[c], sqrt(1 / tau[c])) / summation)
      }
    } else {
      prob <- runif(n = C, min = 0, max = 1)
    }
    z[i, ] <- rmultinom(1, 1, prob)
  }
  return(z)
}


# Gibbs Sampler
gibbs <- function(x, niter, C, alpha_0, mu_0, tau_0, a_0, b_0) {
  # Initialize the Markov Chain
  w <- gtools::rdirichlet(1, alpha_0)
  tau <- rgamma(C, shape = a_0, rate = b_0)
  mu <- rnorm(C, mu_0, sqrt(1 / tau_0))
  z <- matrix(, nrow = length(x), ncol = length(w))
  for (i in 1:N_SAMPLES) {
    z[i, ] <- rmultinom(1, 1, w)
  }

  # Save the Markov Chain of mu
  mu_GS <- matrix(, nrow = length(x), ncol = length(mu))
  mu_GS[1, ] <- mu

  cat("\nGibbs Sampling\n")
  cat(0, "/", niter, "\n")
  for (i in 1:niter) {
    if (i %% 1 == 0) {
      cat("\r", i, "/", niter)
    }

    N <- c()
    for (c in 1:C) {
      N <- c(N, sum(z[, c]))
    }

    w <- sample_w(alpha_0, N)
    tau <- sample_tau(mu, z, x, a_0, b_0, N)
    mu <- sample_mu(tau, z, x, tau_0, mu_0, N)
    z <- sample_z(mu, tau, w, x)

    mu_GS[i, ] <- mu
  }

  return(mu_GS)
}

mu_GS <- gibbs(x, niter, C, alpha_0, mu_0, tau_0, a_0, b_0)

x11()
matplot(mu_GS, main="Markov Chain for mu", type = 'l', xlim = c(0, niter), lty = 1, lwd = 2)

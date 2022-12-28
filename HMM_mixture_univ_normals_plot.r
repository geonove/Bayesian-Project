# This code is plotting an Hidden Markov Model with each observation generated from a mixture of univariate gaussians


# PARAMETERS
H <- 4 # number of hidden states
C <- 3 # number of gaussians for each mixture
alpha_0 <- runif(n = C, min = 1, max = 4)
mu_0 <- 0
tau_0 <- 0.02 # precision
a_0 <- 5
b_0 <- 5

LENGTH <- 100
N <- 100 # number of samples
niter <- 200 # number of Gibbs Sampling iterations


# Defining the unknown mixtures
W_REAL <- matrix(, nrow = H, ncol = C)
MU_REAL <- matrix(, nrow = H, ncol = C)
TAU_REAL <- matrix(, nrow = H, ncol = C)

for (i in 1:H) {
  W_REAL[i, ] <- gtools::rdirichlet(1, alpha_0)
  MU_REAL[i, ] <- rnorm(C, mu_0, sqrt(1 / tau_0))
  TAU_REAL[i, ] <- rgamma(C, shape = a_0, rate = b_0)
  cat("W_REAL[", i, "]  :", W_REAL[i, ], "\n")
  cat("MU_REAL[", i, "] :", MU_REAL[i, ], "\n")
  cat("TAU_REAL[", i, "]:", TAU_REAL[i, ], "\n")
}

# Transition probability matrix 
Q <- rbind(c(0.8, 0.05, 0.05, 0.05), 
          c(0.05, 0.8, 0.05, 0.05), 
          c(0.05, 0.05, 0.8, 0.05), 
          c(0.05, 0.05, 0.05, 0.8))
cat("\n Transition matrix\n")
print(Q)


# Generating the sequence of mixtures
seq_mixtures <- function(n, Q, H) {
  h <- c(1)
  for (i in 1:n) {
    h <- c(h, sample(1:H, prob = Q[h[length(h)], ], size = 1))
  }
  return(h)
}
h_real <- seq_mixtures(n = LENGTH, Q, H)
print(h_real)

# Generating the mixtures
x <- seq(-20, 20, by = 0.001)
MIXTURES <- matrix(, nrow = H, ncol = length(x))
norms <- data.frame(x)
col_names <- c("x")
for (i in 1:H) {
  norms <- data.frame(x)
  for (c in 1:C) {
    norms <- cbind(norms, dnorm(x, MU_REAL[i, c], sqrt(1 / TAU_REAL[i, c])) * W_REAL[i, c])
  }
  norms_matrix <- data.matrix(norms[2:length(norms)])
  MIXTURES[i, ] <- rowSums(norms_matrix)
}
print(dim(MIXTURES))


# Plotting
library(rgl)
z_plot <- c()
y_plot <- c()
x_plot <- c()
for (i in 1:LENGTH) {
  z_plot <- c(z_plot, rep(i, length(x)))
  x_plot <- c(x_plot, x)
  y_plot <- c(y_plot, MIXTURES[h_real[i], ])

  plot3d(x, MIXTURES[h_real[i], ], rep(i, length(x)), type = "l", lwd = 4, col = h_real[i], zlim = c(1, LENGTH))
}
#grid3d(c("x", "y+", "z"))



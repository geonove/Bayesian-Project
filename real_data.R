setwd("C:/Users/andre/Documents/Università/Magistrale/2° Anno/1° Semestre/Bayesian statistics/Progetto/codes")
#setwd("/Users/andrespasinetti/Bayesian-Project-7")

library("readxl")
library(corrplot)
library(caret)
library(gRbase)
library(GLSME)
library(BCDAG)
library(gss)
library(graph)
library(Rcpp)
library(rgl)
library(car)
library('salso')

source("dags/propose_DAG.R")
source("dags/operation.R")
source("dags/acceptreject_DAG.R")
source("dags/new_bcdag.R")
source("dags/update_DAG.R")

set.seed(55)

graphics.off()

data<- read_excel("dataset.xlsx")
data<- data.frame(data)
head(data)

for(i in 1:1110) {
  for (j in 3:44)
  data[i, j] <- (data[i+1, j] - data[i, j]) / data[i, j]
}
data <- data[1:1110, ]
data[3:dim(data)[2]] <- scale(data[3:dim(data)[2]])
print(head(data))
data <- data[, c(which( !is.na(data[1,]), arr.ind=TRUE)[,2])]
print(head(data))

df <- cor(data[, 3:dim(data)[2]])
x11()
corrplot(df, method="circle")

risk <- data[data$Y == 1, ]
nonrisk <- data[data$Y == 0, ]

M1<-cor(risk[, 3:dim(data)[2]])
x11()
corrplot(M1, method="circle")

M2<-cor(nonrisk[, 3:dim(data)[2]])
x11()
corrplot(M2, method="circle")

M <- abs(M1 - M2)
x11()
corrplot(M, method="circle")

df <- data[3:dim(data)[2]]
print(dim(df))


########################################################################

# PARAMETERS
C <- 2 # number of classes
alpha_0 <- rep(1, C) 
q <- dim(df)[2]

niter <- 2500 # number of Gibbs Sampling iterations
burnin <- 300

N_SAMPLES <- 1110
x <- df[1:N_SAMPLES, ]

# Generating an Hidden Markov Chain
HMC <- function(Q) {
  z <- c(1)
  for (i in 2:N_SAMPLES) {
    z <- c(z, sample(1:C, prob = Q[z[length(z)], ], size = 1))
  }
  return(z)
}


# Full conditionals
sample_Q <- function(alpha_0, z) {
  Q <- matrix(, nrow = C, ncol = C)
  NN <- matrix(0, nrow = C, ncol = C)
  for (i in 2:N_SAMPLES) {
    NN[z[i - 1], z[i]] <- NN[z[i - 1], z[i]] + 1
  }
  
  for (i in 1:C) {
    Q[i, ] <- gtools::rdirichlet(1, alpha_0 + NN[i, ])
  }
  
  return(Q)
}

# Forward-Backward 
sample_z <- function(d, Q, mu, Sigma) {
  z <- numeric(N_SAMPLES)
  # Forward recursion
  P <- array(0, dim = c(C, C, N_SAMPLES))
  pi <- matrix(0, nrow = N_SAMPLES, ncol = C)
  pi[1, 1] <- 1
  for (t in 2:N_SAMPLES) {
    for (r in 1:C) {
      for (s in 1:C) {
        P[r, s, t] <- exp(log(pi[t - 1, r]) + log(Q[r, s]) + log(dmvnorm(d[t, ], mu[s,], (Sigma[, , s]))))  
      }
    }
    summation <- sum(P[, , t])
    
    # to reconcile proportionality
    # to avoid division by 0
    if (summation == 0) {
      P[, , t] <- matrix(1 / C^2, nrow = C, ncol = C)
    } else {
      P[, , t] <- P[, , t] / summation
    }
    
    for (s in 1:C) {
      pi[t, s] <-  sum(P[, s, t])
    }
  }
  
  # Backward recursion
  z[N_SAMPLES] <- sample(1:C, prob = pi[N_SAMPLES, ], size = 1)
  for (i in (N_SAMPLES - 1):1) {
    # to avoid division by 0
    if (sum(P[, z[i + 1], i + 1]) == 0) {
      prob <- rep(1 / C, C)
    } else {
      prob <- P[, z[i + 1], i + 1]
    }
    z[i] <- sample(1:C, prob = prob, size = 1)
  }
  return(z)
}


# Gibbs Sampler
gibbs <- function(x, niter, C, alpha_0) {
  x <- data.matrix(x)
  mu <- matrix(0, nrow = C, ncol = q)
  # Hyperparameters
  w_dags <- 0.5 # probability of adjacency in dag
  
  # Initialize DAG
  DAG <- array(dim = c(q, q, C))
  Sigma <- array(dim = c(q, q, C))
  for (c in 1:C) {    
    DAG[, , c] <- rDAG(q = q, w = w_dags)
    L <- matrix(runif(n = q*(q), min = 0, max = 1), q, q)    
    L <- L * DAG[, , c] 
    diag(L) <- 1
    D <- diag(1, q)
    Sigma[, , c] <- solve(t(L))%*%solve(D)%*%solve(L)       
  }
  
  # Initialize the Markov Chain
  w <- gtools::rdirichlet(1, alpha_0) 
  Q <- c()
  for (h in 1:C) {
    Q <- rbind(Q, w)
  }
  # Initialize cluster allocation
  z <- HMC(Q) # Z | Q distr MC
  
  N <- rep(0, C) # Samples per cluster
  
  # Save the Markov Chain 
  Q_GS <- array(dim = c(C, C, niter) )
  Q_GS[, , 1] <- Q
  
  z_GS <- array(dim = c(N_SAMPLES , niter))
  z_GS[, 1] <- z
  
  DAG_GS <- array(dim = c(q, q, C, niter))
  DAG_GS[, , , 1] <- DAG
  
  Sigma_GS <- array(dim = c(q, q, C, niter))
  Sigma_GS[, , , 1] <- Sigma
  
  
  cat("\nGibbs Sampling\n")
  cat(0, "/", niter, "\n")
  for (i in 1:niter) {
    if (i %% 10 == 0) {
      cat(i, "/", niter, "\n")
    }
    
    for (c in 1:C) {
      N[c] <- sum(z == c)
    }
    
    Q <- sample_Q(alpha_0, z)
    z <- sample_z(x, Q, mu, Sigma)
    out <- update_DAGS(DAG = DAG, data = x, z = z, a = q, U=diag(1, q), w = w_dags)
    DAG <- out$DAG
    Sigma <- out$Omega
    
    Q_GS[, , i] <- Q
    z_GS[, i] <- z
    DAG_GS[, , , i] <- DAG
    Sigma_GS[, , , i] <- Sigma
  }
  
  return(list("Q_GS" = Q_GS, "z_GS" = z_GS, "DAG_GS" = DAG_GS, "Sigma_GS" = Sigma_GS))
}
# Running Gibbs
#mix <- gibbs(x, niter, C, alpha_0)


#saveRDS(mix, file = "perc_diff_3000_noscale.rds")

mix <- readRDS("diff_perc_10000.RDS")

z_mcmc <- t(mix$z_GS[,burnin:niter])
z_binder <- salso(z_mcmc, "binder", maxNClusters = C)
z_binder <- abs(z_binder -1)

confusionMatrix(data=as.factor(z_binder), reference=as.factor(data[1:N_SAMPLES, 1]))



z_real <- abs(data[1:N_SAMPLES, 1]+1)


### convergence
z_error <- array(, dim = niter)
for (i in 1:(niter)) {
  z_error[i] <- sum(z_real != mix$z_GS[, i])
  }


png(filename="z_error.png", width = 4, height = 4, units = 'in',res=300)
#x11()
matplot(z_error, type='l', lty = 1, lwd = 2, xlab="iteration", ylab="# misclassified hidden states")
dev.off()
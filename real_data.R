setwd("C:/Users/andre/Documents/Università/Magistrale/2° Anno/1° Semestre/Bayesian statistics/Progetto/codes")
library("readxl")
library(corrplot)

graphics.off()

data<- read_excel("dataset.xlsx")
data<- data.frame(data)
#data[, 3:44] <- scale(data[, 3:44])
head(data)

for(i in 1:1110) {
  for (j in 3:44)
  data[i, j] <- (data[i+1, j] - data[i, j]) 
}
data <- data[1:1110, ]
#data <- scale(data)

# remove higly correlated features
#M <- cor(data[, 3:44])
#cols <- findCorrelation(M, cutoff = 0.7)
#data <- data[, -cols]

risk <- data[data$Y == 1, ]
nonrisk <- data[data$Y == 0, ]

M1<-cor(risk[, 3:44])
head(round(M1,2))
x11()
corrplot(M1, method="circle")

M2<-cor(nonrisk[, 3:44])
head(round(M2,2))
x11()
corrplot(M2, method="circle")


M <- abs(M1 - M2)
head(round(M,2))
x11()
corrplot(M, method="circle")


cols <- c()
data3 <- data[, c(1,2, cols)]
head(data3)

summary(data3)

for(i in 1:1110) {
  data3[i, 3:5] <- (data3[i+1, 3:5] - data3[i, 3:5])
}

head(data3)
summary(data3)
cols <- colnames(data)[c(12, 14, 32)]
df <- data[, c(7, 8, 9)]

open3d()
scatter3d(x = df[, 1], y = df[, 2], z = df[, 3], point.col= data3[1:1110, 1]+2, surface = FALSE)

print(cols)

x11()
plot(df[, 1:2], type="p", col=data3[1:1110, 1]+2, pch=20)


########################################################################

library(gRbase)
library(GLSME)
library(BCDAG)
library(gss)
library(graph)
library(Rcpp)
library(rgl)
library(car)


source("dags/propose_DAG.R")
source("dags/operation.R")
source("dags/acceptreject_DAG.R")
source("dags/new_bcdag.R")
source("dags/update_DAG.R")

### Plot 3d
#install.packages(c("rgl", "car"))
graphics.off()
set.seed(55)

# PARAMETERS
C <- 2 # number of classes
alpha_0 <- runif(n = C, min = 1, max = 4) ## <- rep(1, C)
q <- 3

niter <- 100 # number of Gibbs Sampling iterations


# Generating an Hidden Markov Chain
HMC <- function(Q) {
  z <- c(1)
  for (i in 2:N_SAMPLES) {
    z <- c(z, sample(1:C, prob = Q[z[length(z)], ], size = 1))
  }
  return(z)
}



N_SAMPLES <- 1110
x <- df[1:N_SAMPLES, ]


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
sample_z <- function(d, Q, mu, tau) {
  z <- numeric(N_SAMPLES)
  # Forward recursion
  P <- array(0, dim = c(C, C, N_SAMPLES))
  pi <- matrix(0, nrow = N_SAMPLES, ncol = C)
  pi[1, 1] <- 1
  for (t in 2:N_SAMPLES) {
    for (r in 1:C) {
      for (s in 1:C) {
        P[r, s, t] <- exp(log(pi[t - 1, r]) + log(Q[r, s]) + log(dmvnorm(d[t, ], mu[s,], solve(tau[, , s]))))   # dnorm -> dmvtnorm (libreria mvtfast)
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
  Omega <- array(dim = c(q, q, C))
  for (c in 1:C) {    
    DAG[, , c] <- rDAG(q = q, w = w_dags)
    L <- matrix(runif(n = q*(q), min = -10, max = 10), q, q)     ### va bene mettere questi min e max? 
    L <- L * DAG[, , c] 
    diag(L) <- 1
    D <- diag(1, q)
    Omega[, , c] <- L%*%solve(D)%*%t(L)       
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
  
  Omega_GS <- array(dim = c(q, q, C, niter))
  Omega_GS[, , , 1] <- Omega
  
  
  cat("\nGibbs Sampling\n")
  cat(0, "/", niter, "\n")
  for (i in 1:niter) {
    if (i %% 5 == 0) {
      cat(i, "/", niter, "\n")
    }
    
    for (c in 1:C) {
      N[c] <- sum(z == c)
    }
    
    Q <- sample_Q(alpha_0, z)
    z <- sample_z(x, Q, mu, Omega)
    out <- update_DAGS(DAG = DAG, data = x, z = z, a = q, U=diag(1, q), w = w_dags)
    DAG <- out$DAG
    Omega <- out$Omega
    
    Q_GS[, , i] <- Q
    z_GS[, i] <- z
    DAG_GS[, , , i] <- DAG
    Omega_GS[, , , i] <- Omega
  }
  
  return(list("Q_GS" = Q_GS, "z_GS" = z_GS, "DAG_GS" = DAG_GS, "Omega_GS" = Omega_GS))
}
# Running Gibbs
mix <- gibbs(x, niter, C, alpha_0)



# Plot of true values
if (q==2) {
  x11()
  plot(x, type="p", col=z_real, pch=20)
} else if (q==3) {
  open3d()
  scatter3d(x = x[, 1], y = x[, 2], z = x[, 3], point.col= data3[1:N_SAMPLES, 1] + 2, surface = FALSE)
}
# Plot of sampled values
if (q==2) {
  x11()
  plot(x, type="p", col=mix$z_GS[, niter], pch=20)
} else if (q==3) {
  open3d()
  x <- data.frame(x)
  scatter3d(x = x[, 1], y = x[, 2], z = x[, 3], point.col= mix$z_GS[, niter] + 1, surface = FALSE)
}


cat("\nDAG_real:\n")
print(DAG_real)
cat("\nDAG_GS:\n")
print(mix$DAG_GS[, , , niter])
cat("\nOmega_real:\n")
print(Omega_real)
cat("\nOmega:\n")
print(mix$Omega_GS[, , , niter])

### DOMANDE ###
# - matrice DAG non-triangolare-bassa/ triangolare-alta, va bene? 
# - range valori matrice L va bene? O deve per forza essere [0.1, 1] ? (vedi riga 34)
# - per Q e z ha senso parlare di prior? O semplicemente stiamo inizializzando qualcosa. 
#   (Probabilmente per Q ha senso ed è un vettore di Dirichlet)
# - prior di Q: va bene fatta riga per riga o c'è una prior per Q intera
# - plot della catena: come possiamo fare per la matrice della covarianza che ha q x q componenti?

### TO DO ###
# - aggiungere plot catene di markov
# - sistemare grafici (nome assi e magari legenda)
# - cambiare tutte le h in z
# - sistemare tutto one_hot o no
# - inizializzare a priori tutte le liste e gli array alla dimensione corretta
# - variabili globali
# - rigurdare/riordinare i codici in generale
# - readme

x11()
par(mfrow=c(4,1))
matplot(mix$z_GS[100:160, 1], type='l', lty = 1, lwd = 2)
matplot(mix$z_GS[100:160, niter/2], type='l', lty = 1, lwd = 2)
matplot(mix$z_GS[100:160, niter], type='l', lty = 1, lwd = 2)
matplot(z_real[100:160], type='l', lty = 1, lwd = 2)

sum(z_real == mix$z_GS[, niter])

# get_map (pacchetto BCDAG) -> maximum a posteriori DAG (dag modale)
# get_mpm -> median probability model
# calcolcare distanza (vi.dist o binder loss) tra cluster per ogni iterazione e plottarli
# plot DAG: library "network", "pcalg"
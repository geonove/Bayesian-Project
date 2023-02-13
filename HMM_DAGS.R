library(gRbase)
library(GLSME)
library(BCDAG)
library(gss)
library(graph)
library(Rcpp)
library(rgl)
library(car)

rm(list=ls())
graphics.off()
setwd("C:/Users/andre/Documents/Università/Magistrale/2° Anno/1° Semestre/Bayesian statistics/Progetto/codes")
#setwd("/Users/andrespasinetti/Bayesian-Project-4")

source("dags/propose_DAG.R")
source("dags/operation.R")
source("dags/acceptreject_DAG.R")
source("dags/new_bcdag.R")
source("dags/update_DAG.R")

### Plot 3d
#install.packages(c("rgl", "car"))

set.seed(55)

# PARAMETERS
C <- 3 # number of classes
alpha_0 <- rep(1, C)
q <- 3
N_SAMPLES <- 700 # number of samples
burnin <- 300
niter <- 1000 # number of Gibbs Sampling iterations

w_real <- gtools::rdirichlet(1, alpha_0)
DAG_real <- array(dim = c(q, q, C))
Sigma_real <- array(dim = c(q, q, C))
for (c in 1:C) {    
  DAG_real[, , c] <- rDAG(q = q, w = 0.5)
  L <- matrix(runif(n = q*(q), min = -10, max = 10), q, q)     ### va bene mettere questi min e max? 
  L <- L * DAG_real[, , c] 
  diag(L) <- 1
  D <- diag(1, q)
  Sigma_real[, , c] <- solve(t(L))%*%(D)%*%solve(L)       
}

# Transition probability matrix 
# For finite mixture models, all rows are equal
Q_real <- c()
for (i in 1:C) {
  Q_real <- rbind(Q_real, gtools::rdirichlet(1, alpha_0))
}
cat("\nTransition matrix\n")
print(Q_real)

# Generating an Hidden Markov Chain
HMC <- function(Q) {
  z <- c(1)
  for (i in 2:N_SAMPLES) {
    z <- c(z, sample(1:C, prob = Q[z[length(z)], ], size = 1))
  }
  return(z)
}
z_real <- HMC(Q_real)


rmix <- function(w_real, Omega_real, z_real) {
  x <- matrix(, nrow = N_SAMPLES, ncol = q)
  for (i in 1:N_SAMPLES) {
    x[i, ] <- mvtnorm::rmvnorm(1, sigma=(Sigma_real[, ,z_real[i]]))
  }
  
  return (x)
}
x <- rmix(w_real, Sigma_real, z_real)
x <- data.frame(x)

z_real <- unlist(z_real)


# Full conditionals
sample_Q <- function(alpha_0, z) {
  Q <- matrix(0, nrow = C, ncol = C)
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
    #for (r in 1:C) {
      for (s in 1:C) {
        P[, s, t] <- exp(log(pi[t - 1, ]) + log(Q[, s]) + log(dmvnorm(d[t, ], mu[s,], (Sigma[, , s]))))   # dnorm -> dmvtnorm (libreria mvtfast)
      }
    #}
    
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
  w_dags <- 0.2 # probability of adjacency in dag
  
  # Initialize DAG
  DAG <- array(dim = c(q, q, C))
  Sigma <- array(dim = c(q, q, C))
  for (c in 1:C) {    
    DAG[, , c] <- rDAG(q = q, w = w_dags)
    L <- matrix(runif(n = q*(q), min = -1, max = 1), q, q)     ### va bene mettere questi min e max? 
    L <- L * DAG[, , c] 
    diag(L) <- 1
    D <- diag(1, q)
    Sigma[, , c] <- solve(t(L))%*%(D)%*%solve(L)       
  }
  
  
  # Initialize the Markov Chain
  w <- gtools::rdirichlet(1, alpha_0) 
  Q <- t(matrix(w, nrow = C, ncol = length(w)))
  
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
    if (i %% 5 == 0) {
      cat(i, "/", niter, "\n")
    }
    
    #for (c in 1:C) {
    #  N[c] <- sum(z == c)
    #}
    
    Q <- sample_Q(alpha_0, z)
    z <- sample_z(x, Q, mu, Sigma)
    out <- update_DAGS(DAG = DAG, data = x, z = z, a = q, U=diag(1, q), w = w_dags)
    DAG <- out$DAG
    Sigma <- out$Sigma
    
    Q_GS[, , i] <- Q
    z_GS[, i] <- z
    DAG_GS[, , , i] <- DAG
    Sigma_GS[, , , i] <- Sigma
  }
  
  return(list("Q_GS" = Q_GS, "z_GS" = z_GS, "DAG_GS" = DAG_GS, "Sigma_GS" = Sigma_GS))
}
# Running Gibbs
#mix <- gibbs(x, niter, C, alpha_0)
#saveRDS(mix, "output1000MCMC.RDS")
mix <- readRDS("output1000MCMC.RDS")


# Masking z_real so that
getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

library(hash)

mappa <- array(, dim=c(C))
mask <- matrix(, nrow = N_SAMPLES, ncol = C)
for (c in 1:C) {
  idx <- (z_real == c)
  moda <- getmode(mix$z_GS[idx, niter])
  mappa[c] <- moda
  mask[, moda] <- idx
}




for (c in 1:C) {
  z_real[mask[, c]] <- c 
}

# Plot of true values
if (q==2) {
  x11()
  plot(x, type="p", col=z_real, pch=20)
} else if (q==3) {
  open3d()
  scatter3d(x = x[, 1], y = x[, 2], z = x[, 3], point.col= z_real, surface = FALSE)
  title3d("True values")
  rgl.snapshot('trueValues.png', fmt = 'png')
}
# Plot of sampled values
if (q==2) {
  x11()
  plot(x, type="p", col=mix$z_GS[, niter], pch=20)
} else if (q==3) {
  open3d()
  x <- data.frame(x)
  scatter3d(x = x[, 1], y = x[, 2], z = x[, 3], point.col= mix$z_GS[, niter], surface = FALSE)
  title3d("Sampled values")
}



cat("\nDAG_real:\n")
print(DAG_real)
cat("\nDAG_GS:\n")
print(mix$DAG_GS[, , , niter])
cat("\nSigma_real:\n")
print(Sigma_real)
cat("\nSigma:\n")
print(mix$Sigma_GS[, , , niter])

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
matplot(mix$z_GS[200:260, 1], type='l', lty = 1, lwd = 2)
matplot(mix$z_GS[200:260, niter/2], type='l', lty = 1, lwd = 2)
matplot(mix$z_GS[200:260, niter], type='l', lty = 1, lwd = 2)
matplot(z_real[200:260], type='l', lty = 1, lwd = 2)

sum(z_real == mix$z_GS[, niter])

# get_map (pacchetto BCDAG) -> maximum a posteriori DAG (dag modale)
# get_mpm -> median probability model
# calcolcare distanza (vi.dist o binder loss) tra cluster per ogni iterazione e plottarli
# plot DAG: library "network", "pcalg"


### Infernce about DAG struture
out <- new_bcdag(list(Graphs = mix$DAG_GS[, , 1, burnin:niter]), input = x, type = "collapsed")
get_MAPdag(out)
get_MPMdag(out)
DAG_real[, , mappa[1]]
out <- new_bcdag(list(Graphs = mix$DAG_GS[, , 2, burnin:niter]), input = x, type = "collapsed")
get_MAPdag(out)
get_MPMdag(out)
DAG_real[, , mappa[2]]
out <- new_bcdag(list(Graphs = mix$DAG_GS[, , 3, burnin:niter]), input = x, type = "collapsed")
get_MAPdag(out)
get_MPMdag(out)
DAG_real[, , mappa[3]]


#################################
## Posterior similarity matrix ##
#################################
n <- dim(mix$z_GS)[1]
simil_mat <- matrix(0, nrow = n, ncol = n)
for (t in (burnin + 1):niter){
  
  simil_mat = simil_mat + (matrix(mix$z_GS[,t], nrow = n, ncol = n) == t(matrix(mix$z_GS[,t], nrow = n, ncol = n)))*1
  
}

simil_probs = simil_mat/(niter-burnin)


x11()
par(mfrow = c(1,1))

library(fields)
colori = colorRampPalette(c('white','black'))
par(mar = c(4,4,1,2), oma = c(0.5,0.5,0.5,0.5), cex = 1, mgp = c(2,0.5,0))
image.plot(t(simil_probs), col = colori(100), zlim = c(0,1), cex.sub = 1, xlab = "i'", ylab = "i", axes = F,
           horizontal = F, legend.shrink = 1, axis.args = list(cex.axis = 1), cex.lab = 1)
title("Similarity matrix")
axis(1, at = seq(1/n*10,1,by=1/n*10), lab = seq(10,n, by = 10), las = 1)
axis(2, at = seq(1/n*10,1,by=1/n*10), lab = seq(10,n, by = 10), las = 1)


#######################
## Cluster estimates ##
#######################

library(mcclust)
library(caret)

# Estimated clustering (i and i' in the same cluster iff simil.probs > 0.5)

from_simil_to_clust = function(simil_probs){
  simil_mat = round(simil_probs)
  clust_ind = c()

  for(i in n:1){
    clust_ind[simil_mat[i,] == 1] = i
  }
  
  clust_ind = as.factor(clust_ind)
  print(clust_ind)
  levels(clust_ind) = 1:(length(levels(clust_ind)))
  
  return(clust_ind)
}

z_estimated <- from_simil_to_clust(simil_probs)

sum(z_estimated == z_real)
sum(mix$z_GS[, niter] == z_real)
# Variation of Information (VI) between true and estimated clustering

vi.dist(z_real, from_simil_to_clust(simil_probs))

x11()
par(mfrow=c(1,2))
plot(density(z_real), main = "true clusters", xlab = "clusters", ylab = "density", xlim = c(0, 4))
lines(1:length(table(z_real)), (table(as.numeric(z_real)))/N_SAMPLES, type = 'h', lwd=3, col = "blue")
plot(density(as.numeric(z_estimated)), main = "estimated clusters",  xlab = "clusters", ylab = "density")
lines(1:length(table(as.numeric(z_estimated))), (table(as.numeric(z_estimated)))/N_SAMPLES, type = 'h', lwd=3, col = "blue")

z_estimated[as.numeric(z_estimated) > C] = C 
z_estimated <- as.numeric(z_estimated)
confusionMatrix(data = as.factor(z_estimated), reference = as.factor(z_real))


### SALSO

library('salso')
?salso

dim(mix$z_GS)
z_mcmc <- matrix(t(mix$z_GS[,burnin:niter]), nrow = niter-burnin+1, ncol = N_SAMPLES)

z_binder <- salso(z_mcmc, "binder", maxNClusters = C)
sum(z_binder == z_real)

table(z_binder, z_real)
psm <- salso::psm(z_mcmc, nCores = 4)
x11()
heatmap(psm, main = "posterior similarity matrix")


x11()
par(mfrow=c(1,2))
plot(density(z_real), main = "true clusters", xlab = "clusters", ylab = "density", xlim = c(0, 4))
lines(1:length(table(z_real)), (table(as.numeric(z_real)))/N_SAMPLES, type = 'h', lwd=3, col = "blue")
plot(density(as.numeric(z_binder)), main = "estimated clusters",  xlab = "clusters", ylab = "density", xlim = c(0, 4))
lines(1:length(table(as.numeric(z_binder))), (table(as.numeric(z_binder)))/N_SAMPLES, type = 'h', lwd=3, col = "blue")

confusionMatrix(data = as.factor(z_binder), reference = as.factor(z_real))

###################################
# estimated Sigma

Sigma_hat <- array(0, dim = c(q, q, C))


for (i in (burnin+1):niter) {
  for (j in 1:C) {
    Sigma_hat[, , j] <- Sigma_hat[, , j] + mix$Sigma_GS[, , j, i]
  }
  
}
Sigma_hat <- Sigma_hat / (niter-burnin)
Sigma_hat[, , 1]
Sigma_real[, , 1]

library('SMFilter')
FDist2(scale(Sigma_hat[, , 1]), scale(Sigma_real[, , 3]))


### label switching
x11()
par(mfrow = c(2,2))
for (k in 1:(q-1)) {
  for (l in 1:(q-1)) {
    var_hat <- matrix(, nrow = niter, ncol = C)

    for (i in 1:niter) {
      for (j in 1:C) {
        var_hat[i, j] <- mix$Sigma_GS[k, l, j, i]
      }
    }

    var_real <- matrix(, nrow = 1, ncol = C)

    for (i in 1:C) {
      var_real[i] <- Sigma_real[k, l, i]
    }

    if (k == 1 && l == 1) {
      matplot(var_hat, main=paste("MC Sigma", k,l), type = 'l', xlim = c(0, niter), lty = 1, lwd = 1, 
              ylab = paste("Sigma",k,l), xlab = "niter", ylim = c(0,150))
    }
    else {
    matplot(var_hat, main=paste("MC Sigma", k,l), type = 'l', xlim = c(0, niter), lty = 1, lwd = 1, 
            ylab = paste("Sigma",k,l), xlab = "niter")
    }
    lines(1:niter, rep(var_real[, 1], niter), col = mappa[1])
    lines(1:niter, rep(var_real[, 2], niter), col = mappa[2])
    lines(1:niter, rep(var_real[, 3], niter), col = mappa[3])
    legend("topright", legend=c("Cluster 1", "Cluster 2", "Cluster 3"),
           col=c(mappa[1], mappa[2], mappa[3]), lty=1, cex=0.8)
  }
}



### convergence

z_error <- array(, dim = niter)
for (i in 1:(niter)) {
  z_error[i] <- sum(z_real != mix$z_GS[, i])
  }
print(z_error)
x11()
matplot(z_error, type='l', lty = 1, lwd = 2)


##################

for (i in 1:3) {
  print(norm(Sigma_real[, , i]))
}

norms <- matrix(0, nrow = niter, ncol = C)
for (i in 1:niter) {
  for (j in 1:c) {
    norms[i, j] <- norm(mix$Sigma_GS[, , j, i])
  }
}


#######################
open3d()
#scatter3d(x = x[z_binder = z_real, 1], y = x[z_binder != z_real, 2], z = x[z_binder != z_real, 3], point.col= mix$z_GS[z_binder != z_real, niter], surface = FALSE)
scatter3d(x = x[, 1], y = x[, 2], z = x[, 3], point.col= as.numeric(z_binder == z_real)+1, surface = FALSE)
###############################

x11()
matplot(norms, main="Markov Chain for determinant(Omega)", type = 'l', xlim = c(0, niter), lty = 1, lwd = 2)
lines(1:niter, rep(norm(Sigma_real[, , 1]), niter), col = mappa[1])
lines(1:niter, rep(norm(Sigma_real[, , 2]), niter), col = mappa[2])
lines(1:niter, rep(norm(Sigma_real[, , 3]), niter), col = mappa[3])

x11()
plot(x[,1], x[,2], main="Covariance Matrix estimate, components 1 and 2", lwd = 1, pch = 1)
mixtools::ellipse(c(0,0), Sigma_hat[1:2, 1:2, 1], col = 2, lwd = 3)
mixtools::ellipse(c(0,0), Sigma_real[1:2, 1:2, mappa[1]], col = 1, lwd = 2, lty = 2)
mixtools::ellipse(c(0,0), Sigma_hat[1:2, 1:2, 2], col = 3, lwd = 2)
mixtools::ellipse(c(0,0), Sigma_real[1:2, 1:2, mappa[2]], col = 1, lwd = 2, lty = 3)
mixtools::ellipse(c(0,0), Sigma_hat[1:2, 1:2, 3], col = 4, lwd = 2)
mixtools::ellipse(c(0,0), Sigma_real[1:2, 1:2, mappa[3]], col = 1, lwd = 2, lty = 4)
legend("topleft", 
       legend = c("Sigma_hat, cluster 1", "Sigma_real, cluster 1", 
                  "Sigma_hat, cluster 2", "Sigma_real, cluster 2",
                  "Sigma_hat, cluster 3", "Sigma_real, cluster 3"),
       col = c(2, 1, 3, 1, 4, 1), lty=c(1, 2, 1, 3, 1, 4), cex=0.5)
x11()
plot(x[,2], x[,3], main="Covariance Matrix estimate, components 2 and 3", lwd = 1, pch = 1)
mixtools::ellipse(c(0,0), Sigma_hat[2:3, 2:3, 1], col = 2, lwd = 3)
mixtools::ellipse(c(0,0), Sigma_real[2:3, 2:3, mappa[1]], col = 1, lwd = 2, lty = 2)
mixtools::ellipse(c(0,0), Sigma_hat[2:3, 2:3, 2], col = 3, lwd = 2)
mixtools::ellipse(c(0,0), Sigma_real[2:3, 2:3, mappa[2]], col = 1, lwd = 2, lty = 3)
mixtools::ellipse(c(0,0), Sigma_hat[2:3, 2:3, 3], col = 4, lwd = 2)
mixtools::ellipse(c(0,0), Sigma_real[2:3, 2:3, mappa[3]], col = 1, lwd = 2, lty = 4)
legend("topleft", 
       legend = c("Sigma_hat, cluster 1", "Sigma_real, cluster 1", 
                  "Sigma_hat, cluster 2", "Sigma_real, cluster 2",
                  "Sigma_hat, cluster 3", "Sigma_real, cluster 3"),
       col = c(2, 1, 3, 1, 4, 1), lty=c(1, 2, 1, 3, 1, 4), cex=0.5)


#############################

library('network')
x11()
plot(network(mix$DAG_GS[, , 1, 1]))

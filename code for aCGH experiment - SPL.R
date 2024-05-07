#aCGH - SPL

library("wavethresh")
library('waveslim')
library('changepoint')
library('ggplot2')
library('latex2exp')
library('coda')
library('ggtext')
library('bfw')

# ================================= Data =======================================
data(Lai2005fig4)
y <- Lai2005fig4$GBM29
data_size <- length(y)
wavelet_adapted_size <- 2^ceiling(log2(data_size )) - data_size
serie <- seq_len(data_size)/data_size

# ============================ Global variables ================================
#Gibbs sampler: Setting down the size of the chain, burn-in and lag

burn <- 1000
lags <- 5
nchain <- 1000
BB <- burn + lags*nchain
n.verbose <- 1000
verbose <- TRUE

#Wavelet transform specifications
#Family choice
family_choice <- "Coiflets" #"DaubLeAsymm", "DaubExPhase", "Coiflets"

#Number of vanishing moments
number_vm <- 5

#BayesThresh arguments
#Let's define its arguments
type_Bayesthreshold <- "soft" #"soft" or "hard"

by_level_Bayesthreshold <- FALSE # F:global threshold computed.
#V: a threshold is computed and applied separately to each scale level.

#See Wavelet thresholding via Bayesian approach (1998) for details on specifying 
#hyperparameters.
alpha_Bayesthreshold <- 0.5 
beta_Bayesthreshold <- 1

# ==============================================================================
# ============================== Wavelet Function ==============================
# ==============================================================================

#To simplify DWT and the denoising, we create a function to be used in the MCMC

wavelettransform <- function(data){
  w_transformed <- wd(data, filter.number = number_vm, family = family_choice) 
  Bayes_Thresh<- threshold(w_transformed,levels= 3:(nlevelsWT(w_transformed)-1),
                           type = type_Bayesthreshold, policy = "BayesThresh",
                           by.level = by_level_Bayesthreshold, dev = madmad,
                           boundary = FALSE, alpha = alpha_Bayesthreshold,
                           beta = beta_Bayesthreshold, C1 = NA, C2 = NA, 
                           C1.start = 100)
  thresholded_data <- wr(Bayes_Thresh)
  thresholded_data[thresholded_data<0]=0
  thresholded_data[thresholded_data>1]=1
  return(thresholded_data)
}

# ==============================================================================
# ============================== Gibbs Sampler =================================
# ==============================================================================

#Creating the matrices
mu <- matrix(nrow=BB,ncol = 2)
phi <- matrix(nrow=BB,ncol = 2)
alpha <- matrix(nrow=BB,ncol=data_size) 

#Consider the following priors to \mu_k and \tau_k^{-2}
prior_mean1 <- quantile(y, probs = 0.25)
prior_mean2 <- quantile(y, probs = 0.75)
prior_variance1 <- var(y)
prior_variance2 <- var(y)
prior_alpha_gamma <- 0.01
prior_beta_gamma <- 0.01

#Defining start values
mu[1,]<-c(prior_mean1,prior_mean2)
phi[1,]<- c(1/prior_variance1, 1/prior_variance2)
yaux <- c(y,rev(y)[1:wavelet_adapted_size])
w <- (yaux-mu[1,1])/(mu[1,2]-mu[1,1])
alphaaux <- wavelettransform(w)
alpha[1,] <- alphaaux[1:data_size]
prob0start <- (1-alpha[1,])*dnorm(y, mean = mu[1,1], 
                                  sd = 1/sqrt(phi[1,1]))
prob1start <- alpha[1,]*dnorm(y, mean = mu[1,2], 
                              sd = 1/sqrt(phi[1,2]))
probvalstart <- prob1start/(prob0start + prob1start)

pred_label <- rbinom(data_size, 1, prob = probvalstart)

# ================================== Gibbs =====================================

for(ii in 2:BB){
  #Auxiliary variables for calculating the conditional posterior 
  n1 <- sum(pred_label)
  n0 <- data_size - n1
  y1 <- y[pred_label==1]
  y0 <- y[pred_label==0]
  sy1 <- sum(y1)
  sy0 <- sum(y0)
  
# ====================== Generating mean.lower =================================
  var_mu_lower <- 1/(n0*phi[ii-1,1]+ 1/prior_variance1)
  mean_mu_lower <- var_mu_lower *(sy0*phi[ii-1,1] +
                                    prior_mean1/prior_variance1)
  mu[ii,1] <- rnorm(n = 1, mean = mean_mu_lower, sd = sqrt(var_mu_lower))
  
# ====================== Generating prec.lower =================================  
  alpha_gamma_lower <- prior_alpha_gamma + n0/2
  beta_gamma_lower <- prior_beta_gamma + sum((y0 - mu[ii,1])^2)/2
  phi[ii,1] <- rgamma(1, shape = alpha_gamma_lower, rate = beta_gamma_lower)
  
# ====================== Generating mean.higher ================================  
  var_mu_higher <- 1/(n1*phi[ii-1,2]+ 1/prior_variance2)
  mean_mu_higher <- var_mu_higher *(sy1*phi[ii-1,2] +
                                      prior_mean2/prior_variance2)
  mu[ii,2] <- rnorm(n = 1, mean = mean_mu_higher, sd = sqrt(var_mu_higher))
  
# ====================== Generating prec.higher ================================ 
  alpha_gamma_higher <- prior_alpha_gamma + n1/2
  beta_gamma_higher <- prior_beta_gamma + sum((y1 - mu[ii,2])^2)/2
  phi[ii,2] <- rgamma(1, shape = alpha_gamma_higher, rate = beta_gamma_higher)
  
  
# ============================== Constraint ====================================
  #To dealing with label switching, we stablish that the pairs (\mu_k, \phi_k) 
  #are ordered under the constraint \mu_1 < \mu_2 and \phi_1 < \phi_2
  
  if(mu[ii,2] < mu[ii,1]){
    mu[ii,]<- mu[ii,2:1]
    phi[ii,]<- phi[ii,2:1]
  }
  
# ============================ Generating W ====================================  
  yaux <- c(y,rev(y)[1:wavelet_adapted_size])
  w <- (yaux-mu[ii,1])/(mu[ii,2]-mu[ii,1])
  
  
# ========================== Generating alpha ==================================  
  alphaaux <- wavelettransform(w)
  alpha[ii,] <- alphaaux[1:data_size]
  
  
# ====================== Generating pred_label =================================
  prob0 <- (1-alpha[ii,])*dnorm(x = y, mean = mu[ii,1], 
                                sd = 1/sqrt(phi[ii,1]))
  prob1 <- alpha[ii,]*dnorm(x = y, mean = mu[ii,2], 
                            sd = 1/sqrt(phi[ii,2]))
  probval <- prob1/(prob0 + prob1)
  
  pred_label <- rbinom(n = data_size, size = 1, prob = probval)
  
  if(!ii%%n.verbose & verbose)
    print(paste0(ii, " iterations of ", BB, "."))
  
}

# =========================== Burn-in and lag ==================================
#Creating the matrices with the burn-in and lag

N <- seq(burn+1, BB, lags)

matrix_mu <- mu[N,]
matrix_phi <- phi[N,]
matrix_alpha <- alpha[N,]

#estimated values by median
estimated_mu_median <- apply(matrix_mu,2,median)
estimated_phi_median <- apply(matrix_phi,2,median)
estimated_alpha_median <- apply(matrix_alpha,2,median)

#estimated values by mean
estimated_mu_mean <- apply(matrix_mu,2,mean)
estimated_phi_mean <- apply(matrix_phi,2,mean)
estimated_alpha_mean <- apply(matrix_alpha,2,mean)

# =========================== Checking convergence =============================

#Trace mu1
plot.ts(matrix_mu[,1], main = "", xlab="", 
        ylab = expression(mu[1]))
abline(h = estimated_mu_median[1], lty=2, col = "orange", lwd = 3)
abline(h = estimated_mu_mean[1], lty=3, col = "blue", lwd = 3)

#Correlogram mu1
acf(matrix_mu[,1], main = expression(mu[1]))

#Density mu1
plot(density(matrix_mu[,1]), main = "", col="black", lwd = 2, xlab="", 
     ylab = expression(paste("Density ",mu[1])))
abline(v = estimated_mu_median[1], lty=2, col = "orange", lwd = 3)
abline(v = estimated_mu_mean[1], lty=3, col = "blue", lwd = 3)


# ==============================================================================

#Trace mu2
plot.ts(matrix_mu[,2], main = "", xlab="", 
        ylab = expression(mu[2]))
abline(h = estimated_mu_median[2], lty=2, col = "orange", lwd = 3)
abline(h = estimated_mu_mean[2], lty=3, col = "blue", lwd = 3)


#Correlogram mu2
acf(matrix_mu[,2], main = expression(mu[2]))


#Density mu2
plot(density(matrix_mu[,2]), main = "", col="black", lwd = 2, xlab="", 
     ylab = expression(paste("Density ",mu[2])))
abline(v = estimated_mu_median[2], lty=2, col = "orange", lwd = 3)
abline(v = estimated_mu_mean[2], lty=3, col = "blue", lwd = 3)

# ==============================================================================
#Trace phi1
plot.ts(matrix_phi[,1], main = "", xlab="", 
        ylab = expression(phi[1]))
abline(h = estimated_phi_median[1], lty=2, col = "orange", lwd = 3)
abline(h = estimated_phi_mean[1], lty=3, col = "blue", lwd = 3)

#Correlogram phi1
acf(matrix_phi[,1], main = expression(phi[1]))

#Density phi1
plot(density(matrix_phi[,1]), main = "", col="black", lwd = 2, xlab="", 
     ylab = expression(paste("Density ",phi[1])))
abline(v = estimated_phi_median[1], lty=2, col = "orange", lwd = 3)
abline(v = estimated_phi_mean[1], lty=3, col = "blue", lwd = 3)

# ==============================================================================
#Trace phi2
plot.ts(matrix_phi[,2], main = "", xlab="", 
        ylab = expression(phi[2]))
abline(h = estimated_phi_median[2], lty=2, col = "orange", lwd = 3)
abline(h = estimated_phi_mean[2], lty=3, col = "blue", lwd = 3)


#Correlogram phi2
acf(matrix_phi[,2], main = expression(phi[2]))


#Density phi2
plot(density(matrix_phi[,2]), main = "", col="black", lwd = 2, xlab="", 
     ylab = expression(paste("Density ",phi[2])))
abline(v = estimated_phi_median[2], lty=2, col = "orange", lwd = 3)
abline(v = estimated_phi_mean[2], lty=3, col = "blue", lwd = 3)


plot(estimated_alpha_median,type="l",col="black", lwd = 1, 
      ylab=expression(alpha[t]))
lines(seq(1,data_size,1), estimated_alpha_mean, lty=3, lwd = 3,col = "blue")

#plot(y, ylab="", pch = 1, col = "blue")



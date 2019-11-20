dmouseModel0 <- function()
{
  for ( i in 1 : N_mouse) { 
    for (k in 1:6 ) { 
      R[i,k] ~ dbin(p[i, k],N[i,k]) 
      y_rep[i,k] ~ dbin(p[i, k],N[i,k]) 
      b[i, k] ~ dnorm(0.0,tau);  
      logit(p[i, k]) <- alpha0 + b0[i] + b[i,k] 
      eta_hat[i,k] <- alpha0 
      
      chi[i,k] <- pow(R[i,k]-p[i,k]*N[i,k],2)/p[i,k]/N[i,k]
      chirep[i,k] <- pow(y_rep[i,k]-p[i,k]*N[i,k],2)/p[i,k]/N[i,k]
    }
    b0[i] ~ dnorm(0.0, tau0); 
  }
  
  # priors
  alpha0 ~ dnorm(0.0,1.0E-6)
  tau0 <- pow(sigma0, -2)
  tau <- pow(sigma, -2)
  sigma0 ~ dunif(0,100)
  sigma ~ dunif(0,100)
  
}
mouseModel_0a <- function()
{
  for ( i in 1 : N_mouse) { 
    for (k in 1:6 ) { 
      R[i,k] ~ dbin(p[i, k],N[i,k]) 
      y_rep[i,k] ~ dbin(p[i, k],N[i,k]) 
      b[i, k] ~ dnorm(0.0,tau);  
      logit(p[i, k]) <- alpha0 + alpha1 * female[i] +
        alpha2 * PLT[i] +
        alpha3 * Pound[i] +
        alpha4 * time_3[k] +
        alpha5 * time_10[k] +
        alpha6 * time_17[k] +
        alpha7 * time_24[k] +
        alpha8 * time_31[k] +
        
        b0[i] + 
        b[i,k] 
      eta_hat[i,k] <- alpha0 + alpha1 * female[i] +
        alpha2 * PLT[i] +
        alpha3 * Pound[i] +
        alpha4 * time_3[k] +
        alpha5 * time_10[k] +
        alpha6 * time_17[k] +
        alpha7 * time_24[k] +
        alpha8 * time_31[k] 
      
      chi[i,k] <- pow(R[i,k]-p[i,k]*N[i,k],2)/p[i,k]/N[i,k]
      chirep[i,k] <- pow(y_rep[i,k]-p[i,k]*N[i,k],2)/p[i,k]/N[i,k]
    }
    b0[i] ~ dnorm(0.0, tau0); 
  }
  Tobs <- sum(chi)
  Trep <- sum(chirep)
  # priors
  alpha0 ~ dnorm(0.0,1.0E-6)
  alpha1 ~ dnorm(0.0,1.0E-6)
  alpha2 ~ dnorm(0.0,1.0E-6)
  alpha3 ~ dnorm(0.0,1.0E-6)
  alpha4 ~ dnorm(0.0,1.0E-6)
  alpha5 ~ dnorm(0.0,1.0E-6)
  alpha6 ~ dnorm(0.0,1.0E-6)
  alpha7 ~ dnorm(0.0,1.0E-6)
  alpha8 ~ dnorm(0.0,1.0E-6) 
  beta1 ~ dnorm(0.0,1.0E-6) 
  beta2 ~ dnorm(0.0,1.0E-6)  
  gamma1 ~ dnorm(0.0,1.0E-6)  
  gamma2 ~ dnorm(0.0,1.0E-6)  
  gamma3 ~ dnorm(0.0,1.0E-6)  
  gamma4 ~ dnorm(0.0,1.0E-6)  
  gamma5 ~ dnorm(0.0,1.0E-6)  
  delta1 ~ dnorm(0.0,1.0E-6)  
  delta2 ~ dnorm(0.0,1.0E-6)  
  delta3 ~ dnorm(0.0,1.0E-6)  
  delta4 ~ dnorm(0.0,1.0E-6)  
  delta5 ~ dnorm(0.0,1.0E-6)  
  delta6 ~ dnorm(0.0,1.0E-6)  
  delta7 ~ dnorm(0.0,1.0E-6)  
  delta8 ~ dnorm(0.0,1.0E-6)  
  delta9 ~ dnorm(0.0,1.0E-6)  
  delta10 ~ dnorm(0.0,1.0E-6)  
  tau0 <- pow(sigma0, -2)
  tau <- pow(sigma, -2)
  sigma0 ~ dunif(0,100)
  sigma ~ dunif(0,100)
  
}
mouseModel_0b <- function()
{
  for ( i in 1 : N_mouse) { 
    for (k in 1:6 ) { 
      R[i,k] ~ dbin(p[i, k],N[i,k]) 
      
      b[i, k] ~ dnorm(0.0,tau);  
      
      logit(p[i, k]) <- alpha0 + alpha1 * female[i] +
        alpha2 * PLT[i] +
        alpha3 * Pound[i] +
        b0[i] + 
        b[i,k] 
      eta_hat[i,k] <- alpha0 + alpha1 * female[i] +
        alpha2 * PLT[i] +
        alpha3 * Pound[i] 
      
      
    }
    b0[i] ~ dnorm(0.0, tau0); 
  }
  
  # priors
  alpha0 ~ dnorm(0.0,1.0E-6)
  alpha1 ~ dnorm(0.0,1.0E-6)
  alpha2 ~ dnorm(0.0,1.0E-6)
  alpha3 ~ dnorm(0.0,1.0E-6)
  
  tau0 <- pow(sigma0, -2)
  tau <- pow(sigma, -2)
  sigma0 ~ dunif(0,100)
  sigma ~ dunif(0,100)
  
}
mouseModelNorm0 <- function()
{
  for ( i in 1 : N_mouse) { 
    for (k in 1:6 ) { 
      R[i,k] ~ dnorm(mu[i,k],tau2) 
      y_rep[i,k] ~ dnorm(mu[i,k],tau2) 
      
      mu[i, k] <- alpha0 +  b0[i]
      
      eta_hat[i,k] <- alpha0 
      chirep[i,k] <- pow(y_rep[i,k]-mu[i,k],2) 
      chi[i,k] <- pow(R[i,k]-mu[i,k],2) 
    }
    b0[i] ~ dnorm(0.0, tau0); 
  }
  Tobs <- sum(chi)
  Trep <- sum(chirep)
  # priors
  alpha0 ~ dnorm(0.0,1.0E-6)
  alpha1 ~ dnorm(0.0,1.0E-6)
  alpha2 ~ dnorm(0.0,1.0E-6)
  alpha3 ~ dnorm(0.0,1.0E-6)
  alpha4 ~ dnorm(0.0,1.0E-6)
  alpha5 ~ dnorm(0.0,1.0E-6)
  alpha6 ~ dnorm(0.0,1.0E-6)
  alpha7 ~ dnorm(0.0,1.0E-6)
  alpha8 ~ dnorm(0.0,1.0E-6) 
  beta1 ~ dnorm(0.0,1.0E-6) 
  beta2 ~ dnorm(0.0,1.0E-6)  
  gamma1 ~ dnorm(0.0,1.0E-6)  
  gamma2 ~ dnorm(0.0,1.0E-6)  
  gamma3 ~ dnorm(0.0,1.0E-6)  
  gamma4 ~ dnorm(0.0,1.0E-6)  
  gamma5 ~ dnorm(0.0,1.0E-6)  
  delta1 ~ dnorm(0.0,1.0E-6)  
  delta2 ~ dnorm(0.0,1.0E-6)  
  delta3 ~ dnorm(0.0,1.0E-6)  
  delta4 ~ dnorm(0.0,1.0E-6)  
  delta5 ~ dnorm(0.0,1.0E-6)  
  delta6 ~ dnorm(0.0,1.0E-6)  
  delta7 ~ dnorm(0.0,1.0E-6)  
  delta8 ~ dnorm(0.0,1.0E-6)  
  delta9 ~ dnorm(0.0,1.0E-6)  
  delta10 ~ dnorm(0.0,1.0E-6) 
  
  tau0 <- pow(sigma0, -2)
  sigma0 ~ dunif(0,100)
  
  tau2 <- pow(sigma2, -2)
  sigma2 ~ dunif(0,100)
  
}

mouseModelNorm <- function(){
  for ( i in 1 : 35) { 
    for (k in 1:6 ) { 
      R[i,k] ~ dnorm(mu[i,k],tau2) 
      y_rep[i,k] ~ dnorm(mu[i,k],tau2) 
      
      mu[i, k] <- alpha0 + alpha1 * female[i] +
        alpha2 * PLT[i] +
        alpha3 * Pound[i] +
        alpha4 * time_3[k] +
        alpha5 * time_10[k] +
        alpha6 * time_17[k] +
        alpha7 * time_24[k] +
        alpha8 * time_31[k] + 
        alpha9 * stroke_size[i] +
        beta1 * PLT[i] * female[i] +
        beta2 * Pound[i] * female[i] +
        gamma1 * time_3[k]*female[i] +
        gamma2 * time_10[k]*female[i] +
        gamma3 * time_17[k]*female[i] +
        gamma4 * time_24[k]*female[i] +
        gamma5 * time_31[k]*female[i] +
        delta1 * time_3[k] * PLT[i] +
        delta2 * time_10[k] * PLT[i] +
        delta3 * time_17[k] * PLT[i] +
        delta4 * time_24[k] * PLT[i] +
        delta5 * time_31[k] * PLT[i] +
        delta6 * time_3[k] * Pound[i] +
        delta7 * time_10[k] * Pound[i] +
        delta8 * time_17[k] * Pound[i] +
        delta9 * time_24[k] * Pound[i] +
        delta10 * time_31[k] * Pound[i] +
        b0[i] 
      
      eta_hat[i,k] <- alpha0 + alpha1 * female[i] +
        alpha2 * PLT[i] +
        alpha3 * Pound[i] +
        alpha4 * time_3[k] +
        alpha5 * time_10[k] +
        alpha6 * time_17[k] +
        alpha7 * time_24[k] +
        alpha9 * stroke_size[i] +
        beta1 * PLT[i] * female[i] +
        beta2 * Pound[i] * female[i] +
        gamma1 * time_3[k]*female[i] +
        gamma2 * time_10[k]*female[i] +
        gamma3 * time_17[k]*female[i] +
        gamma4 * time_24[k]*female[i] +
        gamma5 * time_31[k]*female[i] +
        delta1 * time_3[k] * PLT[i]   +
        delta2 * time_10[k] * PLT[i]  +
        delta3 * time_17[k] * PLT[i]  +
        delta4 * time_24[k] * PLT[i]  +
        delta6 * time_3[k] * Pound[i] +
        delta7 * time_10[k] * Pound[i] +
        delta8 * time_17[k] * Pound[i] +
        delta9 * time_24[k] * Pound[i]  
      
      
      chirep[i,k] <- pow(y_rep[i,k]-mu[i,k],2) 
      chi[i,k] <- pow(R[i,k]-mu[i,k],2) 
    }
    b0[i] ~ dnorm(0.0, tau0); 
  }
  Tobs <- sum(chi)
  Trep <- sum(chirep)
  # priors
  alpha0 ~ dnorm(0.0,1.0E-6)
  alpha1 ~ dnorm(0.0,1.0E-6)
  alpha2 ~ dnorm(0.0,1.0E-6)
  alpha3 ~ dnorm(0.0,1.0E-6)
  alpha4 ~ dnorm(0.0,1.0E-6)
  alpha5 ~ dnorm(0.0,1.0E-6)
  alpha6 ~ dnorm(0.0,1.0E-6)
  alpha7 ~ dnorm(0.0,1.0E-6)
  alpha8 ~ dnorm(0.0,1.0E-6) 
  alpha9 ~ dnorm(0.0,1.0E-6) 
  beta1 ~ dnorm(0.0,1.0E-6) 
  beta2 ~ dnorm(0.0,1.0E-6)  
  gamma1 ~ dnorm(0.0,1.0E-6)  
  gamma2 ~ dnorm(0.0,1.0E-6)  
  gamma3 ~ dnorm(0.0,1.0E-6)  
  gamma4 ~ dnorm(0.0,1.0E-6)  
  gamma5 ~ dnorm(0.0,1.0E-6)  
  delta1 ~ dnorm(0.0,1.0E-6)  
  delta2 ~ dnorm(0.0,1.0E-6)  
  delta3 ~ dnorm(0.0,1.0E-6)  
  delta4 ~ dnorm(0.0,1.0E-6)  
  delta5 ~ dnorm(0.0,1.0E-6)  
  delta6 ~ dnorm(0.0,1.0E-6)  
  delta7 ~ dnorm(0.0,1.0E-6)  
  delta8 ~ dnorm(0.0,1.0E-6)  
  delta9 ~ dnorm(0.0,1.0E-6)  
  delta10 ~ dnorm(0.0,1.0E-6) 
  
  tau0 <- pow(sigma0, -2)
  sigma0 ~ dunif(0,100)
  
  tau2 <- pow(sigma2, -2)
  sigma2 ~ dunif(0,100)
  
}



mouseModel <- function()
{
  for ( i in 1 : N_mouse) { 
    for (k in 1:5 ) { 
      R[i,k] ~ dbin(p[i, k],N[i,k]) 
      y_rep[i,k] ~ dbin(p[i, k],N[i,k]) 
      b[i, k] ~ dnorm(0.0,tau);  
      logit(p[i, k]) <- alpha0 + alpha1 * female[i] +
        alpha2 * PLT[i] +
        alpha3 * Pound[i] +
        alpha4 * time_3[k] +
        alpha5 * time_10[k] +
        alpha6 * time_17[k] +
        alpha7 * time_24[k] +
        beta1 * PLT[i]      * female[i] +
        beta2 * Pound[i]    * female[i] +
        gamma1 * time_3[k]  * female[i] +
        gamma2 * time_10[k] * female[i] +
        gamma3 * time_17[k] * female[i] +
        gamma4 * time_24[k] * female[i] +
        delta1 * time_3[k]  * PLT[i]    +
        delta2 * time_10[k] * PLT[i]    +
        delta3 * time_17[k] * PLT[i]    +
        delta4 * time_24[k] * PLT[i]    +
        delta6 * time_3[k]  * Pound[i]  +
        delta7 * time_10[k] * Pound[i]  +
        delta8 * time_17[k] * Pound[i]  +
        delta9 * time_24[k] * Pound[i]  +
        b0[i] + 
        b[i,k] 
      eta_hat[i,k] <- alpha0 +    alpha1 * female[i] +
        alpha2 * PLT[i] +
        alpha3 * Pound[i] +
        alpha4 * time_3[k] +
        alpha5 * time_10[k] +
        alpha6 * time_17[k] +
        alpha7 * time_24[k] + 
        beta1 * PLT[i] * female[i] +
        beta2 * Pound[i] * female[i] +
        gamma1 * time_3[k]*female[i] +
        gamma2 * time_10[k]*female[i] +
        gamma3 * time_17[k]*female[i] +
        gamma4 * time_24[k]*female[i] + 
        delta1 * time_3[k] * PLT[i] +
        delta2 * time_10[k] * PLT[i] +
        delta3 * time_17[k] * PLT[i] +
        delta4 * time_24[k] * PLT[i] +
        delta6 * time_3[k] * Pound[i] +
        delta7 * time_10[k] * Pound[i] +
        delta8 * time_17[k] * Pound[i] +
        delta9 * time_24[k] * Pound[i]  
         
      
      chi[i,k] <- pow(R[i,k]-p[i,k]*N[i,k],2)/p[i,k]/N[i,k]
      chirep[i,k] <- pow(y_rep[i,k]-p[i,k]*N[i,k],2)/p[i,k]/N[i,k]
    }
    b0[i] ~ dnorm(0.0, tau0); 
  }
  Tobs <- sum(chi)
  Trep <- sum(chirep)
  # priors
  phi    ~ dnorm(0.0,1.0E-6)
  alpha0 ~ dnorm(0.0,1.0E-6)
  alpha1 ~ dnorm(0.0,1.0E-6)
  alpha2 ~ dnorm(0.0,1.0E-6)
  alpha3 ~ dnorm(0.0,1.0E-6)
  alpha4 ~ dnorm(0.0,1.0E-6)
  alpha5 ~ dnorm(0.0,1.0E-6)
  alpha6 ~ dnorm(0.0,1.0E-6)
  alpha7 ~ dnorm(0.0,1.0E-6)
  alpha8 ~ dnorm(0.0,1.0E-6) 
  beta1  ~ dnorm(0.0,1.0E-6) 
  beta2  ~ dnorm(0.0,1.0E-6)  
  gamma1 ~ dnorm(0.0,1.0E-6)  
  gamma2 ~ dnorm(0.0,1.0E-6)  
  gamma3 ~ dnorm(0.0,1.0E-6)  
  gamma4 ~ dnorm(0.0,1.0E-6)  
  gamma5 ~ dnorm(0.0,1.0E-6)  
  delta1 ~ dnorm(0.0,1.0E-6)  
  delta2 ~ dnorm(0.0,1.0E-6)  
  delta3 ~ dnorm(0.0,1.0E-6)  
  delta4 ~ dnorm(0.0,1.0E-6)  
  delta5 ~ dnorm(0.0,1.0E-6)  
  delta6 ~ dnorm(0.0,1.0E-6)  
  delta7 ~ dnorm(0.0,1.0E-6)  
  delta8 ~ dnorm(0.0,1.0E-6)  
  delta9 ~ dnorm(0.0,1.0E-6)  
  delta10 ~ dnorm(0.0,1.0E-6)  
  tau0 <- pow(sigma0, -2)
  tau <- pow(sigma, -2)
  sigma0 ~ dunif(0,100)
  sigma ~ dunif(0,100)
  
}


# same as mouseModel, but no sex
mouseModelA <- function()
{
  for ( i in 1 : N_mouse) { 
    for (k in 1:5 ) { 
      R[i,k] ~ dbin(p[i, k],N[i,k]) 
      y_rep[i,k] ~ dbin(p[i, k],N[i,k]) 
      b[i, k] ~ dnorm(0.0,tau);  
      logit(p[i, k]) <- alpha0 + phi * stroke_size[i] + 
        alpha2 * PLT[i] +
        alpha3 * Pound[i] +
        alpha4 * time_3[k] +
        alpha5 * time_10[k] +
        alpha6 * time_17[k] +
        alpha7 * time_24[k] +
        delta1 * time_3[k]  * PLT[i]    +
        delta2 * time_10[k] * PLT[i]    +
        delta3 * time_17[k] * PLT[i]    +
        delta4 * time_24[k] * PLT[i]    +
        delta6 * time_3[k]  * Pound[i]  +
        delta7 * time_10[k] * Pound[i]  +
        delta8 * time_17[k] * Pound[i]  +
        delta9 * time_24[k] * Pound[i]  +
        b0[i] + 
        b[i,k] 
      eta_hat[i,k] <- alpha0 +  
        alpha2 * PLT[i] +
        alpha3 * Pound[i] +
        alpha4 * time_3[k] +
        alpha5 * time_10[k] +
        alpha6 * time_17[k] +
        alpha7 * time_24[k] + 
        delta1 * time_3[k] * PLT[i] +
        delta2 * time_10[k] * PLT[i] +
        delta3 * time_17[k] * PLT[i] +
        delta4 * time_24[k] * PLT[i] +
        delta6 * time_3[k] * Pound[i] +
        delta7 * time_10[k] * Pound[i] +
        delta8 * time_17[k] * Pound[i] +
        delta9 * time_24[k] * Pound[i]  
      
      
      chi[i,k] <- pow(R[i,k]-p[i,k]*N[i,k],2)/p[i,k]/N[i,k]
      chirep[i,k] <- pow(y_rep[i,k]-p[i,k]*N[i,k],2)/p[i,k]/N[i,k]
    }
    b0[i] ~ dnorm(0.0, tau0); 
  }
  Tobs <- sum(chi)
  Trep <- sum(chirep)
  # priors
  phi    ~ dnorm(0.0,1.0E-6)
  alpha0 ~ dnorm(0.0,1.0E-6)
  
  alpha2 ~ dnorm(0.0,1.0E-6)
  alpha3 ~ dnorm(0.0,1.0E-6)
  alpha4 ~ dnorm(0.0,1.0E-6)
  alpha5 ~ dnorm(0.0,1.0E-6)
  alpha6 ~ dnorm(0.0,1.0E-6)
  alpha7 ~ dnorm(0.0,1.0E-6)
  alpha8 ~ dnorm(0.0,1.0E-6) 
   
  delta1 ~ dnorm(0.0,1.0E-6)  
  delta2 ~ dnorm(0.0,1.0E-6)  
  delta3 ~ dnorm(0.0,1.0E-6)  
  delta4 ~ dnorm(0.0,1.0E-6)  
  delta5 ~ dnorm(0.0,1.0E-6)  
  delta6 ~ dnorm(0.0,1.0E-6)  
  delta7 ~ dnorm(0.0,1.0E-6)  
  delta8 ~ dnorm(0.0,1.0E-6)  
  delta9 ~ dnorm(0.0,1.0E-6)  
  delta10 ~ dnorm(0.0,1.0E-6)  
  tau0 <- pow(sigma0, -2)
  tau <- pow(sigma, -2)
  sigma0 ~ dunif(0,100)
  sigma ~ dunif(0,100)
  
}





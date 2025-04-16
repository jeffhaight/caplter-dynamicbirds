# Need to specify:
  # psi.cov: a matrix (or by-season 3D array) of occupancy covariates, including '1' for the intercept
  # rho.cov: a matrix (or by-season 3D array) of detection covariates, including '1' for the intercept
  # y[s,j,tt] : an array of detections per site/species/season, summed across survey occasions
  # J[s,tt]: a matrix of # of survey occasions per site/season
  # ncov.occ : # of occupancy covariates, including the intercept
  # ncov.det : # of detection covariates, including the intercept

# to add:
  # delta.cov: species trait covariates
  # ncov.traits

# Cauchy/Student-t distributions 'dt()' are used as a generalization of the normal distribution 
  # that does not assume known standard deviation
  # This specifies priors that are can be considered *weakly* informative
  

model{
### Priors ###
  # Priors for Occupancy (beta, random intercepts and slopes)
    for(l in 1:ncov.occ){
      beta.comm[l] ~ dt(0, 2.5, 1)  # Community-mean intercept and slopes
      tau.beta.comm[l] ~ dgamma(1, 1)    # Dispersion of community-mean parameters
      # sd.beta.comm[l] <- 1/sqrt(pow(tau.beta.comm[l])
      
      # Priors for the trait effects
      # Include this if allowing environmental responses to vary among species according to traits:
      for(m in 1:n.traits){
        delta.traits[l,m] ~ dnorm(0, 0.1) # priors for species covariate effects
      }
      
      # Priors for the phylogeny effects (similar to trait effects)
      # Include this if allowing environmental responses to vary among species according to taxonomic group:
      # tau.class.psi ~ dgamma(0.1, 0.1)
      # tau.class.phi ~ dgamma(0.1, 0.1)
      # tau.class.gamma ~ dgamma(0.1, 0.1)
      # tau.order ~ dgamma(0.1, 0.1)
      # for(k in 1:n.orders){
      #   eff.order[l,k] ~ dnorm(0, 0.1)
      # }
      
      
      # structure variation of species-specific intercept and slope parameters from community-mean
      for(j in 1:n.species){
        # If species vary unpredictably from the community mean:
          # beta.species[l,j] ~ dnorm(beta.comm[l], tau.beta.comm[l])  # Species-specific
        # If species parameters are allowed to partially vary according to traits:
          beta.species[l,j] ~ dnorm(beta.comm[l] 
                                    + inprod(delta.traits[l,1:n.traits], traits[j,1:n.traits])
                                    # + eff.order[l,order_vec[j]]
                                    , tau.beta.comm[l])
      }
      
      
    }
  
  # Priors for Detection (rho, random intercepts & slopes)
    for(l in 1:ncov.det){
      rho.comm[l] ~ dt(0, 2.5, 1)
      tau.rho.comm[l] ~ dgamma(1, 1)
      # sd.rho.comm[l] <- 1/sqrt(pow(tau.rho.comm[l])
      for(j in 1:n.species){
        rho.species[l,j] ~ dnorm(rho.comm[l], tau.rho.comm[l])
      }
    }
    
  # Priors for Autologistic term
      theta.comm ~ dt(0, 2.5, 1)        # Community-mean autologistic term
      tau.theta.comm ~ dgamma(1, 1)        # Dispersion of community-mean autologistic term
      for(j in 1:n.species){
        theta.species[j] ~ dnorm(theta.comm, tau.theta.comm)   # Species-specific autologistic term
      }
    
  # Priors for the 'Year' effect
  #mu.year ~ dnorm(0, 0.1)
  tau.year ~ dgamma(0.1, 0.1)
  for(t in 1:n.season){ #
    for(j in 1:n.species){
      year[j,t] ~ dnorm(0, tau.year)
      
    }
  }
  
### Likelihood/Regressions ###
  # State Process model
      # Use vector multiplication ('inner product') for greater efficiency
      
      for(j in 1:n.species){
        for(s in 1:n.site){
        # Regression for the 'first' surveys
          logit(psi[s,j,1]) <- inprod(
            beta.species[,j],     
            psi.cov[s,1, 1:ncov.occ]    # a vector of covariates at that site, in that season. Needs to be 3D for seasonally-varying covariates
          ) 
          + year[j,year_vec[1]]
          z[s,j,1] ~ dbern(psi[s,j,1])
        # Regression for subsequent survey 'seasons'
          for(t in 2:n.season){
            logit(psi[s,j,t]) <- inprod(
              beta.species[, j],
              psi.cov[s, t, 1:ncov.occ]     
            ) 
            + theta.species[j] * z[s,j,t-1]
            + year[j,year_vec[t]]
            z[s,j,t] ~ dbern(psi[s,j,t])
          }
        }
      }
  
  # Observation process model
  for(j in 1:n.species){
    for(s in 1:n.site){
      for(t in 1:n.season){
        logit(rho[s,j,t]) <- inprod(
          rho.species[1:ncov.det,j],
          rho.cov[s,t,1:ncov.det]
        )
        y[s,j,t] ~ dbinom(rho[s,j,t] * z[s,j,t], K[s,t])
      }
    }
  }
    
    
    ### Derived Parameters ###
    # Number of sites occupied by each species in each season
    # for(j in 1:n.species){
    #   for(t in 1:n.season){  # Loop over years
    #     n.occ[j,t] <- sum(z[,j,t])
    #   }
    # }
    
    # Alpha diversity metrics   
    # for(s in 1:n.site){
    #   for(t in 1:n.season){
    #     # Species richness (i.e., Hill Number 0)
    #     rich[s,t] <- sum(z[s,,t]) # Number of species occurring at each site, based on the occurrence matrix
    #     #rich[i,t] <- sum(psi[s,,t]) # Number of species occurring at each site, based on the occupancy probabilities
    #     
    #     # Shannon Diversity and Evenness
    #     # for this, we will essentially have to calculate Hill #1, log() that to get Shannon Diversity
    #     # see code and equations from Haight et al 2023, Turrini and Knop 2015, or Boron et al 2019 to see how that works
    #     
    #     # calculating relative "abundances"
    #     # sum.psi[s,t] <- sum(psi[s,,t])  # sum of occupancy probabilities across species
    #     # sum.psi.checked[s,t] <- ifelse(sum.psi[s,t] == 0, 1E6, sum.psi[s,t]) # avoids dividing by 0 when calculating relative psi
    #     # # 
    #     # for(j in 1:n.species){
    #     #   # relative.psi = relative occupancy: occupancy of each species divided by the across-species sum of probabilities
    #     #   relative.psi[s,j,t] <- psi[s,j,t]/sum.psi[s,t]
    #     #   log.relative.psi[s,j,t] <- ifelse(relative.psi[s,j,t]== 0,
    #     #                                     log(relative.psi[s,j,t]+1E-6),
    #     #                                     log(relative.psi[s,j,t]))
    #     # }
    #   }
    # }    
    
    
} # END


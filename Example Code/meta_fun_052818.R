###################################################################################
###### R Code for Exact Inference for the Beta-Binomial Random Effects Model ######
################# Lu Tian, Nie Li, Ying Lu, and Jessica Gronsbell #################
###################################################################################


#######################
## helper functions ##
######################

library(parallel) 
logit = function(xx){log(xx/(1-xx))}
expit = function(xx){1/(1+exp(-xx))}


###################################
## function for data generation ##
##################################

## K_org: number of studies
## N1: sample size in group 1
## N2 sample size in group 2
## alpha, beta: parameterize the beta distribution for 
##    lambda_1/(lambda_1+lambda_2) ~ B(alpha, beta)
## base_rate: baseline event rate in group 1 

gen.data <- function(K_org, N1, N2, alpha, beta, base_rate = 0.01){
  
  base_rate <- rgamma(K_org, 1.44, 1.44/base_rate)
  
  # events in group 1
  y2 <- rpois(K_org, N2*base_rate)
  
  # events in group 2
  rr <- rbeta(K_org, alpha, beta)
  y1 <- rpois(K_org, N1*base_rate*rr/(1-rr))
  
  # data for analysis
  data <- cbind(1:K_org, N1, y1, N2, y2)
  colnames(data) <- c('Study', 'Grp1_Sz', 'Grp1_Ev', 'Grp2_Sz', 'Grp2_Ev')
  
  # remove double zero studies and get new data size
  data =  data[(data[,'Grp1_Ev']+ data[,'Grp2_Ev']) != 0, ]
  K = nrow(data)
  
  # return data, # of non-DZ studies, proportion of DZ studies
  return(list(data = data, K = K, prop_DZ = (K_org - K)/K))
}


#######################################################################
## function to mimic 1-1 randomization via hypergeometric resampling ##
#######################################################################

## data: formatted as above in gen.data function
## i.e. a K x 5 matrix with the columns: 'Study', 'Grp1_Sz', 'Grp1_Ev', 'Grp2_Sz', 'Grp2_Ev'

get.data <- function(data){
  
  new.data = NULL;
  
  for(i in 1:K){
    
    ## More patients in Grp 1 ## 
    if(data[i, 'Grp1_Sz'] >= data[i, 'Grp2_Sz']){
      
      if(data[i, 'Grp1_Ev'] == 0 & data[i, 'Grp2_Ev'] != 0){
        # only one way to sample zero events in Grp1
        p.t = 1
        d.t = data[i, 'Grp1_Ev']  
        n.t = data[i, 'Grp1_Ev'] + data[i, 'Grp2_Ev']
      }else{
        # resample Grp2_Sz number of patients and record weight associated with each value of
        # number of events observed according to the hypergeometric distribution
        p.t <- dhyper(0:data[i, 'Grp1_Ev'], data[i, 'Grp1_Ev'],  
                      data[i, 'Grp1_Sz'] - data[i, 'Grp1_Ev'], data[i, 'Grp2_Sz'])
        d.t <- 0:data[i, 'Grp1_Ev'] 
        n.t <- data[i, 'Grp2_Ev'] + d.t
      }
      
      ## More patients in Grp 2 ## 
    }else{
      
      if(data[i, 'Grp1_Ev'] != 0 & data[i, 'Grp2_Ev'] == 0){
        # only one way to sample zero events in Grp2
        p.t = 1
        d.t = data[i, 'Grp1_Ev']
        n.t = data[i, 'Grp1_Ev'] + data[i, 'Grp2_Ev']
      }else{
        # resample Grp1_Sz number of patients and record weight associated with each value of
        # number of events observed according to the hypergeometric distribution
        p.t <- dhyper(0:data[i, 'Grp2_Ev'], data[i, 'Grp2_Ev'],  
                      data[i, 'Grp2_Sz'] - data[i, 'Grp2_Ev'], data[i, 'Grp1_Sz'])
        d.t <- data[i, 'Grp1_Ev']
        n.t <- 0:data[i, 'Grp2_Ev'] + d.t
      }
      
      
    }
    
    tmp.data = as.matrix(cbind(i, d.t, n.t, p.t, data[i, 'Grp1_Sz'], data[i, 'Grp2_Sz']))
    
    ## if we have zero in the smaller group we employ the following 
    ## rescaling of the probabilities as a continuity correction
    
    ind = which(tmp.data[,'n.t'] == 0)
    if(length(ind) > 0){
      dz = tmp.data[ind, 'p.t']
      tmp.data = matrix(tmp.data[-ind,], ncol = ncol(tmp.data))
      tmp.data[,4] = tmp.data[,4]/sum(tmp.data[,4])
    }
    
    
    ## return the new data ##
    new.data <- rbind(new.data, tmp.data)
    colnames(new.data) <- c('K', 'N_Ev', 'N', 'Wgt', 'N1', 'N2')
  }
  
  return(new.data)
}


#####################################
## function for moment estimation ##
####################################

## data: formatted as above in gen.data function
## i.e. a K x 5 matrix with the columns: 'Study', 'Grp1_Sz', 'Grp1_Ev', 'Grp2_Sz', 'Grp2_Ev'
## K: number of non-DZ studies
## mu.type: indicates if we have 1-1 randomization or not

MOM_est = function(data.org, K = NULL, mu.type = 'unbal'){
  
  if(is.null(K)){K = nrow(data.org)}
  
  # formats data to mimic 1-1 randomization via hypergeometric resampling if not 1-1
  if(mu.type == 'unbal'){
    data = get.data(data.org)
  }else{
    data = cbind(1:K, data.org[, 'Grp1_Ev'], data.org[, 'Grp1_Ev'] + data.org[, 'Grp2_Ev'], 1,
                 data.org[, 'Grp1_Sz'], data.org[, 'Grp2_Sz'])
    colnames(data) = c('K', 'N_Ev', 'N', 'Wgt', 'N1', 'N2')
  }
  
  # estimate of the first moment
  mu_0 = sum(data[,'Wgt']*data[,'N_Ev']/data[,'N'])/K 

  ## continuity correction for the first moment for second moment estimation
  data.var  = data[,]; 
  data.var[, 'N'] = data.var[, 'N']  + 1
    #(1/data.var[,'N1'] + 1/data.var[,'N2'])
  #+ 1
  #+ (1/data.var[,'N1'] + 1/data.var[,'N2'])
  #+ 1
    #(1/data.var[,'N1'] + 1/data.var[,'N2'])
  data.var[, 'N_Ev'] = data.var[, 'N_Ev'] + 0.5
  #+1/data.var[,'N2']
  #+ 0.5
    #data.var[, 'N_Ev'] + 1/data.var[,'N2'] 
  #+ 0.5
    #1/data.var[,'N2']
  
  # continuity corrected first moment
  mu.cc = sum(data.var[,'Wgt']*(data.var[,'N_Ev'])/(data.var[,'N']))/K
 
  ## tau2 estimator using first and second moments
  num = (sum(data.var[,'Wgt']*(data.var[,'N_Ev']/data.var[,'N'])^2)/K - 
           sum(data.var[,'Wgt']*mu.cc/data.var[,'N'])/K)
  
  denum = (sum(data.var[,'Wgt'] *(1 - 1/(data.var[,'N']))))/K 
  
  mu2 = num/denum
  
  tau2 = max(0, mu2 - mu.cc^2)
  
  ## estimated standard error of first moment estimator
  stdE = (sum(data.var[,'Wgt']* (mu.cc*(1-mu.cc)/data.var[,'N']) + data.var[,'Wgt']*(1 -1/data.var[,'N'])*tau2 )/K^2)
  
  ## return results
  return((list(muhat.cc = mu.cc, mu2 = mu2, tau2hat = tau2, 
               muhat = mu_0, v_muhat = stdE)))
}


#########################################
## function for Monte Carlo Simulation ##
#########################################

## bb is number of replications
## mom.o are the moment estimators based on the observed data
## d is (K \times 1) vector with number of events in group 1
## n is the (K \times 1) total number of events 
## N1 and N2 are the (K \times 1) vectors of study specific samples for group 1 and 2 
## a.t and b.t are the alpha and beta of interest

run_MC <- function(bb, mom.o, d, n, N1, N2, a.t, b.t, mu.type = 'unbal'){
  
  # true location and variance 
  mu.t <- a.t/(a.t + b.t)
  tau.t <- mu.t*(1-mu.t)/(a.t + b.t + 1)
  
  # number of studies
  K <- length(n)
  
  # observed test statistics
  T.stat.norm.obs <- (mom.o$muhat - mu.t)^2/(mom.o$v_muhat) # normalized
  T.stat.obs <- (mom.o$muhat - mu.t)^2 # un-normalized
  
  # generate data and get moment estimators
  results = NULL
  for(i in 1:bb){
    
    exp.theta <- rbeta(K, a.t, b.t)
    theta = -logit(exp.theta)
    pi.t <- 1/(1 + exp(theta + log(N2/N1)))
    d.t <- rbinom(K, n, pi.t)
    
    # data under H_0
    data.t = data.t = cbind(1:K, N1, apply(cbind(N1, d.t),  1, min) , N2, apply(cbind(N2, n-d.t),  1, min))
    colnames(data.t) <- c('Study', 'Grp1_Sz', 'Grp1_Ev', 'Grp2_Sz', 'Grp2_Ev')
    
    # moment estimates
    mom.t <- MOM_est(data.t, K, mu.type = mu.type)
    
    # test statistics
    T.stat.norm <- (mom.t$muhat - mu.t)^2/(mom.t$v_muhat) 
    T.stat <- (mom.t$muhat - mu.t)^2
    
    # store result
    results <- rbind(results, c(T.stat.norm, T.stat))
  }
  
  # compute the p-value based on empirical distribution
  p.val.norm <- mean(results[,1] >= T.stat.norm.obs)
  p.val <- mean(results[,2] >= T.stat.obs)
  
  # return result
  return(c(a.t, b.t, mu.t, tau.t, p.val.norm, p.val)) 
  
}


##################################################################################################
## function to compute the alpha and beta corresponding to a given mu and nu along the boundary ##
##################################################################################################

## mu is the mean of interest
## tol is the tolerance (i.e closeness to 1)

compute_ab <- function(mu, tol = 1e-6){
  
  # nu along the boundary for given mu
  v1 <- mu^2*(1-mu)/(1+mu)
  v2 <- mu*(1-mu)^2/(2-mu)
  v.t <- min(v1, v2) 
  
  # corresponding alpha and beta + some tolerance so neither equals 1 so a.t, b.t > 1
  a.t = mu*((mu*(1-mu)-v.t)/v.t)*(1+tol)
  b.t = (1-mu)*((mu*(1-mu)-v.t)/v.t)*(1+tol)
  
  # return result
  return(c(a.t, b.t))
}


############################################################################################
## functions for fast computation of the exact CI via the boundary of the parameter space ##
###########################################################################################

## bb is number of replications
## mom.o are the moment estimators based on the observed data
## d is the (K \times 1) vector with number of events in group 1
## n is the (K \times 1) total number of events 
## N1 and N2 are the (K \times 1) vectors of study specific samples for group 1 and 2 
## a.t and b.t are the alpha and beta of interest
## CI.dsl is the approximate CI for mu based on the moment estimators
## mesh size is distance between each iteration
## numb.ex is number of mesh points to compute beyond failure of rejection
## numb.tau is number of tau mesh points when computing past boundary


## computes upper bound of the CI
outer_bound <- function(bb, mom.o, d, n, N1, N2, CI.dsl, mesh.size = 0.0002, 
                        numb.ex = 20, numb.tau = 10){
  
  res.all = NULL # to store the results
  
  # settings for whether or not have 1-1 randomization in all the studies
  if(sum(N1 == N2) == length(N1)){mu.type = 'one2one'; ind = 5}else{mu.type = 'unbal'; ind = 5}
  
  #######################
  #### UPPPER BOUND ####
  ######################
  
  # Check the upper bound of the approx CI to see if you can start there
  m = min(CI.dsl[2], 1-mesh.size)
  ab.tmp <- compute_ab(m)
  a.t = ab.tmp[1]; b.t = ab.tmp[2]
  
  # run the Monte Carlo
  res.t = run_MC(bb, mom.o, d, n, N1, N2, a.t, b.t, mu.type)
  
  if(res.t[ind] >= 0.05){
    
    # if the bound is in the CI then start there, if not then just start at the moments estimator
    res.all = rbind(res.all, res.t);
    m.val.curr = m
    p.val.curr = res.t[ind]
    
  }else{
    
    p.val.curr = 1 
    m.val.curr = mom.o$muhat
    
  }
  
  # decrease mesh if necessary
  if(m.val.curr >= (1-mesh.size)){mesh.size = 0.5*mesh.size}
  
  # iterate until we reach mu_UB
  while(p.val.curr >= 0.05){
    
    m = m.val.curr + mesh.size
    if(m >= 1){break}
    
    # get alpha and beta along the boundary
    ab.tmp <- compute_ab(m)
    a.t = ab.tmp[1]; b.t = ab.tmp[2]
    
    # run the Monte Carlo
    res.t = run_MC(bb, mom.o, d, n, N1, N2, a.t, b.t, mu.type)
    
    # save the result
    res.all = rbind(res.all, res.t)
    
    # update the while loop
    m.val.curr = m
    p.val.curr = res.t[ind]
  }
  
  
  ## for a few values at and exceeding current value, iterate along all the tau
  m.all = c(m.val.curr, m.val.curr + c((1:numb.ex)*mesh.size))
  m.all = m.all[m.all < 1]; nn1 = length(m.all); nn2 = numb.tau; tol = 1e-6
  
  res.eps.upper = NULL
  for(i in 1:nn1){
    m = m.all[i];
    v1 <- m^2*(1-m)/(1+m)
    v2 <- m*(1-m)^2/(2-m)
    v <- apply(cbind(v1,v2), 1, min)
    v.t  <- seq(1e-6, v, length.out = nn2)
    for(j in 1:nn2){
      v = v.t[j]
      a.t = m*((m*(1-m)-v)/v)*(1+tol)
      b.t = (1-m)*((m*(1-m)-v)/v)*(1+tol)
      if(isTRUE(a.t <= 1 | b.t <= 1)){j = j+1} else{
        
        result.tmp = run_MC(bb, mom.o, d, n, N1, N2, a.t, b.t, mu.type)
      }
      res.eps.upper = rbind(res.eps.upper, result.tmp)
    }
  }
  
  colnames(res.eps.upper) <- c('alpha', 'beta', 'mu', 'tau2', 'pval.norm', 'pval') 
  keep.eps.upper = matrix(res.eps.upper[which(res.eps.upper[,ind] >= 0.05), ], 
                          ncol = ncol(res.eps.upper))
  colnames(keep.eps.upper) <- c('alpha', 'beta', 'mu', 'tau2', 'pval.norm', 'pval') 
  res.all = rbind(res.all, keep.eps.upper)
  
  # record whether or not the given neighborhood was large enough
  if(length(keep.eps.upper) > 0){
    tot = I(round(max(keep.eps.upper[,'mu']),6) == round(max(m.all),6))*1
  }else{
    tot = 0
  }
  
  res.all = cbind(res.all, tot, 0)
  
  # sort and return result
  colnames(res.all) <- c('alpha', 'beta', 'mu', 'tau2', 'pval.norm', 'pval',
                         'epsN_upper', 'epsN_lower') 
  res.all = res.all[order(res.all[, 'mu']), ]
  
  return(res.all)
  
}


## computes lower bound of the CI
lower_bound <- function(bb, mom.o, d, n, N1, N2, CI.dsl, mesh.size = 0.0002, 
                        numb.ex = 20, numb.tau = 10){
  
  res.all = NULL # to store the results
  
  # settings for whether or not have 1-1 randomization in all the studies
  if(sum(N1 == N2) == length(N1)){mu.type = 'one2one'; ind = 5}else{mu.type = 'unbal'; ind = 5}
  
  #######################
  #### UPPPER BOUND ####
  ######################
  
  # Check the upper bound of the approx CI to see if you can start there
  m = max(CI.dsl[1], mesh.size)
  ab.tmp <- compute_ab(m)
  a.t = ab.tmp[1]; b.t = ab.tmp[2]
  
  # run the Monte Carlo
  res.t = run_MC(bb, mom.o, d, n, N1, N2, a.t, b.t, mu.type)
  
  if(res.t[ind] >= 0.05){
    
    # if the bound is in the CI then start there, if not then just start at the moments estimator
    res.all = rbind(res.all, res.t);
    m.val.curr = m
    p.val.curr = res.t[ind]
    
  }else{
    
    p.val.curr = 1 
    m.val.curr = mom.o$muhat
    
  }
  
  # decrease mesh if necessary
  if(m.val.curr <= mesh.size){mesh.size = 0.5*mesh.size}
  
  # iterate until we reach mu_UB
  while(p.val.curr >= 0.05){
    
    m = m.val.curr - mesh.size
    if(m <= 0){break}
    
    # get alpha and beta along the boundary
    ab.tmp <- compute_ab(m)
    a.t = ab.tmp[1]; b.t = ab.tmp[2]
    
    # run the Monte Carlo
    res.t = run_MC(bb, mom.o, d, n, N1, N2, a.t, b.t, mu.type)
    
    # save the result
    res.all = rbind(res.all, res.t)
    
    # update the while loop
    m.val.curr = m
    p.val.curr = res.t[ind]
  }
  
  
  ## for a few values at and exceeding current value, iterate along all the tau
  m.all = c(m.val.curr, m.val.curr - c((1:numb.ex)*mesh.size))
  m.all = m.all[m.all > 0]; nn1 = length(m.all); nn2 = numb.tau; tol = 1e-6
  
  res.eps.lower = NULL
  for(i in 1:nn1){
    m = m.all[i];
    v1 <- m^2*(1-m)/(1+m)
    v2 <- m*(1-m)^2/(2-m)
    v <- apply(cbind(v1,v2), 1, min)
    v.t  <- seq(1e-6, v, length.out = nn2)
    for(j in 1:nn2){
      v = v.t[j]
      a.t = m*((m*(1-m)-v)/v)*(1+tol)
      b.t = (1-m)*((m*(1-m)-v)/v)*(1+tol)
      if(isTRUE(a.t <= 1 | b.t <= 1)){j = j+1} else{
        
        result.tmp = run_MC(bb, mom.o, d, n, N1, N2, a.t, b.t, mu.type)
      }
      res.eps.lower = rbind(res.eps.lower, result.tmp)
    }
  }
  
  colnames(res.eps.lower) <- c('alpha', 'beta', 'mu', 'tau2', 'pval.norm', 'pval') 
  keep.eps.lower = matrix(res.eps.lower[which(res.eps.lower[,ind] >= 0.05), ], 
                          ncol = ncol(res.eps.lower))
  colnames(keep.eps.lower) <- c('alpha', 'beta', 'mu', 'tau2', 'pval.norm', 'pval') 
  res.all = rbind(res.all, keep.eps.lower)
  
  # record whether or not the given neighborhood was large enough
  if(length(keep.eps.lower) > 0){
    tot2 = I(round(min(keep.eps.lower[,'mu']),6) == round(min(m.all),6))*1
  }else{
    tot2 = 0
  }
  res.all = cbind(res.all, 0, tot2)
  
  
  
  # sort and return result
  colnames(res.all) <- c('alpha', 'beta', 'mu', 'tau2', 'pval.norm', 'pval',
                         'epsN_upper', 'epsN_lower') 
  res.all = res.all[order(res.all[, 'mu']), ]
  
  return(res.all)
  
}
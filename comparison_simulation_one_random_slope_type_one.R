library(MASS)
library(nlme)
library(parallel)
library(lme4)
library(stats)
library(numDeriv)
library(mvtnorm)
library(Matrix)
library(optimx)
library(LearnBayes)


set.seed(12345)

no_cores <- detectCores() - 1
cl <- makeCluster(no_cores)

setwd('C:/Study/MSC_Project/comparison simulations one random slope/sample 200')

log_multivariate_t_1_reparm <- function(response,fixed_design,rand_design,parm,sample_size,repeated_measure)
{
  beta1 <- parm[1]
  beta2 <- parm[2]
  lambda0 <- parm[3]
  lambda1 <- parm[4]
  gamma12 <- parm[5]
  
  beta <- matrix(c(beta1,beta2),2,1)
  
  m <- sample_size * repeated_measure
  n <- repeated_measure
  mu <- as.vector(fixed_design %*% beta)
  q <- dim(rand_design)[2]
  
  mat_list <- lapply(c(1:sample_size),function(i){rand_design[((i-1)*n+1):(i*n),]})
  
  Z <- as.matrix(bdiag(mat_list))
  ##reparametization####################################
  lambda_mat <- c(exp(lambda0),exp(lambda1)) * diag(2)
  ######################################################
  tri_mat <- matrix(c(1,gamma12,0,1),2,2)  
  lambda_tri_mat <-  lambda_mat %*% tri_mat
  identity_mat <- diag(sample_size)
  
  W <- kronecker(identity_mat,lambda_tri_mat)
  
  cov_mat <- diag(m) + Z %*% W %*% t(W) %*% t(Z)
  dmvt(response,mu,cov_mat,df=2,log=TRUE)
}
## fixed_design m * p
## rand_design m * q
log_multivariate_t_0_reparm <- function(response,fixed_design,rand_design,parm,sample_size,repeated_measure)
{
  beta1 <- parm[1]
  beta2 <- parm[2]
  lambda0 <- parm[3]
  
  beta <- matrix(c(beta1,beta2),2,1)
  
  m <- sample_size * repeated_measure
  n <- repeated_measure
  mu <- as.vector(fixed_design %*% beta)
  q <- dim(rand_design)[2]
  
  mat_list <- lapply(c(1:sample_size),function(i){rand_design[((i-1)*n+1):(i*n),]})
  
  Z <- as.matrix(bdiag(mat_list))
  ## reparametization ####################
  lambda_mat <- matrix(exp(lambda0),1,1)
  ########################################
  tri_mat <- diag(1)
  lambda_tri_mat <-  lambda_mat %*% tri_mat
  identity_mat <- diag(sample_size)
  
  W <- kronecker(identity_mat,lambda_tri_mat)
  
  cov_mat <- diag(m) + Z %*% W %*% t(W) %*% t(Z)
  
  dmvt(x=response,mu,cov_mat,df=2,log=TRUE)
}
log_prior_1_reparm <- function(parm)
{
  beta1 <- parm[1]
  beta2 <- parm[2]
  lambda0 <- parm[3]
  lambda1 <- parm[4]
  gamma12 <- parm[5]
  ## alternative model has independent lambda and mu
  value <- log(dmvnorm(c(beta1,beta2),c(0,0),10*diag(2)) * dnorm(lambda0,mean = log(0.3),sd = sqrt(2)) * dnorm(lambda1,mean = log(0.3),sd = sqrt(2)) * dnorm(gamma12,0,1))
  value
}
log_prior_0_reparm <- function(parm)
{
  beta1 <- parm[1]
  beta2 <- parm[2]
  lambda0 <- parm[3]
  
  log(dmvnorm(c(beta1,beta2),c(0,0),10*diag(2)) * dnorm(lambda0,mean = log(0.3),sd = sqrt(2)))
}
log_posterior_1_reparm <- function(response,fixed_design,rand_design,parm,sample_size,repeated_measure)
{
  log_post<-log_multivariate_t_1_reparm(response,fixed_design,rand_design,parm,sample_size,repeated_measure) + log_prior_1_reparm(parm)
  log_post
}
log_posterior_0_reparm <- function(response,fixed_design,rand_design,parm,sample_size,repeated_measure)
{
  log_post <- log_multivariate_t_0_reparm(response,fixed_design,rand_design,parm,sample_size,repeated_measure) + log_prior_0_reparm(parm)
  log_post
}
generate_LME_sample <- function(N,n,beta_0,beta_1,rand_sigma)
{
  ## generate dataframe of N samples and n repeated measurements
  y <- matrix(0,N*n,7)
  
  ## error <- rnorm(N*n,0,sigma_square)
  
  for (i in 1:N)
  {
    #b1i and b2i
    rand_mat_i <- mvrnorm(n=1,c(0,0),rand_sigma)
    b_0_i <- rand_mat_i[1]
    b_1_i <- rand_mat_i[2]
    
    times <- seq(1,n)
    standev <- sd(times)
    times_star <- as.vector(scale(times,center = TRUE,scale = FALSE))
    times_star <- times_star/(2*standev)
    
    for(j in 1:n)
    {
      ## mutiply fixed effect coeeficient with explanatory matrix
      error_i_j <- rnorm(1,0,1)
      y[(i-1)*n+j,1] <- beta_0 + b_0_i + (beta_1+b_1_i)*times_star[j] + error_i_j
      y[(i-1)*n+j,2] <- i
      y[(i-1)*n+j,3] <- j
      y[(i-1)*n+j,4] <- 1
      y[(i-1)*n+j,5] <- times_star[j]
      y[(i-1)*n+j,6] <- 1
      y[(i-1)*n+j,7] <- times_star[j]
    }
  }
  colnames(y) <- c('Response','Id','Time','Fixed_intercept','Fixed_ScaledTime','Rand_intercept','Rand_ScaledTime')
  data.frame(y)
}
permute_LME<-function(samples,removed_effect,B,n)
{
  permute_loglik <- numeric(B)
  for(i in 1:B)
  {
    permute_samples <- samples
    for (j in 1:n)
    {
      ## random permute the index
      temp <- as.vector(sample(permute_samples[permute_samples[,'Time'] == j,'Id'],replace =FALSE))
      ## taking data out for jth repeated measure
      Time_sample <- permute_samples[permute_samples[,'Time']==j,]
      ## permute indecies i for jth repeated measurement and assign the shuffle back to the permute sample and continue
      permute_samples[permute_samples[,'Time']==j,'Response'] <- Time_sample[match(temp,Time_sample$Id),'Response']
    }
    ## adding back the fixed effect to make sure that the calculated likelihood are based on the same design matrix. 
    permute_samples[,'Response'] <- permute_samples[,'Response'] + removed_effect
    ## fit null model
    nullmod <- lmer(Response~1 + Fixed_ScaledTime + (1 | Id),data=permute_samples,REML=FALSE)
    ## fit alternative model
    altermod <- lmer(Response~1 + Fixed_ScaledTime + (1 + Fixed_ScaledTime| Id),data=permute_samples,REML =FALSE)
    ## likelihood ratio for observation
    logLik_permute <- as.numeric(-2*(logLik(nullmod) - logLik(altermod)))
    permute_loglik[i] <- logLik_permute
  }
  permute_loglik
}


#############################################################Sigma 1###########################################
func1<-function(samples)
{
  sample_size <- length(unique(samples[,'Id']))
  rep_measure <- 5
  
  B <- 200
  significance <- 0.05
  
  ############################################################ Bayesian #####################################################
  response <- samples[,1]
  fixed_design <- as.matrix(samples[,c('Fixed_intercept','Fixed_ScaledTime')])
  rand_design_0 <- as.matrix(samples[,c('Rand_intercept')])
  rand_design_1 <- as.matrix(samples[,c('Rand_intercept','Fixed_ScaledTime')])
  
  inits <- c(1,1,-2,-2,0.01)
  fit_1 <- laplace(logpost=log_posterior_1_reparm,mode = inits,response = response,fixed_design = fixed_design,rand_design = rand_design_1,sample_size=sample_size,repeated_measure=rep_measure)
  
  #fit_1 <- optim(par=inits, fn=log_posterior_1_reparm, gr = NULL, method = "L-BFGS-B",lower=c(-2,-2,-10,-10,-2),upper=c(3,3,1,1,2),response = response,fixed_design = fixed_design,rand_design = rand_design_1,sample_size=sample_size,repeated_measure=rep_measure, hessian=TRUE,control=list(fnscale=-1))
  #fit_1 <- optim(par=inits, fn=log_posterior_1_reparm, gr = NULL, method = "L-BFGS-B",lower=c(-10,-10,-10,-10,-2),upper=c(10,10,10,10,2),response = response,fixed_design = fixed_design,rand_design = rand_design_1,sample_size=sample_size,repeated_measure=rep_measure, hessian=TRUE,control=list(fnscale=-1))
  #param_mode_1 <- fit_1$par
  #h_1 <- -solve(fit_1$hessian)
  #p_1 <- length(mode)
  #int_1 <- p_1/2 * log(2 * pi) + 0.5 * log(det(h_1)) + log_posterior_1_reparm(response = response,fixed_design = fixed_design,rand_design = rand_design_1,parm=param_mode_1,sample_size=sample_size,repeated_measure=rep_measure)
  int_1 <- fit_1$int
  
  inits <- c(1,1,-2)
  fit_0 <- laplace(logpost=log_posterior_0_reparm,mode = inits,response = response,fixed_design = fixed_design,rand_design = rand_design_0,sample_size=sample_size,repeated_measure=rep_measure)
  #fit_0 <- optim(par=inits, fn=log_posterior_0_reparm, gr = NULL, method = "L-BFGS-B",lower=c(-2,-2,-10),upper=c(3,3,1),response = response,fixed_design = fixed_design,rand_design = rand_design_0,sample_size=sample_size,repeated_measure=rep_measure, hessian=TRUE,control=list(fnscale=-1))
  #fit_0 <- optim(par=inits, fn=log_posterior_0_reparm, gr = NULL, method = "L-BFGS-B",lower=c(-10,-10,-10),upper=c(10,10,10),response = response,fixed_design = fixed_design,rand_design = rand_design_0,sample_size=sample_size,repeated_measure=rep_measure, hessian=TRUE,control=list(fnscale=-1))
  #param_mode_0 <- fit_0$par
  #h_0 <- -solve(fit_0$hessian)
  #p_0 <- length(param_mode_0)
  
  #int_0 <- p_0/2 * log(2 * pi) + 0.5 * log(det(h_0)) + log_posterior_0_reparm(response = response,fixed_design = fixed_design,rand_design = rand_design_0,parm=param_mode_0,sample_size=sample_size,repeated_measure=rep_measure)
  int_0 <- fit_0$int
  
  log.bayes.factor <- int_1 - int_0
  
  compare_results <- exp(log.bayes.factor) >= 1
  
  #bayes_factor_1[itr] <- compare_results
  
  
  ############################################################ Permutation #####################################################
  samples[,'Id'] <- factor(samples[,'Id'])
  nullmod <- lmer(Response~1 + Fixed_ScaledTime + (1 | Id),data=samples,REML=FALSE)
  ## fit alternative model
  altermod <- lmer(Response~1 + Fixed_ScaledTime + (1 + Fixed_ScaledTime| Id),data=samples,REML=FALSE)
  ## likelihood ratio for observation
  LME_logLik_obs <- as.numeric(-2*(logLik(nullmod) - logLik(altermod))) 
  ## response - effect_in_NUll to get ystar
  fixed_design <- as.matrix(samples[,c('Fixed_intercept','Fixed_ScaledTime')])
  raneff <- ranef(altermod)$Id
  fixeff <- matrix(fixef(altermod),nrow=2,ncol=1)
  
  rand_effect <- sapply(c(1:sample_size),function(i){samples[samples[,'Id'] == i,'Rand_intercept'] * raneff[i,1]})
  samples[,'Response'] <- as.vector(sapply(c(1:sample_size),function(i){samples[samples[,'Id'] == i,'Response'] - samples[samples[,'Id'] == i,'Rand_intercept'] * raneff[i,1]}))
  
  # ## remove fixed effect second
  samples[,'Response'] <- samples[,'Response'] - fixed_design %*% fixeff
  ## store total removed effects to be parsed to permutation
  removed_effect <- fixed_design %*% fixeff + matrix(as.vector(rand_effect),nrow=sample_size*rep_measure,ncol=1)
  LME_loglike_permute_1 <- permute_LME(samples,removed_effect,B,rep_measure)
  permutation <- as.numeric(mean(LME_loglike_permute_1 >= LME_logLik_obs) <= significance)
  
  c(Bayes=compare_results,Permutation=permutation)
}

func2 <- function(itr,settings)
{
  sample_size <- settings[1]
  rep_measure <- settings[2]
  beta_0 <- settings[3] 
  beta_1 <- settings[4]
  sqrt_phi_22 <- settings[5]
  rand_sigma <- matrix(c(1,-0.3*sqrt_phi_22,-0.3*sqrt_phi_22,sqrt_phi_22^2),nrow=2,ncol=2)
  
  generate_LME_sample(sample_size,rep_measure,beta_0,beta_1,rand_sigma)
}

#################################### set sample size, lambda and iterations
settings1 <- c(sample_size=200,rep=5,beta0=1,beta1=2,sqrt_phi_22=0)


iterations <- rep(1:300)
###########################################################################

results1 <- NULL


samples1 <- lapply(X=iterations,func2,settings=settings1)




clusterExport(cl, "settings1")


clusterExport(cl, "results1")


clusterExport(cl, "samples1")


clusterExport(cl, "iterations")

clusterExport(cl, "generate_LME_sample")
clusterExport(cl, "permute_LME")
clusterExport(cl, "log_multivariate_t_1_reparm")
clusterExport(cl, "log_multivariate_t_0_reparm")
clusterExport(cl, "log_prior_1_reparm")
clusterExport(cl, "log_prior_0_reparm")
clusterExport(cl, "log_posterior_1_reparm")
clusterExport(cl, "log_posterior_0_reparm")
clusterExport(cl, "func1")
clusterExport(cl, "func2")


clusterEvalQ(cl, library(MASS))
clusterEvalQ(cl, library(nlme))
clusterEvalQ(cl, library(lme4))
clusterEvalQ(cl, library(parallel))
clusterEvalQ(cl, library(stats))
clusterEvalQ(cl, library(numDeriv))
clusterEvalQ(cl, library(mvtnorm))
clusterEvalQ(cl, library(Matrix))
clusterEvalQ(cl, library(optimx))
clusterEvalQ(cl, library(LearnBayes))


results1 <- parLapply(cl,samples1,func1)




stopCluster(cl)


results_mat_1 <- as.matrix(unlist(results1))
bayes_factor_1 <- mean(results_mat_1[which(rownames(results_mat_1) == "Bayes"),],na.rm=TRUE)
permutation_test_1 <- mean(results_mat_1[which(rownames(results_mat_1) == "Permutation"),])


bayes_F <- c(bayes_factor_1)
permute_T <- c(permutation_test_1)

testResults <- data.frame(bayes_F,permute_T)



write.csv(testResults,"comparison_one_slope_sample200_type_error.csv")



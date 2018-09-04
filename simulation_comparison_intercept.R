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

setwd('C:/Study/MSC_Project/comparison simulations intercept/sample 25')


log_multivariate_t_1_reparm <- function(response,parm,sample_size,repeated_measure)
{
  mu<-parm[1]
  phi <- parm[2]
  m <- sample_size * repeated_measure
  
  l_m<- t(t(rep(1,m)))
  l_n<- t(t(rep(1,repeated_measure)))
  w <- as.matrix(bdiag(replicate(sample_size,list(l_n))))
  ##reparametization####################################
  sigma_mat<- diag(m) + exp(2*phi) * w %*% t(w)
  ######################################################
  delta<- as.vector(mu*l_m)
  dmvt(response,delta,sigma_mat, 2, log=TRUE)
}


log_multivariate_t_0_reparm <- function(response,parm,sample_size,repeated_measure)
{
  mu <- parm[1]
  m <- sample_size * repeated_measure
  l_m<- t(t(rep(1,m)))
  delta<- as.vector(mu*l_m)
  dmvt(response,delta,diag(m),2, log=TRUE)
  
}


log_posterior_1_reparm <- function(response,parm,sample_size,repeated_measure)
{
  
  mu<-parm[1]
  phi<-parm[2]
  
  like_m1 <- log_multivariate_t_1_reparm(response,parm,sample_size,repeated_measure)
  prior_mu <- dnorm(mu,0,1)
  prior_lambda <- dnorm(phi,log(0.3),sqrt(2))
  prior <- log(prior_mu*prior_lambda)
  like_m1 + prior
}


log_posterior_0_reparm <- function(response,parm,sample_size,repeated_measure)
{
  mu <- parm[1]
  like_m0 <- log_multivariate_t_0_reparm(response,parm,sample_size,repeated_measure)
  prior_m0 <- log(dnorm(mu,0,1))
  like_m0 + prior_m0
}

generate_ANOVA_sample <- function(N,n,lambda)
{
  ## generate dataframe of N samples and n repeated measurements
  y <- matrix(0,N*n,4)
  
  ## error <- rnorm(N*n,0,sigma_square)
  
  for (i in 1:N)
  {
    #b1i and b2i
    rand_mat_i <- rnorm(1,0,lambda^2)
    for(j in 1:n)
    {
      ## mutiply fixed effect coeeficient with explanatory matrix
      ## 
      error_i_j <- rnorm(1,0,1)
      y[(i-1)*n+j,1] <- rand_mat_i + error_i_j
      y[(i-1)*n+j,2] <- i
      y[(i-1)*n+j,3] <- j
      y[(i-1)*n+j,4] <- 1
    }
  }
  colnames(y) <- c('Response','Id','Time','Fixed_intercept')
  data.frame(y)
}


permute<-function(samples,fixed_effect,B,n)
{
  permute_loglik <- numeric(B)
  for(i in 1:B)
  {
    permute_samples <- samples
    for (j in 1:n)
    {
      ## random permute the index
      temp <- as.vector(sample(permute_samples[permute_samples[,'Time'] == j,'Id'],replace =FALSE))
      ## taking response out for jth repeated measure
      Time_sample <- permute_samples[permute_samples[,'Time']==j,]
      ## permute indecies i for jth repeated measurement and assign the shuffle back to the permute sample and continue
      permute_samples[permute_samples[,'Time']==j,'Response'] <- Time_sample[match(temp,Time_sample$Id),'Response']
    }
    ## adding back the fixed effect to make sure that the calculated likelihood are based on the same design matrix. 
    permute_samples[,'Response'] <- permute_samples[,'Response'] + fixed_effect
    ## fit null model
    nullmod <- lm(Response~1,data=permute_samples)
    ## fit alternative model
    altermod <- lmer(Response~1 + (1 | Id), data=permute_samples,REML=FALSE)
    ## likelihood ratio for observation
    logLik_permute <- as.numeric(-2*(logLik(nullmod) - logLik(altermod)))
    permute_loglik[i] <- logLik_permute
  }
  permute_loglik
}


#####################################################Sample 100###########################################
func1<-function(samples)
{
  sample_size <- length(unique(samples[,'Id']))
  rep_measure <- 5
  B <- 200
  significance <- 0.05
  
  
  ############################################################ Bayesian #####################################################
  response <- samples[,1]
  
  inits <- c(1,-2)
  fit_1 <- laplace(logpost=log_posterior_1_reparm,mode = inits,response = response,sample_size=sample_size,repeated_measure=rep_measure)
  int_1 <- fit_1$int
  
  inits <- c(1)
  fit_0 <- laplace(logpost=log_posterior_0_reparm,mode = inits,response = response,sample_size=sample_size,repeated_measure=rep_measure)
  int_0 <- fit_0$int
  
  log.bayes.factor <- int_1 - int_0
  
  compare_results <- exp(log.bayes.factor) >= 1
  
  ############################################################ Permutation #####################################################
  samples[,'Id'] <- factor(samples[,'Id'])
  nullmod <- lm(Response~1,data=samples)
  ## fit alternative model
  altermod <- lmer(Response~1 + (1 | Id), data=samples,REML=FALSE)
  ## likelihood ratio for observation
  logLik_obs <- as.numeric(-2*(logLik(nullmod) - logLik(altermod))) 
  ## remove fixed effect to get ysar
  
  fixeff <- matrix(nullmod$coefficients,nrow=1,ncol=1)
  ## take the intercept out
  fixed_design <- as.matrix(samples[,'Fixed_intercept']) %*% fixeff
  samples[,'Response'] <- samples[,'Response'] - fixed_design
  
  loglike_permute_1 <- permute(samples,fixed_design,B,rep_measure)
  
  permutation <- as.numeric(mean(loglike_permute_1 >= logLik_obs) <= significance)
  c(Bayes=compare_results,Permutation=permutation)
}


func2 <- function(itr,settings)
{
  sample_size <- settings[1]
  rep_measure <- settings[2]
  lambda <- settings[3]
  generate_ANOVA_sample(sample_size,rep_measure,lambda)
}
#################################### set sample size, lambda and iterations
settings1 <- c(sample_size=25,rep=5,lambda=0.15)
settings2 <- c(sample_size=25,rep=5,lambda=0.3)
settings3 <- c(sample_size=25,rep=5,lambda=0.45)
settings4 <- c(sample_size=25,rep=5,lambda=0.6)

iterations <- rep(1:100)
###########################################################################
results1 <- NULL
results2 <- NULL
results3 <- NULL
results4 <- NULL

samples1 <- lapply(X=iterations,func2,settings=settings1)
samples2 <- lapply(X=iterations,func2,settings=settings2)
samples3 <- lapply(X=iterations,func2,settings=settings3)
samples4 <- lapply(X=iterations,func2,settings=settings4)

clusterExport(cl, "results1")
clusterExport(cl, "results2")
clusterExport(cl, "results3")
clusterExport(cl, "results4")
clusterExport(cl, "samples1")
clusterExport(cl, "samples2")
clusterExport(cl, "samples3")
clusterExport(cl, "samples4")
clusterExport(cl, "settings1")
clusterExport(cl, "settings2")
clusterExport(cl, "settings3")
clusterExport(cl, "settings4")
clusterExport(cl, "iterations")
clusterExport(cl, "generate_ANOVA_sample")
clusterExport(cl, "log_multivariate_t_1_reparm")
clusterExport(cl, "log_multivariate_t_0_reparm")
clusterExport(cl, "log_posterior_1_reparm")
clusterExport(cl, "log_posterior_0_reparm")
clusterExport(cl, "permute")
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


results1 <- parLapply(cl,X=samples1,func1)
results2 <- parLapply(cl,X=samples2,func1)
results3 <- parLapply(cl,X=samples3,func1)
results4 <- parLapply(cl,X=samples4,func1)

stopCluster(cl)


results_mat_1 <- as.matrix(unlist(results1))
bayes_factor_1 <- mean(results_mat_1[which(rownames(results_mat_1) == "Bayes"),],na.rm=TRUE)
permutation_test_1 <- mean(results_mat_1[which(rownames(results_mat_1) == "Permutation"),])

results_mat_2 <- as.matrix(unlist(results2))
bayes_factor_2 <- mean(results_mat_2[which(rownames(results_mat_2) == "Bayes"),],na.rm=TRUE)
permutation_test_2 <- mean(results_mat_2[which(rownames(results_mat_2) == "Permutation"),])

results_mat_3 <- as.matrix(unlist(results3))
bayes_factor_3 <- mean(results_mat_3[which(rownames(results_mat_3) == "Bayes"),],na.rm=TRUE)
permutation_test_3 <- mean(results_mat_3[which(rownames(results_mat_3) == "Permutation"),])

results_mat_4 <- as.matrix(unlist(results4))
bayes_factor_4 <- mean(results_mat_4[which(rownames(results_mat_4) == "Bayes"),],na.rm=TRUE)
permutation_test_4 <- mean(results_mat_4[which(rownames(results_mat_4) == "Permutation"),])

bayes_F <- c(bayes_factor_1,bayes_factor_2,bayes_factor_3,bayes_factor_4)
permute_T <- c(permutation_test_1,permutation_test_2,permutation_test_3,permutation_test_4)

testResults <- data.frame(bayes_F,permute_T)
write.csv(testResults,"comparison_intercept_sample25.csv")






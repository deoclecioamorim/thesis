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
library(ggplot2)
library(data.table)

set.seed(12345)

no_cores <- detectCores() - 1
cl <- makeCluster(no_cores)

setwd('C:/Study/MSC_Project/simulation lme section 3_3')

# settings1 <- c(sample_size=20,rep=5,beta0=1,beta1=2,sqrt_phi_22=0.15)
# 
# 
# 
# plot_samples <- func2(0,settings1)
# 
# 
# 
# 
#  plot_samples[,'Id'] <- factor(plot_samples[,'Id'])
#  sample_table <- as.data.table(plot_samples)
#  ggplot(sample_table, aes(x = Time, y = Response, color = Id)) + 
#    geom_point() + 
#    stat_smooth(method = "lm", se = FALSE)



generate_LME_sample <- function(N,n,beta_0,beta_1,rand_sigma)
{
  ## generate dataframe of N samples and n repeated measurements
  y <- matrix(0,N*n,7)
  
  
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

func2<-function(itr,settings)
{
  sample_size <- settings[1]
  rep_measure <- settings[2]
  beta_0 <- settings[3] 
  beta_1 <- settings[4]
  rand_mat <- settings[5] + 1
  rand_sigma<- list(matrix(c(0.05,0,0,0),nrow=2,ncol=2),matrix(c(0.05,0.01,0.01,0.05),nrow=2,ncol=2),matrix(c(0.07,0.01,0.01,0.07),nrow=2,ncol=2),matrix(c(0.1,0.05,0.05,0.1),nrow=2,ncol=2),matrix(c(0.5,0.08,0.08,0.5),nrow=2,ncol=2))
  generate_LME_sample(sample_size,rep_measure,beta_0,beta_1,rand_sigma[[rand_mat]])
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


func1<-function(samples)
{
  sample_size <- length(unique(samples[,'Id']))
  rep_measure <- 5
  
  B <- 200
  significance <- 0.05
  
  ############################################################ Permutation #####################################################
  samples[,'Id'] <- factor(samples[,'Id'])
  nullmod <- lmer(Response~1 + Fixed_ScaledTime + (1 | Id),data=samples,REML=FALSE)
  ## fit alternative model
  altermod <- lmer(Response~1 + Fixed_ScaledTime + (1 + Fixed_ScaledTime| Id),data=samples,REML=FALSE)
  ## likelihood ratio for observation
  LME_logLik_obs <- as.numeric(-2*(logLik(nullmod) - logLik(altermod))) 
  ## response - effect_in_NUll to get ystar
  fixed_design <- as.matrix(samples[,c('Fixed_intercept','Fixed_ScaledTime')])
  raneff <- ranef(nullmod)$Id
  fixeff <- matrix(fixef(nullmod),nrow=2,ncol=1)
  
  rand_effect <- sapply(c(1:sample_size),function(i){samples[samples[,'Id'] == i,'Rand_intercept'] * raneff[i,1]})
  samples[,'Response'] <- as.vector(sapply(c(1:sample_size),function(i){samples[samples[,'Id'] == i,'Response'] - samples[samples[,'Id'] == i,'Rand_intercept'] * raneff[i,1]}))
  
  # ## remove fixed effect second
  samples[,'Response'] <- samples[,'Response'] - fixed_design %*% fixeff
  ## store total removed effects to be parsed to permutation
  removed_effect <- fixed_design %*% fixeff + matrix(as.vector(rand_effect),nrow=sample_size*rep_measure,ncol=1)
  LME_loglike_permute_1 <- permute_LME(samples,removed_effect,B,rep_measure)
  permutation <- as.numeric(mean(LME_loglike_permute_1 >= LME_logLik_obs) <= significance)
  
  c(Permutation=permutation)
}




#################################### set sample size, lambda and iterations
iterations <- rep(1:1000)
settings0 <- c(sample_size=50,rep=5,beta0=1,beta1=2,rand_mat=0)
samples0 <- lapply(X=iterations,func2,settings=settings0)


settings1 <- c(sample_size=50,rep=5,beta0=1,beta1=2,rand_mat=1)
settings2 <- c(sample_size=50,rep=5,beta0=1,beta1=2,rand_mat=2)
settings3 <- c(sample_size=50,rep=5,beta0=1,beta1=2,rand_mat=3)
settings4 <- c(sample_size=50,rep=5,beta0=1,beta1=2,rand_mat=4)

iterations <- rep(1:100)
###########################################################################

results0 <- NULL
results1 <- NULL
results2 <- NULL
results3 <- NULL
results4 <- NULL

samples1 <- lapply(X=iterations,func2,settings=settings1)
samples2 <- lapply(X=iterations,func2,settings=settings2)
samples3 <- lapply(X=iterations,func2,settings=settings3)
samples4 <- lapply(X=iterations,func2,settings=settings4)


clusterExport(cl, "settings0")
clusterExport(cl, "settings1")
clusterExport(cl, "settings2")
clusterExport(cl, "settings3")
clusterExport(cl, "settings4")

clusterExport(cl, "results0")
clusterExport(cl, "results1")
clusterExport(cl, "results2")
clusterExport(cl, "results3")
clusterExport(cl, "results4")

clusterExport(cl, "samples0")
clusterExport(cl, "samples1")
clusterExport(cl, "samples2")
clusterExport(cl, "samples3")
clusterExport(cl, "samples4")

clusterExport(cl, "iterations")

clusterExport(cl, "generate_LME_sample")
clusterExport(cl, "permute_LME")
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

results0 <- parLapply(cl,samples0,func1)
results1 <- parLapply(cl,samples1,func1)
results2 <- parLapply(cl,samples2,func1)
results3 <- parLapply(cl,samples3,func1)
results4 <- parLapply(cl,samples4,func1)



stopCluster(cl)

results_mat_0 <- as.matrix(unlist(results0))
permutation_test_0 <- mean(results_mat_0[which(rownames(results_mat_0) == "Permutation"),])

results_mat_1 <- as.matrix(unlist(results1))
permutation_test_1 <- mean(results_mat_1[which(rownames(results_mat_1) == "Permutation"),])

results_mat_2 <- as.matrix(unlist(results2))
permutation_test_2 <- mean(results_mat_2[which(rownames(results_mat_2) == "Permutation"),])

results_mat_3 <- as.matrix(unlist(results3))
permutation_test_3 <- mean(results_mat_3[which(rownames(results_mat_3) == "Permutation"),])

results_mat_4 <- as.matrix(unlist(results4))
permutation_test_4 <- mean(results_mat_4[which(rownames(results_mat_4) == "Permutation"),])

permute_T <- c(permutation_test_0,permutation_test_1,permutation_test_2,permutation_test_3,permutation_test_4)

testResults <- data.frame(permute_T)



write.csv(testResults,"simulation_one_slope_sample50_rep_5.csv")




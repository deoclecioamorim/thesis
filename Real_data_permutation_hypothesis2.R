library(ggplot2)
library(data.table)
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

setwd('C:/Study/MSC_Project/Real Data')

set.seed(12345)

plasma <- read.csv("plasma.csv")


plasma_control <- plasma[plasma[,'group'] == 'a',]
plasma_non_hyper <- plasma[plasma[,'group'] == 'b',]
plasma_hyper <- plasma[plasma[,'group'] == 'c',]

#
# settings0 <- c(sample_size=15,rep=5,beta0=1,beta1=2,rand_mat=0)
# plot_samples <- func2(0,settings0)
# 

# plasma_non_hyper[,'subject'] <- factor(plasma_non_hyper[,'subject'])
# sample_table <- as.data.table(plasma_control)
#    ggplot(sample_table, aes(x = time, y = plasma, color = subject)) + 
#    geom_point() + 
#    stat_smooth(method = "lm", se = FALSE)
#  



plasma_control_length <- dim(plasma_control)[1]



group_A_plasma <- plasma_control[,'plasma']
group_A <- rep(1,plasma_control_length)
group_B<- rep(0,plasma_control_length)
group_C<- rep(0,plasma_control_length)
Id <- plasma_control[,'subject']
time_A <- plasma_control[,'time']
time_B<- rep(0,plasma_control_length)
time_c<- rep(0,plasma_control_length)
time_squared_A <- time_A * time_A
time_squared_B<- rep(0,plasma_control_length)
time_squared_c<- rep(0,plasma_control_length)
rep_measure <- length(unique(plasma_control[,'subject']))
Time <- rep(c(1,2,3,4,5,6,7,8),rep_measure)
rand_intercept_control <- rep(1,plasma_control_length)
rand_time_squared_control <- plasma_control[,'time'] * plasma_control[,'time']
plasma_A <- data.frame(plasma=group_A_plasma,group_A=group_A,group_B=group_B,group_C=group_C,Id=Id,Time=Time,time_A=time_A,time_B=time_B,time_C=time_c,time_squared_A=time_squared_A,time_squared_B=time_squared_B,time_squared_C=time_squared_c,rand_intercept = rand_intercept_control,rand_time=plasma_control[,'time'],rand_time_squared=rand_time_squared_control)




plasma_non_hyper_length <- dim(plasma_non_hyper)[1]

group_B_plasma <- plasma_non_hyper[,'plasma']
group_A <- rep(0,plasma_non_hyper_length)
group_B <- rep(1,plasma_non_hyper_length)
group_C <- rep(0,plasma_non_hyper_length)
Id <- plasma_non_hyper[,'subject']
time_A <- rep(0,plasma_non_hyper_length)
time_B<- plasma_non_hyper[,'time']
time_c<- rep(0,plasma_non_hyper_length)
time_squared_A <- rep(0,plasma_non_hyper_length)
time_squared_B<- time_B * time_B
time_squared_c<- rep(0,plasma_non_hyper_length)
rep_measure <- length(unique(plasma_non_hyper[,'subject']))
Time <- rep(c(1,2,3,4,5,6,7,8),rep_measure)
rand_intercept_non_hyper <- rep(1,plasma_non_hyper_length)
rand_time_squared_non_hyper <- plasma_non_hyper[,'time'] * plasma_non_hyper[,'time']
plasma_B <- data.frame(plasma=group_B_plasma,group_A=group_A,group_B=group_B,group_C=group_C,Id=Id,Time=Time,time_A=time_A,time_B=time_B,time_C=time_c,time_squared_A=time_squared_A,time_squared_B=time_squared_B,time_squared_C=time_squared_c,rand_intercept = rand_intercept_non_hyper,rand_time=plasma_non_hyper[,'time'],rand_time_squared=rand_time_squared_non_hyper)





plasma_hyper_length <- dim(plasma_hyper)[1]

group_C_plasma <- plasma_hyper[,'plasma']
group_A <- rep(0,plasma_hyper_length)
group_B <- rep(0,plasma_hyper_length)
group_C <- rep(1,plasma_hyper_length)
Id <- plasma_hyper[,'subject']
time_A <- rep(0,plasma_hyper_length)
time_B <- rep(0,plasma_hyper_length)
time_c <- plasma_hyper[,'time']
time_squared_A <- rep(0,plasma_hyper_length)
time_squared_B <- rep(0,plasma_hyper_length)
time_squared_c <- time_c * time_c
rep_measure <- length(unique(plasma_hyper[,'subject']))
Time <- rep(c(1,2,3,4,5,6,7,8),rep_measure)
rand_intercept_hyper <- rep(1,plasma_hyper_length)
rand_time_squared_hyper <- plasma_hyper[,'time'] * plasma_hyper[,'time']
plasma_C <- data.frame(plasma=group_C_plasma,group_A=group_A,group_B=group_B,group_C=group_C,Id=Id,Time=Time,time_A=time_A,time_B=time_B,time_C=time_c,time_squared_A=time_squared_A,time_squared_B=time_squared_B,time_squared_C=time_squared_c,rand_intercept = rand_intercept_hyper,rand_time=plasma_hyper[,'time'],rand_time_squared=rand_time_squared_hyper)

plasma_design <- rbind(plasma_A,plasma_B,plasma_C)


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
      permute_samples[permute_samples[,'Time']==j,'plasma'] <- Time_sample[match(temp,Time_sample$Id),'plasma']
    }
    ## adding back the fixed effect to make sure that the calculated likelihood are based on the same design matrix. 
    permute_samples[,'plasma'] <- permute_samples[,'plasma'] + fixed_effect
    ## fit null model
    nullmod_hypo_2  <- lmer(plasma ~ group_A + group_B + group_C + time_A + time_B + time_C + time_squared_A + time_squared_B + time_squared_C - 1 + (1 | Id),data=permute_samples,REML=FALSE)
    ## fit alternative model
    altermod_hypo_2 <- lmer(plasma ~ group_A + group_B + group_C + time_A + time_B + time_C + time_squared_A + time_squared_B + time_squared_C - 1 + (1 + rand_time + rand_time_squared| Id),data=permute_samples,REML=FALSE)
    ## likelihood ratio for observation
    logLik_permute <- as.numeric(-2*(logLik(nullmod_hypo_2) - logLik(altermod_hypo_2)))
    permute_loglik[i] <- logLik_permute
  }
  permute_loglik
}



func1<-function(plasma_design)
{
  sample_size <- length(unique(plasma_design[,'Id']))
  rep_measure <- 8
  B <- 1000
  significance <- 0.05
  
  ############################################################ Permutation #####################################################
  plasma_design[,'Id'] <- factor(plasma_design[,'Id'])
  
  nullmod_hypo_2  <- lmer(plasma ~ group_A + group_B + group_C + time_A + time_B + time_C + time_squared_A + time_squared_B + time_squared_C - 1 + (1 | Id),data=plasma_design,REML=FALSE)
  altermod_hypo_2 <- lmer(plasma ~ group_A + group_B + group_C + time_A + time_B + time_C + time_squared_A + time_squared_B + time_squared_C - 1 + (1 + rand_time + rand_time_squared| Id),data=plasma_design,REML=FALSE)
  LME_logLik_obs <- as.numeric(-2*(logLik(nullmod_hypo_2) - logLik(altermod_hypo_2)))
  
  fixed_design <- as.matrix(plasma_design[,c('group_A','group_B','group_C','time_A','time_B','time_C','time_squared_A','time_squared_B','time_squared_C')])
  fixeff <- matrix(fixef(altermod_hypo_2),nrow=9,ncol=1)
  fixed_effects <- fixed_design %*% fixeff
  
  #remove rand effects
  raneff <- ranef(altermod_hypo_2)$Id
  
  rand_effect <- sapply(c(1:sample_size),function(i){plasma_design[plasma_design[,'Id'] == i,'rand_intercept'] * raneff[i,1]})
  
  plasma_design[,'plasma'] <- as.vector(sapply(c(1:sample_size),function(i){plasma_design[plasma_design[,'Id'] == i,'plasma'] - plasma_design[plasma_design[,'Id'] == i,'rand_intercept'] * raneff[i,1]}))
  #remove fixed effects
  plasma_design[,'plasma'] <- plasma_design[,'plasma'] - fixed_effects
  
  
  removed_effect <- fixed_effects + matrix(as.vector(rand_effect),nrow=sample_size*rep_measure,ncol=1)
  
  
  loglike_permute_1 <- permute(plasma_design,removed_effect,B,rep_measure)
  permutation <- as.numeric(mean(loglike_permute_1 >= LME_logLik_obs) <= significance)
  
  list(loglik_obs=LME_logLik_obs,permute_logLik = loglike_permute_1,p_value=mean(loglike_permute_1 >= LME_logLik_obs),Permutation=permutation)
}


func1(plasma_design)

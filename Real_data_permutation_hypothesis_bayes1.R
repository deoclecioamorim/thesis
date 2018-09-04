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
#library(LearnBayes)

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


subjects <- length(unique(plasma_control[,'subject']))

times <- c(0,0.5,1,1.5,2,3,4,5)
times_squared <- times^2
times_standev <- sd(times)
times_star <- as.vector(scale(times,center = TRUE,scale = FALSE))
times_star <- times_star/(2*times_standev)
times_squared_standev <- sd(times_squared)
times_squared_star <- as.vector(scale(times_squared,center = TRUE,scale = FALSE))
times_squared_star <- times_squared_star/(2*times_squared_standev)


rand_time_scaled <- as.numeric(sapply(c(1:subjects),function(i){times_star}))
rand_time_squared_scaled <- as.numeric(sapply(c(1:subjects),function(i){times_squared_star}))

plasma_A <- data.frame(plasma=group_A_plasma,group_A=group_A,group_B=group_B,group_C=group_C,Id=Id,time_A_scaled=rand_time_scaled,time_B_scaled=rep(0,plasma_control_length),time_C_scaled=rep(0,plasma_control_length),time_squared_A_scaled=rand_time_squared_scaled,time_squared_B_scaled=rep(0,plasma_control_length),time_squared_C_scaled=rep(0,plasma_control_length),rand_intercept = rep(1,plasma_control_length),rand_time_scaled=rand_time_scaled,rand_time_squared_scaled=rand_time_squared_scaled)




plasma_non_hyper_length <- dim(plasma_non_hyper)[1]

group_A_plasma <- plasma_non_hyper[,'plasma']
group_A <- rep(0,plasma_non_hyper_length)
group_B<- rep(1,plasma_non_hyper_length)
group_C<- rep(0,plasma_non_hyper_length)
Id <- plasma_non_hyper[,'subject']


subjects <- length(unique(plasma_non_hyper[,'subject']))

rand_time_scaled <- as.numeric(sapply(c(1:subjects),function(i){times_star}))
rand_time_squared_scaled <- as.numeric(sapply(c(1:subjects),function(i){times_squared_star}))
plasma_B <- data.frame(plasma=group_A_plasma,group_A=group_A,group_B=group_B,group_C=group_C,Id=Id,time_A_scaled=rep(0,plasma_non_hyper_length),time_B_scaled=rand_time_scaled,time_C_scaled=rep(0,plasma_non_hyper_length),time_squared_A_scaled=rep(0,plasma_non_hyper_length),time_squared_B_scaled=rand_time_squared_scaled,time_squared_C_scaled=rep(0,plasma_non_hyper_length),rand_intercept = rep(1,plasma_non_hyper_length),rand_time_scaled=rand_time_scaled,rand_time_squared_scaled=rand_time_squared_scaled)





plasma_hyper_length <- dim(plasma_hyper)[1]

group_A_plasma <- plasma_hyper[,'plasma']
group_A <- rep(0,plasma_hyper_length)
group_B<- rep(0,plasma_hyper_length)
group_C<- rep(1,plasma_hyper_length)
Id <- plasma_hyper[,'subject']


subjects <- length(unique(plasma_hyper[,'subject']))

rand_time_scaled <- as.numeric(sapply(c(1:subjects),function(i){times_star}))
rand_time_squared_scaled <- as.numeric(sapply(c(1:subjects),function(i){times_squared_star}))
plasma_C <- data.frame(plasma=group_A_plasma,group_A=group_A,group_B=group_B,group_C=group_C,Id=Id,time_A_scaled=rep(0,plasma_hyper_length),time_B_scaled=rep(0,plasma_hyper_length),time_C_scaled=rand_time_scaled,time_squared_A_scaled=rep(0,plasma_hyper_length),time_squared_B_scaled=rep(0,plasma_hyper_length),time_squared_C_scaled=rand_time_squared_scaled,rand_intercept = rep(1,plasma_hyper_length),rand_time_scaled=rand_time_scaled,rand_time_squared_scaled=rand_time_squared_scaled)

plasma_design <- rbind(plasma_A,plasma_B,plasma_C)


fixed_design <- as.matrix(plasma_design[,c('group_A','group_B','group_C','time_A_scaled','time_B_scaled','time_C_scaled','time_squared_A_scaled','time_squared_B_scaled','time_squared_C_scaled')])

rand_design <- as.matrix(plasma_design[,c('rand_intercept','rand_time_scaled','rand_time_squared_scaled')])

response <- plasma_design[,'plasma']

log_multivariate_t_1_reparm <- function(response,fixed_design,rand_design,parm,sample_size,repeated_measure)
{
  beta1 <- parm[1]
  beta2 <- parm[2]
  beta3 <- parm[3]
  beta4 <- parm[4]
  beta5 <- parm[5]
  beta6 <- parm[6]
  beta7 <- parm[7]
  beta8 <- parm[8]
  beta9 <- parm[9]
  phi0 <- parm[10]
  phi1 <- parm[11]
  phi2 <- parm[12]
  gamma12 <- parm[13]
  gamma13 <- parm[14]
  gamma23 <- parm[15]
  
  
  beta <- matrix(c(beta1,beta4,beta7,beta2,beta5,beta8,beta3,beta6,beta9),9,1)
  
  m <- sample_size * repeated_measure
  n <- repeated_measure
  mu <- as.vector(fixed_design %*% beta)
  q <- dim(rand_design)[2]
  
  mat_list <- lapply(c(1:sample_size),function(i){rand_design[((i-1)*n+1):(i*n),]})
  
  Z <- as.matrix(bdiag(mat_list))
  ##reparametization####################################
  lambda_mat <- c(exp(phi0),exp(phi1),exp(phi2)) * diag(3)
  ######################################################
  tri_mat <- matrix(c(1,gamma12,gamma13,0,1,gamma23,0,0,1),3,3)  
  lambda_tri_mat <-  lambda_mat %*% tri_mat
  identity_mat <- diag(sample_size)
  
  W <- kronecker(identity_mat,lambda_tri_mat)
  
  cov_mat <- diag(m) + Z %*% W %*% t(W) %*% t(Z)
  dmvt(response,mu,cov_mat,df=2,log=TRUE)
}

log_multivariate_t_0_reparm <- function(response,fixed_design,parm,sample_size,repeated_measure)
{
  beta1 <- parm[1]
  beta2 <- parm[2]
  beta3 <- parm[3]
  beta4 <- parm[4]
  beta5 <- parm[5]
  beta6 <- parm[6]
  beta7 <- parm[7]
  beta8 <- parm[8]
  beta9 <- parm[9]
  m <- sample_size * repeated_measure
  beta <- matrix(c(beta1,beta4,beta7,beta2,beta5,beta8,beta3,beta6,beta9),9,1)
  delta <- as.vector(fixed_design %*% beta)
  dmvt(response,delta,diag(m),df=2, log=TRUE)
}

log_prior_1_reparm <- function(parm)
{
  beta1 <- parm[1]
  beta2 <- parm[2]
  beta3 <- parm[3]
  beta4 <- parm[4]
  beta5 <- parm[5]
  beta6 <- parm[6]
  beta7 <- parm[7]
  beta8 <- parm[8]
  beta9 <- parm[9]
  phi0 <- parm[10]
  phi1 <- parm[11]
  phi2 <- parm[12]
  gamma12 <- parm[13]
  gamma13 <- parm[14]
  gamma23 <- parm[15]
  ## alternative model has independent lambda and mu
  value <- log(dmvnorm(c(beta1,beta4,beta7,beta2,beta5,beta8,beta3,beta6,beta9),c(0,0,0,0,0,0,0,0,0),10*diag(9)) * dmvnorm(c(phi0,phi1,phi2),c(log(0.3),log(0.3),log(0.3)),2*diag(3)) * dnorm(gamma12,0,1) * dnorm(gamma13,0,1) * dnorm(gamma23,0,1))
  value
}

log_prior_0_reparm <- function(parm)
{
  beta1 <- parm[1]
  beta2 <- parm[2]
  beta3 <- parm[3]
  beta4 <- parm[4]
  beta5 <- parm[5]
  beta6 <- parm[6]
  beta7 <- parm[7]
  beta8 <- parm[8]
  beta9 <- parm[9]
  
  log(dmvnorm(c(beta1,beta4,beta7,beta2,beta5,beta8,beta3,beta6,beta9),c(0,0,0,0,0,0,0,0,0),10*diag(9)))
}

log_posterior_1_reparm <- function(response,fixed_design,rand_design,parm,sample_size,repeated_measure)
{
  log_post<-log_multivariate_t_1_reparm(response,fixed_design,rand_design,parm,sample_size,repeated_measure) + log_prior_1_reparm(parm)
  log_post
}
log_posterior_0_reparm <- function(response,fixed_design,parm,sample_size,repeated_measure)
{
  log_post <- log_multivariate_t_0_reparm(response,fixed_design,parm,sample_size,repeated_measure) + log_prior_0_reparm(parm)
  log_post
}


laplace <- function (logpost, mode, ...) {
  
  options(warn=-1)
  
  fit=optim(mode, logpost, gr = NULL, ..., hessian=TRUE,
            
            control=list(fnscale=-1,maxit = 20000))
  
  options(warn=0)
  
  mode=fit$par
  
  h=-solve(fit$hessian)
  
  p=length(mode)
  
  int = p/2 * log(2 * pi) + 0.5 * log(det(h)) +
    
    logpost(mode, ...)
  
  output = list(mode = mode, var = h, int = int, 
                
                converge=fit$convergence,
                n_inputs=nargs() - 2)
  
  output
}





sample_size <- 33

rep_measure <- 8

inits <- c(4,-0.72,0.16,4.3,-0.86,0.16,4.7,-0.94,0.16,0.3487837,-0.4881762,-2.446011,-0.4,-0.97,0.35)
fit_1 <- laplace(logpost=log_posterior_1_reparm,mode = inits,response = response,fixed_design = fixed_design,rand_design = rand_design,sample_size=sample_size,repeated_measure=rep_measure)
int_1 <- fit_1$int

inits <- c(4,-0.72,0.16,4.3,-0.86,0.16,4.7,-0.94,0.16)
fit_0 <- laplace(logpost=log_posterior_0_reparm,mode = inits,response = response,fixed_design = fixed_design,sample_size=sample_size,repeated_measure=rep_measure)

int_0 <- fit_0$int

log.bayes.factor <- int_1 - int_0

compare_results <- exp(log.bayes.factor) >= 1

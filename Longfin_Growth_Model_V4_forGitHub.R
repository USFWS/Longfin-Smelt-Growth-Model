#-----------------------------------------------------------------------------
# Author: Bryan Matthias
# Email: bryan_matthias@fws.gov 
# Date: 9/24/2021
# Version: 4
# File: 1 of 3
# Purpose: This R script is designed to estimate the temporal patterns in 
# longfin smelt growth using a modified Lester et al. (2004) growth model
# that separates growth into a pre-reproductive (linear) and reproductive
# (follows von Bertalanffy model) phases. We used a mixed-effects framework
# to predict the effects of both observed (fixed effects) and unobserved 
# (random effects) environmental conditions on growth patterns. The methods
# are described in detail in the longfin smelt technical note V1.0, submitted
# for review with the Longfin Smelt Technical Team. This code uses template
# model builder (TMB) to run the model and I have been having difficulty
# with TMB on certain computers. I have a compiled version of the code
# that I can share so you can get around this step. File 2 is the TMB code
# Longfin_TMB_FracAge_V3.cpp (also .dll and .o in the folder, but these
# files are created by TMB). Note, this code is dynamic and predictor
# variables can be added/changed as needed without needing to recompile
# the TMB code. As it stands, V3 of the TMB model is complete and future
# versions of this analysis will focus on adding/changing the covariates
# used to predict growth patterns (e.g., abundance, prey, additional
# interactions, etc.) and not changes in the underlying growth model. File 3
# is source code for creating filled contour plots. 
#
# Difference with Version 3:
# - This version is intended to add more predictor variables (abundance
#   and prey) and more interaction terms
# - Also rearrange code so that functions are defined at the beginning
#
# Other notes: 
# I don't understand what has happened, but I was having difficulty with Rtools. Traditionally
# Rtools gets installed in the correct location, but mine was not. Therefore, I needed to put
# the following two lines of code:
# Sys.setenv(PATH = paste("C:/RBuildTools/3.5/bin", Sys.getenv("PATH"), sep=";"))
# Sys.setenv(BINPREF = "C:/RBuildTools/3.5/mingw_$(WIN)/bin/")
# into a .Rprofile file located in C:/Users/.../Documents. I had to create this file
# by creating a new txt file and renaming it to ".Rprofile" and make sure you do not include
# ".txt" at the end of the file name (i.e., it is not a text file). This error will be evident
# when compiling an "xxx.cpp" TMB model and there is an error stating the the file cannot be
# compiled. If you don't get this error, ignore above instructions. 
# note: for R 4.0.x, need rtools40 installed and below in document
# Sys.setenv(PATH = paste("C:/rtools40", Sys.getenv("PATH"), sep=";"))
# Sys.setenv(BINPREF = "C:/rtools40/mingw_$(WIN)/bin/")
#-----------------------------------------------------------------------------




#------------------------------------------------------------------------------------
#
# Step 1: Setup and define functions
#
#------------------------------------------------------------------------------------

rm(list=ls(all=T))

require(TMB)
#set.seed(1)

# load source code for creating multi-panel contour plots
source('./filled.contour3.r')

# define functions

# basic Lester growth model to generate starting values for the TMB model
lester_fun<-function(theta,inpt)
{
  h<-exp(theta[1])
  g<-exp(theta[2])
  t1<-theta[3]
  amat<-exp(theta[4])
  sig<-exp(theta[5])
  
  linf<-3*h/g
  k<-log(1+g/3)
  
  ek<- -k*(inpt$FracAgeMonth-amat)					# 
  phase2<-ifelse(inpt$FracAgeMonth>amat,1,0)		# If age>T, mature=1, else 0
  
  # length up to min(age,transition)
  lpred<-h*					# Linear growth rate (depending on age)
    (pmin(inpt$FracAgeMonth,amat)-t1)*	# determines amount of time fish grew according to h[X]
    exp(ek*phase2)+				# exp(-kT) if they are mature, else it equals 1
    linf*(1-exp(ek))*phase2	# growth after maturity if they are mature
  
  nll<- sum(-dnorm(inpt$Length,lpred,sig,log=TRUE),na.rm=TRUE)
  return(nll)
}

# plotting function 
lester_fun_p<-function(theta,inpt)
{
  h<-exp(theta[1])
  g<-exp(theta[2])
  t1<-theta[3]
  amat<-exp(theta[4])
  sig<-exp(theta[5])
  
  age<-seq(0,max(dat$FracAgeMonth),length=50)
  
  linf<-3*h/g
  k<-log(1+g/3)
  
  ek<- -k*(age-amat)					# 
  phase2<-ifelse(age>amat,1,0)		# If age>T, mature=1, else 0
  
  # length up to min(age,transition)
  lpred<-h*					# Linear growth rate (depending on age)
    (pmin(age,amat)-t1)*	# determines amount of time fish grew according to h[X]
    exp(ek*phase2)+				# exp(-kT) if they are mature, else it equals 1
    linf*(1-exp(ek))*phase2	# growth after maturity if they are mature
  
  plot(inpt$Length~inpt$FracAgeMonth,pch=19,xlim=c(0,max(dat$FracAgeMonth)),ylim=c(0,max(dat$Length)))
  lines(lpred~age,col='steelblue',lwd=3)
}

# function generates predicted length-at-age matrix
lpred_fun<-function(sdr,map)
{
  lng<-summary(sdr)[which(rownames(summary(sdr))=='lng'),1]
  if(length(lng)!=12)
  {
    gm<-exp(lng)[map$lng]
  }else{
    gm<-exp(lng)
  }
  
  if(sum(is.na(map$b_ha))==tmbdat$n_beta)
  {
    b_ha<-rep(0,tmbdat$n_beta)
  }else{
    b_ha<-summary(sdr)[which(rownames(summary(sdr))=='b_ha'),1]
    b_ha[is.na(b_ha)]<-0
  }
  if(sum(is.na(map$b_hj))==tmbdat$n_beta)
  {
    b_hj<-rep(0,tmbdat$n_beta)
  }else{
    b_hj<-summary(sdr)[which(rownames(summary(sdr))=='b_hj'),1]
    b_hj[is.na(b_hj)]<-0
  }
  if(sum(is.na(map$alpha))==tmbdat$n_alpha)
  {
    alpha<-rep(0,tmbdat$n_alpha)
  }else{
    alpha<-summary(sdr)[which(rownames(summary(sdr))=='alpha'),1]
    alpha[is.na(alpha)]<-0
  }
  
  lo<-summary(sdr)[which(rownames(summary(sdr))=='lo'),1][map$lo]
  
  hc<-summary(sdr)[which(rownames(summary(sdr))=='hc'),1][map$hc]
  if(sum(is.na(map$hy))==n_cohorts)
  {
    hy<-rep(0,length=n_years)
  }else{
    hy<-summary(sdr)[which(rownames(summary(sdr))=='hy'),1][map$hy]
  }
  hy[which(is.na(hy))]<-0
  a50<-exp(summary(sdr)[which(rownames(summary(sdr))=='lna50'),1][map$lna50])
  tcrit<-summary(sdr)[which(rownames(summary(sdr))=='tcrit'),1]
  
  
  len_pred<-array(NA,c(n_cohorts,amax+1))	# array for predicted length
  ltran<-NA			# array for length-at-transition
  
  for(crt in 1:n_cohorts)
  {
    lo_a<-0
    for(j in 1:tmbdat$n_alpha)
    {
      lo_a<-lo_a+alpha[j]*tmbdat$env_lo_dat[crt,j]
    }
    len_pred[crt,1]<-lo[crt]+lo_a							# set age-0
    
    for(a in 2:(amax+1))		# loop over ages
    {
      # a is index
      # age_temp is age
      age_temp=a-1
      
      h<-exp(hc[crt])
      g<-gm[monvec[a-1]+1]
      k<-log(1+g/3)
      #linf<-3*h/g		# asymptotic length at age infinity
      linf<-3*h/exp(mean(lng))
      
      h_a<-b_ha[1]*(tmbdat$env_dat[tmbdat$env_imat[crt,a-1]+1,1]-tcrit[1])^2   # effect of temperature^2
      h_j<-(b_ha[1]+b_hj[1])*(tmbdat$env_dat[tmbdat$env_imat[crt,a-1]+1,1]-tcrit[2])^2   # effect of temperature^2
      for(j in 2:tmbdat$n_beta)
      {  
        h_a<-h_a+b_ha[j]*tmbdat$env_dat[tmbdat$env_imat[crt,a-1]+1,j]#; // other effects
        h_j<-h_j+(b_ha[j]+b_hj[j])*tmbdat$env_dat[tmbdat$env_imat[crt,a-1]+1,j]#; // other effects
      }
      
      juv<-h*exp(h_j)				# juvenile growth
      maturing<-h*exp(h_j)*(1-(age_temp-a50[crt]))+	# juvenile growth component
        ((linf-(len_pred[crt,a-1]+h*exp(h_j)*(1-(age_temp-a50[crt]))))*	# juvenile growth component
            (1-exp(-(k)*(age_temp-a50[crt]))))*					# adult growth component
        exp(h_a)
      mature<-(linf-len_pred[crt,a-1])*
        (1-exp(-k))*
        exp(h_a)        # adult growth 
      
      pmature<-1/(1+exp(-((age_temp-1)-a50[crt])/0.1))	# probability of being mature
      pmaturing<-1/(1+exp(-(age_temp-a50[crt])/0.1))	# probability of maturing in a given year
      
      inc<-juv*(1-pmaturing)+				# add all growth parts up
        maturing*(pmaturing-pmature)+	
        mature*pmature
      
      len_pred[crt,a]<-len_pred[crt,a-1]+inc*exp(hy[yrmat[crt,a-1]+1])			# length-at-age
      
      if((age_temp-1)==floor(a50[crt]))
      {	
        tme<-1-(age_temp-a50[crt])		# time spent as juvenile
        ltran[crt]<-len_pred[crt,a-1]+inc*tme*exp(hy[yrmat[crt,a-1]+1])		# length-at-transition
      }
    }
  }
  return(list(len_pred,ltran))
}

# function generates predicted length-at-age matrix with variation to generate CIs
lpred_var_fun<-function(sdr,map)
{
  lng<-rnorm(1,summary(sdr)[which(rownames(summary(sdr))=='lng'),1],
    summary(sdr)[which(rownames(summary(sdr))=='lng'),2])
  if(length(lng)!=12)
  {
    gm<-exp(lng)[map$lng]
  }else{
    gm<-exp(lng)
  }
  
  if(sum(is.na(map$b_ha))==tmbdat$n_beta)
  {
    b_ha<-rep(0,tmbdat$n_beta)
  }else{
    b_ha<-rnorm(tmbdat$n_beta,summary(sdr)[which(rownames(summary(sdr))=='b_ha'),1],
      summary(sdr)[which(rownames(summary(sdr))=='b_ha'),2])
    b_ha[is.na(b_ha)]<-0
  }
  if(sum(is.na(map$b_hj))==tmbdat$n_beta)
  {
    b_hj<-rep(0,tmbdat$n_beta)
  }else{
    b_hj<-rnorm(tmbdat$n_beta,summary(sdr)[which(rownames(summary(sdr))=='b_hj'),1],
      summary(sdr)[which(rownames(summary(sdr))=='b_hj'),2])
    b_hj[is.na(b_hj)]<-0
  }
  if(sum(is.na(map$alpha))==tmbdat$n_alpha)
  {
    alpha<-rep(0,tmbdat$n_alpha)
  }else{
    alpha<-rnorm(tmbdat$n_alpha,summary(sdr)[which(rownames(summary(sdr))=='alpha'),1],
      summary(sdr)[which(rownames(summary(sdr))=='alpha'),2])
    alpha[is.na(alpha)]<-0
  }
  
  lo<-rnorm(tmbdat$ncrt,summary(sdr)[which(rownames(summary(sdr))=='lo'),1][map$lo],
    summary(sdr)[which(rownames(summary(sdr))=='lo'),2][map$lo])
  
  hc<-rnorm(tmbdat$ncrt,summary(sdr)[which(rownames(summary(sdr))=='hc'),1][map$hc],
    summary(sdr)[which(rownames(summary(sdr))=='hc'),2][map$hc])
  if(sum(is.na(map$hy))==n_cohorts)
  {
    hy<-rep(0,length=n_years)
  }else{
    hy<-rnorm(max(tmbdat$yrmat+1),summary(sdr)[which(rownames(summary(sdr))=='hy'),1][map$hy],
      summary(sdr)[which(rownames(summary(sdr))=='hy'),2][map$hy])
  }
  hy[which(is.na(hy))]<-0
  a50<-exp(rnorm(tmbdat$ncrt,summary(sdr)[which(rownames(summary(sdr))=='lna50'),1][map$lna50],
    summary(sdr)[which(rownames(summary(sdr))=='lna50'),2][map$lna50]))
  tcrit<-rnorm(2,summary(sdr)[which(rownames(summary(sdr))=='tcrit'),1],
    summary(sdr)[which(rownames(summary(sdr))=='tcrit'),2])
  
  len_pred<-array(NA,c(n_cohorts,amax+1))	# array for predicted length
  ltran<-NA			# array for length-at-transition
  
  for(crt in 1:n_cohorts)
  {
    lo_a<-0
    for(j in 1:tmbdat$n_alpha)
    {
      lo_a<-lo_a+alpha[j]*tmbdat$env_lo_dat[crt,j]
    }
    len_pred[crt,1]<-lo[crt]+lo_a							# set age-0
    
    for(a in 2:(amax+1))		# loop over ages
    {
      # a is index
      # age_temp is age
      age_temp=a-1
      
      h<-exp(hc[crt])
      g<-gm[monvec[a-1]+1]
      k<-log(1+g/3)
      linf<-3*h/exp(mean(lng))
      
      h_a<-b_ha[1]*(tmbdat$env_dat[tmbdat$env_imat[crt,a-1]+1,1]-tcrit[1])^2   # effect of temperature^2
      h_j<-(b_ha[1]+b_hj[1])*(tmbdat$env_dat[tmbdat$env_imat[crt,a-1]+1,1]-tcrit[2])^2   # effect of temperature^2
      
      
      for(j in 2:tmbdat$n_beta)
      {  
        h_a<-h_a+b_ha[j]*tmbdat$env_dat[tmbdat$env_imat[crt,a-1]+1,j]#; // other effects
        h_j<-h_j+(b_ha[j]+b_hj[j])*tmbdat$env_dat[tmbdat$env_imat[crt,a-1]+1,j]#; // other effects
      }
      
      juv<-h*exp(h_j)				# juvenile growth
      maturing<-h*exp(h_j)*(1-(age_temp-a50[crt]))+	# juvenile growth component
        ((linf-(len_pred[crt,a-1]+h*exp(h_j)*(1-(age_temp-a50[crt]))))*	# juvenile growth component
            (1-exp(-(k)*(age_temp-a50[crt]))))*					# adult growth component
        exp(h_a)
      mature<-(linf-len_pred[crt,a-1])*
        (1-exp(-k))*
        exp(h_a)        # adult growth 
      
      pmature<-1/(1+exp(-((age_temp-1)-a50[crt])/0.1))	# probability of being mature
      pmaturing<-1/(1+exp(-(age_temp-a50[crt])/0.1))	# probability of maturing in a given year
      
      inc<-juv*(1-pmaturing)+				# add all growth parts up
        maturing*(pmaturing-pmature)+	
        mature*pmature
      
      len_pred[crt,a]<-len_pred[crt,a-1]+inc*exp(hy[yrmat[crt,a-1]+1])			# length-at-age
      
      if((age_temp-1)==floor(a50[crt]))
      {	
        tme<-1-(age_temp-a50[crt])		# time spent as juvenile
        ltran[crt]<-len_pred[crt,a-1]+inc*tme*exp(hy[yrmat[crt,a-1]+1])		# length-at-transition
      }
    }
  }
  return(list(len_pred,ltran))
}

# function generates predicted growth increments 
l_inc_var_fun<-function(sdr,map)
{
  if(sum(is.na(map$b_ha))==tmbdat$n_beta)
  {
    b_ha<-rep(0,tmbdat$n_beta)
  }else{
    b_ha<-rnorm(tmbdat$n_beta,summary(sdr)[which(rownames(summary(sdr))=='b_ha'),1],
      summary(sdr)[which(rownames(summary(sdr))=='b_ha'),2])
    b_ha[is.na(b_ha)]<-0
  }
  if(sum(is.na(map$b_hj))==tmbdat$n_beta)
  {
    b_hj<-rep(0,tmbdat$n_beta)
  }else{
    b_hj<-rnorm(tmbdat$n_beta,summary(sdr)[which(rownames(summary(sdr))=='b_hj'),1],
      summary(sdr)[which(rownames(summary(sdr))=='b_hj'),2])
    b_hj[is.na(b_hj)]<-0
  }
  
  hc<-rnorm(tmbdat$ncrt,summary(sdr)[which(rownames(summary(sdr))=='hc'),1][map$hc],
    summary(sdr)[which(rownames(summary(sdr))=='hc'),2][map$hc])
  if(sum(is.na(map$hy))==n_cohorts)
  {
    hy<-rep(0,length=n_years)
  }else{
    hy<-rnorm(max(tmbdat$yrmat+1),summary(sdr)[which(rownames(summary(sdr))=='hy'),1][map$hy],
      summary(sdr)[which(rownames(summary(sdr))=='hy'),2][map$hy])
  }
  hy[which(is.na(hy))]<-0
  tcrit<-rnorm(2,summary(sdr)[which(rownames(summary(sdr))=='tcrit'),1],
    summary(sdr)[which(rownames(summary(sdr))=='tcrit'),2])
  
  g_inc<-array(NA,c(n_cohorts,12,2))
  
  for(crt in 1:n_cohorts)
  {
    for(m in 1:12)		# loop over ages
    {
      h<-exp(hc[crt])
      
      h_a<-b_ha[1]*(tmbdat$env_dat[tmbdat$env_imat[crt,m]+1,1]-tcrit[1])^2   # effect of temperature^2
      h_j<-(b_ha[1]+b_hj[1])*(tmbdat$env_dat[tmbdat$env_imat[crt,m]+1,1]-tcrit[2])^2   # effect of temperature^2
      for(j in 2:tmbdat$n_beta)
      {  
        h_a<-h_a+b_ha[j]*tmbdat$env_dat[tmbdat$env_imat[crt,m]+1,j]#; // other effects
        h_j<-h_j+(b_ha[j]+b_hj[j])*tmbdat$env_dat[tmbdat$env_imat[crt,m]+1,j]#; // other effects
      }
      
      juv<-h*exp(h_j)				# juvenile growth
      g_inc[yrmat[crt,m]+1,m,1]<-h*exp(h_j)*exp(hy[yrmat[crt,m]+1])
      g_inc[yrmat[crt,m]+1,m,2]<-h*exp(h_a)*exp(hy[yrmat[crt,m]+1])
    }
  }
  return(g_inc)
}

# function for predicting the length-at-age for an individual observation
lpred_i_fun<-function(sdr,map,tmbdat,lpred,i)
{
  lng<-summary(sdr)[which(rownames(summary(sdr))=='lng'),1]
  if(length(lng)!=12)
  {
    gm<-exp(lng)[map$lng]
  }else{
    gm<-exp(lng)
  }
  
  if(sum(is.na(map$b_ha))==tmbdat$n_beta)
  {
    b_ha<-rep(0,tmbdat$n_beta)
  }else{
    b_ha<-summary(sdr)[which(rownames(summary(sdr))=='b_ha'),1]
    b_ha[is.na(b_ha)]<-0
  }
  if(sum(is.na(map$b_hj))==tmbdat$n_beta)
  {
    b_hj<-rep(0,tmbdat$n_beta)
  }else{
    b_hj<-summary(sdr)[which(rownames(summary(sdr))=='b_hj'),1]
    b_hj[is.na(b_hj)]<-0
  }
  if(sum(is.na(map$alpha))==tmbdat$n_alpha)
  {
    alpha<-rep(0,tmbdat$n_alpha)
  }else{
    alpha<-summary(sdr)[which(rownames(summary(sdr))=='alpha'),1]
    alpha[is.na(alpha)]<-0
  }
  
  lo<-summary(sdr)[which(rownames(summary(sdr))=='lo'),1][map$lo]
  
  hc<-summary(sdr)[which(rownames(summary(sdr))=='hc'),1][map$hc]
  if(sum(is.na(map$hy))==n_cohorts)
  {
    hy<-rep(0,length=n_years)
  }else{
    hy<-summary(sdr)[which(rownames(summary(sdr))=='hy'),1][map$hy]
  }
  hy[which(is.na(hy))]<-0
  a50<-exp(summary(sdr)[which(rownames(summary(sdr))=='lna50'),1][map$lna50])
  tcrit<-summary(sdr)[which(rownames(summary(sdr))=='tcrit'),1]
  
  age_temp=tmbdat$wage[i]+tmbdat$fage[i]
  
  h<-exp(hc[tmbdat$cohort[i]+1])
  g<-gm[monvec[tmbdat$wage[i]+1]+1]
  k<-log(1+g/3)
  #linf<-3*h/g		# asymptotic length at age infinity
  linf<-3*h/exp(mean(lng))
  
  
  h_a<-b_ha[1]*(tmbdat$env_dat[tmbdat$env_imat[tmbdat$cohort[i]+1,tmbdat$wage[i]+1]+1,1]-tcrit[1])^2   # effect of temperature^2
  h_j<-(b_ha[1]+b_hj[1])*(tmbdat$env_dat[tmbdat$env_imat[tmbdat$cohort[i]+1,tmbdat$wage[i]+1]+1,1]-tcrit[2])^2   # effect of temperature^2
  for(j in 2:tmbdat$n_beta)
  {  
    h_a<-h_a+b_ha[j]*tmbdat$env_dat[tmbdat$env_imat[tmbdat$cohort[i]+1,tmbdat$wage[i]+1]+1,j]#; // other effects
    h_j<-h_j+(b_ha[j]+b_hj[j])*tmbdat$env_dat[tmbdat$env_imat[tmbdat$cohort[i]+1,tmbdat$wage[i]+1]+1,j]#; // other effects
  }
  
  juv<-h*tmbdat$fage[i]*exp(h_j)				# juvenile growth
  maturing<-h*tmbdat$fage[i]*exp(h_j)*(1-(age_temp-a50[tmbdat$cohort[i]+1]))+	# juvenile growth component
    ((linf-(lpred[tmbdat$cohort[i]+1,tmbdat$wage[i]+1]+h*tmbdat$fage[i]*exp(h_j)*(1-(age_temp-a50[tmbdat$cohort[i]+1]))))*	# juvenile growth component
        (1-exp(-(k*tmbdat$fage[i])*(age_temp-a50[tmbdat$cohort[i]+1]))))*					# adult growth component
    exp(h_a)
  mature<-(linf-lpred[tmbdat$cohort[i]+1,tmbdat$wage[i]+1])*
    (1-exp(-k*tmbdat$fage[i]))*
    exp(h_a)        # adult growth 
  
  pmature<-1/(1+exp(-((age_temp-1)-a50[tmbdat$cohort[i]])/0.1))	# probability of being mature
  pmaturing<-1/(1+exp(-(age_temp-a50[tmbdat$cohort[i]])/0.1))	# probability of maturing in a given year
  
  inc<-juv*(1-pmaturing)+				# add all growth parts up
    maturing*(pmaturing-pmature)+	
    mature*pmature
  
  len_pred<-lpred[tmbdat$cohort[i]+1,tmbdat$wage[i]+1]+inc*exp(hy[tmbdat$year[i]+1])			# length-at-age
  
  return(len_pred)
}

# function to fill in missing values from abundance index
fill_N_fun<-function(lnN)
{
  # loop over each cohort present in the dataset (i=cohort) to estimate
  # the slope over time for an approximate estimate of Z for each
  # cohort i. Then take the cohort mean between SFBS midwater trawl
  # and otter trawl index (cohrot mean), followed by the mean across
  # cohorts for a global mean to fill in missing data
  slp_SFBSmwt<-NULL
  slp_SFBSot<-NULL
  for(i in 1:(length(lnN[,1])-2))
  {
    temp<-c(lnN$SFBSmwt0[i],lnN$SFBSmwt1[i+1],lnN$SFBSmwt2[i+2])
    if(sum(is.na(temp))==0)
    {
      slp_SFBSmwt[i]<-lm(temp~c(0:2))$coef[2]
    }
    temp<-c(lnN$SFBSot0[i],lnN$SFBSot1[i+1],lnN$SFBSot2[i+2])
    if(sum(is.na(temp))==0)
    {
      slp_SFBSot[i]<-lm(temp~c(0:2))$coef[2]
    }
  }
  mu_slp<-mean(rowMeans(cbind(slp_SFBSmwt,slp_SFBSot),na.rm=TRUE),na.rm=TRUE)
  
  # fill in cohorts with all missing data (1994, 2018 cohorts in midwater trawl,
  # and 2016 from otter trawl)
  lm_mwtxfmwt<-lm(lnN$SFBSmwt0~lnN$FMWT,na.action=na.exclude)
  lm_otxfmwt<-lm(lnN$SFBSot0~lnN$FMWT)
  lm_mwtxot<-lm(lnN$SFBSmwt0~lnN$SFBSot0)
  lm_otxmwt<-lm(lnN$SFBSot0~lnN$SFBSmwt0)
  
  pred_mwtxfmwt<-lm_mwtxfmwt$coef[1]+lm_mwtxfmwt$coef[2]*lnN$FMWT
  pred_otxfmwt<-lm_otxfmwt$coef[1]+lm_otxfmwt$coef[2]*lnN$FMWT
  pred_mwtxot<-lm_mwtxot$coef[1]+lm_mwtxot$coef[2]*lnN$SFBSot0
  pred_otxmwt<-lm_otxmwt$coef[1]+lm_otxmwt$coef[2]*lnN$SFBSmwt0
  
  take<-which(lnN$Year<=1979|lnN$Year==1994|lnN$Year==2016|lnN$Year==2018)
  lnN$SFBSmwt0[take]<-rowMeans(cbind(pred_mwtxfmwt,pred_mwtxot),na.rm=TRUE)[take]
  take<-which(lnN$Year<=1979|lnN$Year==2016)
  lnN$SFBSot0[take]<-rowMeans(cbind(pred_otxfmwt,pred_otxmwt),na.rm=TRUE)[take]
  
  # manually fill in missing values where there is some information available
  take<-which(lnN$Year==1974)
  lnN$SFBSmwt0[take]<-mean(lnN$SFBSmwt0,na.rm=TRUE)
  lnN$SFBSot0[take]<-mean(lnN$SFBSot0,na.rm=TRUE)
  
  # filles in 2019-2020 using only observed data from post 2000
  take<-which(lnN$Year>=2019)
  lnN$SFBSmwt0[take]<-mean((lnN$SFBSmwt0[which(lnN$Year>=2000)])[-c(17,19:21)],na.rm=TRUE)
  lnN$SFBSot0[take]<-mean((lnN$SFBSot0[which(lnN$Year>=2000)])[-c(17,20:21)],na.rm=TRUE)
  
  lnN$SFBSmwt1[1]<-mean(lnN$SFBSmwt1,na.rm=TRUE)
  lnN$SFBSmwt2[c(1:2)]<-mean(lnN$SFBSmwt2,na.rm=TRUE)
  
  lnN$SFBSot1[1]<-mean(lnN$SFBSot1,na.rm=TRUE)
  lnN$SFBSot2[c(1:2)]<-mean(lnN$SFBSot2,na.rm=TRUE)
  
  take<-which(lnN$Year<=1977)
  lnN$SFBSmwt1[take+1]<-lnN$SFBSmwt0[take]+mu_slp
  lnN$SFBSot1[take+1]<-lnN$SFBSot0[take]+mu_slp
  take<-which(lnN$Year<=1976)
  lnN$SFBSmwt2[take+2]<-lnN$SFBSmwt0[take]+mu_slp*2
  lnN$SFBSot2[take+2]<-lnN$SFBSot0[take]+mu_slp*2
  
  take<-which(lnN$Year==1979)
  lnN$SFBSmwt0[take]<-mean(c(lnN$SFBSmwt1[take+1]-mu_slp*1,
    lnN$SFBSmwt2[take+2]-mu_slp*2))
  lnN$SFBSmwt1[take]<-mean(c(lnN$SFBSmwt0[take-1]+mu_slp*1,
    lnN$SFBSmwt2[take+1]-mu_slp*1))
  lnN$SFBSmwt2[take]<-mean(c(lnN$SFBSmwt0[take-2]+mu_slp*2,
    lnN$SFBSmwt1[take-1]+mu_slp*1))
  lnN$SFBSot0[take]<-mean(c(lnN$SFBSot1[take+1]-mu_slp*1,
    lnN$SFBSot2[take+2]-mu_slp*2))
  lnN$SFBSot1[take]<-mean(c(lnN$SFBSot0[take-1]+mu_slp*1,
    lnN$SFBSot2[take+1]-mu_slp*1))
  lnN$SFBSot2[take]<-mean(c(lnN$SFBSot0[take-2]+mu_slp*2,
    lnN$SFBSot1[take-1]+mu_slp*1))
  
  take<-which(lnN$Year==1993)
  lnN$SFBSot2[take]<-mean(c(lnN$SFBSot0[take-2]+mu_slp*2,
    lnN$SFBSot1[take-1]+mu_slp*1))
  
  take<-which(lnN$Year==1995)
  lnN$SFBSmwt1[take]<-lnN$SFBSmwt0[take-1]+mu_slp*1
  lnN$SFBSmwt2[take]<-mean(c(lnN$SFBSmwt0[take-2]+mu_slp*2,
    lnN$SFBSmwt1[take-1]+mu_slp*1))
  
  take<-which(lnN$Year==1996)
  lnN$SFBSmwt1[take]<-mean(c(lnN$SFBSmwt0[take-1]+mu_slp*1,
    lnN$SFBSmwt2[take+1]-mu_slp*1))
  lnN$SFBSmwt2[take]<-mean(c(lnN$SFBSmwt0[take-2]+mu_slp*2,
    lnN$SFBSmwt1[take-1]+mu_slp*1))
  
  take<-which(lnN$Year==2014)
  lnN$SFBSmwt2[take]<-mean(c(lnN$SFBSmwt0[take-2]+mu_slp*2,
    lnN$SFBSmwt1[take-1]+mu_slp*1))
  
  take<-which(lnN$Year==2016)
  lnN$SFBSmwt0[take]<-lnN$SFBSmwt1[take+1]-mu_slp*1
  lnN$SFBSmwt1[take]<-mean(c(lnN$SFBSmwt0[take-1]+mu_slp*1,
    lnN$SFBSmwt2[take+1]-mu_slp*1))
  lnN$SFBSmwt2[take]<-mean(c(lnN$SFBSmwt0[take-2]+mu_slp*2,
    lnN$SFBSmwt1[take-1]+mu_slp*1))
  lnN$SFBSot1[take]<-mean(c(lnN$SFBSot0[take-1]+mu_slp*1,
    lnN$SFBSot2[take+1]-mu_slp*1))
  lnN$SFBSot2[take]<-mean(c(lnN$SFBSot0[take-2]+mu_slp*2,
    lnN$SFBSot1[take-1]+mu_slp*1))
  
  take<-which(lnN$Year==2017)
  lnN$SFBSot1[take]<-lnN$SFBSot0[take-1]+mu_slp*1
  
  for (i in 0:2)
  {
    take<-which(lnN$Year==2018+i)
    lnN$SFBSmwt1[take]<-lnN$SFBSmwt0[take-1]+mu_slp*1
    lnN$SFBSmwt2[take]<-mean(c(lnN$SFBSmwt0[take-2]+mu_slp*2,
      lnN$SFBSmwt1[take-1]+mu_slp*1))
    lnN$SFBSot1[take]<-lnN$SFBSot0[take-1]+mu_slp*1
    lnN$SFBSot2[take]<-mean(c(lnN$SFBSot0[take-2]+mu_slp*2,
      lnN$SFBSot1[take-1]+mu_slp*1))
  }
  return(list('lnN'=lnN,'mu_slp'=mu_slp))
}

#------------------------------------------------------------------------------------
#
# Step 2: Organize data
#
#------------------------------------------------------------------------------------

# load data
dat_raw<-read.csv('./Data/Length_Clusters_20210803.csv',header=TRUE)
colnames(dat_raw)[which(colnames(dat_raw)=='month')]<-"MonthCap"
#head(dat_raw)

# remove extreme outlier in the dataset and NA observations
rmv<-which(dat_raw$Length>=205|is.na(dat_raw$ClusterAge))
dat<-dat_raw[-rmv,]           # there are still some oddities, but it should be okay

# assume that hatch date is Jan 1 (so age-0 fish caught on jan 1 is 0.03 months old)
# However, additional hatching info from: 
# https://www.federalregister.gov/documents/2009/04/09/E9-8087/endangered-and-threatened-wildlife-and-plants-12-month-finding-on-a-petition-to-list-the-san 
# - spawning Nov - June (Moyle 2000 p. 236)
# - peak spawn Feb-April (Moyle 2000 p. 236)
# - Hatching ~40 days depending on temp (CDFG 2001 p. 477)

# extract month and year of capture
dat$MonthCap<-as.numeric(format(as.Date(dat$Date,format='%Y-%m-%d'),'%m'))
dat$YearCap<-as.numeric(format(as.Date(dat$Date,format='%Y-%m-%d'),'%Y'))

mon<-as.Date(paste(dat$MonthCap,"1",dat$YearCap,sep='/'),format='%m/%d/%Y')
temp<-dat$MonthCap+ifelse(dat$MonthCap==12,-11,1)
mon_p1<-as.Date(paste(temp,"1",dat$YearCap+ifelse(temp==1,1,0),sep='/'),format='%m/%d/%Y')

# fraction time from beginning of month to capture 
dat$F_AgeMonth<-as.numeric(difftime(as.Date(dat$Date,format='%Y-%m-%d'),mon,units='days'))/
  as.numeric(difftime(mon_p1,mon,units='days'))

# assume that hatch date is Jan 1 (so age-0 fish caught on jan 1 is 0.03 months old)
dat$W_AgeMonth<-dat$ClusterAge*12+dat$MonthCap-1
dat$FracAgeMonth<-dat$W_AgeMonth+dat$F_AgeMonth

cohorts<-sort(unique(dat$ClusterCohort))
n_cohorts<-length(cohorts)

n_gears<-length(unique(dat$Program))
gear_names<-levels(as.factor(dat$Program))
dat$Gear<-NA
min_size<-NA
for(g in 1:n_gears)
{
  take<-which(dat$Program==gear_names[g])
  dat$Gear[take]<-g
  min_size[g]<-min(dat$Length[take])
}
gears<-sort(unique(dat$Gear))

length_stats_fun<-function(X)
{
  take<-which(dat$Gear==unq[X,1]&dat$Year==unq[X,2])
  return(quantile(dat$Length[take],c(0,0.025,0.25,0.5,0.75,0.975,1)))
}
unq<-unique(cbind(dat$Gear,dat$Year))
out<-t(sapply(c(1:length(unq[,1])),length_stats_fun))
for(i in 1:5)print(cbind(unq[which(unq[,1]==i),2],out[which(unq[,1]==i),]))

# thin out n_rmv observations
# Gear 1 = 20mm survey, has 180k observations
# Gear 5 = SFBS, has 102k observations
# others collectively account for 26k
# so focusing on removing these from these surveys and 
# years in wich we had high numbers of observations (>10000)
#n_rmv<-100000          # 50000            # 
#ctoff<-c(6000,7000)    # c(9500,12000)    # 
#g_rmv<-NULL
#g_rmv[1]<-round(n_rmv*table(dat$Gear)[1]/sum(table(dat$Gear)[c(1,5)]))
#g_rmv[2]<-n_rmv-g_rmv[1]
#rmv<-NA

# remove 2/3 from this
#for(j in c(1,2))
#{
#  take_yrs<-as.numeric(names(which(table(dat$Year[which(dat$Gear==c(1,5)[j])])>ctoff[j])))
#  take<-which(table(dat$Year[which(dat$Gear==c(1,5)[j])])>ctoff[j])
#  tab<-table(dat$Year[which(dat$Gear==c(1,5)[j])])[take]
#  nj_rmv<-rmultinom(1,g_rmv[j],(tab-ctoff[j])/sum(tab-ctoff[j]))
#  for(i in 1:length(take_yrs))
#  {
#    take_obs<-which(dat$Year==take_yrs[i]&dat$Gear==c(1,5)[j])
#    qnt<-quantile(dat$Length[take_obs],c(0.05,0.95))
#    take_obs<-which(dat$Year==take_yrs[i]&dat$Gear==c(1,5)[j]&dat$Length>=qnt[1]&dat$Length<=qnt[2])
#    rmv<-c(rmv,take_obs[sample(1:length(take_obs),nj_rmv[i])])
#  }
#}

#table(dat$Year[rmv])
#table(dat$Gear[rmv])
#take_yrs<-which(table(dat$Year[which(dat$Gear==1)])>ctoff[1])
#tab<-table(dat$Year[which(dat$Gear==1)])[take_yrs]
#tab
#tab-table(dat$Year[rmv[which(dat$Gear[rmv]==1)]])
#table(dat$Year[which(dat$Gear==1)])

#take_yrs<-which(table(dat$Year[which(dat$Gear==5)])>ctoff[2])
#tab<-table(dat$Year[which(dat$Gear==5)])[take_yrs]
#tab
#tab-table(dat$Year[rmv[which(dat$Gear[rmv]==5)]])
#table(dat$Year[which(dat$Gear==5)])

#hist(dat$Length[rmv])
#dat<-dat[-rmv[-1],]

#------------------------------------------------------------------------------------
#
# Step 3: generate starting values
#
#------------------------------------------------------------------------------------

# the code below is used to generate startign values for the TMB model and plot
# the basic growth mdoel to assess model fit (i.e., does it look okay)
theta<-log(c(5,0.07,1,15,50))
take<-1:length(dat$Year)
inpt<-dat[take,]
lester_fun(theta,inpt)

fit2<-optim(theta,lester_fun,control=list(maxit=5000),'inpt'=inpt)
fit2
for(i in 1:5)fit2<-optim(fit2$par,lester_fun,control=list(maxit=2000),'inpt'=inpt)
fit2<-optim(fit2$par,lester_fun,control=list(maxit=2000),hessian=TRUE,'inpt'=inpt)
pars2<-c(exp(fit2$par[1:4]),fit2$par[5])
pars2
startpars<-fit2$par
stdev<-sqrt(diag(solve(fit2$hessian)))

lester_fun_p(fit2$par,'inpt'=inpt)

#------------------------------------------------------------------------------------
#
# Step 3: Process Environmental Data
#
#------------------------------------------------------------------------------------
 
# environmental data is organized by year (row) and columns represent
# monthly estimates for various metrics. We need to convert these to 
# a matrix organized by year-month (row). 

dflow_dat<-read.csv('./Data/Dayflow_Data_Summary.csv',head=TRUE)
n_env_year<-length(dflow_dat$WY)

take<-grep('maxout',names(dflow_dat))

df<-data.frame('WY'=as.vector(as.matrix(dflow_dat[,rep(1,length(take))])),'Year'=NA)
df$month<-rep(match(substr(names(dflow_dat[,take]),1,3),tolower(month.abb)),each=n_env_year)
df$Year<-df$WY-ifelse(df$month>=10,1,0)
df$maxout<-as.vector(as.matrix(dflow_dat[,grep('maxout',names(dflow_dat))]))
df$meanout<-as.vector(as.matrix(dflow_dat[,grep('meanout',names(dflow_dat))]))
df$medout<-as.vector(as.matrix(dflow_dat[,grep('medout',names(dflow_dat))]))
df$minout<-as.vector(as.matrix(dflow_dat[,grep('minout',names(dflow_dat))]))

flow_dat<-df[order(df$Year,df$month),]

rts_dat<-read.csv('./Data/Temperature_and_Salinity_Summary.csv',head=TRUE)
n_env_year<-length(rts_dat$WY)

site<-c('ANH','MAL','MRZ')
msmt<-c('cond','salinity','temp')
metric<-c('min','mean','max')

take1<-grep('ANH_cond',names(rts_dat))
take2<-grep('max',names(rts_dat))
take3<-which(nchar(names(rts_dat))==16)
take<-take1[take1%in%take2&take1%in%take3]

names(rts_dat[,take])

df<-data.frame('WY'=as.vector(as.matrix(rts_dat[,rep(1,length(take))])),Year=NA,
  'month'=rep(match(substr(names(rts_dat[,take]),10,12),month.abb),each=n_env_year),
  NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA)
df$Year<-df$WY-ifelse(df$month>=10,1,0)

iter<-4
for(i in 1:3)
{
  for(j in 1:3)
  {
    for(k in 1:3)
    {
      take1<-grep(paste(site[i],msmt[j],sep='_'),names(rts_dat))
      take2<-grep(metric[k],names(rts_dat))
      take3<-which(nchar(names(rts_dat))==(4+nchar(paste(site[i],msmt[j],metric[k],sep='_'))))
      take<-take1[take1%in%take2&take1%in%take3]
      names(df)[iter]<-paste(site[i],msmt[j],metric[k],sep='_')
      df[,iter]<-as.vector(as.matrix(rts_dat[,take]))
      iter<-iter+1
    }
  }
}

ts_dat<-df[order(df$Year,df$month),]

all_env_dat<-merge(flow_dat,ts_dat,all.x=TRUE,all.y=TRUE)

# assuming Jan 1 hatch date, need to take all env. data after min(cohorts)
env_dat<-all_env_dat[which(all_env_dat$Year>=min(cohorts)),]
env_dat$MAL_salinity_min[which(env_dat$MAL_salinity_min<0)]<-0.0001
env_dat$minout[which(env_dat$minout< -20000)]<-NA
env_dat$minout[which(env_dat$minout<0)]<-0.01
colSums(is.na(env_dat))

# fill in missing env. data with mean from that given month observed across all years
for(i in 7:length(env_dat[1,]))
{
  for(j in 1:12)
  {
    take<-which(is.na(env_dat[,i])&env_dat$month==j)
    if(length(take)>0)
    {
      mtake<-which(env_dat$month==j)
      env_dat[take,i]<-exp(mean(log(env_dat[mtake,i]),na.rm=TRUE))
    }
  }
}
#need env_dat at end because it is used throughout code later
env_dat_temp<-env_dat
env_dat<-NULL

abundance<-read.csv('./Data/countPVAdataforCraig.csv',head=TRUE)
lnN<-rbind(log(abundance),NA,NA)
lnN[,1]<-c(abundance[,1],2019,2020)

for(i in 2:length(lnN[1,]))
{
  lnN[which(!is.finite(lnN[,i])),i]<-NA
}

# just fills in the SFBS midwater trawl and otter trawl surveys
out<-fill_N_fun(lnN)
lnN_na.rm<-out$lnN[,c(1,3:8)]
mu_slp<-out$mu_slp
N_df<-exp(lnN_na.rm)
N_df$Year<-lnN_na.rm$Year

N_df$SFBSmwt1p<-N_df$SFBSmwt1+N_df$SFBSmwt2
N_df$SFBSot1p<-N_df$SFBSot1+N_df$SFBSot2

N_df2<-N_df[,c(1,2,8,5,9)]

# merge environmental and abundance datasets
env_dat_temp1<-merge(env_dat_temp,N_df2,all.x=TRUE)

# assume annual abundance declines over year, this uses the mean M (slope) 
# estimated from the progression of cohorts over time from the abundance indices
# to describe the monthly decline in abundance over time. 
env_dat_temp1$SFBSmwt0_t<-env_dat_temp1$SFBSmwt0*exp(mu_slp*(env_dat_temp1$month-1)/12)
env_dat_temp1$SFBSmwt1p_t<-env_dat_temp1$SFBSmwt1p*exp(mu_slp*(env_dat_temp1$month-1)/12)
env_dat_temp1$SFBSot0_t<-env_dat_temp1$SFBSot0*exp(mu_slp*(env_dat_temp1$month-1)/12)
env_dat_temp1$SFBSot1p_t<-env_dat_temp1$SFBSot1p*exp(mu_slp*(env_dat_temp1$month-1)/12)


# work up E. affiinis data
Eaffinis_dat<-read.csv('./Data/Eurytemora_affinis_Summary.csv',head=TRUE)
n_Ea_year<-length(Eaffinis_dat$WY)

take<-grep('NZ028',names(Eaffinis_dat))

df<-data.frame('WY'=as.vector(as.matrix(Eaffinis_dat[,rep(1,length(take))])),'Year'=NA)
df$month<-rep(match(substr(names(Eaffinis_dat[,take]),7,9),(month.abb)),each=n_Ea_year)
df$Year<-df$WY-ifelse(df$month>=10,1,0)
df$NZ028<-as.vector(as.matrix(Eaffinis_dat[,grep('NZ028',names(Eaffinis_dat))]))
df$NZ048<-as.vector(as.matrix(Eaffinis_dat[,grep('NZ048',names(Eaffinis_dat))]))
df$NZ054<-as.vector(as.matrix(Eaffinis_dat[,grep('NZ054',names(Eaffinis_dat))]))
df$NZ060<-as.vector(as.matrix(Eaffinis_dat[,grep('NZ060',names(Eaffinis_dat))]))
df$NZ064<-as.vector(as.matrix(Eaffinis_dat[,grep('NZ064',names(Eaffinis_dat))]))
df$NZ086<-as.vector(as.matrix(Eaffinis_dat[,grep('NZ086',names(Eaffinis_dat))]))
df$NZ325<-as.vector(as.matrix(Eaffinis_dat[,grep('NZ325',names(Eaffinis_dat))]))

Ea_dat<-df[order(df$Year,df$month),]
# calculate geometric mean across sampling stations
Ea_dat$EaGeoMean<-exp(rowMeans(log(Ea_dat[,-c(1:3)]+1),na.rm=TRUE))-1

# work up Mysis data
Mysis_dat<-read.csv('./Data/Mysida_Summary.csv',head=TRUE)
n_mys_year<-length(Mysis_dat$WY)

take<-grep('NZ028',names(Mysis_dat))

df<-data.frame('WY'=as.vector(as.matrix(Mysis_dat[,rep(1,length(take))])),'Year'=NA)
df$month<-rep(match(substr(names(Mysis_dat[,take]),7,9),(month.abb)),each=n_mys_year)
df$Year<-df$WY-ifelse(df$month>=10,1,0)
df$NZ028<-as.vector(as.matrix(Mysis_dat[,grep('NZ028',names(Mysis_dat))]))
df$NZ048<-as.vector(as.matrix(Mysis_dat[,grep('NZ048',names(Mysis_dat))]))
df$NZ054<-as.vector(as.matrix(Mysis_dat[,grep('NZ054',names(Mysis_dat))]))
df$NZ060<-as.vector(as.matrix(Mysis_dat[,grep('NZ060',names(Mysis_dat))]))
df$NZ064<-as.vector(as.matrix(Mysis_dat[,grep('NZ064',names(Mysis_dat))]))
df$NZ086<-as.vector(as.matrix(Mysis_dat[,grep('NZ086',names(Mysis_dat))]))
df$NZ325<-as.vector(as.matrix(Mysis_dat[,grep('NZ325',names(Mysis_dat))]))

mys_dat<-df[order(df$Year,df$month),]
# calculate geometric mean across sampling stations
mys_dat$MysGeoMean<-exp(rowMeans(log(mys_dat[,-c(1:3)]+1),na.rm=TRUE))-1

# merge geometric means from E. affinis and mysis datasets
prey_dat<-merge(Ea_dat[,-c(4:10)],mys_dat[,-c(4:10)],all.x=TRUE,all.y=TRUE)

# fill in first 3 months where observations are missing (use geometric mean of given month)
for(i in 1:3)
{
  take<-which(prey_dat$month==prey_dat$month[i])
  prey_dat$EaGeoMean[i]<-exp(mean(log(prey_dat$EaGeoMean[take[-1]]+1),na.rm=TRUE))-1
  prey_dat$MysGeoMean[i]<-exp(mean(log(prey_dat$MysGeoMean[take[-1]]+1),na.rm=TRUE))-1
}

# fill in last several months of missing data (only fill in consecuative missing months of data)
# and use geometric mean from a given month
take<-which(is.na(prey_dat$EaGeoMean))
i=0
run<-TRUE
while(run==TRUE)
{
  if(sum(take==(take[length(take)]-i))==1)
  {
    #prey_dat[(take[length(take)]-i),c(4,5)]<-0
    take1<-which(prey_dat$month==prey_dat$month[(take[length(take)]-i)])
    prey_dat$EaGeoMean[(take[length(take)]-i)]<-exp(mean(log(prey_dat$EaGeoMean[take1]+1),na.rm=TRUE))-1
    prey_dat$MysGeoMean[(take[length(take)]-i)]<-exp(mean(log(prey_dat$MysGeoMean[take1]+1),na.rm=TRUE))-1
  }else{
    run<-FALSE
  }
  i<-i+1
}

# note, if missing values from E.affinis, also missing mysis data and vis versa
# use linear interpolation to fill in other missing value 
run<-TRUE
catch<-0
while(run==TRUE)
{
  catch<-catch+1
  if(sum(is.na(prey_dat$EaGeoMean))==0|catch==1000)
  {
    run<-FALSE
  }else{
    take<-which(is.na(prey_dat$EaGeoMean))[1]
    take_plus<-sum(is.na(prey_dat$EaGeoMean[take+c(1,2)]))
    for(i in 4:5)
    {
      s<-seq((prey_dat[take-1,i]),(prey_dat[(take+take_plus+1),i]),
        length=(3+take_plus))
      prey_dat[take:(take+take_plus),i]<-s[-c(1,length(s))]
    }
  }
}
rowSums(is.na(prey_dat))

env_dat<-merge(env_dat_temp1,prey_dat,all.x=TRUE)

#------------------------------------------------------------------------------------
#
# Step 4: organize data to pass to TMB
#
#------------------------------------------------------------------------------------

amax<-ceiling(max(dat$FracAgeMonth))             # goes from 0 to max observed age
temp<-array(NA,c(n_cohorts,amax+1))
for(crt in 1:n_cohorts)
{
  temp[crt,]<-cohorts[crt]+c(rep((1:(amax/12))-1,each=12),amax/12)
}

temp[which(temp>max(dat$Year))]<-max(dat$Year)
yrmat<-temp-min(temp)

years<-sort(unique(as.vector(temp)))
n_years<-length(years)

monvec<-c(rep(0:11,(amax)/12),0)

rmv<-which(gear_names=='SFBS')
min_size2<-min_size
min_size2[-rmv]<-min(min_size[-rmv])


std_env_dat<-env_dat
ln_env_dat<-log(env_dat[,-c(1:3)])
std_env_dat[,-c(1:3)]<-t(t(ln_env_dat)-colMeans(ln_env_dat))#/apply(ln_env_dat,2,sd)
# use 'raw' temperature
take<-grep('temp',names(env_dat))
std_env_dat[,take]<-t(t(env_dat[,take])-colMeans(env_dat[,take]))#/apply(env_dat[,take],2,sd)
# E.affinis and mysis, add 1+
take<-grep('GeoMean',names(env_dat))
temp<-log(env_dat[,take]+1)
std_env_dat[,take]<-t(t(temp)-colMeans(temp))#/apply(env_dat[,take],2,sd)
# center month
std_env_dat$c_month<-std_env_dat$month-mean(c(1:12))



# create vectors for interaction effects 
# temperature * outflow
std_env_dat$temp_out<-std_env_dat$meanout*std_env_dat$MRZ_temp_mean

# otter trawl abundance index * time
std_env_dat$ot_N0_t<-std_env_dat$SFBSot0*std_env_dat$c_month
std_env_dat$ot_N1p_t<-std_env_dat$SFBSot1p*std_env_dat$c_month

# otter trawl abundance index * E. affinis
std_env_dat$ot_N0_Ea<-std_env_dat$SFBSot0*std_env_dat$EaGeoMean
std_env_dat$ot_N1p_Ea<-std_env_dat$SFBSot1p*std_env_dat$EaGeoMean

# otter trawl abundance index * mysis
std_env_dat$ot_N0_mys<-std_env_dat$SFBSot0*std_env_dat$MysGeoMean
std_env_dat$ot_N1p_mys<-std_env_dat$SFBSot1p*std_env_dat$MysGeoMean

# Per capita E. affinis - otter trawl abundance index
std_env_dat$pcap_Ea_ot_N0<-std_env_dat$EaGeoMean-std_env_dat$SFBSot0
std_env_dat$pcap_Ea_ot_N1p<-std_env_dat$EaGeoMean-std_env_dat$SFBSot1p

# Per capita mysis - otter trawl abundance index
std_env_dat$pcap_mys_ot_N0<-std_env_dat$MysGeoMean-std_env_dat$SFBSot0
std_env_dat$pcap_mys_ot_N1p<-std_env_dat$MysGeoMean-std_env_dat$SFBSot1p

# Per capita E. affinis * Time - otter trawl abundance index
std_env_dat$pcap_Ea_ot_N0_t<-(std_env_dat$EaGeoMean-std_env_dat$SFBSot0)*std_env_dat$c_month
std_env_dat$pcap_Ea_ot_N1p_t<-(std_env_dat$EaGeoMean-std_env_dat$SFBSot1p)*std_env_dat$c_month

# Per capita mysis * Time - otter trawl abundance index
std_env_dat$pcap_mys_ot_N0_t<-(std_env_dat$MysGeoMean-std_env_dat$SFBSot0)*std_env_dat$c_month
std_env_dat$pcap_mys_ot_N1p_t<-(std_env_dat$MysGeoMean-std_env_dat$SFBSot1p)*std_env_dat$c_month


# add a zero (i.e., mean) to last row to account for unobserved future growth
# this is an accounting necessity for the TBM model
std_env_dat0<-rbind(std_env_dat,0)
env_dat0<-rbind(env_dat,colMeans(env_dat))

# note, we need the c(which(temperature),which(others)) to ensure
# that the temperature effect is first
take_env_dat<-c(which(names(std_env_dat0)=='MRZ_temp_mean'), # Temperature
  which(names(std_env_dat0)=='SFBSot0'|                            # age-0 abudnance index
  names(std_env_dat0)=='SFBSot1p'|                           # age-1+ abundance index
#  names(std_env_dat0)=='c_month'|                             # month (e.g., time of year)
  names(std_env_dat0)=='EaGeoMean'|                           # E. affinis
  names(std_env_dat0)=='MysGeoMean'|                          # Mysis
  names(std_env_dat0)=='meanout'|                       # outfow
#  names(std_env_dat0)=='MRZ_salinity_mean'|                       # conductivity
  names(std_env_dat0)=='MRZ_cond_mean'))                       # conductivity
  #names(std_env_dat0)=='MRZ_mean_temp_outflow'|               # temperature*outflow
  #names(std_env_dat0)=='MRZ_mean_outflow_cond'|               # outflow*conductivity
#  names(std_env_dat0)=='mwt_N0_Ea'|                           # age-0 * time
#  names(std_env_dat0)=='mwt_N1p_Ea'|                          # Age1+ * time
#  names(std_env_dat0)=='mwt_N0_mys'|                          # age-0 * time
#  names(std_env_dat0)=='mwt_N1p_mys'))                        # Age1+ * time
#  names(std_env_dat0)=='pcap_Ea_mwt_N0'|                           # Age-0 Per capita E. affinis
#  names(std_env_dat0)=='pcap_Ea_mwt_N1p'|                          # Age1+ Per capita E. affinis
#  names(std_env_dat0)=='pcap_mys_mwt_N0'|                          # age-0 Per capita mysis
#  names(std_env_dat0)=='pcap_mys_mwt_N1p'))                        # Age1+ Per capita mysis
      

# index matrix for environmental data
env_imat<-matrix(1:length(monvec)-1,ncol=length(monvec),nrow=n_cohorts,byrow=TRUE)+
  matrix(seq(0,by=12,length=n_cohorts),ncol=length(monvec),nrow=n_cohorts)
env_imat[which(env_imat>(dim(std_env_dat)[1]-1))]<-(dim(std_env_dat0)[1])-1

take_dat<-which(as.Date(dat$Date)<as.Date(paste(env_dat$Year[length(env_dat$Year)],
  env_dat$month[length(env_dat$Year)]+1,'01',sep='/'),format='%Y/%m/%d'))


# compile data to inform lo
take1<-c(1,which(names(dflow_dat)=='spawn_days_over_20k'):length(names(dflow_dat)))
take2<-which(dflow_dat$WY>=min(cohorts))
#dflow_dat[take2,take1]

take3<-c(1,grep('season',names(rts_dat)),grep('day',names(rts_dat)))
#rts_dat[,take3]
env_lo_dat<-merge(dflow_dat[take2,take1],rts_dat[,take3],all.x=TRUE,all.y=TRUE)

# data from 1983 are incomplete and spawn season metrics are biased
take<-grep('Spawn_season',colnames(env_lo_dat))
env_lo_dat[1:5,take]<-NA

# 
std_env_lo_dat<-env_lo_dat

# various metrics associated with number of days under/over a certain
# value (e.g., 9 dev C) often have 0, which prevent taking the log
# when standardizing. Added 1 to all these metrics to allow for log
# transformation and consistency across all metrics. 
take<-grep('days_under',names(std_env_lo_dat))
std_env_lo_dat[,take]<-env_lo_dat[,take]+1
take<-grep('days_over',names(std_env_lo_dat))
std_env_lo_dat[,take]<-env_lo_dat[,take]+1


temp<-std_env_lo_dat[,-1]
temp[temp<0]<-0.00001
#temp[temp<=0]<-0.00001

std_env_lo_dat[,-1]<-t(t(log(temp))-colMeans(log(temp),na.rm=TRUE))
std_env_lo_dat[is.na(std_env_lo_dat)]<-0

std_env_lo_dat$out_temp<-std_env_lo_dat$spawn_days_over_20k*std_env_lo_dat$MRZ_temp_Spawn_season_mean
std_env_lo_dat$out_under9<-std_env_lo_dat$spawn_days_over_20k*std_env_lo_dat$MRZ_days_under_9deg
std_env_lo_dat$cond_temp<-std_env_lo_dat$MRZ_cond_Spawn_season_mean*std_env_lo_dat$MRZ_temp_Spawn_season_mean
std_env_lo_dat$cond_under9<-std_env_lo_dat$MRZ_cond_Spawn_season_mean*std_env_lo_dat$MRZ_days_under_9deg
std_env_lo_dat$sal_temp<-std_env_lo_dat$MRZ_salinity_Spawn_season_mean*std_env_lo_dat$MRZ_temp_Spawn_season_mean
std_env_lo_dat$sal_under9<-std_env_lo_dat$MRZ_salinity_Spawn_season_mean*std_env_lo_dat$MRZ_days_under_9deg

take_env_lo_dat<-which(
  names(std_env_lo_dat)=='spawn_days_over_20k'|         # spawning days over 20k cfs
  names(std_env_lo_dat)=='MRZ_cond_Spawn_season_mean'|  # mean cond over spawning season
  names(std_env_lo_dat)=='MRZ_temp_Spawn_season_mean')  # mean temperature over spawning season

#------------------------------------------------------------------------------------
#
# Step 5: Set up TMB model
#
#------------------------------------------------------------------------------------

tmbdat<-list('amax'=amax+1,					# maximum whole age (in month)
  #'nyr'=n_years,					# number of years (same as number of cohorts in this model)
  'ncrt'=n_cohorts,					# number of years (same as number of cohorts in this model)
  #'crtmat'=crtmat-min(crtmat),	# matrix of pointers to each cohort by year*age
  'yrmat'=yrmat,
  'monvec'=monvec,
  # catch data
  'nobs'=length(take_dat),					# number of catch observations with age & length
  'year'=dat$Year[take_dat]-min(dat$Year[take_dat]),			# year caught
  'MonCap'=dat$MonthCap[take_dat]-1,
  'cohort'=dat$ClusterCohort[take_dat]-min(dat$ClusterCohort[take_dat]),				# cohort
  'len'=dat$Length[take_dat],						# length
  'wage'=dat$ClusterMonth[take_dat]-1,				# whole age
  'fage'=dat$F_AgeMonth[take_dat],
  # selectivity
  'gear'=dat$Gear[take_dat]-1,					# gear fish was caught with
  'min_size'=min_size2,
  'n_beta'=length(take_env_dat)-1,
  'n_alpha'=length(take_env_lo_dat),
  'env_dat'=as.matrix(std_env_dat0[,take_env_dat]),
  'env_imat'=env_imat,   # pointers to environmental variables
  'env_lo_dat'=as.matrix(std_env_lo_dat[,take_env_lo_dat]))

hypsig<-0.5
dunf<-c(-0.05,0.05)

# update to make this a function
start<-function(){list('lnmuh'=startpars[1]+runif(1,dunf[1],dunf[2]),	# hyperprior mean log(h)
  'lnsighy'=log(hypsig)+runif(1,dunf[1],dunf[2]),					# hyperprior sigma log(h)
  'hy'=rep(0,length=n_years),		# RE for log(h)
  'lnsighc'=log(hypsig/10)+runif(1,dunf[1],dunf[2]),					# hyperprior sigma log(h)
  'hc'=rep(startpars[1],n_cohorts),		# RE for log(h)
  'lnsigg'=log(hypsig)+runif(1,dunf[1],dunf[2]),								# hyperprior mean log(g)
  'lng'=rep(startpars[2],12),								# hyperprior mean log(g)
  'mulo'=runif(1,dunf[1],dunf[2]),										# lo - length-at-age 0
  'lnsiglo'=log(hypsig)+runif(1,dunf[1],dunf[2]),										# lo - length-at-age 0
  'lo'=runif(n_cohorts,dunf[1],dunf[2]),										# lo - length-at-age 0
  'lnsig'=log(5)+runif(1,dunf[1],dunf[2]),							# log(cv)
  'lnmua50'=startpars[4]+runif(1,dunf[1],dunf[2]),				# hyperprior mean log(a50) - age @ 50% maturity
  'lnsiga50'=log(hypsig)+runif(1,dunf[1],dunf[2]),						# hyperprior sigma log(a50)
  'lna50'=rep(startpars[4],length=n_cohorts),	# RE for log(a50)
  'b_ha'=rep(0,tmbdat$n_beta),
  'b_hj'=rep(0,tmbdat$n_beta),
  'tcrit'=log(c(15,15)),
  'alpha'=rep(0,tmbdat$n_alpha))}

# the map list is used in TMB to turn on and off estimation of parameter and vectors
# for instance, with REs, we can set the RE_vec=factor(rep(1,n_RE)) so that TMB will  
# only estimate a single value for this vector and keep it constant over the RE. THis
# effectively removes the RE from the calculation. We can also set the RE hyperparameters
# to myp_mu = factor(NA) to prevent TMB from estimating a value for that. 
map<-list(#'hy'=factor(rep(NA,n_years)),'lnsighy'=factor(NA),
  'hy'=factor(c(NA,NA,NA,seq(1,n_years-3,by=1))),
  #'hy'=factor(c(1,1,seq(1,n_years-2,by=1))),
  'hc'=factor(rep(1,n_cohorts)),'lnsighc'=factor(NA),'lnmuh'=factor(NA),
  #'hc'=factor(c(1,1,1,1:37,37,37)),
  'lo'=factor(rep(1,n_cohorts)),'lnsiglo'=factor(NA),'mulo'=factor(NA),
  #'lo'=factor(c(1,1,1:39,39)),
  'lng'=factor(rep(1,12)),'lnsigg'=factor(NA),
  'lna50'=factor(rep(1,n_cohorts)),'lnmua50'=factor(NA),'lnsiga50'=factor(NA)
  #'lna50'=factor(c(1,1,1:37,37))
  #'b_hj'=factor(rep(NA,tmbdat$n_beta)),'tcrit'=factor(c(1,1)),
)

# assign with parameters are random effects
REs<-c('hy','hc','lng','lna50','lo','lnmuh','mulo','b_ha','b_hj','alpha','tcrit')


fil<-'Longfin_TMB_FracAge_V3'		# file name
compile(paste(fil,'.cpp',sep=''))	# compile TMB code
dyn.load(dynlib(fil))				# load dynamic library

# ignore the function below, it is used to test the model and changes to the model 
fun_ignore<-function()
{
  # set function/model objects to NULL
  obj<-NULL
  fit<-NULL
  sdr<-NULL
  
  # create model object
  obj<-MakeADFun(tmbdat,start(),DLL=fil,random=REs,map=map,
    inner.control=list(maxit=5000))
  obj$report()
  
  # fit model, if doesn't converge in iterations, run again
  fit<-nlminb(obj$par,obj$fn,obj$gr,contol=list(eval.max=400,iter.max=300))
  fit

  sdr<-sdreport(obj)
  summary(sdr)
}



# water year totals reported in https://cdec.water.ca.gov/reportapp/javareports?name=WSIHIST
WY_78_17<-c(33.57,18.39,31.8,14.32,44.82,52.69,29.48,14.64,35.33,
  11.35,11.71,18.38,11.72,11.64,11.45,30.59,10.35,46.87,29.51,34.93,
  41.83,27.1,24.8,12.99,18.66,24.18,19.85,27.76,42.53,12.79,13.77,
  17.96,22.09,36.2,14.6,15.24,9.18,10.67,22.54,52.66)

# create some diagnostic plots - not needed unless changing model
fun_ignore2<-function()
{
  out<-lpred_fun(sdr,map)
  
  lpred<-out[[1]]
  ltran<-out[[2]]
  
  hist(lpred[,-1]-lpred[,-amax],breaks=100)
  
  hy_out<-summary(sdr)[which(rownames(summary(sdr))=='hy'),1][map$hy]
  hy_out[is.na(hy_out)]<-0
  
  plot(exp(summary(sdr)[which(rownames(summary(sdr))=='hc'),1][map$hc]))
  par(mar=c(5,5,0.1,.1))
  plot(exp(summary(sdr)[which(rownames(summary(sdr))=='hc'),1][map$hc[1]])*
      exp(hy_out)~years,type='l',lwd=3,las=1,xlab="Year",ylab='Growth Rate (h)',cex.lab=1.7,cex.axis=1.5)
  
  
  plot(exp(summary(sdr)[which(rownames(summary(sdr))=='hc'),1][map$hc])~
      years,type='l',lwd=3,las=1,xlab="Cohort",ylab='Growth Rate (h)',cex.lab=1.7,cex.axis=1.5)
  
  plot(lpred[,1]~years,type='l',lwd=3,las=1,xlab="Year",ylab='Length-at-age 0 (l0)',cex.lab=1.7,cex.axis=1.5)
  
  plot(exp(summary(sdr)[which(rownames(summary(sdr))=='hc'),1][map$hc])*exp(hy_out)~lpred[,1],
    pch=20,lwd=3,las=1,xlab="Length-at-age 0 (l0)",ylab='Max. Annual Growth (mm)',cex.lab=1.7,cex.axis=1.5)
  
  plot(exp(hy_out)~lpred[,1],
    pch=20,lwd=3,las=1,xlab="Length-at-age 0 (l0)",ylab='Year-specific growth effect',cex.lab=1.7,cex.axis=1.5)
  
  plot((hy_out)~lpred[,1],
    pch=20,lwd=3,las=1,xlab="Length-at-age 0 (l0)",ylab='Year-specific growth effect',cex.lab=1.7,cex.axis=1.5)
  
  
  s<-(tmbdat$env_dat[,4])
  windows(height=9,width=9)
  temp<-ceiling(sum(!is.na(c(map$b_km,map$b_hm,map$b_hm2)))/2)*2
  layout(matrix(1:temp,nrow=2))
  
  if(sum(is.na(map$b_km))!=5)
  {
    b<-summary(sdr)[which(rownames(summary(sdr))=='b_km'),1][map$b_km]
    take<-which(!is.na(b))
    for(i in take)plot(b[i]*s~s,type='l',lwd=3,main=paste('km, i = ',i,sep=''))
  }
  if(sum(is.na(map$b_hm2))!=5)
  {
    b<-summary(sdr)[which(rownames(summary(sdr))=='b_hm2'),1][map$b_hm2]
    take<-which(!is.na(b))
    for(i in take)plot(b[i]*s~s,type='l',lwd=3,main=paste('hm2, i = ',i,sep=''))
  }
  if(sum(is.na(map$b_hm))!=5)
  {
    b<-summary(sdr)[which(rownames(summary(sdr))=='b_hm'),1][map$b_hm]
    take<-which(!is.na(b))
    for(i in take)plot(b[i]*s~s,type='l',lwd=3,main=paste('hm, i = ',i,sep=''))
  }
  
  
  s<-(tmbdat$env_dat[,1]-summary(sdr)[which(rownames(summary(sdr))=='tcrit'),1][1])^2
  tpred_a<-s*summary(sdr)[which(rownames(summary(sdr))=='b_ha'),1][1]
  s<-(tmbdat$env_dat[,1]-summary(sdr)[which(rownames(summary(sdr))=='tcrit'),1][2])^2
  tpred_j<-s*(summary(sdr)[which(rownames(summary(sdr))=='b_ha'),1][1]+
    summary(sdr)[which(rownames(summary(sdr))=='b_hj'),1][1])
  
  plot(tpred_a~tmbdat$env_dat[,1],ylim=range(c(tpred_a,tpred_j)),pch=19)
  points(tpred_j~tmbdat$env_dat[,1],pch=19,col='blue')
  
  summary(sdr)[which(rownames(summary(sdr))=='b_hm'),1][map$b_hm][5]+mean(env_dat0[,take_env_dat[4]])
  
  #temp<-1/(1+exp(-(s/1)))
  #plot(s*temp*summary(sdr)[which(rownames(summary(sdr))=='b_hm2'),1][map$b_hm2][5]~s)
  
  
  
  
  
  windows(width=20,height=8,record=TRUE)
  layout(matrix(1:10,nrow=2,byrow=TRUE))
  for(crt in 1:n_cohorts)
  {
    take<-which(dat$ClusterCohort==cohorts[crt])
    plot(dat$Length[take]~dat$FracAgeMonth[take],las=1,xlim=c(0,amax),ylim=c(0,max(dat$Length)),pch=20,
      main=cohorts[crt],xlab='Age in months',ylab='Length (mm)')
    
    age<-(0:amax)
    lines(lpred[crt,]~age,lwd=3,col='steelblue')
  }
}


#------------------------------------------------------------------------------------
#
# Step 6: run TMB models for model comparison 
#
#------------------------------------------------------------------------------------

# model fitting and comparison

fit_lst<-NULL
obj_lst<-NULL
sdr_lst<-NULL
nll<-NULL

# first comparison (i=1 and i=2) is looking at cohort vs year RE structure
# NO INTERACTIONS - keeping simple
i=1
# cohort-specific RE
# mean temperature, mean outflow, abundance (otter trawl) for age-0 and age-1+,
# and prey densities (E. afinnis and Mysis)
take_env_dat<-c(which(names(std_env_dat0)=='MRZ_temp_mean'),
  which(names(std_env_dat0)=='meanout'|
  names(std_env_dat0)=='SFBSot0'|
  names(std_env_dat0)=='SFBSot1p'|
  names(std_env_dat0)=='EaGeoMean'|
  names(std_env_dat0)=='MysGeoMean'))
# Spawning Outflow index (days over 20k cfs) and mean spawn season temperature
take_env_lo_dat<-which(names(std_env_lo_dat)=='spawn_days_over_20k'|
  names(std_env_lo_dat)=='MRZ_temp_Spawn_season_mean')
tmbdat$n_beta<-length(take_env_dat)
tmbdat$n_alpha<-length(take_env_lo_dat)
tmbdat$env_dat<-as.matrix(std_env_dat0[,take_env_dat])
tmbdat$env_lo_dat<-as.matrix(std_env_lo_dat[,take_env_lo_dat])

#map$b_hj<-factor(rep(NA,tmbdat$n_beta))
#map$tcrit<-factor(c(1,1))
map1<-map
map1$hy<-factor(rep(NA,n_years))
map1$lnsighy<-factor(NA)
map1$lnsighc<-factor(1)
map1$lnmuh<-factor(1)
map1$hc<-factor(c(1,1,1,1:37,37,37))

map1$b_hj<-factor(seq(1,tmbdat$n_beta,by=1))
map1$tcrit<-factor(c(1:2))

REs<-c('hy','hc','lng','lna50','lo','lnmuh','mulo','b_ha','b_hj','alpha','tcrit')

for(j in 1:5)
{
  # set function/model objects to NULL
  obj<-NULL

  # create model object
  obj<-MakeADFun(tmbdat,start(),DLL=fil,random=REs,map=map1,
    inner.control=list(maxit=5000))
  obj$report()
  
  if(j==1)
  {
    # fit model, if doesn't converge in alloted iterations, run again
    converge<-FALSE
    iter<-0
    while(converge==FALSE&iter<5)
    {
      fit_lst[[i]]<-nlminb(obj$par,obj$fn,obj$gr,contol=list(eval.max=400,iter.max=300))
      obj_lst[[i]]<-obj
      converge<-ifelse(sum(grep('(4)',fit_lst[[i]]$message))==1,TRUE,FALSE)
      iter<-iter+1
      gc()
    }
  }else{
    fit<-NULL
    fit<-nlminb(obj$par,obj$fn,obj$gr,contol=list(eval.max=400,iter.max=300))
    converge<-ifelse(sum(grep('(8)',fit$message))==1,FALSE,TRUE)
    if(converge==TRUE&fit$objective<fit_lst[[i]]$objective)
    {
      fit_lst[[i]]<-fit
      obj_lst[[i]]<-obj
    }
    gc()
  }
}
fit_lst[[i]]
nll[i]<-fit_lst[[i]]$objective
sdr_lst[[i]]<-sdreport(obj_lst[[i]])



#-------------------------------------------------------------------------------------
i=2
# Year-specific RE
# mean temperature, mean outflow, abundance (otter trawl) for age-0 and age-1+,
# and prey densities (E. afinnis and Mysis)
take_env_dat<-c(which(names(std_env_dat0)=='MRZ_temp_mean'),
  which(names(std_env_dat0)=='meanout'|
  names(std_env_dat0)=='SFBSot0'|
  names(std_env_dat0)=='SFBSot1p'|
  names(std_env_dat0)=='EaGeoMean'|
  names(std_env_dat0)=='MysGeoMean'))
# Spawning Outflow index (days over 20k cfs) and mean spawn season temperature
take_env_lo_dat<-which(names(std_env_lo_dat)=='spawn_days_over_20k'|
  names(std_env_lo_dat)=='MRZ_temp_Spawn_season_mean')
tmbdat$n_beta<-length(take_env_dat)
tmbdat$n_alpha<-length(take_env_lo_dat)
tmbdat$env_dat<-as.matrix(std_env_dat0[,take_env_dat])
tmbdat$env_lo_dat<-as.matrix(std_env_lo_dat[,take_env_lo_dat])

map1<-map
map1$b_hj<-factor(seq(1,tmbdat$n_beta,by=1))
map1$tcrit<-factor(c(1:2))

REs<-c('hy','hc','lng','lna50','lo','lnmuh','mulo','b_ha','b_hj','alpha','tcrit')

for(j in 1:5)
{
  # set function/model objects to NULL
  obj<-NULL

  # create model object
  obj<-MakeADFun(tmbdat,start(),DLL=fil,random=REs,map=map1,
    inner.control=list(maxit=5000))
  obj$report()
  
  if(j==1)
  {
    # fit model, if doesn't converge in alloted iterations, run again
    converge<-FALSE
    iter<-0
    while(converge==FALSE&iter<5)
    {
      fit_lst[[i]]<-nlminb(obj$par,obj$fn,obj$gr,contol=list(eval.max=400,iter.max=300))
      obj_lst[[i]]<-obj
      converge<-ifelse(sum(grep('(4)',fit_lst[[i]]$message))==1,TRUE,FALSE)
      iter<-iter+1
      gc()
    }
  }else{
    fit<-NULL
    fit<-nlminb(obj$par,obj$fn,obj$gr,contol=list(eval.max=400,iter.max=300))
    converge<-ifelse(sum(grep('(8)',fit$message))==1,FALSE,TRUE)
    if(converge==TRUE&fit$objective<fit_lst[[i]]$objective)
    {
      fit_lst[[i]]<-fit
      obj_lst[[i]]<-obj
    }
    gc()
  }
}
fit_lst[[i]]
nll[i]<-fit_lst[[i]]$objective
sdr_lst[[i]]<-sdreport(obj_lst[[i]])

#save.image("./Longfin0.RData")

nll
# compare juv/adult growth effects using year RE structure (best)
# hy is best
#-------------------------------------------------------------------------------------
# betas for adult=juvenile
i=3

# mean temperature, mean outflow, abundance (otter trawl) for age-0 and age-1+,
# and prey densities (E. afinnis and Mysis)
take_env_dat<-c(which(names(std_env_dat0)=='MRZ_temp_mean'),
  which(names(std_env_dat0)=='meanout'|
  names(std_env_dat0)=='SFBSot0'|
  names(std_env_dat0)=='SFBSot1p'|
  names(std_env_dat0)=='EaGeoMean'|
  names(std_env_dat0)=='MysGeoMean'))
# Spawning Outflow index (days over 20k cfs) and mean spawn season temperature
take_env_lo_dat<-which(names(std_env_lo_dat)=='spawn_days_over_20k'|
  names(std_env_lo_dat)=='MRZ_temp_Spawn_season_mean')
tmbdat$n_beta<-length(take_env_dat)
tmbdat$n_alpha<-length(take_env_lo_dat)
tmbdat$env_dat<-as.matrix(std_env_dat0[,take_env_dat])
tmbdat$env_lo_dat<-as.matrix(std_env_lo_dat[,take_env_lo_dat])

map1<-map
map1$b_hj<-factor(rep(NA,tmbdat$n_beta))
map1$tcrit<-factor(c(1,1))

REs<-c('hy','hc','lng','lna50','lo','lnmuh','mulo','b_ha','b_hj','alpha','tcrit')

for(j in 1:5)
{
  # set function/model objects to NULL
  obj<-NULL

  # create model object
  obj<-MakeADFun(tmbdat,start(),DLL=fil,random=REs,map=map1,
    inner.control=list(maxit=5000))
  obj$report()
  
  if(j==1)
  {
    # fit model, if doesn't converge in alloted iterations, run again
    converge<-FALSE
    iter<-0
    while(converge==FALSE&iter<5)
    {
      fit_lst[[i]]<-nlminb(obj$par,obj$fn,obj$gr,contol=list(eval.max=400,iter.max=300))
      obj_lst[[i]]<-obj
      converge<-ifelse(sum(grep('(4)',fit_lst[[i]]$message))==1,TRUE,FALSE)
      iter<-iter+1
      gc()
    }
  }else{
    fit<-NULL
    fit<-nlminb(obj$par,obj$fn,obj$gr,contol=list(eval.max=400,iter.max=300))
    converge<-ifelse(sum(grep('(8)',fit$message))==1,FALSE,TRUE)
    if(converge==TRUE&fit$objective<fit_lst[[i]]$objective)
    {
      fit_lst[[i]]<-fit
      obj_lst[[i]]<-obj
    }
    gc()
  }
}
fit_lst[[i]]
nll[i]<-fit_lst[[i]]$objective
sdr_lst[[i]]<-sdreport(obj_lst[[i]])

nll[-1]

save.image("./Longfin1.RData")

#-----------------------------------------------------------------------------------
#
# best is year RE with separate effects for juv and adult growth 
#
#-----------------------------------------------------------------------------------


# compare temperature parameter for lo
# - Model 2 is the best with mean temperature for lo
# - Model 4 is the Ha with days over 9C for lo

#-------------------------------------------------------------------------------------
i=4
# mean temperature, mean outflow, abundance (otter trawl) for age-0 and age-1+,
# and prey densities (E. afinnis and Mysis)
take_env_dat<-c(which(names(std_env_dat0)=='MRZ_temp_mean'),
  which(names(std_env_dat0)=='meanout'|
  names(std_env_dat0)=='SFBSot0'|
  names(std_env_dat0)=='SFBSot1p'|
  names(std_env_dat0)=='EaGeoMean'|
  names(std_env_dat0)=='MysGeoMean'))
# Spawning Outflow index (days over 20k cfs) and number of spawning days under 9C
take_env_lo_dat<-which(names(std_env_lo_dat)=='spawn_days_over_20k'|
  names(std_env_lo_dat)=='MRZ_days_under_9deg')
tmbdat$n_beta<-length(take_env_dat)
tmbdat$n_alpha<-length(take_env_lo_dat)
tmbdat$env_dat<-as.matrix(std_env_dat0[,take_env_dat])
tmbdat$env_lo_dat<-as.matrix(std_env_lo_dat[,take_env_lo_dat])

map1<-map
map1$b_hj<-factor(seq(1,tmbdat$n_beta,by=1))
map1$tcrit<-factor(c(1:2))

REs<-c('hy','hc','lng','lna50','lo','lnmuh','mulo','b_ha','b_hj','alpha','tcrit')

for(j in 1:5)
{
  # set function/model objects to NULL
  obj<-NULL

  # create model object
  obj<-MakeADFun(tmbdat,start(),DLL=fil,random=REs,map=map1,
    inner.control=list(maxit=5000))
  obj$report()
  
  if(j==1)
  {
    # fit model, if doesn't converge in alloted iterations, run again
    converge<-FALSE
    iter<-0
    while(converge==FALSE&iter<5)
    {
      fit_lst[[i]]<-nlminb(obj$par,obj$fn,obj$gr,contol=list(eval.max=400,iter.max=300))
      obj_lst[[i]]<-obj
      converge<-ifelse(sum(grep('(4)',fit_lst[[i]]$message))==1,TRUE,FALSE)
      iter<-iter+1
      gc()
    }
  }else{
    fit<-NULL
    fit<-nlminb(obj$par,obj$fn,obj$gr,contol=list(eval.max=400,iter.max=300))
    converge<-ifelse(sum(grep('(8)',fit$message))==1,FALSE,TRUE)
    if(converge==TRUE&fit$objective<fit_lst[[i]]$objective)
    {
      fit_lst[[i]]<-fit
      obj_lst[[i]]<-obj
    }
    gc()
  }
}
nll[i]<-fit_lst[[i]]$objective
sdr_lst[[i]]<-sdreport(obj_lst[[i]])
fit_lst[[i]]

nll[c(2,4)]-min(nll[c(2,4)])

save.image("./Longfin2.RData")

#-------------------------------------------------------------------------------------
# compare outflow and conductivity for growth increment and lo
# - Model 2 is the best with outflow
# - Model 5 is the Ha with conductivity
#-------------------------------------------------------------------------------------
i=5
# mean temperature, mean Conductivity, abundance (otter trawl) for age-0 and age-1+,
# and prey densities (E. afinnis and Mysis)
take_env_dat<-c(which(names(std_env_dat0)=='MRZ_temp_mean'),
  which(names(std_env_dat0)=='MRZ_cond_mean'|
  names(std_env_dat0)=='SFBSot0'|
  names(std_env_dat0)=='SFBSot1p'|
  names(std_env_dat0)=='EaGeoMean'|
  names(std_env_dat0)=='MysGeoMean'))
# Spawning mean conductivity and number of spawning days under 9C
take_env_lo_dat<-which(names(std_env_lo_dat)=='MRZ_cond_Spawn_season_mean'|
  names(std_env_lo_dat)=='MRZ_days_under_9deg')
tmbdat$n_beta<-length(take_env_dat)
tmbdat$n_alpha<-length(take_env_lo_dat)
tmbdat$env_dat<-as.matrix(std_env_dat0[,take_env_dat])
tmbdat$env_lo_dat<-as.matrix(std_env_lo_dat[,take_env_lo_dat])

map1<-map
map1$b_hj<-factor(seq(1,tmbdat$n_beta,by=1))
map1$tcrit<-factor(c(1:2))

REs<-c('hy','hc','lng','lna50','lo','lnmuh','mulo','b_ha','b_hj','alpha','tcrit')

for(j in 1:5)
{
  # set function/model objects to NULL
  obj<-NULL

  # create model object
  obj<-MakeADFun(tmbdat,start(),DLL=fil,random=REs,map=map1,
    inner.control=list(maxit=5000))
  obj$report()
  
  if(j==1)
  {
    # fit model, if doesn't converge in alloted iterations, run again
    converge<-FALSE
    iter<-0
    while(converge==FALSE&iter<5)
    {
      fit_lst[[i]]<-nlminb(obj$par,obj$fn,obj$gr,contol=list(eval.max=400,iter.max=300))
      obj_lst[[i]]<-obj
      converge<-ifelse(sum(grep('(4)',fit_lst[[i]]$message))==1,TRUE,FALSE)
      iter<-iter+1
      gc()
    }
  }else{
    fit<-NULL
    fit<-nlminb(obj$par,obj$fn,obj$gr,contol=list(eval.max=400,iter.max=300))
    converge<-ifelse(sum(grep('(8)',fit$message))==1,FALSE,TRUE)
    if(converge==TRUE&fit$objective<fit_lst[[i]]$objective)
    {
      fit_lst[[i]]<-fit
      obj_lst[[i]]<-obj
    }
    gc()
  }
}
nll[i]<-fit_lst[[i]]$objective
sdr_lst[[i]]<-sdreport(obj_lst[[i]])
fit_lst[[i]]

#-------------------------------------------------------------------------------------
nll[c(2,5)]-min(nll[c(2,5)])
# outflow is the best predictor

save.image("./Longfin3.RData")

#-------------------------------------------------------------------------------------


# -------SKIP 6 & 7 ---- memory issues 
#-------------------------------------------------------------------------------------
# compare outflow, conductivity and salinity for growth increment and lo
# - Model 2 is the best with mean abiotic conditions
# - Model 6 is the Ha with min abiotic conditions --- memory issues
# - Model 7 is the Ha with max abiotic conditions --- memory issues
# note, keeping the prey resources at mean because these were collected
# much less regularly (monthly) than the abiotic conditions (daily)
# NOTE 2: frequent problems with convergence onto viable estimates - BAIL ON THIS
#-------------------------------------------------------------------------------------
i=6
# min temperature, min outflow, abundance (otter trawl) for age-0 and age-1+,
# and prey densities (E. afinnis and Mysis)
take_env_dat<-c(which(names(std_env_dat0)=='MRZ_temp_min'),
  which(names(std_env_dat0)=='minout'|
  names(std_env_dat0)=='SFBSot0'|
  names(std_env_dat0)=='SFBSot1p'|
  names(std_env_dat0)=='EaGeoMean'|
  names(std_env_dat0)=='MysGeoMean'))
# Spawning Outflow index (days over 20k cfs) and min spawn season temperature
take_env_lo_dat<-which(names(std_env_lo_dat)=='spawn_days_over_20k'|
  names(std_env_lo_dat)=='MRZ_temp_Spawn_season_min')
tmbdat$n_beta<-length(take_env_dat)
tmbdat$n_alpha<-length(take_env_lo_dat)
tmbdat$env_dat<-as.matrix(std_env_dat0[,take_env_dat])
tmbdat$env_lo_dat<-as.matrix(std_env_lo_dat[,take_env_lo_dat])

map1<-map
map1$b_hj<-factor(seq(1,tmbdat$n_beta,by=1))
map1$tcrit<-factor(c(1:2))

REs<-c('hy','hc','lng','lna50','lo','lnmuh','mulo','b_ha','b_hj','alpha','tcrit')

for(j in 1:5)
{
  # set function/model objects to NULL
  obj<-NULL

  # create model object
  obj<-MakeADFun(tmbdat,start(),DLL=fil,random=REs,map=map1,
    inner.control=list(maxit=5000))
  obj$report()
  
  if(j==1)
  {
    # fit model, if doesn't converge in alloted iterations, run again
    converge<-FALSE
    iter<-0
    while(converge==FALSE&iter<5)
    {
      fit_lst[[i]]<-nlminb(obj$par,obj$fn,obj$gr,contol=list(eval.max=400,iter.max=300))
      obj_lst[[i]]<-obj
      converge<-ifelse(sum(grep('(4)',fit_lst[[i]]$message))==1,TRUE,FALSE)
      iter<-iter+1
      gc()
    }
  }else{
    fit<-NULL
    fit<-nlminb(obj$par,obj$fn,obj$gr,contol=list(eval.max=400,iter.max=300))
    converge<-ifelse(sum(grep('(8)',fit$message))==1,FALSE,TRUE)
    if(converge==TRUE&fit$objective<fit_lst[[i]]$objective)
    {
      fit_lst[[i]]<-fit
      obj_lst[[i]]<-obj
    }
    gc()
  }
}
nll[i]<-fit_lst[[i]]$objective
sdr_lst[[i]]<-sdreport(obj_lst[[i]])
fit_lst[[i]]

#-------------------------------------------------------------------------------------
i=7
# min temperature, min outflow, abundance (otter trawl) for age-0 and age-1+,
# and prey densities (E. afinnis and Mysis)
take_env_dat<-c(which(names(std_env_dat0)=='MRZ_temp_max'),
  which(names(std_env_dat0)=='maxout'|
  names(std_env_dat0)=='SFBSot0'|
  names(std_env_dat0)=='SFBSot1p'|
  names(std_env_dat0)=='EaGeoMean'|
  names(std_env_dat0)=='MysGeoMean'))
# Spawning Outflow index (days over 20k cfs) and min spawn season temperature
take_env_lo_dat<-which(names(std_env_lo_dat)=='spawn_days_over_20k'|
  names(std_env_lo_dat)=='MRZ_temp_Spawn_season_max')
tmbdat$n_beta<-length(take_env_dat)
tmbdat$n_alpha<-length(take_env_lo_dat)
tmbdat$env_dat<-as.matrix(std_env_dat0[,take_env_dat])
tmbdat$env_lo_dat<-as.matrix(std_env_lo_dat[,take_env_lo_dat])

map1<-map
map1$b_hj<-factor(seq(1,tmbdat$n_beta,by=1))
map1$tcrit<-factor(c(1:2))

REs<-c('hy','hc','lng','lna50','lo','lnmuh','mulo','b_ha','b_hj','alpha','tcrit')


for(j in 1:5)
{
  # set function/model objects to NULL
  obj<-NULL

  # create model object
  obj<-MakeADFun(tmbdat,start(),DLL=fil,random=REs,map=map1,
    inner.control=list(maxit=5000))
  obj$report()
  
  if(j==1)
  {
    # fit model, if doesn't converge in alloted iterations, run again
    converge<-FALSE
    iter<-0
    while(converge==FALSE&iter<5)
    {
      fit_lst[[i]]<-nlminb(obj$par,obj$fn,obj$gr,contol=list(eval.max=400,iter.max=300))
      obj_lst[[i]]<-obj
      converge<-ifelse(sum(grep('(4)',fit_lst[[i]]$message))==1,TRUE,FALSE)
      iter<-iter+1
      gc()
    }
  }else{
    fit<-NULL
    fit<-nlminb(obj$par,obj$fn,obj$gr,contol=list(eval.max=400,iter.max=300))
    converge<-ifelse(sum(grep('(8)',fit$message))==1,FALSE,TRUE)
    if(converge==TRUE&fit$objective<fit_lst[[i]]$objective)
    {
      fit_lst[[i]]<-fit
      obj_lst[[i]]<-obj
    }
    gc()
  }
}
nll[i]<-fit_lst[[i]]$objective
sdr_lst[[i]]<-sdreport(obj_lst[[i]])
fit_lst[[i]]

nll[c(2,6,7)]-min(nll[c(2,6,7)])
# mean is the best predictor






# interactions based on model #2

# basic interactions between abiotic (increment and lo) plus #s and time
i=8
# Year-specific RE
# mean temperature * mean outflow, abundance (otter trawl) for age-0 and age-1+ * time,
# and prey densities (E. afinnis and Mysis)
take_env_dat<-c(which(names(std_env_dat0)=='MRZ_temp_mean'),  # temperature^2
  which(names(std_env_dat0)=='meanout'|     # mean outflow
  names(std_env_dat0)=='temp_out'|          # temperature * outflow
  names(std_env_dat0)=='SFBSot0'|           # Age-0 Otter trawl abundance index
  names(std_env_dat0)=='SFBSot1p'|          # age-1+ otter trawl abundance index
  names(std_env_dat0)=='c_month'|           # centered month
  names(std_env_dat0)=='ot_N0_t'|           # Age-0 #s * month
  names(std_env_dat0)=='ot_N1p_t'|          # age-1+ #s * month
  names(std_env_dat0)=='EaGeoMean'|         # E.afinnis abundance index
  names(std_env_dat0)=='MysGeoMean'))       # Mysis abundance index
# Spawning Outflow index (days over 20k cfs) * mean spawn season temperature
take_env_lo_dat<-which(names(std_env_lo_dat)=='spawn_days_over_20k'|
  names(std_env_lo_dat)=='MRZ_temp_Spawn_season_mean'|
  names(std_env_lo_dat)=='out_temp')
tmbdat$n_beta<-length(take_env_dat)
tmbdat$n_alpha<-length(take_env_lo_dat)
tmbdat$env_dat<-as.matrix(std_env_dat0[,take_env_dat])
tmbdat$env_lo_dat<-as.matrix(std_env_lo_dat[,take_env_lo_dat])

map1<-map
map1$b_hj<-factor(seq(1,tmbdat$n_beta,by=1))
map1$tcrit<-factor(c(1:2))

REs<-c('hy','hc','lng','lna50','lo','lnmuh','mulo','b_ha','b_hj','alpha','tcrit')

for(j in 1:5)
{
  # set function/model objects to NULL
  obj<-NULL

  # create model object
  obj<-MakeADFun(tmbdat,start(),DLL=fil,random=REs,map=map1,
    inner.control=list(maxit=5000))
  obj$report()
  
  if(j==1)
  {
    # fit model, if doesn't converge in alloted iterations, run again
    converge<-FALSE
    iter<-0
    while(converge==FALSE&iter<5)
    {
      fit_lst[[i]]<-nlminb(obj$par,obj$fn,obj$gr,contol=list(eval.max=400,iter.max=300))
      obj_lst[[i]]<-obj
      converge<-ifelse(sum(grep('(4)',fit_lst[[i]]$message))==1,TRUE,FALSE)
      iter<-iter+1
      gc()
    }
  }else{
    fit<-NULL
    fit<-nlminb(obj$par,obj$fn,obj$gr,contol=list(eval.max=400,iter.max=300))
    converge<-ifelse(sum(grep('(8)',fit$message))==1,FALSE,TRUE)
    if(converge==TRUE&fit$objective<fit_lst[[i]]$objective)
    {
      fit_lst[[i]]<-fit
      obj_lst[[i]]<-obj
    }
    gc()
  }
}
fit_lst[[i]]
nll[i]<-fit_lst[[i]]$objective
sdr_lst[[i]]<-sdreport(obj_lst[[i]])

save.image("./Longfin8.RData")


fit8_lst<-NULL
obj8_lst<-NULL
sdr8_lst<-NULL

j=1
jj=1
while(j<=10&jj<20)
{
  # set function/model objects to NULL
  obj<-NULL

  # create model object
  obj<-MakeADFun(tmbdat,start(),DLL=fil,random=REs,map=map1,
    inner.control=list(maxit=5000))
  obj$report()
  
  fit<-NULL
  fit<-nlminb(obj$par,obj$fn,obj$gr,contol=list(eval.max=400,iter.max=300))
  converge<-ifelse(sum(grep('(8)',fit$message))==1,FALSE,TRUE)
  if(converge==TRUE&fit$objective<fit_lst[[i]]$objective)
  {
    fit8_lst[[j]]<-fit
    obj8_lst[[j]]<-obj
    sdr8_lst[[j]]<-sdreport(obj8_lst[[j]])
    j=j+1
  }
  jj=jj+1
}

save.image("./Longfin_only8.RData")



# set up to change to linear temperature if needed (or remove temp effect for adults)





# interactions between abiotic (increment and lo), prey * #s (no time)
# memory issues
i=9
# Year-specific RE
# mean temperature * mean outflow, abundance (otter trawl) for age-0 and age-1+ * time,
# and prey densities (E. afinnis and Mysis)
take_env_dat<-c(which(names(std_env_dat0)=='MRZ_temp_mean'),  # temperature^2
  which(names(std_env_dat0)=='meanout'|     # mean outflow
  names(std_env_dat0)=='temp_out'|          # temperature * outflow
  names(std_env_dat0)=='SFBSot0'|           # Age-0 Otter trawl abundance index
  names(std_env_dat0)=='SFBSot1p'|          # age-1+ otter trawl abundance index
  names(std_env_dat0)=='EaGeoMean'|         # E.afinnis abundance index
  names(std_env_dat0)=='MysGeoMean'|        # Mysis abundance index
  names(std_env_dat0)=='ot_N0_Ea'|          # Age-0 #s * E.afinnis abundance index
  names(std_env_dat0)=='ot_N1p_Ea'|         # age-1+ #s * E.afinnis abundance index
  names(std_env_dat0)=='ot_N0_mys'|         # Age-0 #s * Mysis abundance index
  names(std_env_dat0)=='ot_N1p_mys'))       # age-1+ #s * Mysis abundance index
# Spawning Outflow index (days over 20k cfs) * mean spawn season temperature
take_env_lo_dat<-which(names(std_env_lo_dat)=='spawn_days_over_20k'|
  names(std_env_lo_dat)=='MRZ_temp_Spawn_season_mean'|
  names(std_env_lo_dat)=='out_temp')
tmbdat$n_beta<-length(take_env_dat)
tmbdat$n_alpha<-length(take_env_lo_dat)
tmbdat$env_dat<-as.matrix(std_env_dat0[,take_env_dat])
tmbdat$env_lo_dat<-as.matrix(std_env_lo_dat[,take_env_lo_dat])

map1<-map
map1$b_hj<-factor(seq(1,tmbdat$n_beta,by=1))
map1$tcrit<-factor(c(1:2))

REs<-c('hy','hc','lng','lna50','lo','lnmuh','mulo','b_ha','b_hj','alpha','tcrit')

for(j in 1:5)
{
  # set function/model objects to NULL
  obj<-NULL

  # create model object
  obj<-MakeADFun(tmbdat,start(),DLL=fil,random=REs,map=map1,
    inner.control=list(maxit=5000))
  obj$report()
  
  if(j==1)
  {
    # fit model, if doesn't converge in alloted iterations, run again
    converge<-FALSE
    iter<-0
    while(converge==FALSE&iter<5)
    {
      fit_lst[[i]]<-nlminb(obj$par,obj$fn,obj$gr,contol=list(eval.max=400,iter.max=300))
      obj_lst[[i]]<-obj
      converge<-ifelse(sum(grep('(4)',fit_lst[[i]]$message))==1,TRUE,FALSE)
      iter<-iter+1
      gc()
    }
  }else{
    fit<-NULL
    fit<-nlminb(obj$par,obj$fn,obj$gr,contol=list(eval.max=400,iter.max=300))
    converge<-ifelse(sum(grep('(8)',fit$message))==1,FALSE,TRUE)
    if(converge==TRUE&fit$objective<fit_lst[[i]]$objective)
    {
      fit_lst[[i]]<-fit
      obj_lst[[i]]<-obj
    }
    gc()
  }
}
fit_lst[[i]]
nll[i]<-fit_lst[[i]]$objective
sdr_lst[[i]]<-sdreport(obj_lst[[i]])

save.image("./Longfin9.RData")

# interactions between abiotic (increment and lo), prey*#s, and time*#s
# doesn't work, takes way too long
i=10
# Year-specific RE
# mean temperature * mean outflow, abundance (otter trawl) for age-0 and age-1+ * time,
# and prey densities (E. afinnis and Mysis)
take_env_dat<-c(which(names(std_env_dat0)=='MRZ_temp_mean'),  # temperature^2
  which(names(std_env_dat0)=='meanout'|     # mean outflow
  names(std_env_dat0)=='temp_out'|          # temperature * outflow
  names(std_env_dat0)=='SFBSot0'|           # Age-0 Otter trawl abundance index
  names(std_env_dat0)=='SFBSot1p'|          # age-1+ otter trawl abundance index
  names(std_env_dat0)=='c_month'|           # centered month
  names(std_env_dat0)=='ot_N0_t'|           # Age-0 #s * month
  names(std_env_dat0)=='ot_N1p_t'|          # age-1+ #s * month
  names(std_env_dat0)=='EaGeoMean'|         # E.afinnis abundance index
  names(std_env_dat0)=='MysGeoMean'|        # Mysis abundance index
  names(std_env_dat0)=='ot_N0_Ea'|          # Age-0 #s * E.afinnis abundance index
  names(std_env_dat0)=='ot_N1p_Ea'|         # age-1+ #s * E.afinnis abundance index
  names(std_env_dat0)=='ot_N0_mys'|         # Age-0 #s * Mysis abundance index
  names(std_env_dat0)=='ot_N1p_mys'))       # age-1+ #s * Mysis abundance index
# Spawning Outflow index (days over 20k cfs) * mean spawn season temperature
take_env_lo_dat<-which(names(std_env_lo_dat)=='spawn_days_over_20k'|
  names(std_env_lo_dat)=='MRZ_temp_Spawn_season_mean'|
  names(std_env_lo_dat)=='out_temp')
tmbdat$n_beta<-length(take_env_dat)
tmbdat$n_alpha<-length(take_env_lo_dat)
tmbdat$env_dat<-as.matrix(std_env_dat0[,take_env_dat])
tmbdat$env_lo_dat<-as.matrix(std_env_lo_dat[,take_env_lo_dat])

map1<-map
map1$b_hj<-factor(seq(1,tmbdat$n_beta,by=1))
map1$tcrit<-factor(c(1:2))

REs<-c('hy','hc','lng','lna50','lo','lnmuh','mulo','b_ha','b_hj','alpha','tcrit')

for(j in 1:5)
{
  # set function/model objects to NULL
  obj<-NULL

  # create model object
  obj<-MakeADFun(tmbdat,start(),DLL=fil,random=REs,map=map1,
    inner.control=list(maxit=5000))
  obj$report()
  
  if(j==1)
  {
    # fit model, if doesn't converge in alloted iterations, run again
    converge<-FALSE
    iter<-0
    while(converge==FALSE&iter<5)
    {
      fit_lst[[i]]<-nlminb(obj$par,obj$fn,obj$gr,contol=list(eval.max=400,iter.max=300))
      obj_lst[[i]]<-obj
      converge<-ifelse(sum(grep('(4)',fit_lst[[i]]$message))==1,TRUE,FALSE)
      iter<-iter+1

      save.image("../OneDrive - DOI/Longfin/R workspaces/Longfin10.1_02.23.24.RData")

      gc()
    }
  }else{
    fit<-NULL
    fit<-nlminb(obj$par,obj$fn,obj$gr,contol=list(eval.max=400,iter.max=300))
    converge<-ifelse(sum(grep('(8)',fit$message))==1,FALSE,TRUE)
    if(converge==TRUE&fit$objective<fit_lst[[i]]$objective)
    {
      fit_lst[[i]]<-fit
      obj_lst[[i]]<-obj
    }
    
    save.image("./Longfin10.2.RData")
    
    gc()
  }
}
fit_lst[[i]]
nll[i]<-fit_lst[[i]]$objective
sdr_lst[[i]]<-sdreport(obj_lst[[i]])

save.image("./Longfin10.RData")


# interactions between abiotic (increment and lo) with per capita prey (no prey * time interaction)
# DOES NOT WORK
i=11
# Year-specific RE
# mean temperature * mean outflow, abundance (otter trawl) for age-0 and age-1+ * time,
# and prey densities (E. afinnis and Mysis)
take_env_dat<-c(which(names(std_env_dat0)=='MRZ_temp_mean'),  # temperature^2
  which(names(std_env_dat0)=='meanout'|       # mean outflow
  names(std_env_dat0)=='temp_out'|            # temperature * outflow
  names(std_env_dat0)=='pcap_Ea_ot_N0'|       # per capita E.afinnis for age-0
  names(std_env_dat0)=='pcap_Ea_ot_N1p'|      # per capita E.afinnis for age-1+
  names(std_env_dat0)=='pcap_mys_ot_N0'|      # per capita mysis for age-0
  names(std_env_dat0)=='pcap_mys_ot_N1p'))     # per capita mysis for age-1+
# Spawning Outflow index (days over 20k cfs) * mean spawn season temperature
take_env_lo_dat<-which(names(std_env_lo_dat)=='spawn_days_over_20k'|
  names(std_env_lo_dat)=='MRZ_temp_Spawn_season_mean'|
  names(std_env_lo_dat)=='out_temp')
tmbdat$n_beta<-length(take_env_dat)
tmbdat$n_alpha<-length(take_env_lo_dat)
tmbdat$env_dat<-as.matrix(std_env_dat0[,take_env_dat])
tmbdat$env_lo_dat<-as.matrix(std_env_lo_dat[,take_env_lo_dat])

map1<-map
map1$b_hj<-factor(seq(1,tmbdat$n_beta,by=1))
map1$tcrit<-factor(c(1:2))

REs<-c('hy','hc','lng','lna50','lo','lnmuh','mulo','b_ha','b_hj','alpha','tcrit')

for(j in 1:5)
{
  # set function/model objects to NULL
  obj<-NULL

  # create model object
  obj<-MakeADFun(tmbdat,start(),DLL=fil,random=REs,map=map1,
    inner.control=list(maxit=5000))
  obj$report()
  
  if(j==1)
  {
    # fit model, if doesn't converge in alloted iterations, run again
    converge<-FALSE
    iter<-0
    while(converge==FALSE&iter<5)
    {
      fit_lst[[i]]<-nlminb(obj$par,obj$fn,obj$gr,contol=list(eval.max=400,iter.max=300))
      obj_lst[[i]]<-obj
      converge<-ifelse(sum(grep('(4)',fit_lst[[i]]$message))==1,TRUE,FALSE)
      iter<-iter+1
      gc()
    }
  }else{
    fit<-NULL
    fit<-nlminb(obj$par,obj$fn,obj$gr,contol=list(eval.max=400,iter.max=300))
    converge<-ifelse(sum(grep('(8)',fit$message))==1,FALSE,TRUE)
    if(converge==TRUE&fit$objective<fit_lst[[i]]$objective)
    {
      fit_lst[[i]]<-fit
      obj_lst[[i]]<-obj
    }
    gc()
  }
}
fit_lst[[i]]
nll[i]<-fit_lst[[i]]$objective
sdr_lst[[i]]<-sdreport(obj_lst[[i]])



# interactions between abiotic (increment and lo) with per capita prey * time
i=12
# Year-specific RE
# mean temperature * mean outflow, abundance (otter trawl) for age-0 and age-1+ * time,
# and prey densities (E. afinnis and Mysis)
take_env_dat<-c(which(names(std_env_dat0)=='MRZ_temp_mean'),  # temperature^2
  which(names(std_env_dat0)=='meanout'|       # mean outflow
  names(std_env_dat0)=='temp_out'|            # temperature * outflow
  names(std_env_dat0)=='pcap_Ea_ot_N0'|       # per capita E.afinnis for age-0
  names(std_env_dat0)=='pcap_Ea_ot_N1p'|      # per capita E.afinnis for age-1+
  names(std_env_dat0)=='pcap_mys_ot_N0'|      # per capita mysis for age-0
  names(std_env_dat0)=='pcap_mys_ot_N1p'|     # per capita mysis for age-1+
  names(std_env_dat0)=='c_month'|             # centered month
  names(std_env_dat0)=='pcap_Ea_ot_N0_t'|     # Age-0 per capita E.afinnis * month
  names(std_env_dat0)=='pcap_Ea_ot_N1p_t'|    # age-1+ per capita E.afinnis * month
  names(std_env_dat0)=='pcap_mys_ot_N0_t'|    # Age-0 per capita mysis * month
  names(std_env_dat0)=='pcap_mys_ot_N1p_t'))  # age-1+ per capita mysis * month
# Spawning Outflow index (days over 20k cfs) * mean spawn season temperature
take_env_lo_dat<-which(names(std_env_lo_dat)=='spawn_days_over_20k'|
  names(std_env_lo_dat)=='MRZ_temp_Spawn_season_mean'|
  names(std_env_lo_dat)=='out_temp')
tmbdat$n_beta<-length(take_env_dat)
tmbdat$n_alpha<-length(take_env_lo_dat)
tmbdat$env_dat<-as.matrix(std_env_dat0[,take_env_dat])
tmbdat$env_lo_dat<-as.matrix(std_env_lo_dat[,take_env_lo_dat])

map1<-map
map1$b_hj<-factor(seq(1,tmbdat$n_beta,by=1))
map1$tcrit<-factor(c(1:2))

REs<-c('hy','hc','lng','lna50','lo','lnmuh','mulo','b_ha','b_hj','alpha','tcrit')

for(j in 1:5)
{
  # set function/model objects to NULL
  obj<-NULL

  # create model object
  obj<-MakeADFun(tmbdat,start(),DLL=fil,random=REs,map=map1,
    inner.control=list(maxit=5000))
  obj$report()
  
  if(j==1)
  {
    # fit model, if doesn't converge in alloted iterations, run again
    converge<-FALSE
    iter<-0
    while(converge==FALSE&iter<5)
    {
      fit_lst[[i]]<-nlminb(obj$par,obj$fn,obj$gr,contol=list(eval.max=400,iter.max=300))
      obj_lst[[i]]<-obj
      converge<-ifelse(sum(grep('(4)',fit_lst[[i]]$message))==1,TRUE,FALSE)
      iter<-iter+1
      gc()
    }
  }else{
    fit<-NULL
    fit<-nlminb(obj$par,obj$fn,obj$gr,contol=list(eval.max=400,iter.max=300))
    converge<-ifelse(sum(grep('(8)',fit$message))==1,FALSE,TRUE)
    if(converge==TRUE&fit$objective<fit_lst[[i]]$objective)
    {
      fit_lst[[i]]<-fit
      obj_lst[[i]]<-obj
    }
    gc()
  }
}
fit_lst[[i]]
nll[i]<-fit_lst[[i]]$objective
sdr_lst[[i]]<-sdreport(obj_lst[[i]])


#nll[1:5]<-c(985213.9,981953.0,987425.5,982101.1,982510.3)
nll[-c(1,3)]

aic<-nll*2+2*2
daic<-aic-min(aic,na.rm=TRUE)
daic
mdl_take<-which(daic==0)
#mdl_take<-8

summary(sdr_lst[[mdl_take]])


#------------------------------------------------------------------------------------
#
# Step 7: create plots from best model
#
#------------------------------------------------------------------------------------

sdr<-sdr_lst[[mdl_take]]

env_dat_ord<-c('MRZ_temp_mean','meanout','temp_out','SFBSot0','SFBSot1p',
  'c_month','ot_N0_t','ot_N1p_t','EaGeoMean','MysGeoMean')      

take_re<-which(rownames(summary(sdr))=='hy')

take_ln<-which(rownames(summary(sdr))=='lnsighy'|
    rownames(summary(sdr))=='lnsig'|
    rownames(summary(sdr))=='hc'|
    rownames(summary(sdr))=='lng'|
    rownames(summary(sdr))=='lna50')
exp(cbind(summary(sdr)[take_ln,1],
  summary(sdr)[take_ln,1]-1.96*summary(sdr)[take_ln,2],
  summary(sdr)[take_ln,1]+1.96*summary(sdr)[take_ln,2]))

take_raw<-which(rownames(summary(sdr))=='lo'|
    rownames(summary(sdr))=='alpha')
cbind(summary(sdr)[take_raw,1],
  summary(sdr)[take_raw,1]-1.96*summary(sdr)[take_raw,2],
  summary(sdr)[take_raw,1]+1.96*summary(sdr)[take_raw,2])

ord<-match(env_dat_ord,colnames(tmbdat$env_dat))
take_b_ha<-which(rownames(summary(sdr))=='b_ha')[ord]
cbind(summary(sdr)[take_b_ha,1],
  summary(sdr)[take_b_ha,1]-1.96*summary(sdr)[take_b_ha,2],
  summary(sdr)[take_b_ha,1]+1.96*summary(sdr)[take_b_ha,2])

take_b_hj<-which(rownames(summary(sdr))=='b_hj')[ord]
cbind(summary(sdr)[take_b_hj,1],
  summary(sdr)[take_b_hj,1]-1.96*summary(sdr)[take_b_hj,2],
  summary(sdr)[take_b_hj,1]+1.96*summary(sdr)[take_b_hj,2])


take_tcrit<-which(rownames(summary(sdr))=='tcrit')



i=which(names(std_env_dat)=='MRZ_temp_mean')
mean(env_dat[,i])+cbind(summary(sdr)[take_tcrit,1],
  summary(sdr)[take_tcrit,1]-1.96*summary(sdr)[take_tcrit,2],
  summary(sdr)[take_tcrit,1]+1.96*summary(sdr)[take_tcrit,2])



take_env_dat<-c(which(names(std_env_dat0)=='MRZ_temp_mean'),  # temperature^2
  which(names(std_env_dat0)=='meanout'|     # mean outflow
  names(std_env_dat0)=='temp_out'|          # temperature * outflow
  names(std_env_dat0)=='SFBSot0'|           # Age-0 Otter trawl abundance index
  names(std_env_dat0)=='SFBSot1p'|          # age-1+ otter trawl abundance index
  names(std_env_dat0)=='c_month'|           # centered month
  names(std_env_dat0)=='ot_N0_t'|           # Age-0 #s * month
  names(std_env_dat0)=='ot_N1p_t'|          # age-1+ #s * month
  names(std_env_dat0)=='EaGeoMean'|         # E.afinnis abundance index
  names(std_env_dat0)=='MysGeoMean'))       # Mysis abundance index
# Spawning Outflow index (days over 20k cfs) * mean spawn season temperature
take_env_lo_dat<-which(names(std_env_lo_dat)=='spawn_days_over_20k'|
  names(std_env_lo_dat)=='MRZ_temp_Spawn_season_mean'|
  names(std_env_lo_dat)=='out_temp')
tmbdat$n_beta<-length(take_env_dat)
tmbdat$n_alpha<-length(take_env_lo_dat)
tmbdat$env_dat<-as.matrix(std_env_dat0[,take_env_dat])
tmbdat$env_lo_dat<-as.matrix(std_env_lo_dat[,take_env_lo_dat])

map1<-map
map1$b_hj<-factor(seq(1,tmbdat$n_beta,by=1))
map1$tcrit<-factor(c(1:2))

REs<-c('hy','hc','lng','lna50','lo','lnmuh','mulo','b_ha','b_hj','alpha','tcrit')


out<-lpred_fun(sdr,map1)

lpred<-out[[1]]
ltran<-out[[2]]


#---------------------------------------------------------------------------------------------------
# figure 1. Environemtnal effects on growth increment
#---------------------------------------------------------------------------------------------------
# showing maximum growth increment 
b_a<-summary(sdr)[which(rownames(summary(sdr))=="b_ha"),]
b_j<-summary(sdr)[which(rownames(summary(sdr))=="b_hj"),]

plt_xlabs<-c(expression(Temp~(""*degree~C)),'ln(outflow) (cfs)','Temp * ln(outflow)',
  'ln(Age-0 N Index)','ln(Age-1+ N Index)','Month',
  'ln(Age-0 N Index) * Month','ln(Age-1+ N Index) * Month',
  'ln(E. affinis Index)','ln(Mysis Index)')
plt_ylabs<-'Max. Monthly Growth Rate (mm)'
cexl<-0.75
cexa<-1.25
ypadj<-2.5
pbty<-'l'
niter<-1000
n<-100



#---------------------------------------------------------------------------------------------------
# Technical Note Figure 2
#---------------------------------------------------------------------------------------------------
# Prey with abiotoc & interactions
windows(height=9,width=7)
layout(matrix(c(1:4,6,5),nrow=3,byrow=TRUE))
par(mar=c(5,5,0.25,0.25))
plt_ylabs<-'Max. Growth Rate (mm)'
cexl<-1.25
cexa<-1.45
ypadj<-2.5
pbty<-'o'
n<-100

h<-exp(summary(sdr)[which(rownames(summary(sdr))=='hc'),1])     # mean growth
tcrit<-summary(sdr)[which(rownames(summary(sdr))=='tcrit'),]

for(i in 9:10)
{
  j=which(names(std_env_dat[,take_env_dat])==env_dat_ord[i])
  beta_take<-which(colnames(tmbdat$env_dat)==env_dat_ord[i])
  s<-seq(min(std_env_dat[,take_env_dat[j]]),max(std_env_dat[,take_env_dat[j]]),length=n)
  sm<-matrix(s,nrow=length(s),ncol=niter)
  r_s<-s+mean(log(env_dat[,take_env_dat[j]]+1))
  
  pred<-cbind(h*exp(b_a[beta_take,1]*s),
    h*exp((b_a[beta_take,1]+b_j[beta_take,1])*s))
  
  pred_a<-h*exp(rnorm(niter,b_a[beta_take,1],b_a[beta_take,2])*t(sm))
  ci_a<-apply(pred_a,2,quantile,c(0.025,0.975))
  pred_j<-h*exp((rnorm(niter,b_a[beta_take,1],b_a[beta_take,2])+rnorm(niter,b_j[beta_take,1],b_j[beta_take,2]))*t(sm))
  ci_j<-apply(pred_j,2,quantile,c(0.025,0.975))
  
  plot(NA,NA,ylim=range(c(ci_a,ci_j)),xlim=range(r_s),xlab='',ylab='',type='l',lwd=3,las=1,
    cex.lab=cexl,cex.axis=cexa,bty=pbty)#,yaxt='n')
  cordx<-c(r_s,r_s[n:1],r_s[1])
  cordy<-c(ci_a[1,1:n],ci_a[2,n:1],ci_a[1,1])
  polygon(cordx,cordy,col='grey85',lty=0)
  
  cordx<-c(r_s,r_s[n:1],r_s[1])
  cordy<-c(ci_j[1,1:n],ci_j[2,n:1],ci_j[1,1])
  polygon(cordx,cordy,col='grey60',lty=0)

  lines(pred[,1]~r_s,lwd=3,lty=2)
  lines(pred[,2]~r_s,lwd=3,lty=1)
  mtext(plt_xlabs[i],1,cex=cexl,padj=2)#ypadj) 
  #mtext(plt_ylabs,2,cex=cexl,padj=-0.75)  
  mtext(plt_ylabs,2,cex=cexl,padj=-ypadj)  
}
legend(x=-0.2,y=6.051,legend=c(NA,NA),lty=1,lwd=12,pch=NA,cex=1.5,bty='n',col=c('grey60','grey85'))
legend('bottomleft',legend=c('Phase 1','Phase 2'),lty=c(1,2),lwd=3,pch=NA,cex=1.5,bty='n')

#temperature
i=which(names(std_env_dat[,take_env_dat])==env_dat_ord[1])
beta_take<-which(colnames(tmbdat$env_dat)==env_dat_ord[1])
s<-seq(min(std_env_dat[,take_env_dat[i]]),max(std_env_dat[,take_env_dat[i]]),length=n)
sm<-matrix(s,nrow=length(s),ncol=niter)
#sm<-matrix(1,nrow=length(s),ncol=niter)
r_s<-s+mean(env_dat[,take_env_dat[i]])
pred<-cbind(h*exp(b_a[beta_take,1]*((s-tcrit[1,1])^2)),
  h*exp((b_a[beta_take,1]+b_j[beta_take,1])*((s-tcrit[2,1])^2)))

pred_a<-h*exp(rnorm(niter,b_a[beta_take,1],b_a[beta_take,2])*(t(sm)-rnorm(niter,tcrit[1,1],tcrit[1,2]))^2)
ci_a<-apply(pred_a,2,quantile,c(0.025,0.975))
pred_j<-h*exp((rnorm(niter,b_a[beta_take,1],b_a[beta_take,2])+rnorm(niter,b_j[beta_take,1],b_j[beta_take,2]))*
    (t(sm)-rnorm(niter,tcrit[2,1],tcrit[2,2]))^2)
ci_j<-apply(pred_j,2,quantile,c(0.025,0.975))

plot(NA,NA,ylim=range(c(ci_a,ci_j)),xlim=range(r_s),xlab='',ylab='',type='l',lwd=3,las=1,
  cex.lab=cexl,cex.axis=cexa,bty=pbty)

cordx<-c(r_s,r_s[n:1],r_s[1])
cordy<-c(ci_a[1,1:n],ci_a[2,n:1],ci_a[1,1])
polygon(cordx,cordy,col='grey85',lty=0)

cordx<-c(r_s,r_s[n:1],r_s[1])
cordy<-c(ci_j[1,1:n],ci_j[2,n:1],ci_j[1,1])
polygon(cordx,cordy,col='grey60',lty=0)
#polygon(cordx,cordy,col='black',density=20,lwd=1.5,lty=1)

lines(pred[,1]~r_s,lwd=3,lty=2)
lines(pred[,2]~r_s,lwd=3,lty=1)
mtext(plt_xlabs[1],1,cex=cexl,padj=2)#ypadj)  
mtext(plt_ylabs,2,cex=cexl,padj=-ypadj)  

for(i in 2)
{
  j=which(names(std_env_dat[,take_env_dat])==env_dat_ord[i])
  beta_take<-which(colnames(tmbdat$env_dat)==env_dat_ord[i])
  s<-seq(min(std_env_dat[,take_env_dat[j]]),max(std_env_dat[,take_env_dat[j]]),length=n)
  sm<-matrix(s,nrow=length(s),ncol=niter)
  if(j!=6)
  {
    r_s<-s+mean(log(env_dat[,take_env_dat[j]]))
  }else{
    r_s<-s+mean(c(1:12))
  }
  
  pred<-cbind(h*exp(b_a[beta_take,1]*s),
    h*exp((b_a[beta_take,1]+b_j[beta_take,1])*s))
  
  pred_a<-h*exp(rnorm(niter,b_a[beta_take,1],b_a[beta_take,2])*t(sm))
  ci_a<-apply(pred_a,2,quantile,c(0.025,0.975))
  pred_j<-h*exp((rnorm(niter,b_a[beta_take,1],b_a[beta_take,2])+rnorm(niter,b_j[beta_take,1],b_j[beta_take,2]))*t(sm))
  ci_j<-apply(pred_j,2,quantile,c(0.025,0.975))
  
  plot(NA,NA,ylim=range(c(ci_a,ci_j)),xlim=range(r_s),xlab='',ylab='',type='l',lwd=3,las=1,
    cex.lab=cexl,cex.axis=cexa,bty=pbty)
  cordx<-c(r_s,r_s[n:1],r_s[1])
  cordy<-c(ci_a[1,1:n],ci_a[2,n:1],ci_a[1,1])
  polygon(cordx,cordy,col='grey85',lty=0)
  
  cordx<-c(r_s,r_s[n:1],r_s[1])
  cordy<-c(ci_j[1,1:n],ci_j[2,n:1],ci_j[1,1])
  polygon(cordx,cordy,col='grey60',lty=0)

  lines(pred[,1]~r_s,lwd=3,lty=2)
  lines(pred[,2]~r_s,lwd=3,lty=1)
  mtext(plt_xlabs[i],1,cex=cexl,padj=2) 
  mtext(plt_ylabs,2,cex=cexl,padj=-ypadj)  
}
# interactions
for(i in 1)
{
  cov1_take<-which(names(std_env_dat[,take_env_dat])==env_dat_ord[1])
  cov2_take<-which(names(std_env_dat[,take_env_dat])==env_dat_ord[2])
  beta_take<-which(names(std_env_dat[,take_env_dat])==env_dat_ord[3])
  # ignoring 95% CIs on this one
  s1<-seq(min(std_env_dat[,take_env_dat[cov1_take]]),max(std_env_dat[,take_env_dat[cov1_take]]),length=n)
  s2<-seq(min(std_env_dat[,take_env_dat[cov2_take]]),max(std_env_dat[,take_env_dat[cov2_take]]),length=n)
  if(cov1_take==1)
  {
    r_s1<-s1+mean(env_dat[,take_env_dat[cov1_take]])
  }else{
    r_s1<-s1+mean(log(env_dat[,take_env_dat[cov1_take]]))
  }
  r_s2<-s2+mean(log(env_dat[,take_env_dat[cov2_take]]))
  
  sm1<-matrix(s1,nrow=n,ncol=n)
  sm2<-t(matrix(s2,nrow=n,ncol=n))
  
  pred_a<-h*exp(
    b_a[cov1_take,1]*(sm1-tcrit[1,1])^2+      # temperature effect
      b_a[cov2_take,1]*sm2+                                # outflow effect
      b_a[beta_take,1]*(sm1*sm2))              # temp*outflow interactioneffect
  pred_j<-h*exp(
    (b_a[cov1_take,1]+b_j[cov1_take,1])*(sm1-tcrit[2,1])^2+   # temperature effect
      (b_a[cov2_take,1]+b_j[cov2_take,1])*sm2+                                # outflow effect
      (b_a[beta_take,1]+b_j[beta_take,1])*(sm1*sm2))          # temp*outflow interactioneffect
  
  for(j in 1:2)
  {
    if(j==2)
    {
      pred<-pred_j
    }else{
      pred<-pred_a
    }
    
    filled.contour3(x=r_s1,y=r_s2,z=pred,xlab='',ylab='',color.palette=gray.colors,las=1,
      cex.lab=cexl,axes=FALSE)
    contour(x=r_s1,y=r_s2,z=pred,add=TRUE)
    box(which='plot')
    axis(1,cex.axis=cexa)
    axis(2,cex.axis=cexa)
    if(i!=3)
    {
      points(log(env_dat[,take_env_dat[cov2_take]])~env_dat[,take_env_dat[cov1_take]],bg='yellow2',cex=1,pch=21)
    }else{
      points(log(env_dat[,take_env_dat[cov2_take]])~log(env_dat[,take_env_dat[cov1_take]]),bg='yellow2',cex=1,pch=21)
    }
    mtext(plt_xlabs[1],1,cex=cexl,padj=2)#ypadj)  
    mtext(plt_xlabs[2],2,cex=cexl,padj=-1.5,las=3)  
  }
}


#---------------------------------------------------------------------------------------------------
# Technical Note Figure 3
#---------------------------------------------------------------------------------------------------
# biotoc & interactions
windows(height=9,width=7)
layout(matrix(c(rep(c(1,2,3),each=2),rep(c(5,4,7,6),each=3)),nrow=3,byrow=TRUE))
par(mar=c(5,5,0.25,0.25))
plt_ylabs<-'Max. Growth Rate (mm)'
cexl<-1.25
cexa<-1.45
ypadj<-2.5
pbty<-'o'
n<-100

for(i in 4:6)
{
  j=which(names(std_env_dat[,take_env_dat])==env_dat_ord[i])
  beta_take<-which(colnames(tmbdat$env_dat)==env_dat_ord[i])
  s<-seq(min(std_env_dat[,take_env_dat[j]]),max(std_env_dat[,take_env_dat[j]]),length=n)
  sm<-matrix(s,nrow=length(s),ncol=niter)
  if(i!=6)
  {
    r_s<-s+mean(log(env_dat[,take_env_dat[j]]))
  }else{
    r_s<-s+mean(c(1:12))
  }
  
  pred<-cbind(h*exp(b_a[beta_take,1]*s),
    h*exp((b_a[beta_take,1]+b_j[beta_take,1])*s))
  
  pred_a<-h*exp(rnorm(niter,b_a[beta_take,1],b_a[beta_take,2])*t(sm))
  ci_a<-apply(pred_a,2,quantile,c(0.025,0.975))
  pred_j<-h*exp((rnorm(niter,b_a[beta_take,1],b_a[beta_take,2])+rnorm(niter,b_j[beta_take,1],b_j[beta_take,2]))*t(sm))
  ci_j<-apply(pred_j,2,quantile,c(0.025,0.975))
  
  plot(NA,NA,ylim=range(c(ci_a,ci_j)),xlim=range(r_s),xlab='',ylab='',type='l',lwd=3,las=1,
    cex.lab=cexl,cex.axis=cexa,bty=pbty)
  cordx<-c(r_s,r_s[n:1],r_s[1])
  cordy<-c(ci_a[1,1:n],ci_a[2,n:1],ci_a[1,1])
  polygon(cordx,cordy,col='grey85',lty=0)
  
  cordx<-c(r_s,r_s[n:1],r_s[1])
  cordy<-c(ci_j[1,1:n],ci_j[2,n:1],ci_j[1,1])
  polygon(cordx,cordy,col='grey60',lty=0)

  lines(pred[,1]~r_s,lwd=3,lty=2)
  lines(pred[,2]~r_s,lwd=3,lty=1)
  mtext(plt_xlabs[i],1,cex=cexl,padj=2)#ypadj) 
  mtext(plt_ylabs,2,cex=cexl,padj=-ypadj)  
  if(i==6)
  {
    legend(x=3.8,y=23.4,legend=c(NA,NA),lty=1,lwd=12,pch=NA,cex=1.5,bty='n',col=c('grey60','grey85'))
    legend('topright',legend=c('Phase 1','Phase 2'),lty=c(1,2),lwd=3,pch=NA,cex=1.5,bty='n')
  }
}
# interactions
for(i in 4:5)
{
  cov1_take<-which(names(std_env_dat[,take_env_dat])==env_dat_ord[i])
  cov2_take<-which(names(std_env_dat[,take_env_dat])==env_dat_ord[6])
  beta_take<-which(names(std_env_dat[,take_env_dat])==env_dat_ord[i+3])
  # ignoring 95% CIs on this one
  s1<-seq(min(std_env_dat[,take_env_dat[cov1_take]]),max(std_env_dat[,take_env_dat[cov1_take]]),length=n)
  s2<-seq(min(std_env_dat[,take_env_dat[cov2_take]]),max(std_env_dat[,take_env_dat[cov2_take]]),length=n)
  if(cov1_take==1)
  {
    r_s1<-s1+mean(env_dat[,take_env_dat[cov1_take]])
  }else{
    r_s1<-s1+mean(log(env_dat[,take_env_dat[cov1_take]]))
  }
  r_s2<-s2+6.5
  
  sm1<-matrix(s1,nrow=n,ncol=n)
  sm2<-t(matrix(s2,nrow=n,ncol=n))

  pred_a<-h*exp(
    b_a[cov1_take,1]*sm1+                    # NOT temperature effect
    b_a[cov2_take,1]*sm2+                                # outflow effect
    b_a[beta_take,1]*(sm1*sm2))              # temp*outflow interactioneffect
  pred_j<-h*exp(
    (b_a[cov1_take,1]+b_j[cov1_take,1])*sm1+                 # NOT temperature effect
    (b_a[cov2_take,1]+b_j[cov2_take,1])*sm2+                                # outflow effect
    (b_a[beta_take,1]+b_j[beta_take,1])*(sm1*sm2))          # temp*outflow interactioneffect
  
  for(j in 1:2)
  {
    if(j==2)
    {
      pred<-pred_j
    }else{
      pred<-pred_a
    }
    
    filled.contour3(x=r_s2,y=r_s1,z=t(pred),xlab='',ylab='',color.palette=gray.colors,las=1,
      cex.lab=cexl,axes=FALSE)
    contour(x=r_s2,y=r_s1,z=t(pred),add=TRUE)
    box(which='plot')
    axis(1,cex.axis=cexa,at=seq(1,12,by=1),labels=c(1,'',3,'',5,'',7,'',9,'',11,''))
    axis(2,cex.axis=cexa)
    mtext(plt_xlabs[6],1,cex=cexl,padj=2)#ypadj)  
    mtext(plt_xlabs[i],2,cex=cexl,padj=-1.7,las=3)  
  }
}





#---------------------------------------------------------------------------------------------------
# Technical Note Figure 4
#---------------------------------------------------------------------------------------------------
# length-at-age-0 main effects & interactions
windows(height=7,width=6)
layout(matrix(c(1,2,3,3),nrow=2,byrow=TRUE),height=c(1,1.35))
par(mar=c(5,5,0.25,0.25))
plt_xlabs2<-c('ln(Days over 20k cfs)',expression(ln(Mean~Temp.~(""*degree~C))))
plt_ylabs<-'Length-at-age 0 (mm)'
cexl<-1.25
cexa<-1.45
ypadj<-2.5
pbty<-'o'
n<-100

lo<-summary(sdr)[which(rownames(summary(sdr))=='lo'),1]
alpha<-summary(sdr)[which(rownames(summary(sdr))=="alpha"),]

for(i in 1:2)
{
  s<-seq(min(std_env_lo_dat[,take_env_lo_dat[i]]),max(std_env_lo_dat[,take_env_lo_dat[i]]),length=n)
  sm<-matrix(s,nrow=length(s),ncol=niter)
  if(i==2)
  {
    r_s<-s+mean(log(env_lo_dat[,take_env_lo_dat[i]]),na.rm=TRUE)
  }else{
    r_s<-s+mean(log(env_lo_dat[,take_env_lo_dat[i]]))
  }

  pred<-lo+alpha[i,1]*s
  
  pred_a<-lo+rnorm(niter,alpha[i,1],alpha[i,2])*t(sm)
  ci<-apply(pred_a,2,quantile,c(0.025,0.975))

  plot(NA,NA,ylim=range(ci),xlim=range(r_s),xlab='',ylab='',type='l',lwd=3,las=1,
    cex.lab=cexl,cex.axis=cexa,bty=pbty)
  cordx<-c(r_s,r_s[n:1],r_s[1])
  cordy<-c(ci[1,1:n],ci[2,n:1],ci_a[1,1])
  polygon(cordx,cordy,col='grey85',lty=0)
  
  lines(pred~r_s,lwd=3,lty=2)
  mtext(plt_xlabs2[i],1,cex=cexl,padj=2)#ypadj) 
  mtext(plt_ylabs,2,cex=cexl,padj=-ypadj)  
}
# interactions
for(i in 1)
{
  cov1_take<-1
  cov2_take<-2
  beta_take<-3
  # ignoring 95% CIs on this one
  s1<-seq(min(std_env_lo_dat[,take_env_lo_dat[cov1_take]]),max(std_env_lo_dat[,take_env_lo_dat[cov1_take]]),length=n)
  s1<-seq(-1.89,1.09,length=n)
  s2<-seq(min(std_env_lo_dat[,take_env_lo_dat[cov2_take]]),max(std_env_lo_dat[,take_env_lo_dat[cov2_take]]),length=n)
  s2<-seq(-0.142,0.142,length=n)
  r_s1<-s1+mean(log(env_lo_dat[,take_env_lo_dat[cov1_take]]+1))
  r_s2<-s2+mean(log(env_lo_dat[,take_env_lo_dat[cov2_take]]),na.rm=TRUE)
  
  sm1<-matrix(s1,nrow=n,ncol=n)
  sm2<-t(matrix(s2,nrow=n,ncol=n))
  
  pred<-lo+alpha[cov1_take,1]*sm1+alpha[cov2_take,1]*sm2++alpha[beta_take,1]*sm1*sm2
  
  filled.contour3(x=r_s1,y=r_s2,z=pred,xlab='',ylab='',color.palette=gray.colors,las=1,
    cex.lab=cexl,axes=FALSE)
  contour(x=r_s1,y=r_s2,z=pred,add=TRUE)
  box(which='plot')
  axis(1,cex.axis=cexa)
  axis(2,at=seq(2.45,2.7,by=0.05),labels=c('',2.5,'',2.6,'',2.7),cex.axis=cexa)
  
  temp<-log(env_lo_dat[,take_env_lo_dat[cov2_take]])
  temp[which(is.na(temp))]<-mean(log(env_lo_dat[,take_env_lo_dat[cov2_take]]),na.rm=TRUE)
  points(temp~log(env_lo_dat[,take_env_lo_dat[cov1_take]]+1),bg='yellow2',cex=1,pch=21)
  mtext(plt_xlabs2[1],1,cex=cexl,padj=2)#ypadj)  
  mtext(plt_xlabs2[2],2,cex=cexl,padj=-2,las=3)  
}







#---------------------------------------------------------------------------------------------------
# Figure 2.5 in SSA
#---------------------------------------------------------------------------------------------------
# dotplot with main regression parameters
niter<-1000
n<-100

ord<-match(env_dat_ord,colnames(tmbdat$env_dat))

windows(width=9,height=6)
layout(matrix(1:2,nrow=1),width=1.8,1)
par(mar=c(5,13,2,0))

clr<-c('grey','dodgerblue4','lightcoral')

bj_temp<-matrix(rnorm(niter*length(b_a[,1]),b_a[,1],b_a[,2]),nrow=length(b_a[,1]))+
  matrix(rnorm(niter*length(b_j[,1]),b_j[,1],b_j[,2]),nrow=length(b_j[,1]))
bj_mu<-rowMeans(bj_temp)
bj_cis<-t(apply(bj_temp,1,quantile,c(0.025,0.975)))
bj_sig<-rep(1,length(b_j[,1]))
bj_sig[which(bj_cis[,2]<0)]<-2
bj_sig[which(bj_cis[,1]>0)]<-3

ba_cis<-b_a[,1]+b_a[,2]*matrix(c(-1.96,1.96),nrow=length(b_a[,1]),ncol=2,byrow=TRUE)
ba_sig<-rep(1,length(b_a[,1]))
ba_sig[which(ba_cis[,2]<0)]<-2
ba_sig[which(ba_cis[,1]>0)]<-3

plot(-seq(1,length(b_j[,1]),by=1)~bj_mu[ord],pch=19,cex=1.5,xlim=range(bj_cis),col=clr[bj_sig[ord]],
  xlab='',ylab='',yaxt='n',cex.axis=1.2,main='Phase-1 Growth',cex.main=1.5)
abline(v=0,lty=2,lwd=2)
axis(2,at=-(1:length(b_j[ord,1])),labels=plt_xlabs,las=1,cex.axis=1.2)

for(i in 1:length(b_a[,1]))
{
  abline(h=-i,lty=3,lwd=0.3,col='grey')
  lines(-c(i,i)~bj_cis[ord[i],],lwd=3,col=clr[bj_sig[ord[i]]])
}

# phase 2 growth (b_a by itself)
par(mar=c(5,0,2,0.25))

plot(-seq(1,length(b_a[,1]),by=1)~b_a[ord,1],pch=19,cex=1.5,xlim=range(ba_cis)+c(0,0.01),
  col=clr[ba_sig[ord]],yaxt='n',ylab='',xlab='',cex.axis=1.2,main='Phase-2 Growth',cex.main=1.5)
abline(v=0,lty=2,lwd=2)
mtext("Regression Parameter Estimates",1,cex=1.5,padj=3,adj=4)
for(i in 1:length(b_a[,1]))
{
  abline(h=-i,lty=3,lwd=0.3,col='grey')
  lines(-c(i,i)~ba_cis[ord[i],],lwd=3,col=clr[ba_sig[ord[i]]])
}




#---------------------------------------------------------------------------------------------------
# Technical Note Figure 6
#---------------------------------------------------------------------------------------------------
# example fits with residuals
lpred_i_fun(sdr,map=map1,tmbdat,lpred,i=niter) # make sure function is loaded

cexl<-1.5
cexa<-1.5
xadj<- -2.25#2
xpadj<- 2#-2.25
yadj<- -1.1
ypadj<-2.5

windows(height=9,width=7)
layout(matrix(c(1:9,9),ncol=2,byrow=TRUE),height=c(1,1,1,1,0.4))
gcol<-c('lightsalmon','darkorange4','steelblue1','dodgerblue3','darkblue','orange3')
take_crt<-c(7,17,28,39)
for(crt in take_crt)
{
  par(mar=c(2.5,6,1,0.25))
  
  take<-which(dat$ClusterCohort==cohorts[crt])
  plot(dat$Length[take]~dat$FracAgeMonth[take],las=1,xlim=c(0,amax),ylim=c(0,max(dat$Length)),pch=20,
    xlab='',ylab='',col=gcol[dat$Gear[take]],cex=1.25,bty='o',cex.axis=cexa)
  text(4,195,cohorts[crt],cex=1.5)
  if(crt==take_crt[3])mtext("Length (mm)",2,cex=cexl,padj=-2.75,adj=5)
  age<-(0:amax)
  lines(lpred[crt,]~age,lwd=3,col='black')
  
  pred<-NULL
  for(j in 1:length(take))pred[j]<-lpred_i_fun(sdr,map=map1,tmbdat,lpred,i=take[j])
  
  resid<-dat$Length[take]-pred
  plot(resid~dat$FracAgeMonth[take],las=1,xlim=c(0,amax),pch=20,xlab='',ylab='',cex.axis=cexa,
    col=gcol[dat$Gear[take]],cex=1.25,bty='o')
  if(crt==take_crt[4])mtext("Age in Months",1,cex=cexl,padj=ypadj,adj=yadj)
  if(crt==take_crt[3])mtext("Residuals",2,cex=cexl,padj=-2.75,adj=2.5)
  abline(h=0,lwd=3)
}
par(mar=c(0,0,0,0))
plot.new()
legend('bottom',legend=gear_names[c(5,1,4,6,3,2)],pch=20,lwd=NA,col=gcol[c(5,1,4,6,3,2)],ncol=3,cex=1.5,bty='n',pt.cex=3)


# figure showing temporal trends in RE increment, l0, lena-1, lena-2
l0<-matrix(NA,nrow=niter,ncol=n_years)
l1<-matrix(NA,nrow=niter,ncol=n_years)
l2<-matrix(NA,nrow=niter,ncol=n_years)
hy_out<-matrix(NA,nrow=niter,ncol=n_years)
hc_out<-matrix(NA,nrow=niter,ncol=n_years)
for(i in 1:niter)
{
  out<-lpred_var_fun(sdr,map1)
  
  l0[i,]<-out[[1]][,1]
  l1[i,]<-out[[1]][,13]
  l2[i,]<-out[[1]][,25]

  hy_out[i,]<-rnorm(n_years,summary(sdr)[which(rownames(summary(sdr))=='hy'),1][map1$hy],
    summary(sdr)[which(rownames(summary(sdr))=='hy'),2][map1$hy])
  if(sum(map1$hc==1)==n_cohorts)
  {
    hc_out[i,]<-rnorm(1,summary(sdr)[which(rownames(summary(sdr))=='hc'),1][map1$hc],
      summary(sdr)[which(rownames(summary(sdr))=='hy'),2][map1$hc])
  }else{
    hc_out[i,]<-rnorm(n_cohorts,summary(sdr)[which(rownames(summary(sdr))=='hc'),1][map1$hc],
      summary(sdr)[which(rownames(summary(sdr))=='hy'),2][map1$hc])
  }
}
# this loop will generate warnings related rnorm and NA's. THIS IS OKAY
ci_l0<-apply(l0,2,quantile,c(0.025,0.975))
ci_l1<-apply(l1,2,quantile,c(0.025,0.975))
ci_l2<-apply(l2,2,quantile,c(0.025,0.975))

hy_out[which(is.na(hy_out))]<-0
h_mat<-exp(hc_out+hy_out)

ci_h<-apply(h_mat,2,quantile,c(0.025,0.975))


#---------------------------------------------------------------------------------------------------
# Technical Note Figure 5
#---------------------------------------------------------------------------------------------------
windows(height=8,width=7)
layout(matrix(c(1,2,3,4),nrow=4,byrow=TRUE),height=c(1,1,1,1.15))
par(mar=c(2.5,6,0.25,0.25))
plt_xlabs<-NA
plt_ylabs<-c('Max. Growth Rate (mm)    ','Length-at-age 0 (mm)','Length-at-age 1 (mm)','Length-at-age 2 (mm)')
cexl<-0.75
cexa<-1.5
xpadj<- -3.5

plot(NA,NA,ylim=range(ci_h),xlim=range(years),type='l',lwd=3,las=1,
  xlab='',ylab='',cex.lab=cexl,cex.axis=cexa)
axis(1,at=seq(1980,2020,by=5),labels=NA)
mtext(plt_ylabs[1],2,cex.lab=cexl,padj=xpadj)

cordx<-c(years,years[n_years:1],years[1])
cordy<-c(ci_h[1,1:n_years],ci_h[2,n_years:1],ci_h[1,1])
polygon(cordx,cordy,col='grey70',lty=0)

hy_out<-summary(sdr)[which(rownames(summary(sdr))=='hy'),1][map1$hy]
hy_out[is.na(hy_out)]<-0
lines(exp(summary(sdr)[which(rownames(summary(sdr))=='hc'),1][map1$hc[1]])*
    exp(hy_out)~years,lty=2,lwd=3)


# predicted length-at-age-0
plot(NA,NA,ylim=range(ci_l0),xlim=range(years),type='l',lwd=3,
  las=1,xlab='',ylab='',cex.lab=cexl,cex.axis=cexa)
axis(1,at=seq(1980,2020,by=5),labels=NA)
mtext(plt_ylabs[2],2,cex.lab=cexl,padj=xpadj)

cordx<-c(years,years[n_years:1],years[1])
cordy<-c(ci_l0[1,1:n_years],ci_l0[2,n_years:1],ci_l0[1,1])
polygon(cordx,cordy,col='grey70',lty=0)
lines(lpred[,1]~years,type='l',lty=2,lwd=3)


# predicted length-at-age-1
plot(NA,NA,ylim=range(ci_l1),xlim=range(years),type='l',lwd=3,
  las=1,xlab='',ylab='',cex.lab=cexl,cex.axis=cexa)
axis(1,at=seq(1980,2020,by=5),labels=NA)
mtext(plt_ylabs[3],2,cex.lab=cexl,padj=xpadj)

cordx<-c(years,years[n_years:1],years[1])
cordy<-c(ci_l1[1,1:n_years],ci_l1[2,n_years:1],ci_l1[1,1])
polygon(cordx,cordy,col='grey70',lty=0)
lines(lpred[,13]~years,type='l',lty=2,lwd=3)


# predicted length-at-age-2
par(mar=c(5,6,0.25,0.25))
plot(NA,NA,ylim=range(ci_l2),xlim=range(years),type='l',lwd=3,
  las=1,xlab='',ylab='',cex.lab=cexl,cex.axis=cexa)
axis(1,at=seq(1980,2020,by=5),labels=NA)
mtext(plt_ylabs[4],2,cex.lab=cexl,padj=xpadj)
mtext('Year',1,cex.lab=cexl,padj=ypadj)

cordx<-c(years,years[n_years:1],years[1])
cordy<-c(ci_l2[1,1:n_years],ci_l2[2,n_years:1],ci_l2[1,1])
polygon(cordx,cordy,col='grey70',lty=0)
lines(lpred[,25]~years,type='l',lty=2,lwd=3)



#---------------------------------------------------------------------------------------------------
# Technical Note Figure 1
#---------------------------------------------------------------------------------------------------
# time-series covariate data
windows(height=10,width=7)
layout(matrix(c(1:5),nrow=5,byrow=TRUE),height=c(1,1,1,1,1.15))
par(mar=c(2.5,6,0.25,0.25))
plt_xlabs<-NA
plt_ylabs<-c('Temperature (?C)','Outflow (cfs*10,000)','Cond. (uS/cm)','Abundance Index',
  'Prey Index')
abn<-c('Age-0','Age-1+')
prey<-c('E. affinis','Mysis')
clr<-c('black','steelblue')
cexl<-1.2
cexa<-1.5
xpadj<- -3#.5

#temperature
y<-env_dat$Year+(env_dat$month-1)/12
plot(env_dat$MRZ_temp_mean~y,
  ylim=range(env_dat$MRZ_temp_mean),xlim=range(years),type='l',lwd=3,las=1,
  xlab='',ylab='',cex.lab=cexl,cex.axis=cexa,col=clr[1])
axis(1,at=seq(1980,2020,by=5),labels=NA)
mtext(plt_ylabs[1],2,cex=cexl,padj=xpadj)

# outflow
plot(env_dat$meanout/10000~y,
  ylim=range(env_dat$meanout)/10000,xlim=range(years),type='l',lwd=3,
  las=1,xlab='',ylab='',cex.lab=cexl,cex.axis=cexa,col=clr[1])
axis(1,at=seq(1980,2020,by=5),labels=NA)
mtext(plt_ylabs[2],2,cex=cexl,padj=xpadj)

# Conductivity
plot(env_dat$MRZ_cond_mean/10000~y,
  ylim=range(env_dat$MRZ_cond_mean)/10000,xlim=range(years),type='l',lwd=3,
  las=1,xlab='',ylab='',cex.lab=cexl,cex.axis=cexa,col=clr[1])
axis(1,at=seq(1980,2020,by=5),labels=NA)
mtext(plt_ylabs[3],2,cex=cexl,padj=xpadj)

# Abundance indices
plot(env_dat$SFBSot0/max(env_dat$SFBSot0)~env_dat$Year,
  ylim=c(0,1),xlim=range(years),type='l',lwd=3,
  las=1,xlab='',ylab='',cex.lab=cexl,cex.axis=cexa,col=clr[1])
lines(env_dat$SFBSot1p/max(env_dat$SFBSot1p)~env_dat$Year,lwd=3,col=clr[2])
axis(1,at=seq(1980,2020,by=5),labels=NA)
mtext(plt_ylabs[4],2,cex=cexl,padj=xpadj)
legend('topright',legend=abn,lty=1,lwd=3,cex=1.5,bty='n',col=clr,horiz=TRUE)

# Prey availability
par(mar=c(5,6,0.25,0.25))
plot(env_dat$EaGeoMean/max(env_dat$EaGeoMean)~y,
  ylim=c(0,1.03),xlim=range(years),type='l',lwd=3,
  las=1,xlab='',ylab='',cex.lab=cexl,cex.axis=cexa,col=clr[1])
lines(env_dat$MysGeoMean/max(env_dat$MysGeoMean)~y,lwd=3,col=clr[2])
axis(1,at=seq(1980,2020,by=5),labels=NA)
mtext(plt_ylabs[5],2,cex=cexl,padj=xpadj)
mtext('Year',1,cex=cexl,padj=ypadj)
legend('topright',legend=prey,lty=1,lwd=3,cex=1.5,bty='n',col=clr,horiz=TRUE)




##  Model Specifications ----
epi.model='SIRS.bd.climate.stepL' 
lambda.t= 1.06 # inflation factor
obsvar_M = 30 # obs variance multiplier
obsvar_B = 10000 # obs variance baseline
num_ens=10 # number of ensemles
 
## Setup environment ----

# library
library(tgp)
library(MASS)
library(plyr)
library(mvtnorm)
library(truncnorm)
library(tidyverse)
library(data.table)
library(docxtractr)

# directory
dir_home_data="./data/"
dir_home_code="./script/"


### Load Data ----

# ILI+ incidence
N=1e5
load(paste0(dir_home_data,"gamma.RData"))
iliiso = read.csv(paste(dir_home_data,"iliiso_HK_1998_2019.csv",sep='')) 
weeks=iliiso$Week.ending
iliiso=iliiso[,"iliiso.tot"]
iliiso=as.matrix(iliiso,dim(iliiso)[1],dim(iliiso)[2]);
iliiso=iliiso*N;    # scaling
iliiso[597:638,]=iliiso[597:638,]*gamma

## Climate data
load(paste0(dir_home_data,"seasonal_trend.RData"))
da_climate= read.csv(paste0(dir_home_data,"HK_climate_Jan1998_Dec2019.csv")) 
mean.temp.daily = da_climate$temp
sh1.daily = da_climate$sh

## Initial parameter ranges
param_ranges = read.csv(paste(dir_home_data,"param_ranges.csv",sep='')) 
param.bounds = param_ranges %>% filter(model==epi.model) %>% select(-model) %>% column_to_rownames(var="param")

### Load Functions ----
source(paste0(dir_home_code,"functions.R"),local=T)


### Global variables ----

##  DA check parameter boundaries 
Iexp=0.97;
sf_low=.5; sf_high=1.5;
R0max_lower=1.5; R0max_upper=3.5; 
R0diff_lower=0.6; R0diff_upper=1.3;
D_lower=1; D_upper=7;
L_lower=ifelse(grepl("stepL",epi.model), 366, 180);L_upper=750;

DA_low=c(S=0,I=0,newI=0,R0=1.1,R0max=R0max_lower,R0diff=R0diff_lower,
         D=D_lower,L= L_lower,Lshort=60,percLshort=0,SF=sf_low,Iexp=Iexp);
DA_up=c(S=N,I=N*.2,newI=N*.2,R0=3,R0max=R0max_upper,R0diff=R0diff_upper,
        D=D_upper,L=L_upper,Lshort=365,percLshort=0.6,SF=sf_high,Iexp=Iexp)
DA_low.pdm=c(S=0,I=0,newI=0,R0=1.1,R0max=R0max_lower,R0diff=R0diff_lower,
             D=D_lower, L=L_lower,Lshort=60,percLshort=0,SF=sf_low,Iexp=Iexp);
DA_up.pdm=c(S=N,I=N*.2,newI=N*.2,R0=3,R0max=R0max_upper,R0diff=R0diff_upper,
            D=D_upper,L=L_upper,Lshort=365,percLshort=0.6,SF=sf_high,Iexp=Iexp)

## list of parameters for models
allparams= list("SIR"=c("S","I","newI","R0","D","L","SF"),
                "SIRS.bd"=c("S","I","newI","R0","D","L","Iexp","SF"),
                "SIRS.bd.stepL"=c("S","I","newI","R0","D","Lshort","L","percLshort","Iexp","SF"),			
                "SIRS.bd.climate"=c("S","I","newI","R0max","R0diff","D","L","Iexp","SF"),
                "SIRS.bd.climate.stepL"=c("S","I","newI","R0max","R0diff","D","Lshort","L","percLshort","Iexp","SF"),
                "SIRS.bd.climate.normalized"=c("S","I","newI","R0","D","L","Iexp","SF"),
                "SIRS.bd.climate.stepL.normalized"=c("S","I","newI","R0","D","Lshort","L","percLshort","Iexp","SF")
                )  


# other global variables 
discrete=T;
forecast=T
pdmSinit = c(.7, .85) # ~75% susceptible at the start of pandemic; smaller interval to reduce uncertainty
birth.rate.HK = 9.18 / 1000 / 365 # birth rate in HK during 1998-2018, per day
dt=1; 
tmstep=7; # wkly data
nfc=26; # number of weeks to forecast
lambda_range=lambda_range1=lambda_range2=c(0.98,1.02) # applied when lambda is free
lambda_var=.5;
param.names=allparams[[epi.model]]
tm.ini=1;
p=1
obs_i=iliiso*p

# variance of ILI+ data
tmp=rep(0,length(obs_i))
for (i in 4:length(obs_i)){
  tmp[i]=mean(obs_i[(i-3):(i-1)]);
}
obs_vars = obsvar_B+tmp*obsvar_M

# padding for newI
newIpad=iliiso[floor((1:365-1)/7)+1]


# Generate Forecast ----
train.proj = EAKF_AI_contFC_sf_freeR_SR(num_ens, tmstep, epi.model,param.bounds,param.names, obs_i=obs_i, nfc,
                    obs_vars,tm.ini)

      

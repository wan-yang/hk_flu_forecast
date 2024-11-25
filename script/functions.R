# MODE <- function(x) {
#   md=plyr::count(x)
#   mode=md$x[which.max(md$freq)];
#   count=max(md$freq);
#   c(mode,count)
# }

calc_R0_M1_new <- function(in.parms,num.ens, sh1, mean.temp) {
  
  # Create matrices for storing/returning results:
  res.temp = res.temp.red = matrix(0, length(sh1), num.ens)
  
  # Loop through all ensemble members:
  for (ix in 1:num.ens) {
    
    # Assign parameters:
    q.mn <- in.parms[ix, 'qmin']/1000; q.mx <- in.parms[ix, 'qmax']/1000; q.md <- in.parms[ix, 'qmid']/1000
    R0.max <- in.parms[ix, 'R0max']; R0.diff <- in.parms[ix, 'R0diff']
    Tc <- in.parms[ix, 'Tc']; Tdiff <- in.parms[ix, "Tdiff"]; t.exp <- in.parms[ix, 'Texp']
    
    
    if(F){
      if (dim(in.parms)[2] == 9) {
        q.mn.cut <- in.parms[ix, 9]
      } else {
        q.mn.cut <- q.mn
      }
    }
    
    q.mn.cut <- q.mn
    
    # Calculate and correct R0.min
    R0.min <- R0.max - R0.diff
    if (R0.min < 0) {
      R0.min <- 0.1
    }
    Tmin <- Tc - Tdiff
    
    # Calculate parabola params:
    if(F){
      b <- ((R0.max - R0.min) * (q.mx + q.mn)) / ((q.mx - q.md) * (q.mn - q.md))
      a <- (-1 * b) / (q.mx + q.mn)
      c <- R0.min - a * q.md ** 2 - b * q.md
    }
    
    # given the symmetry:
    # for those with values < q.md, use q.mn to determine the corresponding 
    q.mx.left = 2 * q.md - q.mn; 
    b.left <- ((R0.max - R0.min) * (q.mx.left + q.mn)) / ((q.mx.left - q.md) * (q.mn - q.md))
    a.left <- (-1 * b.left) / (q.mx.left + q.mn)
    c.left <- R0.min - a.left * q.md ** 2 - b.left * q.md
    
    # for those with values > q.md, use q.mx to determine the corresponding 
    q.mn.right = 2 * q.md - q.mx
    b.right <- ((R0.max - R0.min) * (q.mx + q.mn.right)) / ((q.mx - q.md) * (q.mn.right - q.md))
    a.right <- (-1 * b.right) / (q.mx + q.mn.right)
    c.right <- R0.min - a.right * q.md ** 2 - b.right * q.md
    
    
    fit1 = fit2 =numeric(length(sh1))
    
    # split the data into two sets (those >=q.md, and those <q.md)
    idx.left = which(sh1 < q.md); idx.right = which(sh1 >= q.md)
    
    # Full model:
    q1.left <- sh1[idx.left]; q1.right = sh1[idx.right]
    t1.left <- mean.temp[idx.left]; t1.right = mean.temp[idx.right]
    fit1[idx.left] <- (a.left * q1.left ** 2 + b.left * q1.left + c.left) * (Tc / t1.left) ** t.exp
    fit1[idx.right] <- (a.right * q1.right ** 2 + b.right * q1.right + c.right) * (Tc / t1.right) ** t.exp
    
    # Reduced model:
    q1 <- sh1
    q1[q1 < q.mn.cut] <- q.mn.cut; q1[q1 > q.mx] <- q.mx
    t1 <- mean.temp; t1[t1 < Tmin] <- Tmin
    q1.left <- q1[idx.left]; q1.right = q1[idx.right]
    t1.left <- t1[idx.left]; t1.right = t1[idx.right]
    fit2[idx.left] <- (a.left * q1.left ** 2 + b.left * q1.left + c.left) * (Tc / t1.left) ** t.exp
    fit2[idx.right] <- (a.right * q1.right ** 2 + b.right * q1.right + c.right) * (Tc / t1.right) ** t.exp
    
    
    # Store results in matrices:
    res.temp[, ix] <- fit1; res.temp.red[, ix] <- fit2
  }
  
  # Return results:
  return(list(res.temp, res.temp.red))
}


SIRS.bd <-function(tm_strt, tm_end, tm_step, states, N, beta, birth.rate = birth.rate.HK,newIpad=NULL, realdata=FALSE){
  # function to integrate to the next time step
  # use SIR model, integrate dicretely with Poisson distributions
  # input: tm_strt: starting time; tm_end: ending time; tm_step: time step
  #         S0, I0: initial states; N: population size
  #         D: infection period, day; L: immune period, day; 
  #         alpha: rate from exposed to infectious; beta: transmission matrix
  # output: S, I for all time steps
  S0=states["S",]; I0=states["I",]; D=states["D",]; L=states["L",]; expoI=states["Iexp",]
  cnt=1;
  # beta stores only data during the time used for the truth
  if(!file.exists('tm.range')) tm.range=tm_strt;
  tm_strt=tm_strt-tm.range[1]+1; # adjust the index to match beta
  tm_end=tm_end-tm.range[1]+1;
  tm_vec=seq(tm_strt,tm_end,by=tm_step)
  tm_sz=length(tm_vec)+1; # including the initial conditions and results of length(tm_vec) integration steps
  Np=length(S0); # number of particles
  S=I=newI=matrix(0,Np,tm_sz)
  S[,1]=S0; I[,1]=I0; # R[,1]=N-S0-I0;
  newI[,1]=0;
 
  if(! exists("discrete")) discrete=FALSE; 
  
  #print(discrete)
  
  if (discrete){
    S[,1]=round(S0,0); I[,1]=round(I0,0); 
    for (t in tm_vec){
      cnt=cnt+1;
      
      # toggle pandemic start:
      if (t == 4240) { 
        S[, cnt - 1] <- round(as.vector(lhs(num_ens, pdmSinit*N)), 0) 
      }
      
      Eimmloss=tm_step*(1/L*(N-S[,cnt-1]-I[,cnt-1]))
      # Einf=tm_step*(beta[t,]*I[,cnt-1]*S[,cnt-1]/N)
      Einf=tm_step*(beta[t,]*pmin(I[,cnt-1], I[,cnt-1]^expoI)*S[,cnt-1]/N)
      Erecov=tm_step*(1/D*I[,cnt-1])
      Eimmloss[Eimmloss<0]=0   # adjust, set <0 to 0
      Einf[Einf<0 | is.na(Einf)]=0
      Erecov[Erecov<0]=0
      smcl=rpois(Np,Eimmloss)
      smci=rpois(Np,Einf)
      smcr=rpois(Np,Erecov)
      
      sk1=smcl-smci
      ik1=smci-smcr
      ik1a=smci
      Ts1=S[,cnt-1]+round(sk1/2,0)
      Ti1=I[,cnt-1]+round(ik1/2,0)
      
      Eimmloss=tm_step*(1/L*(N-Ts1-Ti1))
      Einf=tm_step*(beta[t,]*pmin(Ti1, Ti1^expoI)*Ts1/N)
      Erecov=tm_step*(1/D*Ti1)
      Eimmloss[Eimmloss<0]=0  
      Einf[Einf<0 | is.na(Einf)]=0
      Erecov[Erecov<0]=0
      smcl=rpois(Np,Eimmloss)
      smci=rpois(Np,Einf)
      smcr=rpois(Np,Erecov)
      sk2=smcl-smci
      ik2=smci-smcr
      ik2a=smci;
      Ts2=S[,cnt-1]+round(sk2/2,0)
      Ti2=I[,cnt-1]+round(ik2/2,0)
      
      Eimmloss=tm_step*(1/L*(N-Ts2-Ti2))
      Einf=tm_step*(beta[t,]*pmin(Ti2, Ti2^expoI)*Ts2/N)
      Erecov=tm_step*(1/D*Ti2)
      Eimmloss[Eimmloss<0]=0   
      Einf[Einf<0 | is.na(Einf)]=0
      Erecov[Erecov<0]=0                   
      smcl=rpois(Np,Eimmloss)
      smci=rpois(Np,Einf)
      smcr=rpois(Np,Erecov)
      sk3=smcl-smci
      ik3=smci-smcr
      ik3a=smci;
      Ts3=S[,cnt-1]+round(sk3,0)
      Ti3=I[,cnt-1]+round(ik3,0)
      
      Eimmloss=tm_step*(1/L*(N-Ts3-Ti3))
      Einf=tm_step*(beta[t,]*pmin(Ti3, Ti3^expoI)*Ts3/N)
      Erecov=tm_step*(1/D*Ti3)
      Eimmloss[Eimmloss<0]=0  
      Einf[Einf<0 | is.na(Einf)]=0
      Erecov[Erecov<0]=0
      smcl=rpois(Np,Eimmloss)
      smci=rpois(Np,Einf)
      smcr=rpois(Np,Erecov)
      sk4=smcl-smci
      ik4=smci-smcr
      ik4a=smci;
      
      seed=rpois(Np,.1)
      mu = birth.rate
      S[,cnt]=S[,cnt-1]+round(sk1/6+sk2/3+sk3/3+sk4/6 + mu * N - mu * S[, cnt-1],0)-seed
      I[,cnt]=I[,cnt-1]+round(ik1/6+ik2/3+ik3/3+ik4/6 - mu * I[, cnt-1],0)+seed
      newI[,cnt]=round(newI[,cnt-1]+ik1a/6+ik2a/3+ik3a/3+ik4a/6+seed,0);
    }
  } else {
    # run continuously
    for (t in tm_vec){
      cnt=cnt+1;
      
      # toggle pandemic start:
      if (t == 4240) { 
        S[, cnt - 1] <- as.vector(lhs(num_ens, c(0.6 * N, 0.8 * N)))
      }
      
      Eimmloss=tm_step*(1/L*(N-S[,cnt-1]-I[,cnt-1]))
      Einf=tm_step*(beta[t,]*I[,cnt-1]*S[,cnt-1]/N)
      Erecov=tm_step*(1/D*I[,cnt-1])
      Eimmloss[Eimmloss<0]=0  
      Einf[Einf<0 | is.na(Einf)]=0
      Erecov[Erecov<0]=0
      smcl=Eimmloss
      smci=Einf
      smcr=Erecov
      
      sk1=smcl-smci
      ik1=smci-smcr
      ik1a=smci
      Ts1=S[,cnt-1]+sk1/2
      Ti1=I[,cnt-1]+ik1/2
      
      Eimmloss=tm_step*(1/L*(N-Ts1-Ti1))
      Einf=tm_step*(beta[t,]*Ti1*Ts1/N)
      Erecov=tm_step*(1/D*Ti1)
      Eimmloss[Eimmloss<0]=0
      Einf[Einf<0 | is.na(Einf)]=0
      Erecov[Erecov<0]=0
      smcl=Eimmloss
      smci=Einf
      smcr=Erecov
      sk2=smcl-smci
      ik2=smci-smcr
      ik2a=smci;
      Ts2=S[,cnt-1]+sk2/2
      Ti2=I[,cnt-1]+ik2/2
      
      Eimmloss=tm_step*(1/L*(N-Ts2-Ti2))
      Einf=tm_step*(beta[t,]*Ti2*Ts2/N)
      Erecov=tm_step*(1/D*Ti2)
      Eimmloss[Eimmloss<0]=0
      Einf[Einf<0 | is.na(Einf)]=0
      Erecov[Erecov<0]=0                   
      smcl=Eimmloss
      smci=Einf
      smcr=Erecov
      sk3=smcl-smci
      ik3=smci-smcr
      ik3a=smci;
      Ts3=S[,cnt-1]+sk3
      Ti3=I[,cnt-1]+ik3
      
      Eimmloss=tm_step*(1/L*(N-Ts3-Ti3))
      Einf=tm_step*(beta[t,]*Ti3*Ts3/N)
      Erecov=tm_step*(1/D*Ti3)
      Eimmloss[Eimmloss<0]=0  
      Einf[Einf<0 | is.na(Einf)]=0
      Erecov[Erecov<0]=0
      smcl=Eimmloss
      smci=Einf
      smcr=Erecov
      sk4=smcl-smci
      ik4=smci-smcr
      ik4a=smci;
      
      seed = 1 
      
      mu = birth.rate.HK 
      
      S[,cnt]=S[,cnt-1]+sk1/6+sk2/3+sk3/3+sk4/6-seed + mu * N - mu * S[,cnt-1]
      I[,cnt]=I[,cnt-1]+ik1/6+ik2/3+ik3/3+ik4/6+seed - mu * I[,cnt-1] # natural mortality
      newI[,cnt]=newI[,cnt-1]+ik1a/6+ik2a/3+ik3a/3+ik4a/6+seed;
    }
  }
  
  S=t(S); I=t(I); newI=t(newI);
  if (realdata==FALSE){
    rec=list(S=S,I=I); 
  } else {
    rec=list(S=S,I=I,newI=newI); 
  }
  rec;
}

# stepwise model
SIRS.bd.stepL=function(tm_strt, tm_end, tm_step, states, N, beta,
                       birth.rate = birth.rate.HK, expoI=1, expoS=1, newIpad=newIpad, realdata=FALSE){
  # run stochastically, using random Poisson draws
  # immunity lose modeled as a step function, rather than gradully
  # function to integrate to the next time step
  # use SIR model, integrate dicretely with Poisson distributions
  # input: tm_strt: starting time; tm_end: ending time; tm_step: time step
  #         S0, I0: initial states; N: population size
  #         D: infection period, day; L: immune period, day; 
  #         alpha: rate from exposed to infectious; beta: transmission matrix
  # output: S, I for all time steps
  S0=states["S",]; I0=states["I",]; D=states["D",]; Lshort=states["Lshort",];
  Llong=states["L",];percLshort=states["percLshort",];expoI=states["Iexp",]
  cnt=1;
  # beta stores only data during the time used for the truth
  if(!file.exists('tm.range')) tm.range=tm_strt;
  tm_strt.t=tm_strt-tm.range[1]+1; # adjust the index to match beta
  tm_end.t=tm_end-tm.range[1]+1;
  tm_vec=seq(tm_strt.t,tm_end.t,by=tm_step)
  tm_sz=length(tm_vec)+1; # including the initial conditions and results of length(tm_vec) integration steps
  Np=length(S0); # number of particles
  S=I=newI=newi=matrix(0,Np,tm_sz)
  S[,1]=S0; I[,1]=I0; 
  newI[,1]=0; newi[,1]=0
  if(! is.null(newIpad) & is.null(dim(newIpad))) newIpad=matrix(newIpad,Np,length(newIpad),byrow=T) # if pad newI at the very begining
  # run continuously
  for (t in tm_vec){
    cnt=cnt+1;
    if (F){
      # toggle pandemic start:
      if (t+tm_strt[1]-1 == 4271) {
        # if (t == 4165) {
        if (Np ==1){
          S[, cnt - 1] <- mean(round(as.vector(lhs(num_ens, pdmSinit*N)), 0))
        }else{
          S[, cnt - 1] <- round(as.vector(lhs(num_ens, pdmSinit*N)), 0)
        }
        
      }
    }
    
    Eimmloss=tm_step*(1/Llong*(N-S[,cnt-1]-I[,cnt-1]) ) * (1-percLshort) # long term immunity * (1-percLshort)
    Einf=tm_step*(beta[t,]* pmin(I[,cnt-1], I[,cnt-1]^expoI) * S[,cnt-1]^expoS /N) # raise I to expoI - nonlinear
    Erecov=tm_step*(1/D*I[,cnt-1])
    Eimmloss[Eimmloss<0]=0   
    Einf[Einf<0 | is.na(Einf)]=0
    Erecov[Erecov<0]=0
    smcl=rpois(Np,Eimmloss)
    smci=rpois(Np,Einf)
    smcr=rpois(Np,Erecov)
    
    sk1=smcl-smci
    ik1=smci-smcr
    ik1a=smci;
    Ts1=S[,cnt-1]+sk1/2
    Ti1=I[,cnt-1]+ik1/2
    
    Eimmloss=tm_step*(1/Llong*(N-Ts1-Ti1)) * (1-percLshort)
    Einf=tm_step*(beta[t,]*pmin(Ti1, Ti1^expoI)*Ts1^expoS/N)
    Erecov=tm_step*(1/D*Ti1)
    Eimmloss[Eimmloss<0]=0   
    Einf[Einf<0 | is.na(Einf)]=0
    Erecov[Erecov<0]=0
    smcl=rpois(Np,Eimmloss)
    smci=rpois(Np,Einf)
    smcr=rpois(Np,Erecov)
    sk2=smcl-smci
    ik2=smci-smcr
    ik2a=smci;
    Ts2=S[,cnt-1]+sk2/2
    Ti2=I[,cnt-1]+ik2/2
    
    Eimmloss=tm_step*(1/Llong*(N-Ts2-Ti2)) * (1-percLshort)
    Einf=tm_step*(beta[t,]*pmin(Ti2,Ti2^expoI)*Ts2^expoS/N)
    Erecov=tm_step*(1/D*Ti2)
    Eimmloss[Eimmloss<0]=0   
    Einf[Einf<0 | is.na(Einf)]=0
    Erecov[Erecov<0]=0                   
    smcl=rpois(Np,Eimmloss)
    smci=rpois(Np,Einf)
    smcr=rpois(Np,Erecov)
    sk3=smcl-smci
    ik3=smci-smcr
    ik3a=smci;
    Ts3=S[,cnt-1]+sk3
    Ti3=I[,cnt-1]+ik3
    
    Eimmloss=tm_step*(1/Llong*(N-Ts3-Ti3)) * (1-percLshort)
    Einf=tm_step*(beta[t,]*pmin(Ti3,Ti3^expoI)*Ts3^expoS/N)
    Erecov=tm_step*(1/D*Ti3)
    Eimmloss[Eimmloss<0]=0   
    Einf[Einf<0 | is.na(Einf)]=0
    Erecov[Erecov<0]=0
    smcl=rpois(Np,Eimmloss)
    smci=rpois(Np,Einf)
    smcr=rpois(Np,Erecov)
    sk4=smcl-smci
    ik4=smci-smcr
    ik4a=smci;
    
    seed=rpois(Np,.1); 
    S[,cnt]=round(S[,cnt-1]+sk1/6+sk2/3+sk3/3+sk4/6-seed + birth.rate*(N-S[,cnt-1]),0)
    I[,cnt]=round(I[,cnt-1]+ik1/6+ik2/3+ik3/3+ik4/6+seed - birth.rate*I[,cnt-1],0) # death
    # add immune loss
    Tidx=pmax(0,cnt-round(Lshort,0))
    J=1:Np
    Jnon0=J[Tidx!=0]
    if(any(Jnon0)){
      for(jj in Jnon0){
        S[jj,cnt]=S[jj,cnt] + newi[jj,Tidx[jj]] * percLshort[jj]
      }
    }
    if(!is.null(newIpad)){
      Tidx=cnt-round(Lshort,0)
      J=1:Np
      Jneg=J[Tidx<=0]
      if(any(Jneg)){
        for(jj in Jneg){
          S[jj,cnt]=S[jj,cnt] + newIpad[jj,Tidx[jj]+365] * percLshort[jj]
        }
      }
    }
    S[,cnt]=pmin(S[,cnt],N)
    newI[,cnt]=round(newI[,cnt-1]+ik1a/6+ik2a/3+ik3a/3+ik4a/6+seed,0);
    newi[,cnt]=round(ik1a/6+ik2a/3+ik3a/3+ik4a/6+seed,0); 
  }
  
  S=t(S); I=t(I); newI=t(newI);
  if (realdata==FALSE){
    rec=list(S=S,I=I); 
  } else {
    rec=list(S=S,I=I,newI=newI); 
  }
  rec;
}

Fn_checkDA<-function(xnew,bound.low,bound.up){
  b.low=bound.low;
  b.up=bound.up;
  n.var=nrow(xnew); n.ens=ncol(xnew);
  for(vi in 1:n.var){
    #  Corrects if <b.low
    ug=min(xnew[vi,]);
    if (ug<b.low[vi]){  
      for (jj in 1:n.ens){
        if (xnew[vi,jj]<b.low[vi]){
          xnew[vi,jj]=pmax(b.low[vi],runif(1, min=pmax(b.low[vi],quantile(xnew[vi,],.25)), max=pmax(b.low[vi] * 1.5, quantile(xnew[vi,],.75))));
        }
      }
    }
    ug=max(xnew[vi,]);
    if (ug>b.up[vi]){  
      for (jj in 1:n.ens){
        if (xnew[vi,jj]>b.up[vi]){
          xnew[vi,jj]=runif(1, min=min(b.up[vi]/2,quantile(xnew[vi,],.5)), max=min(quantile(xnew[vi,],.80),b.up[vi]));
        }
      }
    }
  }
  xnew;
}

# function to generate evaluation of the forecast
if(T){
  bins.case = c(seq(0, 1, by = .05) / 100 * N, N)
  bins.hosp =  c(seq(0, 2, by = .1) / 100 /20 * N, N) # ~ 1/20
  bins.death = c(seq(0, 2, by = .1) / 100 / 100 * N, N)
  
  bins.tot.case = c(seq(0, N, 2500)[1:21],N)
  bins.tot.hosp = c(seq(0, 10, by = 2)/100 / 50, seq(15, 50, by = 5)/100 / 50, 1) * N # for the cumulative
  bins.tot.death = c(seq(0, 10, by = 2)/100 / 100, seq(15, 50, by = 5)/100 / 100, 1) * N # for the cumulative
}


fn_getProbDist = function(fcast.case, 
                          fcast.hosp = NULL, 
                          fcast.death = NULL){
  num_wk.fast = fcast.case %>% nrow
  num_ens = fcast.case %>% ncol
  
  meas = c('case')  # , 'death'
  
  if(! is.null(fcast.hosp)){
    meas = c(meas,  'hosp')
  }
  
  if(! is.null(fcast.death)){
    meas = c(meas, 'death')
  }
  
  # find prob. distributions for peak weeks, peak intensities and forecasts for the next n weeks
  ProbDist = NULL
  for(ftype in meas){
    
    fcast.t = get(paste0('fcast.', ftype))
    bins.t = get(paste0('bins.', ftype))
    bins.tot.t = get(paste0('bins.tot.', ftype))
    
    # seasonal target over entire forecast period - in our case 6 months
    peakWeeks.t1= fcast.t %>% apply(2, which.max)
    peakIntensities.t1 = fcast.t %>% apply(2, max)
    totals.t1 = fcast.t %>% apply(2, sum) # sum over the fcast period
    # seasonal target over half forecast period - in our case 3 months
    peakWeeks.t2= as.matrix(fcast.t[1:(nrow(fcast.t)/2),]) %>% apply(2, which.max)
    peakIntensities.t2 = as.matrix(fcast.t[1:(nrow(fcast.t)/2),]) %>% apply(2, max)
    totals.t2 = as.matrix(fcast.t[1:(nrow(fcast.t)/2),]) %>% apply(2, sum) # sum over the fcast period
    
    peakWeeksDist.t1 = matrix(NA, nrow=length(unique(peakWeeks.t1)), ncol=3)
    row=1
    for(i in sort(unique(peakWeeks.t1))){
      peakWeeksDist.t1[row, 1] = i;
      peakWeeksDist.t1[row, 2] = i;
      peakWeeksDist.t1[row, 3] = round(length(peakWeeks.t1[peakWeeks.t1==i])/length(peakWeeks.t1), 4)
      row=row+1
    } 
    
    peakIntensitiesDist.t1  = matrix(NA, nrow=length(bins.t)-1, ncol=3)
    row=1
    for(i in 2:length(bins.t)){
      peakIntensitiesDist.t1 [row, 1] = bins.t[i-1]
      peakIntensitiesDist.t1 [row, 2] = bins.t[i]
      peakIntensitiesDist.t1 [row, 3] = round(length(peakIntensities.t1[peakIntensities.t1 >= bins.t[i-1] & peakIntensities.t1 < bins.t[i]])/length(peakIntensities.t1), 4)
      row=row+1    
    }
    
    totalsDist.t1  = matrix(NA, nrow=length(bins.tot.t)-1, ncol=3)
    row=1
    for(i in 2:length(bins.tot.t)){
      totalsDist.t1 [row, 1] = bins.tot.t[i-1]
      totalsDist.t1 [row, 2] = bins.tot.t[i]
      totalsDist.t1 [row, 3] = round(length(totals.t1[totals.t1 >= bins.tot.t[i-1] & totals.t1 < bins.tot.t[i]])/length(totals.t1), 4)
      row=row+1    
    }
    
    peakWeeksDist.t2 = matrix(NA, nrow=length(unique(peakWeeks.t2)), ncol=3)
    row=1
    for(i in sort(unique(peakWeeks.t2))){
      peakWeeksDist.t2[row, 1] = i;
      peakWeeksDist.t2[row, 2] = i;
      peakWeeksDist.t2[row, 3] = round(length(peakWeeks.t2[peakWeeks.t2==i])/length(peakWeeks.t2), 4)
      row=row+1
    } 
    
    peakIntensitiesDist.t2  = matrix(NA, nrow=length(bins.t)-1, ncol=3)
    row=1
    for(i in 2:length(bins.t)){
      peakIntensitiesDist.t2 [row, 1] = bins.t[i-1]
      peakIntensitiesDist.t2 [row, 2] = bins.t[i]
      peakIntensitiesDist.t2 [row, 3] = round(length(peakIntensities.t2[peakIntensities.t2 >= bins.t[i-1] & peakIntensities.t2 < bins.t[i]])/length(peakIntensities.t2), 4)
      row=row+1    
    }
    
    totalsDist.t2  = matrix(NA, nrow=length(bins.tot.t)-1, ncol=3)
    row=1
    for(i in 2:length(bins.tot.t)){
      totalsDist.t2 [row, 1] = bins.tot.t[i-1]
      totalsDist.t2 [row, 2] = bins.tot.t[i]
      totalsDist.t2 [row, 3] = round(length(totals.t2[totals.t2 >= bins.tot.t[i-1] & totals.t2 < bins.tot.t[i]])/length(totals.t2), 4)
      row=row+1    
    }
    
    # calculate prob. distribution for next n week forecasts
    nextDist.t = matrix(NA, nrow=length(bins.t)-1, ncol=nrow(fcast.t)+2)
    row=1
    for(i in 2:length(bins.t)){
      nextDist.t[row, 1] = bins.t[i-1] # lower bound
      nextDist.t[row, 2] = bins.t[i] # upper bound
      for(j in 1:nrow(fcast.t)){
        values = fcast.t[j,]
        values = values[!is.na(values)]
        nextDist.t[row, j+2] = round(length(values[values >= bins.t[i-1] & values < bins.t[i]]) / length(values), 4)
      }
      row=row+1
    }
    
    colnames(peakWeeksDist.t1) = colnames(peakWeeksDist.t2) = c('lwr','upr', 'prob')
    colnames(peakIntensitiesDist.t1) =colnames(peakIntensitiesDist.t2) = c('lwr','upr', 'prob')
    colnames(totalsDist.t1) =colnames(totalsDist.t2) = c('lwr','upr', 'prob')
    colnames(nextDist.t) = c('lwr','upr', paste0('w', 1:num_wk.fast))
    nextDist.t = melt(nextDist.t %>% data.table(), id.vars = c('lwr','upr')) %>% setnames(c('variable','value'),c('target','prob'))
    nextDist.t$target = factor(nextDist.t$target, levels = paste0('w', 1:num_wk.fast), labels = paste0(1:num_wk.fast,'wk ahead'))
    peakWeeksDist.t1 = peakWeeksDist.t1 %>% data.table()
    peakIntensitiesDist.t1 = peakIntensitiesDist.t1 %>% data.table()
    totalsDist.t1 = totalsDist.t1 %>% data.table()
    peakWeeksDist.t2 = peakWeeksDist.t2 %>% data.table()
    peakIntensitiesDist.t2 = peakIntensitiesDist.t2 %>% data.table()
    totalsDist.t2 = totalsDist.t2 %>% data.table()
    peakWeeksDist.t1$target = 'peak week-6 months'
    peakIntensitiesDist.t1$target = 'peak intensity-6 months'
    totalsDist.t1$target = 'total-6 months'
    peakWeeksDist.t2$target = 'peak week-3 months'
    peakIntensitiesDist.t2$target = 'peak intensity-3 months'
    totalsDist.t2$target = 'total-3 months'
    
    ProbDist = rbind(ProbDist, data.table(data.type = ftype, rbind(peakWeeksDist.t1, peakIntensitiesDist.t1, totalsDist.t1, 
                                                                   peakWeeksDist.t2, peakIntensitiesDist.t2, totalsDist.t2, nextDist.t)))
  }
  
  ProbDist
}

fn_eval = function(fcast.t, Week.fcast.t, Week.start.t){
  # fcast.t is a matrix: time x ens
  # first, get summary stats 
  prob_vec = c(.5, .25, .75, .01, .99, .025, .975,
               .05, .95, .1, .9, .15, .85, 
               .2, .8, .25, .75, .3, .7, 
               .35, .65, .4, .6, .45, .55)
  prob.name_vec = c('median', 'iqr.lwr','iqr.upr','ci98.lwr','ci98.upr','ci95.lwr','ci95.upr',
                    'ci90.lwr','ci90.upr','ci80.lwr','ci80.upr','ci70.lwr','ci70.upr',
                    'ci60.lwr','ci60.upr','ci50.lwr','ci50.upr','ci40.lwr','ci40.upr',
                    'ci30.lwr','ci30.upr','ci20.lwr','ci20.upr','ci10.lwr','ci10.upr')
  tmp.inf = fcast.t %>% apply(1, quantile, prob = prob_vec) %>% t
  colnames(tmp.inf) = prob.name_vec
  
  # seasonal target forecast (13 weeks and 26 weeks)
  periods = c(13,26) 
  fcast_stats_ss = periods %>% lapply(function(p.t){
    p.t = min(p.t,length(Week.start.t))
    Week.start.tt = Week.start.t[1:p.t]
    periods_name = c("13" = "3 months", "26"= "6 months")
    true_peak = obs_i[Week.start.tt] %>% which.max()
    peak_wk = fcast.t[1:p.t,] %>% matrix(ncol=ncol(fcast.t)) %>% apply(2,which.max)
    peak_wk_diff = abs(peak_wk-true_peak) %>% quantile(pro=prob_vec) %>% as.matrix() %>% t()
    peak_ili= 1:ncol(fcast.t) %>% lapply(function(x){
      peak_wk.t=peak_wk[x]
      fcast.t[peak_wk.t,x]
    }) %>% 
      unlist() %>% 
      quantile(prob=prob_vec) %>% as.matrix() %>% t()
    tot_ili = fcast.t[1:p.t,] %>% as.matrix() %>% apply(2,sum) %>% quantile(prob=prob_vec) %>% as.matrix() %>% t()
    colnames(peak_wk_diff) = prob.name_vec;colnames(peak_ili)= prob.name_vec; colnames(tot_ili)= prob.name_vec
    d.peak_wk = data.table(Week.fcast = Week.fcast.t[1],
                           target = paste("peak week",periods_name[as.character(p.t)] ,sep = "-"),
                           peak_wk_diff)
    d.peak_ili = data.table(Week.fcast = Week.fcast.t[1],
                            target = paste("peak intensity",periods_name[as.character(p.t)] ,sep = "-"),
                            peak_ili)
    d.tot_ili = data.table(Week.fcast = Week.fcast.t[1],
                           target = paste("total",periods_name[as.character(p.t)] ,sep = "-"),
                           tot_ili)
    rbind(d.peak_wk,d.tot_ili,d.peak_ili)
  }) %>% rbindlist()
  
  fcast_stats = rbind(data.table(Week.fcast = Week.fcast.t, 
                                 Week.start = Week.start.t,
                                 tmp.inf))
  
  # get the distribution fo different targets
  fcastDist = fn_getProbDist(fcast.t)
  
  return(list(fcast_stats = fcast_stats, fcastDist = fcastDist, fcast_stats_ss = fcast_stats_ss))
}

EAKF_AI_contFC_sf_freeR_SR<-function(num_ens,tmstep,epi.model,param.bounds,param.names,
                                     obs_i=obs_i,nfc=26,obs_vars,tm.ini=1,
                                     fcast.deflat = c(.95,.9) 
){
  ## R Function to run Ensemble Adjustment Kalman Filter to forecast flu time series
  ## The filter is run continuously, i.e., from the begining to the end, across season
  ## developed to forecast flu in Hong Kong, where there are no regular season
  if (grepl("stepL",epi.model) & grepl("SIRS",epi.model)){
    fn_epi=get("SIRS.bd.stepL");
  }else if (grepl("SIRS",epi.model)){
    fn_epi=get("SIRS.bd");
  }
  
  num_times=length(obs_i); num_obs=1;
  num_var= length(param.names); 
  cnt_low=rep(0,2); # counter to track the # weeks each strain has low activity
  reinit_cnt=NULL
  
  
  # initalize parameters
  So=matrix(0,num_var,num_ens,dimnames = list(c(param.names))); 
  xprior=array(0,c(num_var,num_ens,num_times+1),dimnames = list(c(param.names)));
  xpost=array(0,c(num_var,num_ens,num_times),dimnames = list(c(param.names)));
  fcast=fc.metrics=NULL;
  
  # 1/31/23
  fcast_mn_stats = fcast_mn_Dist = fcast_mn_stats_ss = NULL # to save the forecasts using the mean ensemble estimates
  fcast_ens_stats = fcast_ens_Dist = fcast_ens_stats_ss = NULL 
  fcast_mn_stats_sort = fcast_mn_Dist_sort = fcast_mn_stats_ss_sort = NULL
  fcast_ens_stats_sort = fcast_ens_Dist_sort = fcast_ens_stats_ss_sort = NULL
  fcast_ens_all = array(0,c(26,num_ens,num_times))
  
  for (i in 1:length(fcast.deflat)){
    assign(paste0("fcast_ens.df",i,"_stats"), NULL)
    assign(paste0("fcast_ens.df",i,"_stats_ss"), NULL)
    assign(paste0("fcast_ens.df",i,"_Dist"), NULL) 
    assign(paste0("fcast_ens.df",i,"_stats_sort"), NULL)
    assign(paste0("fcast_ens.df",i,"_stats_ss_sort"), NULL)
    assign(paste0("fcast_ens.df",i,"_Dist_sort"), NULL) 
    assign(paste0("fcast_ens.df",i,"_all"), array(0,c(26,num_ens,num_times))) 
  }
  
  
  
  
  paramsEns=lhs(num_ens,as.matrix(param.bounds));
  colnames(paramsEns)=rownames(param.bounds)
  paramsEns[,"SF"]=paramsEns[,"SF"]+.25
  So[param.names[!param.names %in% c("newI")],]=t(paramsEns[,param.names[!param.names %in% c("newI")]]);
  
  if (grepl("climate",epi.model)){
    if (grepl("normalized",epi.model)){
      R0 = paramsEns[,"R0"]
      beta=R0*matrix(sn_trend[1:tmstep],tmstep,num_ens,byrow = F)/matrix(So["D",],tmstep,num_ens,byrow=TRUE);
      
    }else{
      R0 = calc_R0_M1_new(paramsEns, num_ens, sh1.daily[1:tmstep], 
                          mean.temp.daily[1:tmstep])[[2]]
      beta=R0/matrix(So["D",],tmstep,num_ens,byrow=TRUE);
    }
  }else{
    R0 = paramsEns[,"R0"]
    beta = matrix(R0/So["D",],tmstep,num_ens,byrow=TRUE)
  }
  
  tcurrent=tm.ini; 
  
  Sr_tmp=fn_epi(tm_strt=tcurrent+dt, tm_end=tcurrent+tmstep, tm_step=dt,
                states=So, N, beta=beta,newIpad=newIpad, realdata=T)
  xprior["S",,1]=tail(Sr_tmp$S,1);
  xprior["I",,1]=tail(Sr_tmp$I,1);
  xprior["newI",,1]=tail(Sr_tmp$newI,1)*So["SF",]; # times the scaling factor
  xprior[param.names[!param.names %in% c("S","I","newI")],,1]=
    So[param.names[!param.names %in% c("S","I","newI")],];
  
  
  #### Begin looping through observations
  #### Training process
  # prior for lambda
  if (is.numeric(lambda.t)){
    lambda=rep(lambda.t^2,num_var)
  }else{
    lambda=rep(1.03^2,num_var) 
  }
  lambdas=lambda_vars=matrix(0,num_times,num_var); 
  ## maximize this funtion (posterior of lambda):
  fit.lambda=function(xl,rel,lambda_prior,lambda_var,prior_var,obs_var,prior_mean,obs){
    sqrt(2*pi*((1+rel*(sqrt(xl)-1))^2*prior_var+obs_var))^(-1)*
      exp(-(prior_mean-obs)^2/2/((1+rel*(sqrt(xl)-1))^2*prior_var+obs_var))*dnorm(xl,lambda_prior,sqrt(lambda_var));
  }
  ## filter start
  for (tt in 1:num_times){
    print(tt);
    xmn=rowMeans(xprior[,,tt]);
    
    # inflate all states and parameters, except the observed (i.e., newI)
    lambda[3]=1
    xnew0 = xprior[,,tt]
    xnew=diag(sqrt(lambda),num_var,num_var)%*%(xprior[,,tt]-xmn%*%matrix(1,1,num_ens))+
      xmn%*%matrix(1,1,num_ens)
    rownames(xnew) = c(param.names)
    
    ####  Get the variance of the ensemble
    obs_var = obs_vars[tt]
    
    prior_var = var(xnew["newI",]);
    
    post_var = prior_var*obs_var/(prior_var+obs_var);
    
    prior_mean = mean(xnew["newI",]);
    post_mean = post_var*(prior_mean/prior_var + obs_i[tt]/obs_var);
    
    #### Compute alpha and adjust distribution to conform to posterior moments
    alp = sqrt(obs_var/(obs_var+prior_var));
    
    dy = post_mean + alp*((xnew["newI",])-prior_mean)-xnew["newI",];
    
    ### RECORD THE TREND OF ADJUSTMENT TO INFORM VALUE FOR REINITIALIZATION 
    if(!exists('adj.trend')) adj.trend=NULL;
    if(mean(dy)<0) {
      adj.trend=append(adj.trend,-1)
    }else {
      adj.trend=append(adj.trend,1)
    }
    
    ###  Getting the covariance of the prior state space and
    ###  observations  (which could be part of state space, e.g. infections)
    rr=NULL;
    for (j in 1:dim(xnew)[1]){
      C=cov(xnew[j,],xnew["newI",])/prior_var;  # covariance/varance of x.obs
      rr=append(rr,C);
    }
    dx=rr%*%t(dy);
    
    ###  Get the new ensemble and save prior and posterior
    xnew = xnew + dx;
    
    ## updata lambda
    ## per Anderson 2009 Tellus 61A (2009),1 pp76
    if(lambda.t == "freelambda"){
      for(j in 1:(num_var)){
        lambda_prior=lambda[j];
        # find the maximum of the lambda posterior density
        lambda[j]=optimize(f=fit.lambda,rel=rr[j],lambda_prior=lambda_prior,lambda_var=lambda_var,
                           prior_var=prior_var,obs_var=obs_var,prior_mean=prior_mean,
                           obs=obs_i[tt],interval=lambda_range,maximum=TRUE,tol=.0001)$maximum;
      }
      ## save lambda and lambda_var
      lambdas[tt,]=lambda; lambda_vars[tt,]=lambda_var;
    }
    
    
    if(tt %in% 595:620){  # pandemic 
      xnew=Fn_checkDA(xnew,bound.low=DA_low.pdm[param.names],bound.up=DA_up.pdm[param.names]); 
    } else {
      xnew=Fn_checkDA(xnew,bound.low=DA_low[param.names],bound.up=DA_up[param.names]);
    }
    xpost[,,tt]=xnew;
    if (is.null(dim(newIpad))){
      newIpad=matrix(newIpad,num_ens,length(newIpad),byrow=T)
    }
    # update newIpad 
    newIpad[,1:(365-tmstep)] = newIpad[,(tmstep+1):365]
    newIpad[,(365-tmstep+1):365]=matrix(xpost["newI",,tt]/tmstep,num_ens,tmstep,byrow=F)
    
    ## Reinitialize if nesserary BEFORE forecast
    xprior[,,tt+1]=xnew;
    if(obs_i[tt]>quantile(xpost["newI",,tt],.975) | obs_i[tt]<quantile(xpost["newI",,tt],.025)){
      if(cnt_low[1]==0){
        cnt_low[1]=tt;
        cnt_low[2]=1;
      } else {
        if (cnt_low[1]==tt-1){ # continuously diverging
          cnt_low[1]=tt;
          cnt_low[2]=cnt_low[2]+1;
        } else {
          cnt_low[1]=tt;
          cnt_low[2]=1;
        }
      }
    }
    
    ## reinitialization
    if(cnt_low[2]>=2){
      if(T){
        print(sum(xpost[,,tt]==xprior[,,tt+1]))
        idx.t = sample(1:num_ens, size = round(num_ens * .5, 0), replace = F)
        
        paramsEns[,param.names[param.names!="newI"]]=lhs(num_ens,as.matrix(param.bounds[param.names[param.names!="newI"],]));
        
        if(sum(tail(adj.trend,min(length(adj.trend),2)))<(-1)){ # adjust S downward, b/c prior>obs
          
          if(tt %in% 595:620){  # pandemic
            xprior["S",,tt+1]=runif(num_ens,0.5*N,0.6*N); # S
          } else {
            xprior["S",,tt+1]= runif(num_ens,min(xpost["S",,tt])*.9,quantile(xpost["S",,tt],prob=.8))
          }
        } else if (sum(tail(adj.trend,min(length(adj.trend),2)))>1){
          
          if(tt %in% 595:620){  # pandemic
            xprior["S",,tt+1]=runif(num_ens,0.70*N,.85*N); # S
            
          } else {
            xprior["S",,tt+1]=runif(num_ens,quantile(xpost["S",,tt],prob=.1), pmin(N, max(xpost["S",,tt])*1.1));
            
          }
        } else {
          xprior["S",idx.t,tt+1]=t(paramsEns[idx.t,"S"])
          
        }
        
        param.t = allparams[[epi.model]]
        param.t = param.t[!param.t %in% c("S","newI","Iexp")]
        xprior[param.t,idx.t,tt+1]=t(paramsEns[idx.t,param.t])
        xprior["I",idx.t,tt+1]=rep(obs_i[tt],length(idx.t))
        
      }
      
      cnt_low[1]=0; # reset the counter
      cnt_low[2]=0; # reset the counter
      
      reinit_cnt = c(reinit_cnt,tt)
      print(c(tt,' reinitialized'),quote=F);
      reinit=T;
      if(!exists('cnt.reinit')) cnt.reinit=0;
      cnt.reinit=cnt.reinit+1;
    }
    
    ## do a forecast for the following nfc weeks before the SR 
    if (forecast){
      if(tt>=3){
        tcurrent = tm.ini+tmstep*tt; 
        update_params=intersect(colnames(paramsEns),rownames(xprior))
        paramsEns[,update_params]=xprior[update_params,,tt+1] %>% t
        if (grepl("climate",epi.model)){
          if (grepl("normalized",epi.model)){
            R0 = paramsEns[,"R0"]
            beta1 = R0*matrix(sn_trend[tcurrent:(tcurrent+tmstep*nfc)],tmstep*(nfc)+tm.ini,num_ens,byrow = F)/matrix(xprior["D",,tt+1],tmstep*(nfc)+tm.ini,num_ens,byrow=TRUE);
            beta2 = matrix(mean(R0)*sn_trend[tcurrent:(tcurrent+tmstep*nfc)]/mean(xprior["D",,tt+1]),tmstep*(nfc)+1,1)
            
          }else{
            R0 = calc_R0_M1_new(paramsEns, num_ens, sh1.daily[tcurrent:(tcurrent+tmstep*nfc)], 
                                mean.temp.daily[tcurrent:(tcurrent+tmstep*nfc)])[[2]]
            beta1 = R0/matrix(xprior["D",,tt+1],tmstep*(nfc)+tm.ini,num_ens,byrow=TRUE); 
            beta2 = matrix(rowMeans(R0)/mean(xprior["D",,tt+1]),tmstep*(nfc)+1,1)
          }
        }else{
          R0 = paramsEns[,"R0"]
          beta1 = matrix(R0/xprior["D",,tt+1],tmstep*(nfc)+tm.ini,num_ens,byrow=TRUE)
          beta2  = matrix(mean(R0)/mean(xprior["D",,tt+1]),tmstep*(nfc)+1,1)
        }
        Sr_tmp=fn_epi(tcurrent+dt,tcurrent+tmstep*nfc,dt,states=xprior[,,tt+1], N, 
                      beta=beta1,newIpad=newIpad,realdata=T)
        fcast_newI.t=(Sr_tmp$newI[tmstep*(1:nfc)+1,]-Sr_tmp$newI[tmstep*(0:(nfc-1))+1,])%*%diag(xprior["SF",,tt+1],num_ens,num_ens); # 1/31/23 this is the ensemble forecast w/o deflation
        
        # call a fcast.eval function to get the evaluation for the forecast this week
        
        # time stamps: may need to fix these
        Week.fcast.t = rep(tm.ini+tt-1,nfc)
        Week.start.t =seq(tm.ini+tt,tm.ini+tt+nfc-1,by=1)
        
        tmp = fn_eval(fcast_newI.t, Week.fcast.t, Week.start.t)
        tmp2 = fn_eval(t(apply(fcast_newI.t,1,sort)), Week.fcast.t, Week.start.t)
        fcast_stats.t = tmp$fcast_stats
        fcast_stats_ss.t = tmp$fcast_stats_ss
        fcast_Dist.t = tmp$fcastDist
        fcast_Dist.t$Week.fcast = tt
        # save it
        fcast_ens_stats = rbind(fcast_ens_stats, fcast_stats.t)
        fcast_ens_stats_ss = rbind(fcast_ens_stats_ss, fcast_stats_ss.t)
        fcast_ens_Dist = rbind(fcast_ens_Dist, fcast_Dist.t)
        fcast_ens_all[,,tt]=fcast_newI.t
        rm(fcast_stats.t, fcast_stats_ss.t, fcast_Dist.t)
        
        # sorted results
        fcast_stats.t = tmp2$fcast_stats
        fcast_stats_ss.t = tmp2$fcast_stats_ss
        fcast_Dist.t = tmp2$fcastDist
        fcast_Dist.t$Week.fcast = tt
        fcast_ens_stats_sort = rbind(fcast_ens_stats_sort, fcast_stats.t)
        fcast_ens_stats_ss_sort = rbind(fcast_ens_stats_ss_sort, fcast_stats_ss.t)
        fcast_ens_Dist_sort = rbind(fcast_ens_Dist_sort, fcast_Dist.t)
        rm(fcast_stats.t, fcast_stats_ss.t, fcast_Dist.t)
        
        ## change made on Dec 4, 2014
        ## forecasts generated using the mean state variables and parameters
        ## as opposed to using the ensemble and then take the mean
        Sr_tmp=fn_epi(tcurrent+dt,tcurrent+tmstep*nfc,dt,states=t(t(rowMeans(xprior[,,tt+1]))),
                      N,beta=beta2,newIpad=newIpad,realdata=T)
        fc_newI.mn=(Sr_tmp$newI[tmstep*(1:nfc)+1,]-Sr_tmp$newI[tmstep*(0:(nfc-1))+1,])*mean(xprior["SF",,tt+1]);  
        fc_mean=cbind(Week.fcast.t, 
                      Week.start.t,
                      fc_newI.mn);
        fcast=rbind(fcast,fc_mean);
        
        # get evaluation
        tmp = fn_eval(as.matrix(fc_newI.mn), Week.fcast.t, Week.start.t)
        tmp2 = fn_eval(as.matrix(fc_newI.mn), Week.fcast.t, Week.start.t)
        fcast_stats.t = tmp$fcast_stats
        fcast_stats_ss.t = tmp$fcast_stats_ss
        fcast_Dist.t = tmp$fcastDist
        fcast_Dist.t$Week.fcast = tt
        # save it
        fcast_mn_stats = rbind(fcast_mn_stats, fcast_stats.t)
        fcast_mn_stats_ss = rbind(fcast_mn_stats_ss, fcast_stats_ss.t)
        fcast_mn_Dist = rbind(fcast_mn_Dist, fcast_Dist.t)
        rm(fcast_stats.t, fcast_stats_ss.t, fcast_Dist.t)
        
        # sorted results
        fcast_stats.t = tmp2$fcast_stats
        fcast_stats_ss.t = tmp2$fcast_stats_ss
        fcast_Dist.t = tmp2$fcastDist
        fcast_Dist.t$Week.fcast = tt
        fcast_mn_stats_sort = rbind(fcast_mn_stats_sort, fcast_stats.t)
        fcast_mn_stats_ss_sort = rbind(fcast_mn_stats_ss_sort, fcast_stats_ss.t)
        fcast_mn_Dist_sort = rbind(fcast_mn_Dist_sort, fcast_Dist.t)
        rm(fcast_stats.t, fcast_stats_ss.t, fcast_Dist.t)
        
        # 1/31/23
        # generate ensemble forecasts with deflation
        for (i in 1:length(fcast.deflat)){
          v.deflat = c('S','I'); # variables to apply deflation
          state0 = xprior[,,tt+1]
          fcast_newI_deflat.t = NULL
          for(iwk in 1:nfc){ # do it week by week
            if(fcast.deflat[i] !=1){
              xmn=rowMeans(state0); 
              # 10/18/22 - only adjust the state variables
              mat.deflat = matrix(1, nrow = nrow(state0), ncol = ncol(state0))
              rownames(mat.deflat) = rownames(state0)
              mat.deflat[v.deflat,] = fcast.deflat[i] # only adjust the state variables
              
              xnew=mat.deflat*(state0-xmn%*%matrix(1,1,num_ens))+xmn%*%matrix(1,1,num_ens); # same deflation for all
              
              state0 = xnew
            }
            
            beta.t = beta1[(iwk*tmstep-6):(tmstep*iwk+1), ] # crop the beta's for this week
            Sr_tmp=fn_epi(tm_strt=tcurrent+iwk*tmstep-6,tm_end=tcurrent+iwk*tmstep,tm_step=dt,
                          states=state0, N, beta=beta.t,newIpad=newIpad,realdata=T)
            
            fcast_newI_deflat.t = rbind(fcast_newI_deflat.t, 
                                        tail(Sr_tmp$newI, 1) * state0['SF',]) # total number of infections?
            
            # update state0
            state0['S',] = tail(Sr_tmp$S, 1)
            state0['I',] = tail(Sr_tmp$I, 1)
          }
          
          # get evaluation
          tmp = fn_eval(fcast_newI_deflat.t, Week.fcast.t, Week.start.t)
          tmp2 = fn_eval(t(apply(fcast_newI_deflat.t,1,sort)), Week.fcast.t, Week.start.t)
          fcast_stats.t = tmp$fcast_stats
          fcast_stats_ss.t = tmp$fcast_stats_ss
          fcast_Dist.t = tmp$fcastDist
          fcast_Dist.t$Week.fcast = tt
          # save it
          assign(paste0("fcast_ens.df",i,"_stats"), rbind(get(paste0("fcast_ens.df",i,"_stats")), fcast_stats.t))
          assign(paste0("fcast_ens.df",i,"_stats_ss"), rbind(get(paste0("fcast_ens.df",i,"_stats_ss")), fcast_stats_ss.t))
          assign(paste0("fcast_ens.df",i,"_Dist"), rbind(get(paste0("fcast_ens.df",i,"_Dist")), fcast_Dist.t))
          if (i == 1){
            fcast_ens.df1_all[,,tt]=fcast_newI_deflat.t
          }else{
            fcast_ens.df2_all[,,tt]=fcast_newI_deflat.t
          }
          
          rm(fcast_stats.t, fcast_stats_ss.t, fcast_Dist.t)
          
          fcast_stats.t = tmp2$fcast_stats
          fcast_stats_ss.t = tmp2$fcast_stats_ss
          fcast_Dist.t = tmp2$fcastDist
          fcast_Dist.t$Week.fcast = tt
          # save it
          assign(paste0("fcast_ens.df",i,"_stats_sort"), rbind(get(paste0("fcast_ens.df",i,"_stats_sort")), fcast_stats.t))
          assign(paste0("fcast_ens.df",i,"_stats_ss_sort"), rbind(get(paste0("fcast_ens.df",i,"_stats_ss_sort")), fcast_stats_ss.t))
          assign(paste0("fcast_ens.df",i,"_Dist_sort"), rbind(get(paste0("fcast_ens.df",i,"_Dist_sort")), fcast_Dist.t))
          rm(fcast_stats.t, fcast_stats_ss.t, fcast_Dist.t)
          
          
        }
        
        
        
      } # end forecast and evaluation
    }
    #  Integrate forward one time step 
    tcurrent = tm.ini+tmstep*tt; 
    update_params=intersect(colnames(paramsEns),rownames(xprior))
    paramsEns[,update_params]=xprior[update_params,,tt+1] %>% t
    if (grepl("climate",epi.model)){
      if (grepl("normalized",epi.model)){
        R0 = paramsEns[,"R0"]
        beta=R0*matrix(sn_trend[tcurrent:(tcurrent+tmstep-1)],tmstep,num_ens,byrow = F)/matrix(xprior["D",,tt+1],tmstep,num_ens,byrow=TRUE);
        
      }else{
        R0 = calc_R0_M1_new(paramsEns,num_ens,sh1.daily[tcurrent:(tcurrent+tmstep-1)], 
                            mean.temp.daily[tcurrent:(tcurrent+tmstep-1)])[[2]]
        beta=R0/matrix(xprior["D",,tt+1],tmstep,num_ens,byrow=T); 
      }
    }else{
      R0 = paramsEns[,"R0"]
      beta = matrix(R0/xprior["D",,tt+1],tmstep,num_ens,byrow=TRUE)
    }
    #print(newIpad[1:5,361:365])
    Sr_tmp=fn_epi(tm_strt=tcurrent+dt,tm_end=tcurrent+tmstep,tm_step=dt,
                  states=xprior[,,tt+1], 
                  N, beta=beta,newIpad=newIpad,realdata=T)
    xprior["S",,tt+1]=tail(Sr_tmp$S,1);
    xprior["I",,tt+1]=tail(Sr_tmp$I,1);
    xprior["newI",,tt+1]=tail(Sr_tmp$newI,1)*xprior["SF",,tt+1]; 
    
  } # end for-loop
  
  # get the ensemble mean
  xpost_mean=xpost_sd=xprior_mean=xprior_sd=xpost_95CI_upr=xpost_95CI_lwr=xprior_95CI_upr=xprior_95CI_lwr=matrix(0,num_times,num_var); 
  colnames(xpost_mean) = param.names
  xpost_R0 = xpost_Re = numeric(num_times)
  for(ti in 1:num_times){
    xpost_mean[ti,]=rowMeans(xpost[,,ti]);
    xprior_mean[ti,]=rowMeans(xprior[,,ti]);
    xpost_sd[ti,]=apply(xpost[,,ti],1,sd);
    xprior_sd[ti,]=apply(xprior[,,ti],1,sd);
    tcurrent = tm.ini+tmstep*(ti-1); 
    if (grepl("climate",epi.model)){
      common_params = intersect(colnames(xpost_mean),colnames(paramsEns))
      paramsEns.t = paramsEns[1,]
      paramsEns.t[common_params]=xpost_mean[1,common_params] 
      paramsEns.t=as.matrix(paramsEns.t) %>% t
      R0=calc_R0_M1_new(paramsEns.t, 1, sh1.daily[tcurrent:(tcurrent+tmstep-1)], 
                        mean.temp.daily[tcurrent:(tcurrent+tmstep-1)])[[2]] %>% mean
      
    }else{
      R0=xpost_mean[ti,"R0"]
    }
    
    Re=xpost_mean[ti,"S"]*R0/N;
    xpost_R0[ti] = R0
    xpost_Re[ti] = Re
  }
  colnames(xpost_mean)=c(param.names)
  
  ## output data
  tstep=seq(tm.ini,num_times+tm.ini-1,by=1); # time indices
  # metrics for comparison
  Y=xpost_mean[,"newI"];  # newI for each strain
  sf=xpost_mean[,"SF"];
  Y.tot=Y/sf;
  
  
  xpost_mean=cbind(xpost_mean,xpost_Re)
  colnames(xpost_mean)[ncol(xpost_mean)]="Re"
  out1=cbind(tstep,xpost_mean,Y.tot);
  out2=cbind(tstep,xpost_sd);
  out3=cbind(tstep,xprior_mean,Y.tot);
  out4=cbind(tstep,xprior_sd);
  out5=cbind(tstep,lambdas);
  colnames(out1)=c("time",param.names,'Re','newI.tot');
  colnames(out3)=c("time",param.names,'newI.tot');
  
  colnames(out2)=colnames(out4)=colnames(out5)=c("time",param.names);
  if (forecast){
    # save the 6 outputs into two
    fcast_mn_Dist$fcast.type = "mean";fcast_ens_Dist$fcast.type="ensembles"
    fcast_mn_stats$fcast.type = "mean";fcast_ens_stats$fcast.type="ensembles"
    fcast_mn_stats_ss$fcast.type = "mean";fcast_ens_stats_ss$fcast.type="ensembles"
    fcast_mn_Dist$fcast.deflat = 1;fcast_ens_Dist$fcast.deflat=1
    fcast_mn_stats$fcast.deflat = 1;fcast_ens_stats$fcast.deflat=1
    fcast_mn_stats_ss$fcast.deflat = 1;fcast_ens_stats_ss$fcast.deflat=1
    
    fcast_mn_Dist_sort$fcast.type = "mean";fcast_ens_Dist_sort$fcast.type="ensembles"
    fcast_mn_stats_sort$fcast.type = "mean";fcast_ens_stats_sort$fcast.type="ensembles"
    fcast_mn_stats_ss_sort$fcast.type = "mean";fcast_ens_stats_ss_sort$fcast.type="ensembles"
    fcast_mn_Dist_sort$fcast.deflat = 1;fcast_ens_Dist_sort$fcast.deflat=1
    fcast_mn_stats_sort$fcast.deflat = 1;fcast_ens_stats_sort$fcast.deflat=1
    fcast_mn_stats_ss_sort$fcast.deflat = 1;fcast_ens_stats_ss_sort$fcast.deflat=1
    for (i in 1:length(fcast.deflat)){
      tmp_Dist = get(paste0("fcast_ens.df",i,"_Dist"))
      tmp_Dist$fcast.type = "ensembles_deflated"
      tmp_Dist$fcast.deflat = fcast.deflat[i]
      assign(paste0("fcast_ens.df",i,"_Dist"),tmp_Dist)
      
      tmp_stats = get(paste0("fcast_ens.df",i,"_stats"))
      tmp_stats$fcast.type = "ensembles_deflated"
      tmp_stats$fcast.deflat = fcast.deflat[i]
      assign(paste0("fcast_ens.df",i,"_stats"),tmp_stats)
      
      tmp_stats_ss = get(paste0("fcast_ens.df",i,"_stats_ss"))
      tmp_stats_ss$fcast.type = "ensembles_deflated"
      tmp_stats_ss$fcast.deflat = fcast.deflat[i]
      assign(paste0("fcast_ens.df",i,"_stats_ss"),tmp_stats_ss)
    }
    rm(tmp_Dist); rm(tmp_stats); rm(tmp_stats_ss)
    fcastDist = do.call("rbind",mget(ls(pattern = "^fcast_.*Dist$")))
    fcast_stats = do.call("rbind",mget(ls(pattern = "^fcast_.*stats$")))
    fcast_stats_ss = do.call("rbind",mget(ls(pattern = "^fcast_.*stats_ss$")))
    
    for (i in 1:length(fcast.deflat)){
      tmp_Dist = get(paste0("fcast_ens.df",i,"_Dist_sort"))
      tmp_Dist$fcast.type = "ensembles_deflated"
      tmp_Dist$fcast.deflat = fcast.deflat[i]
      assign(paste0("fcast_ens.df",i,"_Dist_sort"),tmp_Dist)
      
      tmp_stats = get(paste0("fcast_ens.df",i,"_stats_sort"))
      tmp_stats$fcast.type = "ensembles_deflated"
      tmp_stats$fcast.deflat = fcast.deflat[i]
      assign(paste0("fcast_ens.df",i,"_stats_sort"),tmp_stats)
      
      tmp_stats_ss = get(paste0("fcast_ens.df",i,"_stats_ss_sort"))
      tmp_stats_ss$fcast.type = "ensembles_deflated"
      tmp_stats_ss$fcast.deflat = fcast.deflat[i]
      assign(paste0("fcast_ens.df",i,"_stats_ss_sort"),tmp_stats_ss)
    }
    rm(tmp_Dist); rm(tmp_stats); rm(tmp_stats_ss)
    fcastDist_sort = do.call("rbind",mget(ls(pattern = "^fcast_.*Dist_sort$")))
    fcast_stats_sort = do.call("rbind",mget(ls(pattern = "^fcast_.*stats_sort$")))
    fcast_stats_ss_sort = do.call("rbind",mget(ls(pattern = "^fcast_.*stats_ss_sort$")))
    
    
    colnames(fcast)=c('fc_start','week',"iliiso.tot")
    
  }
  print("done")
  # get the distribution of states at each timestep
  states_stats=NULL
  for (i in 1:dim(xprior)[3]){
    tmp = xprior[,,i] %>% apply(1,function(x){
      mean = mean(x)
      median = median(x)
      sd = sd(x)
      iqr.lwr = quantile(x,0.25)[[1]]
      iqr.upr = quantile(x,0.75)[[1]]
      for (confidence_level in c(95)){
        assign(paste0("ci",confidence_level,".lwr"),quantile(x,(1-confidence_level/100)/2)[[1]])
        assign(paste0("ci",confidence_level,".upr"),quantile(x,1-(1-confidence_level/100)/2)[[1]])
      }
      df = mget(c("mean", "median", "iqr.lwr", "iqr.upr", ls(pattern = "^ci"))) %>% unlist
      df %>% t() %>% as.data.frame(row.names = T) 
    }) %>% rbindlist()
    tmp$state=param.names
    tmp$Week.start=weeks[i-1]
    states_stats=rbind(states_stats,tmp)
  }
  xprior_stats=states_stats
  states_stats=NULL
  for (i in 1:dim(xpost)[3]){
    tmp = xpost[,,i] %>% apply(1,function(x){
      mean = mean(x)
      median = median(x)
      sd = sd(x)
      iqr.lwr = quantile(x,0.25)[[1]]
      iqr.upr = quantile(x,0.75)[[1]]
      for (confidence_level in c(95)){
        assign(paste0("ci",confidence_level,".lwr"),quantile(x,(1-confidence_level/100)/2)[[1]])
        assign(paste0("ci",confidence_level,".upr"),quantile(x,1-(1-confidence_level/100)/2)[[1]])
      }
      df = mget(c("mean", "median", "iqr.lwr", "iqr.upr", ls(pattern = "^ci"))) %>% unlist
      df %>% t() %>% as.data.frame(row.names = T) 
    }) %>% rbindlist()
    tmp$state=param.names
    tmp$Week.start=weeks[i]
    states_stats=rbind(states_stats,tmp)
  }
  xpost_stats=states_stats
  
  
  # save the last ens
  xlast=t(xprior[,,tt+1]); 
  colnames(xlast)=c(param.names);
  
  out=list(xprior_stats=xprior_stats,xpost_stats=xpost_stats,xpost_mean=out1,xpost_sd=out2,xprior_mean=out3,xprior_sd=out4,fcast_stats=fcast_stats,
           fcast_stats_ss=fcast_stats_ss,fcastDist=fcastDist, xpost=xpost,
           fcast_ens.df1_all=fcast_ens.df1_all,fcast_ens.df2_all=fcast_ens.df2_all,fcast_ens_all=fcast_ens_all,
           loc="HK",variant="iliiso.tot",seasonality=epi.model); 
  return(out)
}


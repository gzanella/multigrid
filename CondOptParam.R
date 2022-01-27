###  R CODE  #######

### SIMULATION SETTING DESCRIBED IN SECTION 3.2 AND FIGURE 5 OF THE PAPER
### "Multilevel Linear Models, Gibbs Samplers and Multigrid Decompositions"
### GIBBS SAMPLER FOR SYMMETRIC THREE LEVEL HIERARCHICAL MODELS (SAMPLING BOTH MEANS AND VARIANCES)
### CODE TO IMPLEMENT THE CENTERED ADN THE CONDITIONALLY-OPTIMAL PARAMETRIZATIONS
###

## SET HYPERPARAMETERS OF PRIORS FOR VARIANCES
shape_a<-shape_b<-shape_e<-0.01;rate_a<-rate_b<-rate_e<-0.01

#### SIMULATE DATA #####
I<-100;J<-100;K<-5;real_tau_a<-0.01;real_tau_b<-10;real_tau_e<-0.01
real_mu<-0
real_a<-rnorm(n=I,mean = 0, sd = (real_tau_a)^(-0.5))
real_b<-matrix(rnorm(n=I*J,mean = 0,sd = (real_tau_b)^(-0.5)),nrow=I,ncol=J)
Y<-matrix(rnorm(n=I*J,mean = real_mu+rep(real_a,J)+real_b,sd = (K*real_tau_e)^(-0.5)),nrow=I,ncol=J)
Ymean<-mean(Y)
Yi<-rowMeans(Y)
Y_var<-matrix(rchisq(n = I*J,df = K-1)/real_tau_e,nrow=I,ncol=J)
sum_Y_var<-sum(Y_var)

#### SET MCMC parameters #####
T <- 10^4;T_plot<-10^4;T_burn<-10^3
means_samp <- matrix(0, nrow=T, ncol=3)
means_samp_list<-list();means_samp_list[[2]]<-NA
prec_samp <- matrix(0, nrow=T, ncol=3)
prec_samp_list<-list();prec_samp_list[[2]]<-NA

### STARTING VALUES FOR THE CHAINS
start_mu<-rnorm(1,0,sqrt(10))
start_a<-rnorm(I,rnorm(1,sqrt(5)),sqrt(5))
start_b<-matrix(rnorm(I*J,rnorm(1,sqrt(5)),sqrt(5)),nrow=I,ncol=J)
start_tau_a<-rgamma(n = 1,shape = 0.1,rate = 0.1)
start_tau_b<-rgamma(n = 1,shape = 0.1,rate = 0.1)
start_tau_e<-rgamma(n = 1,shape = 0.1,rate = 0.1)

######## SAMPLER 1 - CENTERED PARAMETRIZATION ########
mu<-start_mu;gamma<-start_mu+start_a;eta<-start_mu+start_a+start_b
tau_a<-start_tau_a;tau_b<-start_tau_b;tau_e<-start_tau_e
print(system.time(
  # main sampling loop
  for (t in 1:T) {
    mu<-rnorm(1, mean = mean(eta), 
              sd = (I*tau_a)^(-0.5) )
    gamma<-rnorm(I, mean = (tau_a*mu+J*tau_b*rowMeans(eta))/(tau_a+J*tau_b),
                 sd = (tau_a+J*tau_b)^(-0.5) )
    eta<-matrix(
      rnorm(I*J,mean = (tau_b*rep(gamma,J)+K*tau_e*Y)/(tau_b+K*tau_e),
            sd = (tau_b+K*tau_e)^(-0.5)),
      nrow=I,ncol=J)
    tau_a<-rgamma(n=1,shape = shape_a+I/2,
                  rate = rate_a+sum((gamma-mu)^2)/2)
    tau_b<-rgamma(n=1,shape = shape_b+I*J/2, 
                  rate = rate_b+sum((eta-rep(gamma,J))^2)/2)
    tau_e<-rgamma(n=1,shape = shape_e+I*J*K/2, 
                  rate = rate_e+(sum_Y_var+K*sum((Y-eta)^2))/2)
    means_samp[t,1] <- mu
    means_samp[t,2] <- mean(gamma)
    means_samp[t,3] <- mean(eta)
    prec_samp[t,]<-c(tau_a,tau_b,tau_e)
  }
))
means_samp_list[[1]]<-means_samp
prec_samp_list[[1]]<-prec_samp

######## SAMPLER 2 - CONDITIONALLY-OPTIMAL PARAMETRIZATION ########
mu<-start_mu;a<-start_a;gamma<-a+mu;b<-start_b;eta<-mu+rep(a,J)+b
tau_a<-start_tau_a;tau_b<-start_tau_b;tau_e<-start_tau_e
print(system.time(
  # main sampling loop
  for (t in 1:T) {
    if(1/tau_a>1/(J*tau_b)+1/(J*K*tau_e)){#USE GAMMA
      if(1/tau_b>+1/(K*tau_e)){#USE ETA
        mu<-rnorm(1, mean = mean(eta), 
                  sd = (I*tau_a)^(-0.5) )
        gamma<-rnorm(I, mean = (tau_a*mu+J*tau_b*rowMeans(eta))/(tau_a+J*tau_b),
                     sd = (tau_a+J*tau_b)^(-0.5) )
        eta<-matrix(
          rnorm(I*J,mean = (tau_b*rep(gamma,J)+K*tau_e*Y)/(tau_b+K*tau_e),
                sd = (tau_b+K*tau_e)^(-0.5)),
          nrow=I,ncol=J)
        a<-gamma-mu
        b<-eta-rep(gamma,J)
      }else{#USE b
        mu<-rnorm(1, mean = mean(eta), 
                  sd = (I*tau_a)^(-0.5) )
        gamma<-rnorm(I, mean = (tau_a*mu+J*K*tau_e*(Yi-rowMeans(b)))/(tau_a+J*K*tau_e),
                     sd = (tau_a+J*K*tau_e)^(-0.5) )
        b<-matrix(
          rnorm(I*J,mean = (K*tau_e/(tau_b+K*tau_e))*(Y-rep(gamma,J)),
                sd = (tau_b+K*tau_e)^(-0.5)),
          nrow=I,ncol=J)
        a<-gamma-mu
        eta<-rep(gamma,J)+b
      }
    }else{#USE a
      if(1/tau_b>+1/(K*tau_e)){#USE ETA
        mu<-rnorm(1, mean = mean(eta)-mean(a), 
                  sd = (I*J*tau_b)^(-0.5) )
        a<-rnorm(I, mean = (J*tau_b*(rowMeans(eta)-mu))/(tau_a+J*tau_b),
                 sd = (tau_a+J*tau_b)^(-0.5) )
        eta<-matrix(
          rnorm(I*J,mean = (tau_b*(mu+rep(a,J))+K*tau_e*Y)/(tau_b+K*tau_e),
                sd = (tau_b+K*tau_e)^(-0.5)),
          nrow=I,ncol=J)
        gamma<-a+mu
        b<-eta-rep(gamma,J)
      }else{#USE b
        mu<-rnorm(1, mean = Ymean-mean(a)-mean(b) , 
                  sd = (I*J*K*tau_e)^(-0.5) )
        a<-rnorm(I, mean = (J*K*tau_e/(tau_a+J*K*tau_e))*(Yi-mu-rowMeans(b)),
                 sd = (tau_a+J*K*tau_e)^(-0.5) )
        b<-matrix(
          rnorm(I*J,mean = (K*tau_e/(tau_b+K*tau_e))*(Y-mu-rep(a,J)),
                sd = (tau_b+K*tau_e)^(-0.5)),
          nrow=I,ncol=J)
        gamma<-a+mu
        eta<-mu+rep(a,J)+b
      }
    }
    tau_a<-rgamma(n=1,shape = shape_a+I/2,rate = rate_a+sum(a^2)/2)
    tau_b<-rgamma(n=1,shape = shape_b+I*J/2,rate = rate_b+sum(b^2)/2)
    tau_e<-rgamma(n=1,shape = shape_e+I*J*K/2,rate = rate_e+(sum_Y_var+K*sum((Y-mu-rep(a,J)-b)^2))/2)
    means_samp[t,1] <- mu
    means_samp[t,2] <- mean(gamma)
    means_samp[t,3] <- mean(eta)
    prec_samp[t,]<-c(tau_a,tau_b,tau_e)
  }
))
means_samp_list[[2]]<-means_samp
prec_samp_list[[2]]<-prec_samp

##### PLOTTING ACFs#############
###########################
names<-c("Centred","Optimal")
par(mfrow=c(2,6))
for (i in 1:2){
  par(mar=c(0,3,0,0))
  acf(means_samp_list[[i]][T_burn:T,1],lag.max = 50,xaxt='n',ylim=c(-0.05,1.05),ylab=names[i])  
  par(mar=c(0,0,0,0))
  acf(means_samp_list[[i]][T_burn:T,2],lag.max = 50,xaxt='n',yaxt='n',ylim=c(-0.05,1.05))
  acf(means_samp_list[[i]][T_burn:T,3],lag.max = 50,xaxt='n',yaxt='n',ylim=c(-0.05,1.05))
  acf(1/sqrt(prec_samp_list[[i]][T_burn:T,1]),lag.max = 50,xaxt='n',yaxt='n',ylim=c(-0.05,1.05))
  acf(1/sqrt(prec_samp_list[[i]][T_burn:T,2]),lag.max = 50,xaxt='n',yaxt='n',ylim=c(-0.05,1.05))
  acf(1/sqrt(prec_samp_list[[i]][T_burn:T,3]),lag.max = 50,xaxt='n',yaxt='n',ylim=c(-0.05,1.05))
}

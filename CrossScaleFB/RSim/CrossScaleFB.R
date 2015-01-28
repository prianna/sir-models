################################
################################
## the cross scale simulator
## with fixed bottle neck
## created 12/2013 by Sebastian Schreiber
################################
################################

# basic idea of the model is that are k viral strains that compete within any infected host. The within-host dynamics of these viral strains is specified by the within.host function (e.g. a linear model with ceiling below for illustrative purposes). An infected host begins with an initial viral load N. For a host with a given viral load that has made a contact with an infected individual, the function transmission uses multinomial sampling to determines the initial viral load in the contacted individual. Finally, the function between.host runs the full dynamics by tracking each generation of infected hosts. These infected hosts can infect others for T units. During these T units of time, they may contacts at a rate beta and transmit succesfully acording to the prob.transmission function.

# need the following library to compute the matrix exponential 

library(Matrix)

# global parameters

r=c(0.1,0.11) # within host per-capita growth rates
b=c(1,2) # weights for determining transmission
K=10000 # within-host carrying capacity
k=2 # number of viral types
mu=0.001 #mutation rate
N=100 # initial viral load size
T=15 # length of infection
beta=5 # contact rate

##############################
# the within host function
##############################
# input: initial condition v, vector of times 
# output: list with matrix with viral numbers at the desired times 
# this within host model gives exponential growth until the population hits the maximal viral load N at which time the dynamics are normalized.

within.host=function(v,times){  
  M=matrix(0,k,k) # mutation matrix (one steps to the right)
  M[-1,-k]=diag(k-1)*mu
  A=M+diag(r) # the viral dynamic matrix
  
  
  L=length(times) # number of evaluations
  k=length(v) # number of viral types
  V=matrix(0,k,L) # matrix for output; columns are the viral loads at the specified times. 
  
  for(i in 1:L){
    temp=as.matrix(expm(A*times[i])%*%v);
    if(sum(temp)>K)temp=K*temp/sum(temp) # normalize if viral load is greater than N
    V[,i]=temp
  }
  return(list(V=V))
}




############################
# single strain R0 function
############################
# input: nada
# output: R0 (vector)
# computes the R0s for the single strains 
# NOTE: If the function transmission is changed below, it also needs to be changed here. 
R.naught=function(){
  R=numeric(k)
  
  for(i in 1:k){
    S=function(t)1-exp(-b[i]*pmin(N*exp(r[i]*t),K)/N) 
    R[i]=beta*(integrate(S,0,T)$value)
  }
  return(R)	
}





###################################
# the transmission functions
###################################
# input: matrix of viral column vectors
# output: viral type of newly infected

transmission=function(V,n){
  V.weighted=diag(b)%*%V
  if(n>1)V.weighted=V.weighted%*%diag(1/colSums(V.weighted))
  if(n==1)V.weighted=V.weighted/sum(V.weighted)
  W=c()
  for(i in 1:n)W=cbind(W,rmultinom(n=1,size=N,prob=V.weighted[,i]))
  return(W)
}


# input: matrix of viral column distribution vectors
# output: viral type of newly infected

prob.transmission=function(V)1-exp(-colSums(diag(b)%*%V)/N)

#######################################
# between host dynamics (the cross scale model)
########################################
# input: initial infected type v, maximum cases
# output: list with who (who infected this individual, 0 for case 1), when (when was this individual infected), type (matrix of types)

between.host=function(v,cases){
  type=matrix(v,k,1) # holds the intial viral types for all individuals; dynamically increased by adding columns. 
  who=c(0) # vector of who was the parent of the infected individual
  when=c(0) # when the individual was initially infected
  gen=c(0) # generation when individual was infected
  total=1 # total cases
  current=1 # current cases
  gen.count=1
  
  while(total<(cases+1)&&current>0){
    # next line determines the number of contacts each individual in the current generation has with other individuals 
    kids=rpois(n=current,lambda=beta*T)
    # the next loop determines at what times the contact events occur, the transmission of viral particles, and removes any instances where transmitted load is zero
    for(i in 1:current){ # loop for all individuals in current generation
      if(kids[i]>0){ # if current individual has a positive number of contacts do the following 
        v.start=type[,total-current+i] # get the initial viral load of the current individual 
        t=runif(kids[i])*T # find the times at which the contact events occur
        out=within.host(v.start,t) # find the viral load at the contact events 
        v.out=out$V
        success=rbinom(kids[i],1,prob=prob.transmission(v.out)) # which transmissions where succesful
        w=transmission(v.out,kids[i])*success # find the transmitted viral load at these contact events
        viral.counts=colSums(w) # find the total viral load due to the transmission events
        positive=which(viral.counts>0) # find which transmission events included a positive number of viral particles
        if(length(positive)>0){ # if there were some real transmission events do the following
          who=c(who,rep(total-current+i,length(positive))) # make the current individual the parent of all the newly infected individuals 
          mom.when=when[total-current+i] # find the time when mom was infected. 
          when=c(when,t[positive]+mom.when) # add the times at which the newly infected individuals got infected
          
          type=cbind(type,w[,positive]) # add the types for the newly infected individuals 
          kids[i]=length(positive) # reset kids[i] to the number of newly infected individuals	
        }
        if(length(positive)==0)kids[i]=0 # reset kids[i] when there where no infections
        
      }
    }
    gen=c(gen,rep(gen.count,sum(kids))) # add the generation times for all the newly infected individuals for the current generation 
    current=sum(kids) # update the count for the current number of infected individuals 
    total=total+current # update the count for the cummulative number of infected individuals
    gen.count=gen.count+1 # increate the generation count
    
  }
  
  return(list(when=when,who=who,type=type,gen=gen))
  
  
}
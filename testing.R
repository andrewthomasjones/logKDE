library(Conake)
library(VGAM)

pos_claw<-function(x,q){
    y=0.5*dnorm(x, mean = 0, sd=q)
    subzero=0.5*pnorm(0, mean = 0, sd=q)
    for(i in 1:7){
      y=y+(1/14)*dnorm(x, mean=((i/2)+0), sd=0.1)
      subzero=subzero+(1/14)*pnorm(0, mean=((i/2)+0), sd=0.1)
    }

    y[x<=0]=0
    y<-y/(1-subzero)
    return(y)
}

sin_mad_mix<-function(x,q){
  y=0.4*dsinmad(x,2.8,0.193,1.7)+0.6*dsinmad(x,5.8,5.93,q)
  return(y)
}

sin_mad<-function(x,q){
  y=dsinmad(x,2.8,0.193,q)
  return(y)
}

l_normal<-function(x,q){
  y=dlnorm(x, meanlog = 0, sdlog = q)
  return(y)
}

abs_normal<-function(x,q){
  y=2*(dnorm(x, mean =0, sd=q))
  y[x<=0]<-0
  return(y)
}

function_list<-list(pos_claw=pos_claw, sin_mad_mix=sin_mad_mix, sin_mad=sin_mad, l_normal=l_normal, abs_normal=abs_normal)

MISE<-function(dens, x, q, flist, k){
  n=length(dens)
  true_dens<-flist[[k]](x,q)
  (1/n)*(trapz(x,(dens-true_dens)^2))
}

ranorm<-function(n,q){
  abs(rnorm(n, sd=q))
}
rlnorm2<-function(n,q){
  rlnorm(n, meanlog = 0, sdlog = q)
}
rsin_mad<-function(n,q){
  y=rsinmad(n,2.8,0.193,q)
  return(y)
}

rsin_mad_mix<-function(n,q){
  y1<-rsinmad(n,2.8,0.193,1.7)
  y2<-rsinmad(n,5.8,5.93,q)
  y<-y1
  ind<-runif(n)>0.6
  y[ind]<-y2[ind]
  return(y)
}


rpos_claw<-function(n,q){
  y<-array(0,c(n,8))
  y[,1]<-abs(rnorm(n, mean = 0, sd=q))

   for(i in 1:7){
    y[,i+1] = (((rnorm(n, mean=((i/2)+0), sd=0.2))))
   }

  for(i in 1:7){
    for(j in which(y[,i+1]<0)){
      tmp<-1
      while(tmp<0){

        tmp<-(rnorm(1, mean=((i/2)+0), sd=0.2))

      }
      y[j,i+1]<-tmp
    }
  }

  chooser<-rmultinom(n, size=1, prob=c(0.5, rep(1/14, 7)))

  y2<-array(0,n)

  for(i in 1:n){
    j=which.max(chooser[,i])
    y2[i]<-y[i,j]
  }

  return(y2)

}

V<-rgamma(100,1.5,2.6)
h<-cvbw(V, bw = NULL, ker="RIG")
x<-seq.int(from, to, length.out = n.user)
est<-dke(V,"RIG",h$hcv, x=x)
y<-est$f_n

V<-rgamma(100,1.5,2.6)
h<-cvbw(V, bw = NULL, ker="GA")
x<-seq.int(from, to, length.out = n.user)
est<-dke(V,"GA",h$hcv, x=x)
y<-est$f_n




### Sample
XX <- rlnorm(1000)
## Grid size
GRID <- 10
# Maximum value for h
HMAX <- 5*logdensity(XX)$bw
# Min value for h
HMIN <- 0.1*HMAX
### Number of bins
NB <- 100
### Bandwidth
HH <- seq(HMIN,HMAX,length.out = GRID)
## Storage for CV
CVSTORE <- c()

for (hh in 1:GRID) {
  LD <- logdensity(XX,bw=HH[hh],n=NB)
  ### Compute CV for given bandwidth and NB
  CV <- 0
  for (xx in 1:length(XX)) {
    CV <- CV + sum((LD$y - logdensity(XX[-xx],bw=HH[hh],from = min(LD$x),to = max(LD$x),n=NB)$y)^2)/NB
  }
  CVSTORE[hh] <- CV/length(XX)
}
HH[which.min(CVSTORE)]




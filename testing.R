library(Conake)
library(VGAM)
library(pracma)


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
  y=0.4*dsinmad(x,0.193,2.8,1.7)+0.6*dsinmad(x,0.593,5.8,q)
  return(y)
}

sin_mad<-function(x,q){
  y=dsinmad(x,0.193,2.8,q)
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

ranorm<-function(n,q){
  abs(rnorm(n, sd=q))
}
rlnorm2<-function(n,q){
  rlnorm(n, meanlog = 0, sdlog = q)
}
rsin_mad<-function(n,q){
  y=rsinmad(n,0.193,2.8,q)
  return(y)
}

rsin_mad_mix<-function(n,q){
  y1<-rsinmad(n,0.193,2.8,1.7)
  y2<-rsinmad(n,0.593,5.8,q)
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



gamma_test<-function(y, bw_method, from, to, se_len, kern=NULL){
  if(bw_method=="CV") {h<-cvbw(y, bw = NULL, ker="GA")$hcv}
  if(bw_method=="silverman"){h<-bw.nrd0(y)}
  if(bw_method=="lsilverman"){h<-bw.logG(y)}

  x<-seq.int(from, to, length.out = length(y))
  est<-dke(vec_data=as.vector(y),"GA",bw=h, x=x)
  return_list<-list(y=est$f_n,x=x, bw=h)

  return(return_list)
}

rig_test<-function(y, bw_method, from, to, se_len, kern=NULL){
  if(bw_method=="CV") {h<-cvbw(y, bw = NULL, ker="RIG")$hcv}
  if(bw_method=="silverman"){h<-bw.nrd0(y)}
  if(bw_method=="lsilverman"){h<-bw.logG(y)}

  x<-seq.int(from, to, length.out = length(y))
  est<-dke(y,"RIG",h, x=x)

  return_list<-list(y=est$f_n, x=x, bw=h)

  return(return_list)
}

norm_test<-function(y, bw_method, from, to, se_len, kern="gaussian"){
  if(bw_method=="CV") {h<-bw.logCV(y)}
  if(bw_method=="silverman"){h<-bw.nrd0((y))}
  if(bw_method=="lsilverman"){h<-bw.logG(y)}
  est<-density(y, bw=h, from=from, to=to, n=se_len, kernel = kern)

  return_list<-list(y=est$y, x=est$x, bw=h)
  return(return_list)
}

log_test<-function(y, bw_method, from, to, se_len, kern="gaussian"){
  if(bw_method=="CV") {h<-bw.logCV(y)}
  if(bw_method=="silverman"){h<-bw.nrd0(log(y))}
  if(bw_method=="lsilverman"){h<-bw.logG(y)}

  est<-logdensity(y,bw=h, from=from, to=to, n=se_len, kernel = kern)

  return_list<-list(y=est$y, x=est$x,  bw=h)
  return(return_list)
}




function_list<-list(sin_mad_mix=sin_mad_mix, sin_mad=sin_mad, l_normal=l_normal)#, abs_normal=abs_normal,pos_claw=pos_claw)
rfunction_list<-list(rsin_mad_mix=rsin_mad_mix, rsin_mad=rsin_mad, rlnorm2=rlnorm2)#, ranorm=ranorm, rpos_claw=rpos_claw)
q_list<-list(rsin_mad_mix=c(0.7,0.5,0.3), rsin_mad=c(1.45,1.07,0.75), rlnorm2=c(0.5,1,2))# ranorm=c(0.5,1,2), rpos_claw=c(0.5,1,2))

to_list<-list(rsin_mad_mix=c(50,700,12000), rsin_mad=c(30,80,1050), rlnorm2=c(25,250,70000))

method_list<-list(normalkde = norm_test, gamma = gamma_test, rig=rig_test, logKDE= log_test)
kernel_list<-list("gaussian", "epanechnikov", "triangular", "uniform", "laplace", "logistic")
kernel_list2<-list("gaussian", "epanechnikov", "triangular", "rectangular")
kern_list<-list(normalkde = kernel_list2, gamma = "gamma", rig="RIG", logKDE= kernel_list)

# #for making a guess at what a good max is
# for(i in 1:fcount){
#   for(k in 1:length(q_list[[i]])){
#     q<-q_list[[i]][k]
#     y<-rfunction_list[[i]](10^8,q)
#     to_list[[i]][k]<-ceiling(max(y))
#     print(to_list[[i]][k])
#
#   }
# }



bw_list<-list("silverman", "lsilverman",  "CV")

MISE<-function(dens, x, q, flist){
  n=length(dens)
  true_dens<-flist(x,q)
  (trapz(x,(dens-true_dens)^2))
}


M=1000
n=512

fcount<-length(function_list)
count<-0
big_list<-list()
for(i in 1:fcount){
  for(k in 1:length(q_list[[i]])){

    q<-q_list[[i]][k]
    from<-0.01
    to<-to_list[[i]][k]

    se_len<-512

    for(j in 1:M){
      y<-rfunction_list[[i]](n,q)
      fileConn<-file(paste0("/home/andrew/Dropbox/logKDEresults_", Sys.Date(), ".csv"),"a")
      for(meth in 1:length(method_list)){

        for(bw_method in bw_list){

          for(kern in 1:length(kern_list[[meth]])){
            count<-count+1
            kern_type=kern_list[[meth]][[kern]]
            test_result<-method_list[[meth]](y, bw_method, from, to, se_len, kern_type)
            test_result$MISE <- MISE(test_result$y, test_result$x, q, function_list[[i]])
            results<-paste(names(function_list)[i],names(rfunction_list)[i], j, q,  kern_type, bw_method, names(method_list)[meth], test_result$MISE, test_result$bw, sep=",")
            print(results)
            big_list[[count]]<-results
            writeLines(results, fileConn)
          }
        }
      }
      close(fileConn)

    }
  }
}

save(big_list,"/home/andrew/Dropbox/logKDEresults.RData" )


















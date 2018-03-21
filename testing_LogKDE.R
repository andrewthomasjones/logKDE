args <- commandArgs(trailingOnly = TRUE)
print(args[1])

#force install of packages
list.of.packages <- c("Conake","pracma","devtools")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, repos="https://cloud.r-project.org")
library(Conake)
library(pracma)
library(devtools)

list.of.packages2 <- c("VGAM")
new.packages2 <- list.of.packages2[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages2)) devtools::install_version("VGAM", version = "1.0-0", repos = "http://cran.us.r-project.org")

library(VGAM)

#load latest code
install_github("andrewthomasjones/logKDE")

library(logKDE)

#test functions
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


quart_normal<-function(x,q){
  y=1/(1-pnorm(0, mean=1, sd=q))* (dnorm(x, mean = 1, sd=q))
  y[x<=0]<-0
  return(y)
}

#test random generators

rquart_normal<-function(n,q){
  y<-(rnorm(1000*n, mean =1, sd=q))
  y<-y[y>0]
  return(y[1:n])
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

#tests of each Kerneal method

gamma_test<-function(y, bw_method, from, to, se_len, kern=NULL){
  if(bw_method=="CV") {h<-cvbw(y, bw = NULL, ker="GA")$hcv}
  #if(bw_method=="CV"){h<-bw.nrd0(y)}
  if(bw_method=="silverman"){h<-bw.nrd0(y)}
  if(bw_method=="lsilverman"){h<-logKDE::bw.logG(y)}
  if(is.finite(bw_method)){h<-bw_method}

  x<-seq.int(from, to, length.out = length(y))
  est<-dke(vec_data=as.vector(y),"GA",bw=h, x=x)
  return_list<-list(y=est$f_n,x=x, bw=h)

  return(return_list)
}

rig_test<-function(y, bw_method, from, to, se_len, kern=NULL){
  if(bw_method=="CV") {h<-cvbw(y, bw = NULL, ker="RIG")$hcv}
  #if(bw_method=="CV"){h<-bw.nrd0(y)}
  if(bw_method=="silverman"){h<-bw.nrd0(y)}
  if(bw_method=="lsilverman"){h<-logKDE::bw.logG(y)}
  if(is.finite(bw_method)){h<-bw_method}

  x<-seq.int(from, to, length.out = length(y))
  est<-dke(y,"RIG",h, x=x)

  return_list<-list(y=est$f_n, x=x, bw=h)

  return(return_list)
}

norm_test<-function(y, bw_method, from, to, se_len, kern="gaussian"){
  if(bw_method=="CV") {h<-bw.ucv(y)}
  #if(bw_method=="CV") {h<-bw.nrd0(y)}
  if(bw_method=="silverman"){h<-bw.nrd0((y))}
  if(bw_method=="lsilverman"){h<-logKDE::bw.logG(y)}
  if(is.finite(bw_method)){h<-bw_method}

  est<-density(y, bw=h, from=from, to=to, n=se_len, kernel = kern)

  return_list<-list(y=est$y, x=est$x, bw=h)
  return(return_list)
}

log_test<-function(y, bw_method, from, to, se_len, kern="gaussian"){
  #if(bw_method=="CV") {h<-bw.nrd0(log(y))}
  if(bw_method=="CV") {h<-bw.logCV(y)}
  if(bw_method=="silverman"){h<-bw.nrd0(log(y))}
  if(bw_method=="lsilverman"){h<-logKDE::bw.logG(y)}
  if(is.finite(bw_method)){h<-bw_method}

  est<-logdensity(y,bw=h, from=from, to=to, n=se_len, kernel = kern)

  return_list<-list(y=est$y, x=est$x,  bw=h)
  return(return_list)
}

#error functions
MISE<-function(dens, x, q, flist){
  n=length(dens)
  true_dens<-flist(x,q)
  (trapz(x,(dens-true_dens)^2))
}

MIAE<-function(dens, x, q, flist){
  n=length(dens)
  true_dens<-flist(x,q)
  (trapz(x,abs(dens-true_dens)))
}

# test set up ##############################################################################################################################
#lists of each of methods and matching functions
function_list<-list(sin_mad_mix=sin_mad_mix, sin_mad=sin_mad,l_normal=l_normal,  abs_normal=abs_normal,pos_claw=pos_claw, quart_normal=quart_normal)
rfunction_list<-list(rsin_mad_mix=rsin_mad_mix, rsin_mad=rsin_mad,rlnorm2=rlnorm2,  ranorm=ranorm, rpos_claw=rpos_claw,rquart_normal=rquart_normal )

#parameter lists
q_list<-list(rsin_mad_mix=c(0.7,.5,.36), rsin_mad=c(1.145,1.07,.75),rlnorm2=c(0.5,1,2),  ranorm=c(0.5,1,2), rpos_claw=c(0.5,1,2), rquart_normal = c(0.5,1,1.5) )

#sets upper bound of integration region
to_list<-list(rsin_mad_mix=c(5,5,5), rsin_mad=c(5,5,5), rlnorm2=c(5,5,5), ranorm=c(5,5,5),rpos_claw=c(5,5,5), rquart_normal = c(5,5,5))

#methods
method_list<-list(normalkde = norm_test, gamma = gamma_test, rig=rig_test, logKDE= log_test)

#make sure kernel sets go with right methods
kernel_list<-list("gaussian", "epanechnikov", "triangular", "uniform", "laplace", "logistic")
kernel_list2<-list("gaussian", "epanechnikov", "triangular", "rectangular")
kern_list<-list(normalkde = kernel_list2, gamma = "gamma", rig="RIG", logKDE= kernel_list)
#bandwidth methods
bw_list<-list("silverman", "lsilverman",  "CV")
#for speed
#bw_list<-list("silverman", "lsilverman")#  "CV")



#replicates
M=100 #1000

#points per replicate
#n=c(128,256,512,1024)
n=c(20,50,100)
#action
fcount<-length(function_list)
count<-0
big_list<-list()
i <- as.numeric(args[1])
i=2
###for(i in 1:fcount){
  for(k in 1:length(q_list[[i]])){

    q<-q_list[[i]][k]
    from<-0.01
    to<-to_list[[i]][k]

    se_len<-512

    for(j in 1:M){
      tic()
      for(l in 1:length(n)){
      y<-rfunction_list[[i]](n[l],q)

      for(meth in 1:length(method_list)){

        for(bw_method in bw_list){

          for(kern in 1:length(kern_list[[meth]])){
            count<-count+1
            kern_type=kern_list[[meth]][[kern]]
            test_result<-method_list[[meth]](y, bw_method, from, to, se_len, kern_type)
            test_result$MISE <- MISE(test_result$y, test_result$x, q, function_list[[i]])
            test_result$MIAE <- MIAE(test_result$y, test_result$x, q, function_list[[i]])
            results<-paste(names(function_list)[i],names(rfunction_list)[i], j, q, n[l], kern_type, bw_method, names(method_list)[meth], test_result$MISE,test_result$MIAE, test_result$bw, sep=",")

            big_list[[count]]<-results
            gc()
          }
          }
      }

      save(big_list,file=paste0("./logKDEresults_",i,"_", Sys.time(),".RData" ))
      }


      gc()
      time<-toc()
      print(paste0(j, " of ", M, " in ", time, " seconds."))
    }

    #save results
    #save(big_list,file=paste0("./logKDEresults_", Sys.time(),".RData" ))
  }
###}

#save results
save(big_list,file=paste0("./logKDEresults_final_",i,"_", Sys.time(),".RData" ))






#
# # ### post processing ############3##############################################################################################################################
# # #
# zz <- file("res1.csv", "w")
# for(i in 1:length(big_list)){
#   writeLines(big_list[[i]], zz)
# }
# close(zz)
# data1<-read.csv("res1.csv", header=F)
# #names(data1)<-c("Function", "Random", "Replication", "Parameter",  "Kernel", "BWMethod", "Method", "MISE", "MIAE", "BW")
# names(data1)<-c("Function", "Random", "Replication", "Parameter", "Samples", "Kernel", "BWMethod", "Method", "MISE", "MIAE", "BW")


#
#
#
#  MIAE<- aggregate(MIAE~Function + Parameter+Kernel+BWMethod+Method+Samples, data1, mean)


#  MISE<- aggregate(MISE~Function + Samples+Parameter+Kernel+BWMethod+Method, data1, mean)
#  BW<- aggregate(BW~Function + Samples+ Parameter+Kernel+BWMethod+Method, data1, mean)
#  MIAE$MISE<-MISE$MISE
#  MIAE$BW<-BW$BW
#
#  MIAE<-MIAE[order(MIAE$MIAE),]
#
#  #MIAE<-MIAE[order(MIAE$MIAE),]
#
#  dat1<-read.csv("Preliminary.csv")
#  dat2<-read.csv("Preliminary2.csv")
#
#  write.csv(MIAE, file="Prelimary3_samples.csv")


#subset(MIAE, Samples==1024)

#
# x1<-x1[order(x1$V8),]
# x2<-x2[order(x2$V9),]
# x1
# x2
#
#
#
# x2<-read.csv("res1.csv", header=F)
#
# x2$V7<-factor(x2$V7,levels(x2$V7)[c(1,4,2,3)])
#
# x_mad<-subset(x2, x2$V1!="sin_mad_mix")
# x_mad_mix<-subset(x2, x2$V1=="sin_mad_mix")
#
# for(i in 1:3){
# res1<-aggregate(data = split(x_mad,x_mad$V4)[[i]], V8~V6+V7+V5 , FUN=mean)
# res2<-as.data.frame(matrix(res1[,4], byrow=TRUE, ncol = 3))
# res3<-cbind(res1[seq(1, nrow(res1), by=3),2:3],res2)
# res4<-res3[order(res3[,1],res3[,2] ),]
# names(res4)<-c("Package", "Kernel","CV", "Log Silverman", "Silverman" )
# levels(res4$Package) <- c("Conake","Conake", "logKDE", "stats")
# q=levels(factor((x_mad$V4)))[[i]]
# print(xtable(res4, digits=4, caption = paste0("MISE for KDE based on 512 samples drawn from a Singh-Maddala distribution with scale parameter 0.193 and shape parameters 2.8 and ",q, " based on $M=1000$ replications")), include.rownames=FALSE)
# }
#
# for(i in 1:3){
# res1<-aggregate(data = split(x_mad_mix,x_mad_mix$V4)[[i]], V8~V6+V7+V5 , FUN=mean)
# res2<-as.data.frame(matrix(res1[,4], byrow=TRUE, ncol = 3))
# res3<-cbind(res1[seq(1, nrow(res1), by=3),2:3],res2)
# res4<-res3[order(res3[,1],res3[,2] ),]
# names(res4)<-c("Package", "Kernel","CV", "Log Silverman", "Silverman" )
# levels(res4$Package) <- c("Conake","Conake", "logKDE", "stats")
# q=levels(factor((x_mad$V4)))[[i]]
# print(xtable(res4, digits=4, caption = paste0("MISE for KDE based on 512 samples drawn from a mixture of two Singh-Maddala distributions. The mixture consists of a S-M distribution with scale parameter 0.193 and shape parameters 2.8 and 1.7 at a 0.4 mixing proportion, and second with parameters of 0.593, 5.8  and ",q, " based on $M=1000$ replications")), include.rownames=FALSE)
# }
#
#
# # #
# # #
# # #
# # #
# # #
# # #
# # # interaction(dat2$V7, dat2$V5)
# # # logKDE.epanechnikov      logKDE.triangular        logKDE.gaussian         logKDE.uniform        logKDE.logistic         logKDE.laplace     normalkde.gaussian   normalkde.triangular
# # # 0.004220972            0.004276010            0.004432093            0.004582786            0.004755713            0.005600958            0.007091096            0.007239872
# # # normalkde.epanechnikov  normalkde.rectangular            gamma.gamma                rig.RIG
# # # 0.007383662            0.008404307            0.029510654            0.035846674
# #
# #
# #
# #
# #
# #
# #

data<-read.csv("~/Dropbox/Final_set.csv")

levels(data$Kernel)<-  c("Epanechnikov", "gamma"    ,    "Gaussian"    , "Laplace"    ,  "logistic"   ,  "rectangular",  "RIG"   ,       "triangular" , "uniform")
levels(data$BWMethod)<-c("CV", "log-Silverman" ,"Silverman")

library(reshape)

data$Method<- factor(data$Method,levels(data$Method)[c(2,3,4,1)])
levels(data$Method) <- c("logKDE", "stats", "Conake", "Conake")


data_20 <- subset(data, data$Samples==20)
data_50 <- subset(data, data$Samples==50)
data_100 <- subset(data, data$Samples==100)

data_20_b<-melt(data_20, id=c("Method", "Kernel", "BWMethod", "Function", "Parameter"), measure.vars =c("MIAE"))
data_20_c<-melt(data_20, id=c("Method", "Kernel", "BWMethod", "Function", "Parameter"), measure.vars =c("MISE"))


data_20_b2<-cast(data_20_b,  Method + Kernel + BWMethod ~ Function + Parameter, fill=F)
data_20_c2<-cast(data_20_c,  Method + Kernel + BWMethod ~ Function + Parameter, fill=F)


data_100_b<-melt(data_100, id=c("Method", "Kernel", "BWMethod", "Function", "Parameter"), measure.vars =c("MIAE"))
data_100_c<-melt(data_100, id=c("Method", "Kernel", "BWMethod", "Function", "Parameter"), measure.vars =c("MISE"))


data_100_b2<-cast(data_100_b,  Method + Kernel + BWMethod ~ Function + Parameter, fill=F)
data_100_c2<-cast(data_100_c,  Method + Kernel + BWMethod ~ Function + Parameter, fill=F)






b2i<-sapply(data_20_b2[,4:12], function(x) order(x)[1:5])
c2i<-sapply(data_20_c2[,4:12], function(x) order(x)[1:5])
b2i00<-sapply(data_100_b2[,4:12], function(x) order(x)[1:5])
c2i00<-sapply(data_100_c2[,4:12], function(x) order(x)[1:5])


data_20_b2[,4:12] <- format(data_20_b2[,4:12], digits=4)
data_20_c2[,4:12] <- format(data_20_c2[,4:12], digits=4)
data_100_b2[,4:12] <- format(data_100_b2[,4:12], digits=4)
data_100_c2[,4:12] <-format(data_100_c2[,4:12], digits=4)



for(i in 4:12){
  for (j in 1:5){
    data_20_b2[b2i[j, i-3],i] <- paste0("\\textbf{",data_20_b2[b2i[j, i-3], i],"}")
    data_20_c2[c2i[j, i-3],i] <- paste0("\\textbf{",data_20_c2[c2i[j, i-3],i],"}")
    data_100_b2[b2i00[j, i-3],i] <- paste0("\\textbf{",data_100_b2[b2i00[j, i-3],i],"}")
    data_100_c2[c2i00[j, i-3],i] <- paste0("\\textbf{",data_100_c2[c2i00[j, i-3],i],"}")
  }
}

library(xtable)
print(xtable(data_20_b2, digits=4), include.rownames=FALSE, sanitize.text.function = function(x) x)
print(xtable(data_20_c2, digits=4), include.rownames=FALSE, sanitize.text.function = function(x) x)
print(xtable(data_100_b2, digits=4), include.rownames=FALSE, sanitize.text.function = function(x) x)
print(xtable(data_100_c2, digits=4), include.rownames=FALSE, sanitize.text.function = function(x) x)


library(devtools)
install_github("andrewthomasjones/logKDE")
library(logKDE)

chisq10<-rchisq(100,10)

fit1<-logdensity(chisq10)
plot(fit1, ylim=c(0, .1))

print(fit1)

fit2<-logdensity(chisq10, bw ="logCV", kernel = "triangular")
plot(fit2, ylim=c(0, .1))
grid(10,10,2)
x<-(seq(min(chisq10), max(chisq10), 0.01))
lines(x, dchisq(x,10), col =4)




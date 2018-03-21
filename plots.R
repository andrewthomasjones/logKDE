data1<-read.csv("~/Dropbox/baseball.csv", header=F)
x<-data1$V1
library(ggplot2)
library(logKDE)

y1<-logKDE::logdensity(x,  from=0.0001, to=6500)$y
x1<-logKDE::logdensity(x, from=0.0001, to=6500)$x
x2<-stats::density(x, from=0.0001, to=6500)$x
y2<-stats::density(x, from=0.0001, to=6500)$y

plotframe<-data.frame(x1=x1, x2=x2,y1=y1,y2=y2)

p1<-ggplot(data=data1, aes(x=V1))+geom_histogram(aes(y=..density..),binwidth = 100, alpha=0.8)+theme_bw(base_size =18)
p1<-p1+ggtitle("1992 Major League Baseball Salaries")+ylab("Frequency")+xlab("Salary in $1000's")
p1<-p1+geom_line(data=plotframe, aes(x=x2, y=y2, colour="Standard density"), size=1.1)
p1<-p1+geom_line(data=plotframe, aes(x=x1, y=y1, colour="logKDE density"), size=1.1)
p1+scale_color_discrete(name="Density")+theme(legend.position= c(.7,.6))

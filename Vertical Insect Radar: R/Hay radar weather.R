setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(ggplot2)
library(circular)
library(akima)
df1 = read.table("Ha2022Jan26.df",header=T,skip=2)
df1[df1$ft=="M","a0dB"] = df1[df1$ft=="M","a0dB"]-4
df1$hour = (df1$hr*60 + df1$tm)/60
df1$duration_h = with(df1,ifelse(ft=="S",durSmh,durMmh))
df1$flux_contribution = with(
  df1,
  1 / ((1+max(hr)-min(hr)) * 60 * duration_h * ebw)
)
df1$density_contribution = df1$flux_contribution/df1$sp
df1$dca_norm = 2*abs(df1$dca)/df1$ebw
df1[df1==-101] = NA
df1$deltaX = df1$sp * sin(df1$ga*pi/180)
df1$deltaY = df1$sp * cos(df1$ga*pi/180)
wind1 = read.table("metdat-2022-01-26.dat",header=T,skip=1)
wind1$windir = (wind1$windir-180)%%360
wind1$deltaX = wind1$winsp * sin(wind1$windir*pi/180)
wind1$deltaY = wind1$winsp * cos(wind1$windir*pi/180)

iterative_rayleigh = function(angles,deg=T){
  if (!deg==T){angles = angles*180 / pi}
  angles_binned = cut(angles,seq(0,360,5))
  levels(angles_binned) = seq(2.5,357.5,5)
  t1= table(angles_binned)
  degs = circular(as.numeric(rownames(t1)),template="geographics",units="degrees")
  for (i in 1:max(t1)){
    if (rayleigh.test(degs[t1>=i])$p.value<0.05){
      nsig=i
      break
    }
  }
  sum(t1[t1>=nsig] +1 -nsig)/sum(t1)
  return(sum(t1[t1>=nsig] +1 -nsig)/sum(t1))
}

# 2D interpolation of wind vectors
library(metR)
ggplot(wind1,aes(x=ih,y=ht))+
  geom_arrow(aes(mag=winsp,angle=windir),size=0.1) 
# https://stackoverflow.com/questions/47184572/2d-linear-interpolation-in-r
?interpp
sapply(wind1,class)
# Testing 2D interpolation
windInterp = with(wind1,interpp(x=ih,y=ht,z=deltaY,linear=F,extrap=T,
                                xo = rep(seq(min(ih),max(ih),length=100),100),
                                yo = rep(seq(min(ht),max(ht),length=100),each=100)))
ggplot(wind1,aes(x=ih,y=ht,size=deltaY,color=deltaY))+
  geom_point()
fields::quilt.plot(windInterp$x, windInterp$y, windInterp$z)

# Making wind vector interpolation for all data points
df1$windX = with(wind1,interpp(x=ih,y=ht,z=deltaX,xo=df1$hour,yo=df1$zm),linear=F,extrap=T)$z
df1$windY = with(wind1,interpp(x=ih,y=ht,z=deltaY,xo=df1$hour,yo=df1$zm),linear=F,extrap=T)$z
df1$windDir = (atan2(df1$windX,df1$windY)*180/pi)%%360
df1$windSp = sqrt(df1$windX^2+df1$windY^2)
df1$insectX = df1$deltaX - df1$windX
df1$insectY = df1$deltaY - df1$windY
df1$insectDir = (atan2(df1$insectX,df1$insectY)*180/pi)%%360
df1$insectSp = sqrt(df1$insectX^2+df1$insectY^2)
# Calculate how much mean direction disagrees with wind direction (0-1)
df1$ins_wind_diff_lin = with(df1,pmin((ga-insectDir)%%360,(insectDir-ga)%%360))/180
df1$ins_wind_diff_circ = (cos(df1$ins_wind_diff_lin*pi)+1)/2

# Investigating insect-specific directions and how they contrast with the wind and overall
circ = list("ga"=circular(df1$ga,template="geographics",unit="deg"))
circ$wind = circular(df1$windDir,template="geographics",unit="deg")
circ$insect = circular(df1$insectDir,template="geographics",unit="deg")

par(mfrow=c(1,3))
for (ty in names(circ)){rose.diag(circ[ty],bins=36,prop=1,col="lightblue",ylab=ty,
                                  shrink=1)}
dev.off()


library(plotrix)
plot(insectSp~a0dB,data=df1,col=rgb(red=0.5,green=0.5,blue=0.5,alpha=0.2))

windagg = bin.wind.records(df1$insectDir,df1$insectSp,ndir = 12)
oz.windrose(windagg,scale.factor = 30)

windrose(circ$insect,df1$insectSp)

source("WindRose.R")
plot.windrose(data=df1,spd="insectSp",dir="insectDir",spdmin=0,title="Insect",spdmax=30)
plot.windrose(data=df1,spd="windSp",dir="windDir",spdmin=0,title="Wind")
plot.windrose(data=df1,spd="sp",dir="ga",spdmin=0,title="Echo")
# plot.windrose(data=wind1,spd="winsp",dir="windir",spdmin=0,title="Wind")

dev.off()
png(file="RGraphWindroseInsect.png",res=700,width=4800,height=3200,pointsize=14)


# Iterative Rayleigh
iterative_rayleigh(df1$ga)
iterative_rayleigh(df1[df1$sp>6,"ga"])
iterative_rayleigh(df1$insectDir)
iterative_rayleigh(df1[df1$insectSp>6,"insectDir"])

sapply(wind1,class)

# See if insectDir and disagreement varies with season


# Miscellaneous radar data exploration
# Wingbeat frequency
hist(df1$ff)
sum(df1$ff==-101)/nrow(df1)
# Replace -101 values with NA

#####
# Basic visualisation ####
plot(insectSp~al2,data=df1,col=rgb(red=0.5,green=0.5,blue=0.5,alpha=0.2))

plot(al4~al2,data=df1,col=rgb(red=0.5,green=0.5,blue=0.5,alpha=0.2))
ggplot(df1,aes(x=al2,y=al4,color=a0dB))+
  geom_point(shape=1,size=0.8,alpha=0.3)


ggplot(df1,aes(x=al4,y=insectSp,color=al2))+
  geom_point(shape=1,size=0.8,alpha=0.3)

# Statistical analysis ####
sapply(df1,class)
# Insect speed
hist(df1$insectSp)
hist(df1$a0dB)

hist(df1$al2)
hist(df1$al4)

m1 = lm(insectSp ~ a0dB,data=df1)
m1 = lm(insectSp ~ a0dB+al2+al4+zm,data=df1)
summary(m1)
AIC(m1)
plot(m1)
newdat = with(df1,expand.grid(a0dB=seq(min(a0dB),max(a0dB),length=100)))
preds = predict(m1,newdata=newdat,se.fit=T,type="response")
newdat[,c("fit","se")] = preds[1:2]
newdat = with(newdat,data.frame(newdat,insectSp = fit,SEL = fit-se,SEU = fit+se))


m1 = glm(insectSp ~ a0dB+al2+al4+zm,family=Gamma(link=identity),data=df1)
summary(m1)
plot(m1)
with(summary(m1),1-deviance/null.deviance)
newdat = with(df1,expand.grid(a0dB=seq(min(a0dB),max(a0dB),length=100)))
preds = predict(m1,newdata=newdat,se.fit=T,type="link")
newdat[,c("fit","se")] = preds[1:2]
newdat = with(newdat,data.frame(newdat,insectSp = exp(fit),
                                SEL = exp(fit-se),SEU = exp(fit+se)))

ggplot(df1,aes(x=a0dB,y=insectSp,color=ft))+
  geom_point(shape=1,size=0.8,alpha=0.3)+
  #geom_ribbon(data=newdat,linetype=0,alpha=0.9,aes(ymin=SEL,ymax=SEU),color="grey")+
  #geom_line(data=newdat,color="black")+
  labs(y="Insect-specific speed (m/s)")

# Circular lm


ggplot(df1,aes(x=a0dB,y=ff,color=al4))+
  geom_point(shape=1,size=0.8,alpha=0.3)

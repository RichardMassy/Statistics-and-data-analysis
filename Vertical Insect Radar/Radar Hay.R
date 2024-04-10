# 0 Preparation and function load ####
setwd("C:/Users/rm669/OneDrive - University of Exeter/Data and analysis/23 Radar")
library(circular)
library(ggplot2)

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

#####
# 1.0 Visualise Hay 2023 radar data OLD ####
# df1 = read.csv("Ha2023Jan26.csv")
# sapply(df1,class)
# hist(df1$hr)
# range(df1$hr)
# hist(df1$nfh)
# range(df1$nfh)
# table(df1$hr)
# table(df1$nfh)
# plot(nfh~hr,data=df1,ylim=c(0,13))
# df2 = aggregate(nfh~hr,length,data=df1)
# plot(nfh~hr,data=df2,ylim=c(0,max(df2$nfh)))

#####
# 1.1 New open using read.table() ####
df1 = read.table("Ha2023Jan26night.df",header=T,skip=2)
df1 = read.table("Ha2023Jan26.df",header=T,skip=2)
df1 = df1[df1$hr>=0&df1$ft=="S",]
# df1[df1$ft=="M","a0dB"] = df1[df1$ft=="M","a0dB"]-4
df1$hour = (df1$hr*60 + df1$tm)/60
df1$duration_h = with(df1,ifelse(ft=="S",durSmh,durMmh))
df1$flux_contribution = with(
  df1,
  1 / ((1+max(hr)-min(hr)) * 60 * duration_h * ebw)
  )
sum(df1$flux_contribution) # Insects per vertical metre per second
# Old flux contribution:
# (60/duration_h) / ((1+max(hr)-min(hr)) * 3600 * (ebw/mean(ebw)))
df1$density_contribution = df1$flux_contribution/df1$sp
df1$dca_norm = 2*abs(df1$dca)/df1$ebw

# Verification of a few values
sapply(df1,class)
mean(df1$ebw)
mean(df1$duration_h)
# Create small subsample
df1 = df1[seq(1,21000,100),]

# Exploring intermediate minutes - 13 categories
# Missing: 10,13,25,28,31,40,43,55 & 58, extra: 31 & 45
table(df1$tm)


hist(df1$ga) # Direction
hist(df1$zm) # Height
min(df1$zm)
hist(df1$a0dB) # Cross section


# Verifying ebw dca
plot(dca~ebw,data=df1)
abline(0,1)
abline(0,-1)


hist(df1$dca_norm,breaks=18)

plot(zm~abs(dca),data=df1,col=rgb(red=0.5,green=0.5,blue=0.5,alpha=0.2))

png(file="RGraphHistDcaNormSM.png",res=700,width=4800,height=3200,pointsize=14)
ggplot(df1,aes(x=dca_norm,fill=zm<400))+
  geom_histogram(color="black",position="stack",boundary=0)
dev.off()

# Size ~ height
plot(a0dB~zm,data=df1,col=rgb(red=0.5,green=0.5,blue=0.5,alpha=0.2))
plot(zm~a0dB,data=df1,col=rgb(red=0.5,green=0.5,blue=0.5,alpha=0.2),ylim=c(0,max(df1$zm)))
png(file="RGraphHeightSizeShape.png", res=700, width = 4800, height=3200,
    pointsize=14,type="windows", antialias="cleartype")
ggplot(data=df1,aes(x=a0dB,y=zm,color=al2))+
  geom_point(alpha=0.1,stroke=0)
dev.off()

# Comparing S & M (y = height)####
png(file="RGraphHeightSizePulseM-4.png", res=700, width = 4800, height=3200,
    pointsize=14,type="windows", antialias="cleartype")
# Size
# df1[df1$ft=="M","a0dB"] = df1[df1$ft=="M","a0dB"]-4
ggplot(data=df1,aes(x=a0dB,y=zm,color=ft))+
  geom_point(alpha=0.1,stroke=0)
# Time
ggplot(data=df1,aes(x=hour,y=zm,color=ft))+
  geom_point(alpha=0.1,stroke=0)
# Speed
ggplot(data=df1,aes(x=sp,y=zm,color=ft))+
  geom_point(alpha=0.1,stroke=0)
dev.off()

# Circular
ga_circ = circular(df1$ga,template="geographics",unit="deg")
rose.diag(ga_circ,bins=36,prop=2.2,col="lightblue")
mean.circular(ga_circ)
rayleigh.test((ga_circ))
sd(ga_circ)*180/pi

# Directedness
iterative_rayleigh(df1$ga)

# Directedness ~ size
sapply(seq(-20,10,10),function(sizeBand){
  iterative_rayleigh(df1[df1$a0dB>sizeBand&df1$a0dB<sizeBand+10,"ga"])
})
png(file="RGraphSizeDirection.png", res=700, width = 4800, height=3200,
    pointsize=14,type="windows", antialias="cleartype")
ggplot(data=df1,aes(x=a0dB,y=(ga+180)%%360,color=al2))+
  geom_point(alpha=0.1,stroke=0)+
  scale_y_continuous(breaks=seq(0,360,45),
                     labels=c("S","SW","W","NW","N","NE","E","SE","S"))
dev.off()

# Directedness ~ height
dft = as.data.frame(t(sapply(seq(150,1200,150),function(heightUnit){
  subset1 = df1[df1$zm>heightUnit&df1$a0dB<heightUnit+10,"ga"]
  it_ray = iterative_rayleigh(subset1)
  directedness = rayleigh.test(circular(subset1,template="geographics",unit="deg"))$statistic
  c(heightUnit+75,directedness,it_ray)
})))
names(dft) = c("Height","r","Prop_directed")

plot(r~Height,data=dft,ylim=c(0.8,1),type="b",col="red")
lines(Prop_directed~Height,data=dft,col="blue",type="b")
text(x=300,y=c(0.97,0.9),adj=0,
     labels=c("Proportion of directed movement","Directedness (r)"))

# Direction ~ size
plot(ga~a0dB,data=df1,col=rgb(red=0.5,green=0.5,blue=0.5,alpha=0.2))
# Direction ~ height
plot(ga~zm,data=df1,col=rgb(red=0.5,green=0.5,blue=0.5,alpha=0.2))
# Direction ~ time
plot(ga~hour,data=df1,col=rgb(red=0.5,green=0.5,blue=0.5,alpha=0.2))
png(file="RGraphDirectionTime.png", res=700, width = 4800, height=3200,
    pointsize=14,type="windows", antialias="cleartype")
ggplot(data=df1,aes(x=hour,y=(ga+180)%%360))+
  geom_point(alpha=0.1,stroke=0)+
  scale_y_continuous(breaks=seq(0,360,45),
                     labels=c("S","SW","W","NW","N","NE","E","SE","S"))
dev.off()
# Height ~ time
plot(zm~hour,data=df1,col=rgb(red=0.5,green=0.5,blue=0.5,alpha=0.2))
# size ~ time
plot(a0dB~hour,data=df1,col=rgb(red=0.5,green=0.5,blue=0.5,alpha=0.2))
# Speed ~ size
plot(sp~a0dB,data=df1,col=rgb(red=0.5,green=0.5,blue=0.5,alpha=0.2))
# beta ~ gamma alignment
plot(ga~be,data=df1,col=rgb(red=0.5,green=0.5,blue=0.5,alpha=0.2))

# Abundance ~ time
df2 = aggregate(zm~hr,length,data=df1)
colnames(df2) = c("Hour","Echos")
plot(Echos~Hour,data=df2,ylim=c(0,max(df2$Echos)))

df2$Directedness = sapply(df2$Hour,function(hour){
  iterative_rayleigh(df1[df1$hr==hour,"ga"])
})
#####
# Kernel density ####
hist(df1$flux_contribution)
hist(df1$flux_con_norm)

d1 = as.data.frame(density(df1$zm,weights = df1$flux_contribution)[c("x","y")])
plot(d1)

png(file="RGraphHeightFlux.png", res=700, width = 4800, height=3200,
    pointsize=14,type="windows", antialias="cleartype")
ggplot(data=d1,aes(x=y*3600,y=x))+
  geom_polygon(alpha=0.2)+
  geom_path()+
  labs(x=bquote("Insects transiting m"^-2 ~ h^-1),y="Altitude (m)",
       title="Insects detected by short-pulse observations 26/01/2023")
dev.off()


d2 = as.data.frame(density(df1$zm,weights = df1$density_contribution)[c("x","y")])
plot(d2)
png(file="RGraphHeightDensity.png", res=700, width = 4800, height=3200,
    pointsize=14,type="windows", antialias="cleartype")
ggplot(data=d2,aes(x=y*150^3,y=x))+
  geom_polygon(alpha=0.2)+
  geom_path()+
  labs(x=bquote("Insect density 150 m"^-3),y="Altitude (m)",
       title="Insects detected by short-pulse observations 26/01/2023")
dev.off()

# with S + M
d3 = do.call(rbind,lapply(list("S","M"),function(ft){
  df = df1[df1$ft==ft,]
  dens = suppressWarnings(density(df$zm,weights=df$flux_contribution,
                                  from=0,to=2500)[c("x","y")])
  data.frame(dens,"ft"=ft)
}))
png(file="RGraphHeightFluxAll.png", res=700, width = 4800, height=3200,
    pointsize=14,type="windows", antialias="cleartype")
ggplot(data=d3,aes(x=y*3600,y=x,fill=ft))+
  geom_polygon(alpha=0.2)+
  geom_path()+
  labs(x=bquote("Insects transiting m"^-2 ~ h^-1),y="Altitude (m)",
       title="Insects detected on 26/01/2023",fill="Pulse")
dev.off()
#####
# 2.1 Iterative Rayleigh test for calculating proportion of directional movement ####
# Generate directional angles (kappa 2.7 produces ~ 40° circular standard deviation)
angles = list("dir"=rvonmises(1000,mu=circular(pi,template="geographics",
                                               modulo="2pi"),kappa=2.7))
round(sd(angles$dir)*180/pi,digits=3)
hist(as.numeric(angles$dir))
rose.diag(angles$dir,bins=36,prop=2.2,axes=F,ticks=F,col="lightblue")

# Pad with non-directional angles
angles$nondir = circular(runif(1000,min=0,max=2*pi),template="geographics",modulo="2pi")
hist(as.numeric(angles$nondir))
rose.diag(angles$nondir,bins=36,prop=2.2,axes=F,ticks=F,col="lightblue")

# Combine directional and non-directional angles
angles$all = c(angles$dir,angles$nondir)
hist(as.numeric(angles$all))
rose.diag(angles$all,bins=36,prop=2.2,axes=F,ticks=F,col="lightblue")

# Allocate angles into 5° bins (clunky but works)
angles$binned = cut(angles$all,seq(0,2*pi,pi/36))
levels(angles$binned) = seq(2.5,357.5,5) # rename bin factors to corresponding mid angles

t1= table(angles$binned)
degs = circular(as.numeric(rownames(t1)),template="geographics",units="degrees")

# Iterative Rayleigh test. Records n number when p <0.05
for (i in 1:max(t1)){
  if (rayleigh.test(degs[t1>=i])$p.value<0.05){
    nsig=i
    break
  }
}
mov_directional = sum(t1[t1>=nsig] +1 -nsig)
mov_non_directional = sum(t1) - mov_directional
print(paste("Proportion of directional movement =",mov_directional/sum(t1)))


# The above but wrapped into a function to visualise angles with result
iterative_rayleigh_test = function(dir=pi,ndir=1000,nnondir=1000,kapp=2.7,p=0.05,plot=T){
  ang_dir = rvonmises(ndir,mu=circular(dir,template="geographics",modulo="2pi"),kappa=kapp)
  ang_nondir = circular(runif(nnondir,min=0,max=2*pi),template="geographics",modulo="2pi")
  ang_all = c(ang_nondir,ang_dir)
  ang_bin = cut(ang_all,seq(0,2*pi,pi/36))
  levels(ang_bin) = seq(2.5,357.5,5)
  t1= table(ang_bin)
  degs = circular(as.numeric(rownames(t1)),template="geographics",units="degrees")
  nsig = max(t1)
  for (i in 1:max(t1)){
    if (rayleigh.test(degs[t1>=i])$p.value < 0.05){
      nsig=i
      break
    }
  }
  mov_dir_proportion = sum(t1[t1>=nsig] +1 -nsig)/sum(t1)
  
  if (plot==T){
    rose.diag(ang_all,bins=36,prop=3,axes=F,ticks=F,col="lightblue",
              main=paste("Directional movement:",mov_dir_proportion),
              sub=paste("Directional proportion =",ndir/(ndir+nnondir)))
  }
  return(c(mov_dir_proportion, mean.circular(ang_all)))
}

# Run test for 10 levels of directedness, producing a graph for each
for (directedness in seq(0,1,0.1)){
  iterative_rayleigh_test(ndir=directedness*1000,nnondir=(1-directedness)*1000)
}

#####
# 2.2 Plotting the iterative Rayleigh test for a simulated migration season ####

# Generate data for a 180 day season
days = 1:180
# Set directedness to a normal distribution with a peak at day 90
directedness = dnorm(seq(1:180),mean=90,sd=15)
# Add small amount of random jitter
directedness = jitter(directedness,amount=0.005)
# Normalise to 0-1
directedness = (directedness-min(directedness))/(max(directedness)-min(directedness))
plot(directedness)

# Add bonus to directedness for migration season
directedness = directedness + c(rep(0,60),rep(0.4,60),rep(0,60))
# Normalise again
directedness = (directedness-min(directedness))/(max(directedness)-min(directedness))
plot(directedness)

# Generate random variation to mean direction (so pi when directed, random when not)
plot(pi+(1-directedness)*runif(length(days),min=-1,max=1)*pi)

# Apply iterative_rayleigh_test to directedness values with added direction jitter
res = t(sapply(directedness,function(x){
  dir = pi+(1-x)*runif(1,min=-1,max=1)*pi
  iterative_rayleigh_test(ndir=1000*x,nnondir=100*(1-x),plot=F,dir=dir)
  }))

df3 = data.frame(days,directedness,res)
colnames(df3)[3:4] = c("Computed_directedness","Mean_direction")

# Plot results
hist(df3$Computed_directedness)
ggplot(data=df3,aes(x=days,y=Mean_direction,color=Computed_directedness))+
  annotate("rect",xmin=min(days),xmax=max(days),ymin=pi-(pi/8),ymax=pi+(pi/8),alpha=0.2)+
  geom_point(size=2)+
  scale_color_continuous(type="viridis",begin=0,end=0.9,direction=-1,
                         breaks=c(0,0.5,1),labels=c("0%","50%","100%"))+
  theme_test()+
  theme(legend.position="bottom")+
  scale_x_continuous(breaks=seq(0,150,30),
                     labels=c("Jul","Aug","Sep","Oct","Nov","Dec"))+
  scale_y_continuous(breaks=c(0,0.5*pi,0.75*pi,pi,1.25*pi,1.5*pi,2*pi),
                     labels=c("N","E","SE","S","SW","W","N"))+
  labs(y=element_blank(),x=element_blank(),color="Proportion of Directional Movement")

# Export as .png
png(file="RGraphDirectionDirectednessTest.png", res=700, width = 4800, height=3200,
    pointsize=14,type="windows", antialias="cleartype")
# (Run ggplot)
dev.off()

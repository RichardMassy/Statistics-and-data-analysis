# 0: Data preparation - functions and flies ####
library(ggplot2)
library(glmmTMB)
setwd("~/Data S2")
cbPalette = c("#000000","#E69F00","#56B4E9","#228B22","#F0E442","#0072B2",
              "#D55E00","#fc8d59","#8da0cb","#CC79A7","#DAA520","#990000")
AICc = function(model){p = length(model$coefficients)
  return((-2*logLik(model))+((2*p*(p+1))/(nobs(model)-p-1)))}
Model_output = function(model,file=NULL){
  if (is.null(file)) file = "Model output.txt"
  co = summary(model)$coefficients
  if("cond"%in%names(co)) co = co$cond
  if(all(c("count","zero")%in%names(co))) co = rbind(co$count,co$zero)
  p = co["(Intercept)",4]
  estimate = format(round(co["(Intercept)","Estimate"],2),nsmall=2)
  se = format(round(co["(Intercept)","Std. Error"],2),nsmall=2)
  sig=ifelse(p<0.001,"***",ifelse(p<0.01,"**",ifelse(p<0.05,"*",ifelse(p<0.1,"·",""))))
  results = list(paste0(format(round(co["(Intercept)","Estimate"],2),nsmall=2),
                        "±",se,sig))
  for (c in rownames(co)[2:nrow(co)]){
    p = co[c,4]
    estimate = format(round(co[c,"Estimate"],2),nsmall=2)
    se = format(round(co[c,"Std. Error"],2),nsmall=2)
    sig=ifelse(p<0.001,"***",ifelse(p<0.01,"**",ifelse(p<0.05,"*",ifelse(p<0.1,"·",""))))
    results = append(results,paste0(" + ",estimate,"±",se," x ",c,sig))
  }
  output = cat(deparse(substitute(model)),paste(results,collapse=""),
               file=file,append=T,sep="\n")
}
Newdat_output = function(digits=1,df=newdat,file="Newdat output.csv"){
  df[,c("fit","se")]=format(round(df[,c("fit","se")],digits=digits),nsmall=digits)
  df$output = paste(df$fit,df$se,sep="±")
  write.csv(df,file)
}

prep = function(dataSource,flyType=NULL){
  if (dataSource %in% c("Flight_stats","All flights","Transects")){
    df = read.csv(paste0("RD ",dataSource,".csv"))
  } else {df = read.csv(dataSource)}
  if (dataSource == "Transects"){
    df["Speed"] = df[["Events"]]*0.03534292
    df["Flying"] = !is.na(df[["Events"]])
  }#else{df$Direction = factor(df$dom_dir,levels=c("F","R"),
            #                  labels=c("Anti-clockwise","Clockwise"))}
  df$Direction=factor(df$dom_dir,levels=c("F","R"),
                      labels=c("Anti-clockwise","Clockwise"))
  df$Exp_date = as.Date(substr(df$Fly,1,4),format = "%m%d")
  df$Ecl_date = as.Date(gsub("(.*)(.{2})","\\1-\\2",df$Ecl_date),format="%m-%d")
  df$Age = df$Exp_date - df$Ecl_date
  df$Exp_start = ((df$Exp_start-(df$Exp_start%%100))/400) + (df$Exp_start%%100)/240
  df[,c("Sex","Condition","Size")] = do.call(rbind,strsplit(df$SexWeightSize,NULL))
  df$Condition = ordered(df$Condition,levels = c("t","m","f"))
  df$Size = ordered(df$Size,levels = c("s","a","l"))
  df$SourceMorph = paste0(df$Source,df$Morph)
  factor_cols = c("Fly","Source","Morph","SourceMorph","Sex")
  df[factor_cols] = lapply(df[factor_cols],as.factor)

  # Optional: replace factor letter codes with words
  levels(df$Size) = c("Small","Average","Large")
  levels(df$Condition) = c("Thin","Medium","Fat")
  levels(df$Sex) = c("Female","Male")
  levels(df$Morph) = c("Autumn","Summer")
  df$Morph = relevel(df$Morph,"Summer")
  if (!is.null(flyType)){
    df = df[df$Source==flyType,]
  }
  return (df)
}
df1 = prep("Flight_stats",flyType="W")
df2 = prep("All flights",flyType="W")
df3 = prep("Transects",flyType="W")


#####
# 1.0: Overall analysis ####
sapply(df1,class)
colnames(df1)

# Descriptives
hist(log(df1[df1$Morph=="Autumn",]$Distance),breaks=12)
hist(df1[df1$Morph=="Autumn",]$Distance,breaks=12)
hist(df1[df1$Morph=="Autumn",]$Flight_time)
shapiro.test(log(df1$Distance))
hist(log(df1$Mean_speed),breaks=12)
shapiro.test(log(df1$Mean_speed))
hist(df1[df1$Distance<2000,"Distance"])
aggregate(Age~Morph,data=df1,median) # 0 days summer, 4 days autumn
table(df1$Condition,df1$Size)
table(df1$Condition,df1$Morph)
qqnorm(log(df1$Distance))
qqline(log(df1$Distance),col="blue",lwd=2)
aggregate(Distance~Morph+Sex,data=df1,function(x){c(median(x),"N"=length(x))})
aggregate(log(Distance)~Morph+Sex,data=df1,function(x){c(exp(mean(x)),"N"=length(x))})
aggregate(log(Distance)~Condition+Direction,data=df1,mean)
range(df1$Mean_speed)
aggregate(Distance~Morph + Direction + Condition,data=df1,mean)
boxplot(log(Distance)~Sex,data=df1[df1$Morph=="Summer",])
# Test of difference in variance between males and females
# Different result between versions of data extraction
# (error checking more robust = fewer clipped flights)
var.test(log(df1[df1$Morph=="Summer"&df1$Sex=="Male","Distance"]),
         log(df1[df1$Morph=="Summer"&df1$Sex=="Female","Distance"]))
var.test(df1[df1$Morph=="Summer"&df1$Sex=="Male","Distance"],
         df1[df1$Morph=="Summer"&df1$Sex=="Female","Distance"])
var(df1[df1$Morph=="Summer"&df1$Sex=="Male","Distance"])/var(df1[df1$Morph=="Summer"&df1$Sex=="Female","Distance"])

# Distance full model
m1 = glm(Distance ~ Morph + Sex + Direction + Condition + Size
         + Morph:Condition + Morph:Sex,
         family=Gamma(link=log),data=df1)
summary(m1)
# Final model. D.O.: Morph:Condition -3.9, Size -3.6, Morph:Sex -0.1, Sex -1.3 USED
m1 = glm(Distance ~ Morph + Direction + Condition,family=Gamma(link=log),data=df1)
summary(m1)
drop1(m1)

plot(m1)
confint(m1)
plot(resid(m1)~df1$Condition)

# Pseudo R^2
with(summary(m1),1-deviance/null.deviance)
# Predict
newdat = with(df1,expand.grid(Morph=unique(Morph),Condition=unique(Condition),
                               Direction=unique(Direction)))
preds = predict(m1,newdata=newdat,se.fit=T,type="response")
newdat[,c("fit","se")] = preds[1:2]

# Compare model fit to aggregate
agg = merge(newdat,aggregate(Distance~Morph+Direction+Condition,data=df1,function(x){
  c(Mean=mean(x),N=length(x))}))

png(file="RGraphFlyDist4.png", res=900, width = 4800, height=3200,
    pointsize=14,type="windows", antialias="cleartype")#res700 plot 900 model
ggplot(data=agg,aes(x=Distance[,"Mean"],y=fit))+
  geom_abline()+
  labs(x="Group mean distance (m)",y="Model estimate + standard error",
       color="N (mean)",subtitle="Distance ~ Morph + Condition + Direction")+
  geom_pointrange(aes(ymin=fit-se,ymax=fit+se,
                      color=Distance[,"N"]),fatten=1)
dev.off()
# Ordinal variable contrasts are multiplied by ordinal variables
contr.poly(3)

# Female only
m1f = glm(Distance ~ Morph + Direction + Condition,
          family = Gamma(link=log),data=df1[df1$Sex=="Female",])
summary(m1f)

# Mean speed ####
hist(df1$Mean_speed)
range(df1$Mean_speed)
hist(log(df1[,"Mean_speed"]))
qqnorm(log(df1$Mean_speed))
qqline(log(df1$Mean_speed),col="blue",lwd=2)
table(df1$Morph,df1$Sex)
table(df1$Condition,df1$Sex)
table(df1$Morph,df1$Condition)
mean(df1$Mean_speed)
aggregate(Mean_speed ~ Morph + Sex, data=df1,median)
aggregate(Mean_speed ~ Morph, data=df1,median)
# Weighted aggregate
df1$sw = df1$Mean_speed*df1$Flight_time
agg = aggregate(sw~Morph,data=df1,sum)
agg2 = aggregate(Flight_time~Morph,data=df1,sum)
agg$weighted_speed = agg$sw/agg2$Flight_time

# Full model
m11 = glm(Mean_speed ~ Morph + Sex + Direction + Condition + Size
         + Morph:Condition + Sex:Morph,
         family = Gamma(link=log),data=df1)
summary(m11)

# Final model.D.O.: M:C-3.7,Cond-2.8,M:Sex-1.5,Morph-1.6,Sex-0.47 USED
m11 = glm(Mean_speed ~ Direction + Size,
          family = Gamma(link=log),data=df1)
summary(m11)
drop1(m11)
plot(m11)

# Pseudo R^2
with(summary(m11),1-deviance/null.deviance)

# Predict
newdat = with(df1,expand.grid(Size=unique(Size),Direction=unique(Direction)))
preds = predict(m11,newdata=newdat,se.fit=T,type="response")
newdat[,c("fit","se")] = preds[1:2]
# newdat$fitR = exp(newdat$fit)
# newdat$seR = exp(newdat$fit + newdat$se)-newdat$fitR

# Compare model fit to aggregate
agg = merge(newdat,aggregate(Mean_speed~Size+Direction,data=df1,function(x){
  c(Mean=mean(x),N=length(x))}))
ggplot(data=agg[agg$Mean_speed[,2]>1,],aes(x=Mean_speed[,"Mean"],y=fit))+
  geom_abline()+
  labs(x="Group mean speed (m/s)",y="Model estimate + standard error",
       color="N (mean)",subtitle="Mean speed ~ Size + Direction")+
  geom_pointrange(aes(ymin=fit-se,ymax=fit+se,
                      color=Mean_speed[,"N"]),fatten=1)


# N_flights ####
hist(log(df1$N_flights))
hist(df1$N_flights)
shapiro.test(log(df1$N_flights))
aggregate(N_flights~Morph,data=df1,median)
range(df1$N_flights)
nrow(df1[df1$N_flights==1,])
plot(N_flights~Distance,data=df1)

# Flight initiation transformation (alternative to N_flights)
df1$Flight_time = pmin(df1$Flight_time,14390)
df1$Flight_initiation = df1$N_flights / (4-(df1$Flight_time/3600)) # Flights per hour

hist(df1$Flight_initiation,breaks=20)
hist(log(df1$Flight_initiation))
shapiro.test(log(df1$Flight_initiation)) # Just normal
plot(Flight_initiation~N_flights,data=df1)
plot(N_flights~Morph,data=df1)
plot(log(Flight_initiation)~Morph,data=df1)
aggregate(log(Flight_initiation)~Condition,data=df1,median)
# 1 = flights initiated immediately with no downtime
# 0 = whole experiment worth of downtime
# Full model lm
m10 = lm(log(Flight_initiation) ~ Morph + Sex + Direction + Condition + Size
          + Condition:Morph + Sex:Morph,data=df1)
summary(m10)
# Final model. D.O: Size-2.35, Direction-1.97, Morph:Sex-1.33, Sex-0.64
m10 = lm(log(Flight_initiation) ~ Morph + Condition + Condition:Morph,data=df1)
summary(m10)
drop1(m10)



# Pseudo R^2
with(summary(m10),1-deviance/null.deviance)

# Predict
newdat = with(df1,expand.grid(Condition=unique(Condition),Morph=unique(Morph)))
preds = predict(m10,newdata=newdat,se.fit=T,type="response")
newdat[,c("fit","se")] = preds[1:2]
# Normalise by Dir and Sex numbers if necessary
# newdat[,c("fit","se")] = newdat[,c("fit","se")]*
#   (table(df1$Direction)[newdat$Direction]/nrow(df1))*
#   (table(df1$Sex)[newdat$Sex]/nrow(df1))
# newdat = aggregate(newdat[c("fit","se")],newdat[c("Morph","Condition")],sum)
newdat[c("SEL","SEU")] = with(newdat,c(fit-se,fit+se))
# Exp if necessary
newdat[c("fit","SEL","SEU")] = exp(newdat[c("fit","SEL","SEU")])
# Take average of SE intervals for response SE
newdat$se = with(newdat,(SEU-SEL)/2)
# Transform to time to fly if wanted
# newdat[c("fit","SEU","SEL")] = 60/newdat[c("fit","SEL","SEU")]

# Compare model fit to aggregate
agg=merge(newdat,aggregate(log(df1["Flight_initiation"]),
                           df1[c("Condition","Morph")],function(x){
                             c(Mean=mean(x),N=length(x))}))
png(file="RGraphFlyInitiationMorphMod.png", res=900, width = 4800, height=3200,
    pointsize=14,type="windows", antialias="cleartype")#res700 plot 900 model
ggplot(data=agg,aes(x=Flight_initiation[,"Mean"],y=fit))+
  geom_abline()+
  labs(x="Group mean flight initiation /hour",y="Model estimate + standard error",
       color="N (mean)",caption="log link",
       subtitle="Flights per hour ~ Morph + Condition + Morph:Condition")+
  geom_pointrange(aes(ymin=fit-se,ymax=fit+se,
                      color=Flight_initiation[,"N"]),fatten=1)
dev.off()
table(df1$Sex,df1$Morph,df1$Condition)


# Females-only comparison
m10f = lm(log(Flight_initiation) ~ Morph + Condition,data=df1[df1$Sex=="Female",])
summary(m10f)


# 1.1: Plotting ####
# Output
# res700 plot 900 model, width 3200 legend on top
dev.off()
png(file="RGraphFlyInitiationMorphCond.png", res=700, width = 3200, height=3200,
    pointsize=14,type="windows", antialias="cleartype")


# Distance
ggplot(df1,aes(x=Condition,y=Distance,color=Morph,fill=Morph))+
  scale_y_log10()+
  geom_boxplot(notch=F,outlier.shape=NA,alpha=0.1)+
  geom_point(position=position_jitterdodge(),size=1)+
  scale_color_manual(values=cbPalette[c(6,7)],name="Morph:")+
  scale_fill_manual(values=cbPalette[c(6,7)],name="Morph:")+
  theme_classic()+
  theme(legend.position="top",legend.justification="right",text=element_text(size=18),
        panel.grid.minor=element_blank(),panel.grid.major.x=element_blank())+
  labs(y="Distance (m)",x="\nCondition")
ggplot(df1,aes(x=Condition,y=Distance,color=Direction,fill=Direction))+
  scale_y_log10()+
  geom_boxplot(notch=F,outlier.shape=NA,alpha=0.05)+
  geom_point(position=position_jitterdodge(),size=1)+
  scale_color_manual(values=cbPalette[c(1,4)],name="Direction:")+
  scale_fill_manual(values=cbPalette[c(1,4)],name="Direction:")+
  theme_classic()+
  theme(legend.position="top",legend.justification="right",text=element_text(size=18),
        panel.grid.minor=element_blank(),panel.grid.major.x=element_blank())+
  labs(y="Distance (m)",x="\nCondition")
palette.pals()

# Speed
table(df1$Sex,df1$Size)
ggplot(df1,aes(x=Size,y=Mean_speed,fill="Black"))+
  geom_boxplot(notch=F,outlier.shape=NA,alpha=0.05)+
  geom_point(position=position_jitterdodge(),size=1)+
  theme_classic()+
  theme(legend.position="None",text = element_text(size=15),
        panel.grid.minor=element_blank(),panel.grid.major.x=element_blank())+
  labs(y="Mean speed (m/s)",x="\nSize")
ggplot(df1,aes(x=Direction,y=Mean_speed))+
  geom_boxplot()

# Flights per hour
ggplot(data=df1,aes(x=Condition,y=60/Flight_initiation,color=Morph,fill=Morph))+
  scale_y_log10()+
  geom_boxplot(notch=F,outlier.shape=NA,alpha=0.1)+
  geom_point(position=position_jitterdodge(),size=1)+
  scale_color_manual(values=cbPalette[c(6,7)],name="Morph:")+
  scale_fill_manual(values=cbPalette[c(6,7)],name="Morph:")+
  theme_classic()+
  theme(legend.position="top",legend.justification="right",text=element_text(size=15),
        panel.grid.minor=element_blank(),panel.grid.major.x=element_blank())+
  labs(y="Flights initiated per hour of inactivity",x="\nCondition")
# Flights including model
ggplot(data=newdat,aes(x=Condition,y=fit,color=Morph,fill=Morph))+
  geom_boxplot(data=df1,aes(y=Flight_initiation),notch=F,outlier.shape=NA,alpha=0.1)+
  geom_pointrange(aes(ymin=SEL,ymax=SEU),position=position_dodge(width=0.2),
                  fatten=4,linewidth=1)+
  scale_y_log10()+
  scale_color_manual(values=cbPalette[c(6,7)],name="Morph:")+
  scale_fill_manual(values=cbPalette[c(6,7)],name="Morph:")+
  theme_classic()+
  theme(legend.position="top",legend.justification="right",text=element_text(size=15),
        panel.grid.minor=element_blank(),panel.grid.major.x=element_blank())+
  labs(y="Flights initiated per hour of inactivity",x="\nCondition")
  

#####
# 2: Flight-specific analysis ####

# Flight performance over time.
# Variables: flight performance, flight start time, fly, morph
sapply(df2,class)
colnames(df2)
df2$Time = df2$T_start/14400 # 14400 = experimental duration
df2$TimeLog = log(df2$Time)
df2$T_mid_norm = with(df2,(T_start+(0.5*Duration))/14400) # Sample at middle of flight
df2$timeOfDay = df2$Exp_start + df2$Time
df2$lateness = df2$timeOfDay - min(df2$timeOfDay)
df2$Speed_var = df2$Speed_std^2

sapply(df2,function(col){
  any(is.na(col))
})
#Descriptives
table(df2$Condition,df2$Morph)
hist(log(df2[df2$Time>0.25,"Distance"]))
hist(log(df2[df2$Morph=="Autumn","Distance"]))
range(df2$Distance)
qqnorm(log(df2$Distance))
qqline(log(df2$Distance),col="blue",lwd=2)
range(df2$T_start)
hist(df2$T_start)
summary(df2$T_start)
hist(log(df2$Behav_rep))
hist(df2$Behav_add)
aggregate(Distance~Morph,data=df2,mean)
aggregate(Distance~Morph,data=df2,median)
aggregate(Duration~Morph,data=df2,mean)
aggregate(Duration~Morph,data=df2,median)
hist(log(df2$Duration))
qqnorm(log(df2$Duration))
qqline(log(df2$Duration),col="blue",lwd=2)
plot(df2[df2$Duration>3600&df2$T_start<3600,"T_start"])
max(table(df2[df2$T_start<5,"Fly"]))
length(df2[df2$T_start<5,"Fly"])
df2[(df2$T_end>14399),"T_start"]
median(df2[df2$T_start>5,"Distance"])
median(df2[df2$T_start<5,"Distance"])
aggregate(Distance~Morph,data=df2[df2$T_start>5&df2$T_start<13800,],median)
aggregate(Distance~Morph,data=df2[df2$T_start<5,],median)
12754/60


# All flights full model. M:C, M:Sex, Cond:TimeLog, M:TimeLog, Size, Sex
mm2 = glmmTMB(Distance ~ Morph + Condition + Size + Sex + Direction + Time + TimeLog
              + Morph:Condition + Morph:Sex + Morph:Time + Morph:TimeLog
              + Condition:Time + Condition:TimeLog
              + (1|Fly),
              family=Gamma(link="log"),data=df2[df2$T_start<13800,])
summary(mm2)

# Final model. D.O.: M:C-2.2,M:TimeLog-2.0,M:Sex-1.8,Sex-1.5,Cond:TimeLog-0.74 USED
mm2 = glmmTMB(Distance ~ Morph + Condition + Size + Direction + Time + TimeLog
              + Morph:Time + Condition:Time
              + (1|Fly),
              family=Gamma(link="log"),data=df2[df2$T_start<13800,])
summary(mm2)
d=drop1(mm2)
View(d)

AICc(mm2)
AICc(update(mm2,.~.-Direction))
AIC(mm2)
AIC(update(mm2,.~.-Direction))

anova(mm2,glmmTMB(Distance ~ Morph + Time + (1|Fly),
                  family = Gamma(link = "log"),data=df2[df2$T_start>3600,]))

# Pseudo R^2
library(MuMIn)
r.squaredGLMM(mm2)

# Predict like magic glmm
# http://bbolker.github.io/mixedmodels-misc/glmmFAQ.html#glmmtmb
newdat = with(df2,expand.grid(Morph=unique(Morph),Condition=unique(Condition),
                              Direction=unique(Direction),Size=unique(Size),
                              Time=seq(0.0005,0.9583,length.out=200)))
newdat$TimeLog = log(newdat$Time)
mm = model.matrix(delete.response(terms(mm2)),newdat)
newdat$fit = drop(mm %*% fixef(mm2)[["cond"]])
newdat$pred = diag(mm %*% vcov(mm2)[["cond"]] %*% t(mm))
# Normalise variables by size weighted mean
newdat[,c("pred","fit")] = newdat[,c("pred","fit")]*
  (table(df2$Direction)[newdat$Direction]/nrow(df2))*
  (table(df2$Size)[newdat$Size]/nrow(df2))*
  (table(df2$Condition)[newdat$Condition]/nrow(df2))
newdat = aggregate(newdat[,c("pred","fit")],newdat[c("Morph","Time")],sum)
newdat$FE <- sqrt(newdat$pred)
newdat$TE <- sqrt(newdat$pred+sigma(mm2)^2)
newdat$T_start = newdat$Time*3600
newdat = with(newdat,data.frame(newdat,
                                TEL = fit-TE,TEU = fit+TE,
                                Distance = exp(fit),
                                RTEL = exp(fit-TE),RTEU = exp(fit+TE),
                                RFEL = exp(fit-FE),RFEU = exp(fit+FE)))
# Model validation
plot(mm2)
ranef(mm2)
# Normalised residuals
sresid = resid(mm2,type="pearson")
hist(log(sresid))
fits = fitted(mm2)
plot(sresid ~ fits)

# 2.1 Plotting ####

# From newdats earlier USED
ggplot(df2,aes(x=Time*4,y=Distance,color=Morph,fill=Morph))+
  scale_y_log10()+
  geom_point(shape=1,size=0.8,alpha=0.4)+
  geom_ribbon(data=newdat,linetype=0,alpha=0.3,
              aes(ymin=RFEL,ymax=RFEU),color=NA)+
  geom_line(data=newdat)+
  scale_color_manual(values=cbPalette[c(6,7)],name="Morph:")+
  scale_fill_manual(values=cbPalette[c(6,7)],name="Morph:")+
  theme_light()+
  theme(legend.position="top",legend.justification="right",text=element_text(size=12))+
  labs(y="Distance per flight (m)",x="\nFlight start time (hours)")


# OUTPUT
png(file="RGraphFlightsDistance3.png", res=700, width = 4800, height=3200,
    pointsize=14,type="windows", antialias="cleartype")
# 4800/4000 x 3200
dev.off()


#####
# 3: Experiment transect ####
df3 = df3[!df3$Index==0,]
df3$Time = (df3$Index+1)/240
df3$TimeLog = log(df3$Time)
range(df3$TimeLog)
range(df3$Time)

hist(df3$Speed,breaks=20)
hist(log(df3$Speed))
hist(sqrt(df3$Speed))
qqnorm(sqrt(df3$Speed))
qqline(sqrt(df3$Speed),col="blue",lwd=2)
qqnorm(log(df3$Speed))
qqline(log(df3$Speed),col="blue",lwd=2)
for (i in 1:240){
  print(paste(i,mean(df3[df3$Index==i,"Flying"])))
}
mean(mean(df3[df3$Index%in%235:239,"Flying"]))
# Descriptives
table(df1$Morph,df1$Size)
aggregate(Speed~Morph,data=df3,mean)
dft = aggregate(log(Speed)~Morph + Size,data=df3,mean)
dft$Response = exp(dft[3])

# Full model
mm3 = glmmTMB(Speed ~ Morph + Condition + Size + Sex + Direction + Time + TimeLog
              + Morph:Condition + Morph:Sex + Morph:Time + Morph:TimeLog
              + Condition:Time + Condition:TimeLog
              + (1|Fly),family=Gamma(link="log"),data=df3)
summary(mm3)
step(mm3)
sapply(df3,function(col){sum(is.na(col))})

# Final model. D.O.: Morph:Condition-3.9,Morph:Sex-1.3,Sex-1.6 USED
mm3 = glmmTMB(Speed ~ Morph + Condition + Size + Direction + Time + TimeLog
              + Morph:Time + Morph:TimeLog + Condition:Time + Condition:TimeLog
              + (1|Fly),family=Gamma(link="log"),data=df3)
summary(mm3)
drop1(mm3)

plot(resid(mm3,type="pearson")~df3[!is.na(df3$Events),"Time"])
AIC(mm3)

# applied magic preds WORKING
# FE fixed error, TE total error
newdat = with(df3,expand.grid(Morph=unique(Morph),Speed=0,Size=unique(Size),
                              Time=unique(Time),Condition=unique(Condition),
                              Direction=unique(Direction)))
newdat$TimeLog = log(newdat$Time)
mm = model.matrix(terms(mm3),newdat)
newdat$fit = drop(mm %*% fixef(mm3)[["cond"]])
newdat$predvar = diag(mm %*% tcrossprod(vcov(mm3)[["cond"]],mm))
# Normalise variables by size weighted mean
newdat[,c("predvar","fit")] = newdat[,c("predvar","fit")]*
  (table(df3$Size)[newdat$Size]/nrow(df3))*
  (table(df3$Direction)[newdat$Direction]/nrow(df3))*
  (table(df3$Condition)[newdat$Condition]/nrow(df3))
newdat = aggregate(newdat[,c("predvar","fit")],newdat[c("Morph","Time")],sum)
newdat$FE = sqrt(newdat$predvar)
newdat$TE = sqrt(newdat$predvar + VarCorr(mm3)$cond$Fly[1])
newdat = with(newdat,data.frame(newdat,
                                TEL = fit-TE,TEU = fit+TE,Speed = exp(fit),
                                RTEL = exp(fit-TE),RTEU = exp(fit+TE),
                                RFEL = exp(fit-FE),RFEU = exp(fit+FE)))
newdat2=newdat[newdat$Direction=="Clockwise"&newdat$Size=="Average",c(1,4,5,14)]
aggregate(Speed~Time,mean,data=newdat)
# R^2
library(MuMIn)
r.squaredGLMM(mm3)
# Model validation
plot(mm3)
library(DHARMa)
mm3simres = simulateResiduals(mm3)
plot(mm3simres)
car::Anova(mm3,type="III")

# 3.1 Plotting ####

# From newdat earlier
ggplot(df3,aes(x=Time*4,y=Speed,color=Morph,fill=Morph))+
  scale_y_log10()+
  geom_point(shape=1,alpha=0.2,size=0.8)+
  # geom_ribbon(data=newdat,aes(ymin=RTEL,ymax=RTEU),alpha=0.2)+
  #facet_grid(~Condition)+
  geom_ribbon(data=newdat,aes(ymin=RFEL,ymax=RFEU),alpha=0.3,color=NA)+
  geom_line(data=newdat,aes(y=Speed))+
  scale_color_manual(values=cbPalette[c(6,7)],name="Morph:")+
  scale_fill_manual(values=cbPalette[c(6,7)],name="Morph:")+
  theme_light()+
  theme(legend.position="top",legend.justification="right",text=element_text(size=12))+
  labs(y="Speed (m/s)",x="\nTime (hours)")

# OUTPUT
dev.off()
png(file="RGraphTransectSpeed4.png", res=700, width = 4000, height=3200,
    pointsize=14,type="windows", antialias="cleartype")
#facet grid= 6400x3200,normal=4800/4000x3200

# Flying

df3x = aggregate(Flying ~ Time + Morph + Condition,data=df3,mean)

# Facet grid for condition
ggplot(df3x,aes(x=Time*4,y=Flying,color=Morph))+
  geom_point(size=0.8,alpha=0.2)+
  theme_light()+
  facet_grid(~Condition)+
  geom_ribbon(data=newdat,alpha=0.2,linetype=0,
              aes(ymin=FEL,ymax=FEU,fill=Morph))+
  geom_line(data=newdat)+
  scale_color_manual(values=cbPalette[c(6,7)],name="Morph:")+
  scale_fill_manual(values=cbPalette[c(6,7)],name="Morph:")+
  theme(legend.position="top",legend.justification="right",text=element_text(size=15))+
  labs(x="\nTime (hours)",y="Proportion flying")
# Old theme(legend.position=c(0.95,0.85))

# 3.2 Modelling flying ####
table(df3$Morph,df3$Condition)/239
res = resid(mm31)
plot(res~df3$Time)
abline(0,0,col="blue")
m31$deviance/m31$df.residual
sum(residuals(m31,"pearson")^2)/m31$df.residual
library(DHARMa)
testDispersion(m31)
sim = simulateResiduals(m31)
plot(sim,asFactor=T)
plotResiduals(sim,df3$Time,quantreg=T)

# Dir-sex mixed model full model
mm31 = glmmTMB(Flying ~ Morph + Condition + Time + TimeLog
               + Morph:Condition + Morph:Sex + Morph:Time + Morph:TimeLog
               + Condition:Time + Condition:TimeLog
               + (1|Sex) + (1|Direction),
               family = binomial(link="logit"),data = df3)
# Final model. D.O.: Condition:Time+28,Morph:Time+3,Morph:TimeLog+1 USED
mm31 = glmmTMB(Flying ~ Morph + Condition + Time + TimeLog
               + Morph:Condition
               + Condition:TimeLog + Condition:Time
               + (1|Sex) + (1|Direction),
               family = binomial(link="logit"),data = df3)
summary(mm31)
drop1(mm31)


# Test for overdispersion
sum(residuals(mm31,type="deviance")^2)/33685#mm31$df.residual
# Pseudo R^2
library(MuMIn)
r.squaredGLMM(mm31) # need to change df3 from 3 char ending in number (eg to dft)

aggregate(Flying~Time,data=df3,mean)
# predict glmmTMB
library(boot)
newdat = with(df3,expand.grid(Morph=unique(Morph),Condition=unique(Condition),
                              Size=unique(Size),Sex=unique(Sex),
                              Time=unique(Time),Direction=unique(Direction)))
newdat$TimeLog = log(newdat$Time)
#newdat$TimeFun = (1-abs(newdat$Time-0.15))^5
mm = model.matrix(delete.response(terms(mm31)),newdat)
newdat$fit = drop(mm %*% fixef(mm31)[["cond"]])
newdat$pred = diag(mm %*% vcov(mm31)[["cond"]] %*% t(mm))
newdat$FE <- sqrt(newdat$pred)
newdat$TE <- sqrt(newdat$pred+sigma(mm31)^2)
# Weight by Direction
newdat[,c("fit","FE","TE")] = newdat[,c("fit","FE","TE")]*
  (table(df3$Direction)[newdat$Direction]/nrow(df3))*
  (table(df3$Size)[newdat$Size]/nrow(df3))*
  (table(df3$Sex)[newdat$Sex]/nrow(df3))
newdat = aggregate(newdat[c("fit","FE","TE")],newdat[c("Morph","Condition","Time")],sum)
newdat = with(newdat,data.frame(newdat,
                                Flying = inv.logit(fit),
                                CIL = inv.logit(fit-FE*1.96),
                                CIU = inv.logit(fit+FE*1.96),
                                FEL = inv.logit(fit-FE),
                                FEU = inv.logit(fit+FE),
                                TEL = inv.logit(fit-TE),
                                TEU = inv.logit(fit+TE)))

#####
# END OF INCLUDED ANALYSIS
#####

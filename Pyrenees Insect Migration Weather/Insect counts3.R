# Part 0: transform and explore data ####

# Install packages if not already installed
install.packages(setdiff(c("ggplot2","DHARMa","MASS","pscl","lmtest","Hmisc","mvtnorm"),
                         rownames(installed.packages())))
library(ggplot2) # ggplot
library(DHARMa) # simulateResiduals, testZeroInflation
library(pscl) # zeroinfl
library(lmtest) # lrtest
library(marginaleffects) # plot_predictions
library(MuMIn) # dredge, get.models

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
df1 = na.omit(read.csv("Pyrenees Modelling.csv",fileEncoding="UTF-8-BOM"))
names(df1) = substr(names(df1),1,nchar(names(df1))-5)
# Normalise and transform variables for analysis

df1$windLog = log(df1$windspeed)
df1$countsLog = log(df1$counts+1) # log counts starts from 0
df1$rainLog = log(df1$rain+1)
df1$sunNorm = df1$Sun/max(df1$Sun)
df1$Headwind = -cos(df1$Windheading/180*pi)
df1$Headwind_pass = -cos(((df1$Windheading-60)%%360)/180*pi)

#####
# Part 1 Descriptive statistics and descriptive plots ####
sapply(df1,class)
hist(df1$Windheading,breaks=12)
hist(df1$Temp)
hist(df1$Sun)
hist(df1$windspeed)
hist(log(df1$windspeed))
hist(df1$Headwind)
hist(df1$Headwind_pass)
shapiro.test(log(df1$windspeed))
qqnorm(log(df1$windspeed))
hist(df1[df1$rain>0,"rain"])
table(df1$rain>0)
hist(df1$counts)
hist(df1$countsLog)
table(df1$counts>0) # Counts N: 95 = 0, 96 >0

# Single plots
plot(countsLog~sunNorm,data=df1)
plot(countsLog~Temp,data=df1)
plot(countsLog~rain,data=df1)
plot(countsLog~rain,data=df1[df1$counts>0&df1$rain>0,])
plot(countsLog~rainLog,data=df1)
plot(countsLog~rainLog,data=df1[df1$counts>0&df1$rain>0,])
plot(countsLog~windspeed,data=df1)
plot(countsLog~windLog,data=df1)
plot(countsLog~Windheading,data=df1)
plot(countsLog~Headwind,data=df1)
plot(countsLog~Headwind_pass,data=df1)

# Testing linear correlation of variables
cor.test(~ countsLog + Headwind,data=df1[df1$counts>1,])
cor.test(~ countsLog + Headwind_pass,data=df1[df1$counts>1,])# higher correlation - use
cor.test(~ countsLog + sunNorm,data=df1[df1$counts>1,])
cor.test(~ countsLog + rain,data=df1[df1$counts>1&df1$rain>0,])
# Testing the effect of rainfall on days with migration
exp(mean(log(df1[df1$rain==0&df1$counts>0,"counts"]),))
exp(mean(log(df1[df1$rain>0&df1$counts>0,"counts"]),))
cor.test(~ countsLog + windspeed,data=df1[df1$counts>1,])
cor.test(~ countsLog + windLog,data=df1[df1$counts>1,])

# Exploratory multi-variable plots
ggplot(data=df1,aes(x=Headwind_pass,Temp,color=countsLog))+
  geom_point()
ggplot(data=df1,aes(x=windLog,y=sunNorm,color=countsLog))+
  geom_point()
ggplot(data=df1,aes(x=Temp,y=countsLog,color=sunNorm))+
  geom_point()
ggplot(data=df1,aes(x=Headwind_pass,y=Temp,color=countsLog))+
  geom_point()
ggplot(data=df1,aes(x=sunNorm,y=countsLog,color=Temp))+
  geom_point()

# Predictor variable collinearity (sunNorm, Temp, rain, Headwind_pass)
Hmisc::rcorr(as.matrix(df1[,c("sunNorm","Temp","rain","Headwind_pass")]),type="pearson")

ggplot(data=df1,aes(y=Headwind_pass,x=Temp,color=sunNorm))+
  geom_point()
ggplot(data=df1,aes(y=sunNorm,x=Temp,color=rain))+
  geom_point()

# Common sense:
# Sun correlates with temperature (0.45)
# Rain correlates negatively with sun (-0.36) + temperature (-0.22)
# Other:
# Headwind_pass correlates with temperature (0.2)

# None of these correlations are particularly high (sun vs temp @ 0.45 is the highest)
# The models won't have trouble discerning between the effects of different variables


#####
# Part 2 model family selection ####

m1 = list() # Create a list to store different models

# Poisson
m1$p = glm(counts ~ sunNorm + Temp + rain + Headwind,
           family=poisson(link=log),data=df1)
summary(m1$p)
# Validation: test for overdispersion
m1$p$deviance/m1$p$df.residual # Dispersion = 304785! - this is above 2 so neg bin needed

# Negative binomial
m1$nb = MASS::glm.nb(counts ~ sunNorm + Temp + Headwind, data=df1)
summary(m1$nb)
# Validation: test for zero inflation
res = simulateResiduals(m1$nb)
testZeroInflation(res) # Sig = data zero inflated

# Zero inflated poisson
# formula: dependent variable ~ vars affecting count | vars affecting presence/ absence
m1$zip = zeroinfl(counts ~ sunNorm + Temp + Headwind|sunNorm + Temp + Headwind,
                  dist="poisson",link="logit",data=df1)
summary(m1$zip)
# Validation: test for overdispersion
res = residuals(m1$zip,type="pearson")
sum(res^2) / (nrow(df1) - length(coef(m1$zip))) # 11.6 = still overdispersed but much less

# Zero-inflated negative binomial
m1$zinb = zeroinfl(counts ~ sunNorm+Temp+Headwind|sunNorm+Temp+Headwind,
                   dist="negbin",link="logit",data=df1)
summary(m1$zinb)
# Validation: test for overdispersion
res = residuals(m1$zinb,type="pearson")
sum(res^2) / (nrow(df1) - length(coef(m1$zinb))) # 1.21 = no longer overdispersed

# Comparison between zip and zinb
lrtest(m1$zip,m1$zinb) # zinb significantly better than zip

# Validation: examining residuals
preds = predict(m1$zinb)
hist(res)
plot(res~preds)
plot(res~df1$sunNorm)
plot(res~df1$Temp)
plot(res~df1$Headwind)

#####
# Part 3: fitting the best model possible. Variable selection and optimisation ####

# 3.1 Compare link functions (only for zero infl model - count model is always log)
for (link in c("logit","probit","cloglog","cauchit","log")){
  print(link)
  tryCatch({
    mn = zeroinfl(counts ~ sunNorm+Temp+rain+Headwind_pass+windspeed
                  +Headwind_pass:windspeed+sunNorm:Temp+sunNorm:rain
                  |sunNorm+Temp+rain+Headwind_pass+windspeed
                  +Headwind_pass:windspeed+sunNorm:Temp+sunNorm:rain,
                  dist="negbin",link=link,data=df1)
    print(paste("AIC:",AIC(mn),"loglik:",logLik(mn)))
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}
# Cauchit is best: AIC 2582.98 < 2586.05 of logit

# 3.2 Which headwind variable performs best?
# Headwind: wind vector that goes south (1 = perfect southerly, -1 = northerly)
# Headwind_pass: wind vector relative to the pass (1 = headwind, -1 = tailwind)

mo = zeroinfl(counts ~ sunNorm+Temp+rain+Headwind+windLog+Headwind:windLog
              |sunNorm+Temp+rain+Headwind+windLog+Headwind:windLog,
              dist="negbin",link="cauchit",data=df1)
mn = zeroinfl(counts ~ sunNorm+Temp+rain+Headwind_pass+windLog+Headwind_pass:windLog
              |sunNorm+Temp+rain+Headwind_pass+windLog+Headwind_pass:windLog,
              dist="negbin",link="cauchit",data=df1)
lrtest(mo,mn) # Headwind pass is better than headwind (higher LogLik)

# 3.3 Which wind variable (windspeed or windLog) is a better predictor?

mo = zeroinfl(counts ~ sunNorm+Temp+rain+Headwind_pass+windLog+Headwind_pass:windLog
              |sunNorm+Temp+rain+Headwind_pass+windLog+Headwind_pass:windLog,
              dist="negbin",link="cauchit",data=df1)
mn = zeroinfl(counts ~ sunNorm+Temp+rain+Headwind_pass+windspeed+Headwind_pass:windspeed
              |sunNorm+Temp+rain+Headwind_pass+windspeed+Headwind_pass:windspeed,
              dist="negbin",link="cauchit",data=df1)
lrtest(mo,mn) # windspeed marginally better

# 3.4 Which rain variable (rain or rainLog) is a better predictor?

mo = zeroinfl(counts ~ sunNorm+Temp+rain+Headwind_pass+windspeed+Headwind_pass:windspeed
              |sunNorm+Temp+rain+Headwind_pass+windspeed+Headwind_pass:windspeed,
              dist="negbin",link="cauchit",data=df1)
mn = zeroinfl(counts ~ sunNorm+Temp+rainLog+Headwind_pass+windspeed+Headwind_pass:windspeed
              |sunNorm+Temp+rainLog+Headwind_pass+windspeed+Headwind_pass:windspeed,
              dist="negbin",link="cauchit",data=df1)
lrtest(mo,mn) # rain is better


# 3.5 Next step: step-wise variable selection by AIC

# Full model
mo = zeroinfl(counts ~ sunNorm+Temp+rain+Headwind_pass+windspeed
              +Headwind_pass:windspeed+sunNorm:Temp+sunNorm:rain
              |sunNorm+Temp+rain+Headwind_pass+windspeed
              +Headwind_pass:windspeed+sunNorm:Temp+sunNorm:rain,
              dist="negbin",link="logit",data=df1)
summary(mo)
AIC(mo)

# Backwards step-wise variable selection (removed in favour of all subsets approach)

# Start with full model above (copied below)
# summary(model) to see which variable performs worst (highest p), delete from mn and rerun
# If model is improved (lower AIC), delete from old model (mo) too and repeat
# If model worsens (higher AIC), put back in new model (mn) and try a different variable
# Note: treat each variable separately for count model and zero-inflation model
# When the model is no longer improved, move on to validation

# # 1st effort. Starting AIC = 2584.07
# # D.O.: HWP:windspeed(b)-1.94,sun(b)-1.73,windspeed(b)-0.90,HWP:windspeed(c)-0.74
# # windspeed(c)-1.20
# mo = zeroinfl(counts ~ sunNorm+Temp+rain+Headwind_pass
#               |Temp+rain+Headwind_pass,
#               dist="negbin",link="logit",data=df1)
# mn = zeroinfl(counts ~ sunNorm+Temp + Headwind_pass
#               |Temp+rain+Headwind_pass,
#               dist="negbin",link="logit",data=df1)


# Starting with extra interactions - sun:rain and sun:temp
# D.O.: HWP:windspeed(b)-1.96,sun:rain(c)-1.91,sun:temp(b)-1.78, HWD:windspeed(c)-1.63
# windspeed(c)-1.27,sun:rain(b)-0.89,sun(b)-1.74,windspeed(b)-0.95
# mo = zeroinfl(counts ~ sunNorm+Temp+rain+Headwind_pass
#               +sunNorm:Temp
#               |Temp+rain+Headwind_pass+windspeed,
#               dist="negbin",link="cauchit",data=df1)
# mn = zeroinfl(counts ~ sunNorm+Temp+rain+Headwind_pass
#               +sunNorm:Temp
#               |Temp+rain+Headwind_pass,
#               dist="negbin",link="cauchit",data=df1)
# AIC(mo,mn)
# summary(mn)
# AIC(mn)-AIC(mo) # AIC change

# Full model subsets approach
mo = zeroinfl(counts ~ sunNorm+Temp+rain+Headwind_pass+windspeed
              +Headwind_pass:windspeed+sunNorm:Temp+sunNorm:rain
              |sunNorm+Temp+rain+Headwind_pass+windspeed
              +Headwind_pass:windspeed+sunNorm:Temp+sunNorm:rain,
              dist="negbin",link="cauchit",data=df1,na.action="na.fail")
model.set = dredge(mo) # Takes time
top.models = get.models(model.set,subset=delta <2)
sapply(top.models,"[[","call")
sapply(top.models,AIC)


# Final model (top model of top.models - same model as when running stepwise)
m1 = zeroinfl(counts ~ sunNorm+Temp+rain+Headwind_pass
              +sunNorm:Temp
              |Temp+rain+Headwind_pass,
              dist="negbin",link="cauchit",data=df1)
summary(m1)
confint(m1)
AIC(m1)

#####
# Part 4: model validation ####

# Validation: test for overdispersion
res = residuals(m1,type="pearson")
sum(res^2) / (nrow(df1) - length(coef(m1))) # dispersion 0.892 - fine

# Examining residuals
hist(res)
preds = predict(m1)
plot(res~preds) # Variation to be expected but at least the effect is even
plot(residuals(m1,type="response")~preds)
plot(res~df1$sunNorm)
plot(res~df1$Temp) # May correlate with magnitude
plot(res~df1$rain) # 0 inflated so probably okay
plot(res~df1$Headwind_pass)


#####
# Part 5: making predictions ####
# Variables to predict = sunNorm, rain, Temp, Headwind_pass
range(df1$Headwind_pass)
range(df1$Temp)
range(df1$rain)
range(df1$sunNorm)

# Dummy data to predict from (all variables apart from dependent)
# Data to predict all possible variations (although n=30 of each parameter)
newdat = with(df1,expand.grid(sunNorm=seq(0,1,length.out=30),
                              rain=seq(0,max(rain),length.out=30),
                              Headwind_pass=seq(-1,1,length.out=30),
                              Temp=seq(min(Temp),max(Temp),length.out=30)))

# Make the predictions
preds = predict(m1,newdata=newdat,type="response",se=T)
head(preds)
# No standard error :(
# Incorporate into DF anyway
newdat[,c("Counts")] = preds

# Alternative method: bootstrap
# predict.infl bootstrap?
source("predict.zeroinfl.R")

# Predict main data
preds = predictSE_short.zeroinfl(m1,newdata=newdat,type="response",se=T)
saveRDS(preds,"predict.zeroinflPreds") # Save predictions to rerun later
preds = readRDS("predict.zeroinflPreds") # Reload predictions from earlier
newdat = cbind(newdat,preds)
colnames(newdat)[ncol(newdat)-3] = "counts"

# Optional extra: predict target variable combinations
newdat2 = with(df1,expand.grid(sunNorm=seq(0,1,length.out=3),rain=0,Headwind_pass=0,
                               Temp=seq(5,15,length.out=3)))
preds2 = predictSE_short.zeroinfl(m1,newdata=newdat2,type="response",se=T)
tempdat = cbind(newdat2,preds2)
colnames(newdat2)[ncol(newdat2)-3] = "counts"
write.csv(newdat2,"Newdat output.csv")

#####
# Part 6: plotting ####
sapply(df1,median)
colnames(df1)

single_plot = function(Independent,ndat=newdat,df=df1,fill="grey",fly_quantile=F){
  # Find position of median values of other predictors- TRUE when median
  vectors=sapply(setdiff(c("sunNorm","Temp","rain","Headwind_pass"),Independent),
                 function(col){
                   if(fly_quantile){
                     # This outrageous line finds the value that contains the median fly
                     medFly=max(which(cumsum(df[order(df[[col]]),"counts"])
                                      <=sum(df$counts)*fly_quantile))
                     median_fly=df[order(df[[col]])[medFly],col]
                     ndat[,col]==ndat[which.min(abs(median_fly-ndat[,col])),col]
                   } else {
                     # Previous functionality using medians of weather data
                     ndat[,col]==ndat[which.min(abs(median(df[,col])-ndat[,col])),col]
                   }
                 }
  )
  # subset newdat to when all columns are medians
  ndat = ndat[rowSums(vectors)==ncol(vectors),]
  # Hack to change the x variable to "ind" so that it works with ggplot
  colnames(ndat) = sub(Independent,"ind",colnames(ndat))
  colnames(df) = sub(Independent,"ind",colnames(df))
  xlab = list(sunNorm="Daily sunshine proportion",
              Headwind_pass="Daily relative wind direction (1 = headwind)",
              Temp="Average daily temperature (°C)",
              rain="Total daily precipitation (mm)")[[Independent]]
  plt = ggplot(df,aes(x=ind,y=counts+1))+
    scale_y_log10(label=scales::comma)+
    geom_point(alpha=0.3)+
    geom_ribbon(data=ndat,aes(y=counts+1,ymin=lower+1,ymax=upper+1),
                alpha=0.2,fill=fill)+
    geom_line(data=ndat,aes(y=counts+1),linewidth=0.5)+
    theme_classic()+
    theme(text=element_text(size=14))+
    labs(x=xlab,y="Number of insects")
  if(Independent=="rain") plt = plt + scale_x_log10() # Use for only rain
    plt
}
single_plot("Temp",fly_quantile=F)
single_plot("sunNorm",fly_quantile=F)
single_plot("rain")
single_plot("Headwind_pass")


# Alternative wrapper method for plotting higher resolution: only medians
bootstrap_plot = function(Independent,df=df1,fly_quantile=F,fill="grey"){
  if(fly_quantile){
    vars = sapply(c("rain","sunNorm","Temp","Headwind_pass"),function(col){
      medFly=max(which(cumsum(df[order(df[[col]]),"counts"])
                       <=sum(df$counts)*fly_quantile))
      Median_fly = df[order(df[[col]])[medFly],col]
      return(Median_fly)
    },USE.NAMES=T,simplify=F)
  } else {
    vars = with(df,list(rain=median(rain),sunNorm=median(sunNorm),
                        Headwind_pass=median(Headwind_pass),Temp=median(Temp)))
  }
  vars[[Independent]] = seq(min(df[Independent]),max(df[Independent]),length.out=3000)
  newdat = expand.grid(get("vars"))
  if(!exists("predictSE_short.zeroinfl")) source("predict.zeroinfl.R")
  preds = predictSE_short.zeroinfl(m1,newdata=newdat,type="response",se=T)
  newdat = cbind(newdat,preds)
  colnames(newdat)[ncol(newdat)-3] = "counts"
  single_plot(Independent,ndat=newdat,fill=fill)
}
bootstrap_plot("Temp")
bootstrap_plot("rain")
bootstrap_plot("Headwind_pass")
bootstrap_plot("sunNorm")


# Export plots

export_single_plot = function(Independent="sunNorm",bootstrap=F,fly_quantile=F){
  fill=list(sunNorm="yellow",Temp="red",rain="blue",Headwind_pass="grey")[[Independent]]
  # This is how I output. Might need to tweak res.
  if(fly_quantile==F) ID = "medians" else ID = fly_quantile
  png(file=paste0("PyreneesCounts_",Independent,"_",ID,".png"), res=700,width = 4800,
      height=3200, pointsize=14,type="windows", antialias="cleartype")
  if(bootstrap==T){
    print({bootstrap_plot(Independent=Independent,fill=fill,fly_quantile=fly_quantile)})
  } else {
    print({single_plot(Independent = Independent,fill=fill)})
  }
  # Turn off device to stop picture writing and finish
  dev.off()
}

# Make and save all plots
for (independent in c("sunNorm","rain","Temp","Headwind_pass")){
  export_single_plot(independent,bootstrap=T,fly_quantile=0.1)
}

# marginaleffects package to plot zero effects
summary(m1)
plot_predictions(m1, condition = "Temp",type="zero",vcov=T)
plot_predictions(m1, condition = "sunNorm",type="zero",vcov=T)
plot_predictions(m1, condition = "rain",type="zero",vcov=T)
plot_predictions(m1, condition = "Headwind_pass",type="zero",vcov=T)

# Make and save all zero plots
for (independent in c("sunNorm","rain","Temp","Headwind_pass")){
  png(file=paste0("PyreneesZero_",independent,".png"),res=900,width=4800,height=3200,
      pointsize=14,type="windows",antialias="cleartype")
  print(plot_predictions(m1, condition=independent, type="zero",vcov=T)+
          geom_rug(data=df1, aes(x=!!sym(independent)), inherit.aes = FALSE))
  dev.off()
}

# Part 7 #### creating models of various combinations to identify relative AIC weights

# Wrapped in function
model_run = function(tier=1){
  vars=list("sunNorm","windspeed","Temp","rain","Headwind_pass")
  if (tier==1) j = NULL
  non_interacts = list()
  for (interact in list(c(1,3),c(1,4),c(2,5))){
    if (all(interact%in%(1:5)[-c(i,j)])){
      vars = append(vars,paste0(vars[interact[1]],":",vars[interact[2]]))
    } else non_interacts = append(non_interacts,interact[2]+3)
  }
  formula1 = as.formula(paste("counts~",paste(vars[-c(i,j)],collapse="+"),
                              "|",paste(vars[-c(i,j)],collapse="+")))
  m1 = zeroinfl(formula1,dist="negbin",link="cauchit",data=df1)
  r = rep(AIC(m1),8)
  r[c(i,j,unlist(non_interacts))] = NA
  print(paste(paste(vars[-c(i,j)],collapse="+"),AIC(m1)))
  return(r)
}
model_run = function(tier=1,part="both"){
  if (tier==1) j = NULL
  non_interacts = list()
  for (interact in list(c(1,3),c(1,4),c(2,5))){
    if (all(interact%in%(1:5)[-c(i,j)])){
      vars = append(vars,paste0(vars[interact[1]],":",vars[interact[2]]))
    } else non_interacts = append(non_interacts,interact[2]+3)
  }
  countvars = paste(vars[-c(i,j)],collapse="+")
  zerovars = paste(vars[-c(i,j)],collapse="+")
  if(part=="count") zerovars="Temp+rain+Headwind_pass"
  if(part=="zero") countvars="sunNorm+Temp+rain+Headwind_pass+sunNorm:Temp"
  formula1 = as.formula(paste("counts~",countvars,"|",zerovars))
  m1 = zeroinfl(formula1,dist="negbin",link="cauchit",data=df1)
  r = rep(AIC(m1),8)
  r[c(i,j,unlist(non_interacts))] = NA
  print(paste(paste(vars[-c(i,j)],collapse="+"),AIC(m1)))
  return(r)
}
res=list()
p="zero"
for (i in 1:5){
  res[[paste((1:5)[-i],collapse="")]] = model_run(part=p)
  for (j in (i:5)[-1]){
    res[[as.character(paste((1:5)[-c(i,j)],collapse=""))]]=model_run(tier=2,part=p)
  }
}
df2 = do.call(rbind,res)
colnames(df2) = c(vars,"sunNorm:Temp","sunNorm:rain","windspeed:Headwind_pass")
res_mean_N = cbind(colMeans(df2,na.rm=T),colSums(!is.na(df2)))
res_mean_N = res_mean_N[order(res_mean_N[,1]),]
res_mean_N


write.csv(df2,"Round robin models.csv")
write.csv(res_mean_N,"Round robin analysis parameter results.csv")
# END OF SCRIPT ####


# Workings

sum(df1$counts)
median(df1$counts)
mean(df1$counts)
exp(mean(df1$countsLog))

# Find better than median - find median fly!
sapply(df1[,c("rain","sunNorm","Temp","Headwind_pass")],median)
sapply(c("rain","sunNorm","Temp","Headwind_pass"),function(col){
  medFly=max(which(cumsum(df1[order(df1[[col]]),"counts"])<=sum(df1$counts)*0.1))
  Median_fly = df1[order(df1[[col]])[medFly],col]
  return(Median_fly)
})#,USE.NAMES=T,simplify=F)
# Medians: Temp = 9.695°, rain = 0mm, Headwind_pass = 0.7031, sunNorm = 0.3924
# 0.5 fly: Temp = 12.039°, rain = 0mm, Headwind_pass = 0.9543, sunNorm = 0.8522
# 0.2 fly: Temp = 11.21°, rain = 0mm, Headwind_pass = 0.838, sunNorm = 0.548
# 0.1 fly: Temp = 9.821°, rain = 0mm, Headwind_pass = 0.826, sunNorm = 0.376

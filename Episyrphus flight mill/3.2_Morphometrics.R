library(ggplot2)
setwd(paste0("~/Data S2"))
read_sheets = function(file){
  # Base functionality
  require(readxl)
  sheets = excel_sheets(file)
  dfs = setNames(lapply(sheets,function(sheet){
    as.data.frame(read_excel(file,sheet=sheet))
  }),sheets)
  # Specific processing
  lapply(dfs,function(df){
    df = df[df$Species=="Ebal",]
    df[,c("Sex","Condition","Size")] = do.call(rbind,strsplit(df$Criteria,NULL))
    df$Condition = ordered(df$Condition,levels = c("t","m","f"))
    df$Size = ordered(df$Size,levels = c("s","a","l"))
    levels(df$Size) = c("Small","Average","Large")
    levels(df$Condition) = c("Thin","Medium","Fat")
    return(df)
  })
}
dfs = read_sheets("Morphometric information.xlsx")

#####
# Stats ####
# Descriptives
table(dfs$Size$Criteria)
table(dfs$Weights$Criteria)

table(dfs$Weights$Condition)
table(dfs$Weights$Size)
table(dfs$Weights$Size,dfs$Weights$Condition)
table(dfs$Weights$Sex)
# Weights ####
qqnorm(dfs$Weights$Weight_g)
qqline(dfs$Weights$Weight_g,col="blue",lwd=2)
qqnorm(log(dfs$Weights$Weight_g))
qqline(log(dfs$Weights$Weight_g),col="blue",lwd=2)
shapiro.test(dfs$Weights$Weight_g)
shapiro.test(log(dfs$Weights$Weight_g))
# Full model
m1 = lm(Weight_g ~ Size * Condition + Sex, data = dfs$Weights)
summary(m1)
# Final model. DO: Size:Condition, Sex
m1 = lm(Weight_g ~ Size + Condition, data = dfs$Weights)
drop1(m1)
summary(m1)

# Size ####
diff(aggregate(Length_mm~Size,data=dfs$Size,mean)[[2]])
qqnorm(log(dfs$Size$Length_mm))
qqline(log(dfs$Size$Length_mm),col="blue",lwd=2)
qqnorm(log(dfs$Size$Area_mm))
qqline(log(dfs$Size$Area_mm),col="blue",lwd=2)
shapiro.test(dfs$Size$Length_mm)
# Full model
m2 = lm(Length_mm ~ Size * Condition + Sex, data = dfs$Size)
summary(m2)
# Final model. DO: Size:Condition, Sex, Condition
m2 = lm(Length_mm ~ Size, data = dfs$Size)
drop1(m2)
summary(m2)
#####
# Plots ####
ggplot(dfs$Weights,aes(x=Size,y=Weight_g))+
  geom_boxplot(fill="black",alpha=0.05)+
  theme_classic()+
  theme(text=element_text(size=24),axis.text=element_text(size=20))+
  labs(y="Weight (grams)",x="\nSize")
ggplot(dfs$Weights,aes(x=Condition,y=Weight_g))+
  geom_boxplot(fill="black",alpha=0.05)+
  theme_classic()+
  theme(text=element_text(size=24),axis.text=element_text(size=20))+
  labs(y="Weight (grams)",x="\nCondition")
ggplot(dfs$Size,aes(x=Size,y=Length_mm))+
  geom_boxplot(fill="black",alpha=0.05)+
  theme_classic()+
  theme(text=element_text(size=24),axis.text=element_text(size=20))+
  labs(y="Wing length (mm)",x="\nSize")
ggplot(dfs$Size,aes(x=Condition,y=Length_mm))+
  geom_boxplot(fill="black",alpha=0.05)+
  theme_classic()+
  theme(text=element_text(size=24),axis.text=element_text(size=20))+
  labs(y="Wing length (mm)",x="\nCondition")

# Export
dev.off()
png(file="RGraphMMLengthSize.png", res=700, width = 4800, height=3200, pointsize=14,
    type="windows", antialias="cleartype")


#####
# Predict ####
Predict = function(model,predDat,type="response",re.form=NA){
  preds = predict(model,newdata=predDat,type=type,se.fit=T,re.form=re.form)
  predFrame = cbind(predDat,fit=preds$fit,se=preds$se.fit)
  predFrame$upperSE = predFrame$fit+predFrame$se
  predFrame$lowerSE = predFrame$fit-predFrame$se
  return(predFrame)
}
P=list()
P$W = Predict(m1,with(dfs$Weights,expand.grid(Size=unique(Size),
                                              Condition=unique(Condition))))
P$W = merge(P$W,aggregate(Weight_g~Size+Condition,data=dfs$Weights,mean),all=T)
P$W = merge(P$W,aggregate(Criteria~Size+Condition,data=dfs$Weights,length),all=T)
P$S = Predict(m2,with(dfs$Size,expand.grid(Size=unique(Size))))
P$S = merge(P$S,aggregate(Length_mm~Size,data=dfs$Size,mean),all=T)
P$S = merge(P$S,aggregate(Criteria~Size,data=dfs$Size,length),all=T)
View(P$W)
View(P$S)



ggplot(data=P$S,aes(x=Length_mm,y=fit))+
  geom_abline()+
  labs(x="Group mean wing length (mm)",y="Model estimate + standard error",
       color="N (mean)",subtitle="Wing_length_mm ~ Size")+
  geom_pointrange(aes(ymin=lowerSE,ymax=upperSE,
                      color=Criteria),fatten=1)
ggplot(data=P$W,aes(x=Weight_g,y=fit))+
  geom_abline()+
  labs(x="Group mean weight (g)",y="Model estimate + standard error",
       color="N (mean)",subtitle="Weight_g ~ Condition + Size")+
  geom_pointrange(aes(ymin=lowerSE,ymax=upperSE,
                      color=Criteria),fatten=1)


png(file="RGraphMMWeightMod.png", res=900, width = 4800, height=3200, pointsize=14,
    type="windows", antialias="cleartype")
dev.off()
plot(fit~Weight_g,data=P$W)
abline(0,1)
plot(fit~Length_mm,data=P$S)
abline(0,1)
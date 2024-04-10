Predict = function(model,predDat,type="response",re.form=NA){
  preds = predict(model,newdata=predDat,type=type,se.fit=T,re.form=re.form)
  predFrame = cbind(predDat,fit=preds$fit,se=preds$se.fit)
  predFrame$upperSE = predFrame$fit+predFrame$se
  predFrame$lowerSE = predFrame$fit-predFrame$se
  return(predFrame)
}

Lines = function(predFrame,predVar="fit",col="black",lwd=2){
  lines(as.formula(paste("fit ~",predVar)),data=predFrame,col=col,lwd=lwd)
  lines(as.formula(paste("upperSE ~",predVar)),data=predFrame,col=col,lty=2,lwd=lwd)
  lines(as.formula(paste("lowerSE ~",predVar)),data=predFrame,col=col,lty=2,lwd=lwd)
}

# See below for GLM
# http://bbolker.github.io/mixedmodels-misc/glmmFAQ.html#glmmtmb
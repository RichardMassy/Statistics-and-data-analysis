setwd(paste0("~/Data S2"))

studies = c("21 Autumn","22 Summer")

inputs = list()

for (study in studies){
  folder = paste0(study,"/")
  for (file in c("All flights","Flight_stats","Errors","Transects")){
    inputs[[file]][[study]] = read.csv(paste0(folder,file,".csv"),check.names=F)
  }
}

flies = do.call(rbind,inputs[["Flight_stats"]])[,c("Fly","dom_dir")]

inputs[["Transects"]] = lapply(inputs[["Transects"]],function(df){
  df = df[,names(df)%in% c("Index",flies$Fly)]
  df = reshape(df,varying=list(names(df)[-1]),direction="long",v.names="Events",
               timevar="Fly",times=names(df)[-1])[-4]
  df = merge(df,flies)
  return (df)
})

info = read.csv("Fly information.csv")

for (input in names(inputs)){
  df1 = do.call(rbind,Map(cbind,inputs[[input]],"Origin"=names(inputs[[input]])))
  if (input != "Errors"){
    df1 = merge(df1,info,by="Fly")
  }
  write.csv(df1,paste0("RD ",input,".csv"),row.names=F)
}
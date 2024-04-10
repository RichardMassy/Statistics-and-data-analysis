angle.diff = function(angle1,angle2 = NULL){
  if (is.null(angle2)){
    angle2 = group.dir
  }
  min((angle1-angle2)%%360,(angle2-angle1)%%360)
}

Rectification <- function(x){
  if (!is.na(x)){
    if (!exists("experiment")){x <- (x - angle_diff) %% 360}
    else if (experiment %in% c("D","S","U","T","R","RS","RW","P")){
      x <- (x - angle_diff) %% 360
    }
    else if (experiment %in% c("I","J")){x <- (x + angle.diff) %% 360}
    else {cat("Experiment not found in list of experiments")}
  }
  return(x)
}
angle.calc = function(X,Y){
  if (X < 0 && Y > 0){angle = 2.5*pi - atan2(Y,X)
  } else {angle = 0.5*pi - atan2(Y,X)}
  return(angle*180/pi)
}

angle.opp = function(angle){ 
}

vector.angle.diff = function(angle1,angle2,degrees = T){
  if (degrees==T){
    angle1 = angle1/180*pi
    angle2 = angle2/180*pi
  }
  x = angle1 - angle2
  return((ifelse(x>pi,x-2*pi,ifelse(x<= -pi,x+2*pi,x)))/pi*180)
}

Moores.Rayleigh <- function(angles,R,degrees=F,bidirectional=F){
  
  N <- length(angles)
  
  if (!N == length(R)){
    cat("Angle and R vectors must be the same length.")
    return(NA)
  }
  # Check for object formats
  
  # Data transformations
  if (degrees == T){angles <- angles*(pi/180)}
  if (bidirectional == T){
    angles <- sapply(angles,function(x){
      if (x > pi){x <- x - pi}
      x*2})}
  
  # R*
  weights <- rank(R)
  X <- sum(sin(angles)*weights)
  Y <- sum(cos(angles)*weights)
  Rmo <- (sqrt(X^2+Y^2)) / N^1.5
  
  # Weighted mean direction
  if (X < 0 && Y > 0){angle <- 2.5*pi - atan2(Y,X)
  } else {angle <- 0.5*pi - atan2(Y,X)}
  
  if (bidirectional == T){angle <- angle / 2}
  angle.deg <- angle*(180/pi)
  
  # p-value table look-up
  d1 <- matrix(c(
    0,0.354,0.354,0.354,0.356,0.362,0.387,0.791,1.049,1.058,1.06,1.061,1.061,1.061,
    0,0.014,0.04,0.064,0.116,0.18,0.272,0.693,1.039,1.095,1.124,1.143,1.149,1.154,
    0,0.02,0.048,0.072,0.116,0.168,0.243,0.62,1.008,1.09,1.146,1.192,1.212,1.238,
    0,0.023,0.051,0.072,0.115,0.163,0.234,0.588,0.988,1.084,1.152,1.216,1.25,1.298,
    0,0.022,0.05,0.071,0.111,0.158,0.226,0.568,0.972,1.074,1.152,1.23,1.275,1.345,
    0,0.021,0.049,0.069,0.107,0.154,0.22,0.556,0.959,1.066,1.15,1.238,1.291,1.373,
    0,0.02,0.048,0.067,0.105,0.151,0.216,0.546,0.949,1.059,1.148,1.242,1.3,1.397,
    0,0.02,0.047,0.066,0.104,0.148,0.213,0.538,0.94,1.053,1.146,1.245,1.307,1.416,
    0,0.02,0.046,0.065,0.103,0.146,0.21,0.532,0.934,1.048,1.144,1.248,1.313,1.432,
    0,0.019,0.045,0.064,0.101,0.144,0.206,0.523,0.926,1.042,1.14,1.252,1.322,1.456,
    0,0.019,0.044,0.063,0.1,0.142,0.204,0.518,0.92,1.037,1.136,1.252,1.325,1.47,
    0,0.019,0.044,0.063,0.099,0.141,0.202,0.514,0.914,1.031,1.132,1.25,1.327,1.48,
    0,0.019,0.044,0.062,0.099,0.14,0.2,0.51,0.91,1.027,1.129,1.248,1.328,1.487,
    0,0.019,0.043,0.062,0.098,0.139,0.199,0.507,0.906,1.024,1.127,1.247,1.329,1.492,
    0,0.019,0.043,0.061,0.097,0.138,0.198,0.505,0.903,1.022,1.126,1.246,1.33,1.496,
    0,0.019,0.043,0.061,0.097,0.137,0.197,0.503,0.901,1.021,1.125,1.246,1.331,1.499,
    0,0.019,0.043,0.06,0.096,0.137,0.196,0.502,0.899,1.019,1.124,1.246,1.332,1.501,
    0,0.019,0.043,0.06,0.096,0.136,0.196,0.5,0.897,1.018,1.124,1.246,1.333,1.502,
    0,0.018,0.042,0.06,0.095,0.136,0.195,0.499,0.896,1.016,1.123,1.245,1.334,1.502,
    0,0.018,0.042,0.059,0.094,0.134,0.193,0.494,0.891,1.012,1.119,1.243,1.332,1.504,
    0,0.018,0.042,0.059,0.093,0.133,0.191,0.489,0.887,1.007,1.115,1.241,1.329,1.506,
    0,0.018,0.041,0.059,0.093,0.132,0.19,0.487,0.883,1.005,1.113,1.24,1.329,1.508,
    0,0.018,0.041,0.058,0.093,0.132,0.189,0.485,0.881,1.004,1.112,1.24,1.329,1.509,
    0,0.018,0.041,0.058,0.092,0.131,0.187,0.481,0.876,0.999,1.109,1.239,1.329,1.517),
    nrow=24,byrow=T,dimnames = list(
      c(2,3,4,5,6,7,8,9,10,12,14,16,18,20,22,24,26,28,30,40,60,80,100,">100"),
      c(1,0.999,0.995,0.99,0.975,0.95,0.9,0.5,0.1,0.05,0.025,0.01,0.005,0.001)))
  
  if (N > 100){
    Rdist <- d1[">100",]
    
  } else if (N %in% rownames(d1)){
    Rdist <- d1[toString(N),]
    
  } else {
    h <- 10
    for (i in as.numeric(rownames(d1)[10:23])){
      if (N < i){
        Rdist <- d1[toString(h),]+(N-h)/(i-h)*(d1[toString(i),]-d1[toString(h),])
        break
      }
      h <- i
    }
  }
  h = 1
  for (i in Rdist){
    if (Rmo > i){
      if (h == length(Rdist)){
        p <- paste("p <",names(Rdist)[[h]])
        break
      }
      h <- h + 1
      next}
    p <- paste(names(Rdist)[[h]],"< p <",names(Rdist)[[h-1]])
    break
  }
  return(list("Weighted mean direction (radians)" = angle,
              "Weighted mean direction (degrees)" = angle.deg,
              "R*" = Rmo,
              "p" = p))
}
Moores.Rayleigh.bootstrap <- function(angles,R,degrees=F,bidirectional=F,interval=0.95,reps=100000){
  
  N <- length(angles)
  
  if (!N == length(R)){
    cat("Angle and R vectors must be the same length.")
    return(NA)
  }
  
  # Data transformations
  if (degrees == T){angles <- angles*(pi/180)}
  if (bidirectional == T){
    angles <- sapply(angles,function(x){
      if (x > pi){x <- x - pi}
      x*2})}
  
  weights <- rank(R)
  X <- sum(sin(angles)*weights)
  Y <- sum(cos(angles)*weights)
  if (X < 0 && Y > 0){mean.dir <- 2.5*pi - atan2(Y,X)
  } else {mean.dir <- 0.5*pi - atan2(Y,X)}
  
  gen.directions <- vector(length = reps)
  
  # Bootstrap a vector of resampled mean directions
  for (i in 1:reps){
    dist <- sample.int(N,N,replace=T)
    X <- sum(sin(angles[dist])*weights[dist])
    Y <- sum(cos(angles[dist])*weights[dist])
    if (X < 0 && Y > 0){gen.directions[[i]] <- 2.5*pi - atan2(Y,X)
    } else {gen.directions[[i]] <- 0.5*pi - atan2(Y,X)}
  }
  
  # Angle differences of generated directions to mean direction
  angle.diff <- function(x){x <- x - mean.dir
  if(x<=-pi){x <- x+2*pi}
  else if(x>=pi){x <- x %% pi}
  return(x)
  }
  gen.differences <- sapply(gen.directions,angle.diff)
  intervals <- c(interval+(1-interval)/2,(1-interval)/2)
  difference <- quantile(gen.differences,intervals)
  
  clockwise.int <- (mean.dir + difference[[1]]) %% (2*pi)
  a.clockwise.int <- (mean.dir + difference[[2]]) %% (2*pi)
  
  if (bidirectional == T){
    clockwise.int <- clockwise.int / 2
    a.clockwise.int <- a.clockwise.int / 2
  }
  
  if (degrees == T){
    clockwise.int <- clockwise.int *(180/pi)
    a.clockwise.int <- a.clockwise.int *(180/pi)
  }
  
  return(list("clockwise bound" = clockwise.int,
              "anti-clockwise bound" = a.clockwise.int))
}
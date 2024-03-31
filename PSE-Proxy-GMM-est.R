##Estimator 2: GMM estimator
library(numDeriv)

##Functions from Miao et al.(2018)
# GMM function
GMMF <- function (mrf ,para , data ){
  g0 <- mrf ( para =para , data = data )
  g <- apply (g0 ,2, mean )
  gmmf <- sum (g^2)
  return ( gmmf )
}

# Derivative of score equations
G1 <- function (bfun ,para , data ){
  G1 <- apply ( bfun (para , data ),2, mean )
  return (G1)
}
G <- function (bfun ,para , data ){
  G <- jacobian ( func =G1 , bfun =bfun ,x=para , data = data )
  return (G)
}

# Variance estimation
VAREST <- function (bfun ,para , data ){
  bG <- solve (G(bfun ,para , data ))
  bg <- bfun (para , data)
  spsz <- dim (bg )[1]
  Omega <- t(bg)%*%bg/ spsz
  Sigma <- bG%*% Omega %*%t(bG)
  return ( Sigma / spsz )
}

# Confidence interval
CNFINTVl <- function (esti , ci ){
  esti <- as.matrix ( esti )
  dm <- dim ( esti )[2]
  para <- esti [ ,1:( dm/2)]
  dvar <- esti [,( dm/ 2+1): dm]
  z <- -qnorm ((1 - ci)/2)
  dsd <- sqrt ( dvar )
  return ( list (lb=para -z*dsd , ub= para +z*dsd ))
}

# Coverage probability
CVRPRB <- function (esti ,ci , trvlu ){
  esti <- as.matrix ( esti )
  dm <- dim ( esti )[2]
  para <- esti [ ,1:( dm/2)]
  dvar <- esti [,( dm/ 2+1): dm]
  z <- -qnorm ((1 - ci)/2)
  dsd <- sqrt ( dvar )
  lb <- para -z*dsd ; ub <- para +z*dsd
  return (trvlu >= lb&trvlu <= ub)
}

# P- value based on normal approximation
PVALUE <- function ( esti ){
  esti <- as.matrix ( esti )
  dm <- dim ( esti )[2]
  para <- esti [ ,1:( dm/2)]
  dvar <- esti [,( dm/ 2+1): dm]
  dsd <- sqrt ( dvar )
  return ((1 - pnorm (abs ( para ), mean =0, sd=dsd ))*2)
}

# Moment restriction function
NCmrf <- function (para , data1 ){
  if(is.na(data1$C[1])==T) {
    X <- as.matrix ( data1 $X); Y <- as.matrix ( data1 $Y)
    Z <- as.matrix ( data1 $Z); W <- as.matrix ( data1 $W)
    hlink <- cbind (1,X,W) %*% para
    g0 <- cbind (1,X,Z)
  }
  if(is.na(data1$C[1])==F) {
    X <- as.matrix ( data1 $X); Y <- as.matrix ( data1 $Y)
    Z <- as.matrix ( data1 $Z); W <- as.matrix ( data1 $W)
    C <- as.matrix ( data1 $C)
    hlink <- cbind (1,X,C,W) %*% para
    g0 <- cbind (1,X,C,Z)
  }
  g <- (as.vector (Y - hlink)) * g0
  return (g)
}

PSE_Proxy_GMM_est <- function(data1,effect){
  #'@param: data1 a list including vectors X,Y,Z,C and W, where C can be a dataframe, e.g., data1 <- list(X=data$M1,Y=data$M3,Z=data$M2,C=data$X,W=data$W)
  #'@param: effect the coefficients from linear regression of Y~X+C+W as the entered initial value of parameters
  
  inioptim = effect
  hpar <- optim (par = inioptim,
                 fn = GMMF,
                 mrf = NCmrf, data = data1,
                 method = "BFGS", hessian = FALSE)$par
  # This is the NC estimator of the structural parameter
  est <- as.numeric (hpar [2])
  var_est <- diag(VAREST(NCmrf,hpar,data1))[2]
  CILOW <- est-1.96*sqrt(var_est)
  CIUP <- est+1.96*sqrt(var_est)
  result <- data.frame(matrix(c(est,var_est,CILOW,CIUP),1,4))
  colnames(result) <- c("est","var_est","CILOW","CIUP")
  return (result)
}

##example##
##data1 <- list(X=data$M1,Y=data$M3,Z=data$M2,C=data$X,W=data$W)
##effectinp <- as.numeric(lm(M3~M1+X+W,data = data)$coefficients)##Y~X+C+W
##a_m1m3_GMM <- PSE_Proxy_GMM_est(data1,effectinp)

##data1 <- list(X=data$M3,Y=data$M4,Z=data$W,C=cbind(data$X,data$M1),W=data$M2)
##effectinp <- as.numeric(lm(M4~M3+X+M1+M2,data = data)$coefficients)##Y~X+C+W
##a_m3m4_GMM <- PSE_Proxy_GMM_est(data1,effectinp)
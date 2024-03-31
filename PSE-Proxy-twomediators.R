##PSE-Proxy for two non-ordered mediators
##Oservations with missing in any varaibles are removed in this function, in addition, the results are obtained through normalizing all variables

twocausallynonordered <- function(nameuse,data){
  #'@param: data a dataframe including all variables
  #'@param: nameuse a vector of colnames of X,M1,M2,Y,W in turn, e.g., nameuse <- c("BMI","TC","PLT","lbapwv","TB")

  data_x <- data[,c(1,which(colnames(data) %in% nameuse[1]))]
  data_m1 <- data[,c(1,which(colnames(data) %in% nameuse[2]))]
  data_m2 <- data[,c(1,which(colnames(data) %in% nameuse[3]))]
  data_y <- data[,c(1,which(colnames(data) %in% nameuse[4]))]
  data_w <- data[,c(1,which(colnames(data) %in% nameuse[5]))]
  
  data_sum <- merge(data_x,data_m1,all.x= T, by = "id_num")
  data_sum <- merge(data_sum,data_m2,all.x= T, by = "id_num")
  data_sum <- merge(data_sum,data_y,all.x= T, by = "id_num")
  data_sum <- merge(data_sum,data_w,all.x= T, by = "id_num")
  
  colnames(data_sum) <- c("id_num","X","M1","M2","Y","W")
  
  data <- data_sum[-1]
  data <- na.omit(data)
  data <- as.data.frame(scale(data))
  
  ##naive ols
  a_M1XEST_naive <- lm(M1~X,data = data)$coefficients[2]
  a_YM1EST_naive <- lm(Y~M1+X,data = data)$coefficients[2]
  a_M2XEST_naive <- lm(M2~X,data = data)$coefficients[2]
  a_YM2EST_naive <- lm(Y~M2+X,data = data)$coefficients[2]
  
  var_M1XEST_naive <- (summary(lm(M1~X,data = data))$coefficients[2,2])^2
  var_YM1EST_naive <- (summary(lm(Y~M1+X,data = data))$coefficients[2,2])^2
  var_M2XEST_naive <- (summary(lm(M2~X,data = data))$coefficients[2,2])^2
  var_YM2EST_naive <- (summary(lm(Y~M2+X,data = data))$coefficients[2,2])^2  
  
  b_NIE1_naive <- a_M1XEST_naive*a_YM1EST_naive
  b_NIE2_naive <- a_M2XEST_naive*a_YM2EST_naive
  b_NDE_naive <- lm(Y~X+M1+M2,data=data)$coefficients[2]
  
  var_NIE1_naive <- (a_M1XEST_naive^2) * var_YM1EST_naive + (a_YM1EST_naive^2) * var_M1XEST_naive + var_YM1EST_naive*var_M1XEST_naive
  var_NIE2_naive <- (a_M2XEST_naive^2) * var_YM2EST_naive + (a_YM2EST_naive^2) * var_M2XEST_naive + var_YM2EST_naive*var_M2XEST_naive
  var_NDE_naive <- (summary(lm(Y~X+M1+M2,data=data))$coefficients[2,2])^2
  
  #####PSE_Proxy#####
  #step1 X-->M
  a_M1XEST_P2SLS <- PSE_Proxy_P2SLS_est("X","M1","M2","W",NA,data)[1]
  a_M2XEST_P2SLS <- PSE_Proxy_P2SLS_est("X","M2","M1","W",NA,data)[1]
  var_M1XEST_P2SLS <- PSE_Proxy_P2SLS_est("X","M1","M2","W",NA,data)[2]
  var_M2XEST_P2SLS <- PSE_Proxy_P2SLS_est("X","M2","M1","W",NA,data)[2]
  
  
  data1 <- list(X=data$X,Y=data$M1,Z=data$M2,C=NA,W=data$W)
  effectinp <- as.numeric(lm(M1~X+W,data = data)$coefficients)
  a_M1XEST_GMM <- PSE_Proxy_GMM_est(data1,effectinp)[1]
  var_M1XEST_GMM <- PSE_Proxy_GMM_est(data1,effectinp)[2]
  
  data1 <- list(X=data$X,Y=data$M2,Z=data$M1,C=NA,W=data$W)
  effectinp <- as.numeric(lm(M2~X+W,data = data)$coefficients)
  a_M2XEST_GMM <- PSE_Proxy_GMM_est(data1,effectinp)[1]
  var_M2XEST_GMM <- PSE_Proxy_GMM_est(data1,effectinp)[2]
  
  #step1 M-->Y
  a_YM1EST_P2SLS <- PSE_Proxy_P2SLS_est("M1","Y","W","M2","X",data)[1]
  a_YM2EST_P2SLS <- PSE_Proxy_P2SLS_est("M2","Y","W","M1","X",data)[1]
  var_YM1EST_P2SLS <- PSE_Proxy_P2SLS_est("M1","Y","W","M2","X",data)[2]
  var_YM2EST_P2SLS <- PSE_Proxy_P2SLS_est("M2","Y","W","M1","X",data)[2]
  
  data1 <- list(X=data$M1,Y=data$Y,Z=data$W,C=data$X,W=data$M2)
  effectinp <- as.numeric(lm(Y~M1+X+M2,data = data)$coefficients)
  a_YM1EST_GMM <- PSE_Proxy_GMM_est(data1,effectinp)[1]
  var_YM1EST_GMM <- PSE_Proxy_GMM_est(data1,effectinp)[2]
  
  data1 <- list(X=data$M2,Y=data$Y,Z=data$W,C=data$X,W=data$M1)
  effectinp <- as.numeric(lm(Y~M2+X+M1,data = data)$coefficients)
  a_YM2EST_GMM <- PSE_Proxy_GMM_est(data1,effectinp)[1]
  var_YM2EST_GMM <- PSE_Proxy_GMM_est(data1,effectinp)[2]
  
  
  #Y~
  Y_P2SLS <- data$Y - a_YM1EST_P2SLS*data$M1 - a_YM2EST_P2SLS*data$M2
  Y_GMM <- data$Y - a_YM1EST_GMM*data$M1 - a_YM2EST_GMM*data$M2
  data2 <- data.frame(data$W,data$X,data$M1,data$M2,Y_P2SLS,Y_GMM)
  colnames(data2) <- c("W","X","M1","M2","Y_P2SLS","Y_GMM")
  
  b_NDE_P2SLS <- PSE_Proxy_P2SLS_est("X","Y_P2SLS",c("M1","M2"),"W",NA,data2)[1]
  var_NDE_P2SLS <- PSE_Proxy_P2SLS_est("X","Y_P2SLS",c("M1","M2"),"W",NA,data2)[2]
  
  data1 <- list(X=data2$X,Y=data2$Y_GMM,Z=data2$M1,C=NA,W=data2$W)##only choose M1 as Z
  effectinp <- as.numeric(lm(Y_GMM~X+W,data = data2)$coefficients)
  b_NDE_GMM <- PSE_Proxy_GMM_est(data1,effectinp)[1]
  var_NDE_GMM <- PSE_Proxy_GMM_est(data1,effectinp)[2]
  
  b_NIE1_P2SLS <- a_M1XEST_P2SLS*a_YM1EST_P2SLS
  var_NIE1_P2SLS <- (a_M1XEST_P2SLS^2) * var_YM1EST_P2SLS + (a_YM1EST_P2SLS^2) * var_M1XEST_P2SLS + var_YM1EST_P2SLS*var_M1XEST_P2SLS
  b_NIE2_P2SLS <- a_M2XEST_P2SLS*a_YM2EST_P2SLS
  var_NIE2_P2SLS <- (a_M2XEST_P2SLS^2) * var_YM2EST_P2SLS + (a_YM2EST_P2SLS^2) * var_M2XEST_P2SLS + var_YM2EST_P2SLS*var_M2XEST_P2SLS
  b_NIE1_GMM <- a_M1XEST_GMM*a_YM1EST_GMM
  var_NIE1_GMM <- (a_M1XEST_GMM^2) * var_YM1EST_GMM + (a_YM1EST_GMM^2) * var_M1XEST_GMM + var_YM1EST_GMM*var_M1XEST_GMM
  b_NIE2_GMM <- a_M2XEST_GMM*a_YM2EST_GMM
  var_NIE2_GMM <- (a_M2XEST_GMM^2) * var_YM2EST_GMM + (a_YM2EST_GMM^2) * var_M2XEST_GMM + var_YM2EST_GMM*var_M2XEST_GMM
  
  result <- data.frame(a_M1XEST_naive,a_YM1EST_naive,a_M2XEST_naive,a_YM2EST_naive,b_NIE1_naive,b_NIE2_naive,b_NDE_naive,
                       var_M1XEST_naive,var_YM1EST_naive,var_M2XEST_naive,var_YM2EST_naive,var_NIE1_naive,var_NIE2_naive,var_NDE_naive,
                       a_M1XEST_P2SLS,a_YM1EST_P2SLS,a_M2XEST_P2SLS,a_YM2EST_P2SLS,b_NIE1_P2SLS,b_NIE2_P2SLS,b_NDE_P2SLS,
                       var_M1XEST_P2SLS,var_YM1EST_P2SLS,var_M2XEST_P2SLS,var_YM2EST_P2SLS,var_NIE1_P2SLS,var_NIE2_P2SLS,var_NDE_P2SLS,
                       a_M1XEST_GMM,a_YM1EST_GMM,a_M2XEST_GMM,a_YM2EST_GMM,b_NIE1_GMM,b_NIE2_GMM,b_NDE_GMM,
                       var_M1XEST_GMM,var_YM1EST_GMM,var_M2XEST_GMM,var_YM2EST_GMM,var_NIE1_GMM,var_NIE2_GMM,var_NDE_GMM)
  
  
  names(result) <- c("a_M1XEST_naive","a_YM1EST_naive","a_M2XEST_naive","a_YM2EST_naive","b_NIE1_naive","b_NIE2_naive","b_NDE_naive",
                     "var_M1XEST_naive","var_YM1EST_naive","var_M2XEST_naive","var_YM2EST_naive","var_NIE1_naive","var_NIE2_naive","var_NDE_naive",
                     "a_M1XEST_P2SLS","a_YM1EST_P2SLS","a_M2XEST_P2SLS","a_YM2EST_P2SLS","b_NIE1_P2SLS","b_NIE2_P2SLS","b_NDE_P2SLS",
                     "var_M1XEST_P2SLS","var_YM1EST_P2SLS","var_M2XEST_P2SLS","var_YM2EST_P2SLS","var_NIE1_P2SLS","var_NIE2_P2SLS","var_NDE_P2SLS",
                     "a_M1XEST_GMM","a_YM1EST_GMM","a_M2XEST_GMM","a_YM2EST_GMM","b_NIE1_GMM","b_NIE2_GMM","b_NDE_GMM",
                     "var_M1XEST_GMM","var_YM1EST_GMM","var_M2XEST_GMM","var_YM2EST_GMM","var_NIE1_GMM","var_NIE2_GMM","var_NDE_GMM")
  
  matched_columns <- grep(paste(c("NDE","NIE"), collapse = "|"), names(result), value = TRUE)
  result_use <- result[, matched_columns]
  
  matched_columns <- grep("naive", names(result_use), value = TRUE)
  result_use_naive <- result_use[, matched_columns]
  matched_columns <- grep("P2SLS", names(result_use), value = TRUE)
  result_use_P2SLS <- result_use[, matched_columns]
  matched_columns <- grep("GMM", names(result_use), value = TRUE)
  result_use_GMM <- result_use[, matched_columns]
  
  result_use_naive$method <- "Naive"
  result_use_P2SLS$method <- "P2SLS"
  result_use_GMM$method <- "GMM"
  
  colnames(result_use_naive) <- c("b_NIE1","b_NIE2","b_NDE","var_NIE1","var_NIE2","var_NDE","method")
  colnames(result_use_P2SLS) <- c("b_NIE1","b_NIE2","b_NDE","var_NIE1","var_NIE2","var_NDE","method")
  colnames(result_use_GMM) <- c("b_NIE1","b_NIE2","b_NDE","var_NIE1","var_NIE2","var_NDE","method")
  
  result_use <- rbind(result_use_naive,result_use_P2SLS,result_use_GMM)
  result_use$ci_NIE1 <- paste0("(",round(result_use$b_NIE1-1.96*sqrt(result_use$var_NIE1),3),",",round(result_use$b_NIE1+1.96*sqrt(result_use$var_NIE1),3),")")
  result_use$ci_NIE2 <- paste0("(",round(result_use$b_NIE2-1.96*sqrt(result_use$var_NIE2),3),",",round(result_use$b_NIE2+1.96*sqrt(result_use$var_NIE2),3),")")
  result_use$ci_NDE <- paste0("(",round(result_use$b_NDE-1.96*sqrt(result_use$var_NDE),3),",",round(result_use$b_NDE+1.96*sqrt(result_use$var_NDE),3),")")
  
  result_use$b_NIE1 <- round(result_use$b_NIE1,3)
  result_use$b_NIE2 <- round(result_use$b_NIE2,3)
  result_use$b_NDE <- round(result_use$b_NDE,3)
  
  result_use <- result_use[,c("method","b_NDE","ci_NDE","b_NIE1","ci_NIE1","b_NIE2","ci_NIE2")]
  return(result_use)
}

###example:
# nameuse <- c("BMI","TC","PLT","lbapwv","TB")
# cau.hat1=twocausallynonordered(nameuse = nameuse, data = data)
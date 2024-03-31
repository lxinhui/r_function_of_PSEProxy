##Estimator 3: P2SLS
library(ivreg)

PSE_Proxy_P2SLS_est <- function(exposure,outcome,Eproxy,Oproxy,covariants,datause){
  #'@param: datause a dataframe including all variables;
  #'@param: exposure colname of the exposure in datause, e.g., exposure = "X";
  #'@param: outcome colname of the outcome in datause, e.g., outcome = "Y";
  #'@param: Eproxy colname of the Eproxy in datause (can be a vector), e.g., Eproxy = "Z";
  #'@param: Oproxy colname of the Oproxy in datause, e.g., Oproxy = "W";
  #'@param: covariants colname of the covariants in datause (can be a vector), e.g., covariants = c("C1","C2");
  
  Eproxy <- paste(Eproxy,collapse = "+")
  
  #no covariate
  if(is.na(covariants[1])==T){
    func <- paste0(outcome,"~",exposure,"+",Oproxy,"|",exposure,"+",Eproxy)
    estimate <- as.numeric(ivreg(func,data = datause)$coefficients[2])
    var_estimate <- as.numeric(summary(ivreg(func,data = datause))$coefficients[2,2])^2
  }  
  
  if(is.na(covariants[1])==F){ 
    func_cov <- paste(covariants,collapse = "+")
      
    func <- paste0(outcome,"~",exposure,"+",func_cov,"+",Oproxy,"|",exposure,"+",func_cov,"+",Eproxy)
    estimate <- as.numeric(ivreg(func,data = datause)$coefficients[2])
    var_estimate <- (as.numeric(summary(ivreg(func,data = datause))$coefficients[2,2]))^2
  }
  
  results <- data.frame(matrix(c(estimate,var_estimate),1,2))
  colnames(results) <- c("est","var_est")
  return(results)
}

###example:
##PSE_Proxy_P2SLS_est("X","M1","M2","W",NA,data)
##PSE_Proxy_P2SLS_est("M2","M4","W","M3",c("X","M1"),data)
##PSE_Proxy_P2SLS_est("X","M1",c("M22","M32"),"W",NA,data)
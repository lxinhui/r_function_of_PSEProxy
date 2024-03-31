##Estimator 1: Estimator through Identification function

PSE_Proxy_IF_est <- function(exposure,outcome,Eproxy,Oproxy,covariants,datause){
  #'@param: datause a dataframe including all variables;
  #'@param: exposure colname of the exposure in datause, e.g., exposure = "X";
  #'@param: outcome colname of the outcome in datause, e.g., outcome = "Y";
  #'@param: Eproxy colname of the Eproxy in datause, e.g., Eproxy = "Z";
  #'@param: Oproxy colname of the Oproxy in datause, e.g., Oproxy = "W";
  #'@param: covariants colname of the covariants in datause (can be a vector), e.g., covariants = c("C1","C2");
  
  if(is.na(covariants[1])==T){
    func1 <- paste0(outcome,"~",exposure,"+",Eproxy)  
    func2 <- paste0(Oproxy,"~",exposure,"+",Eproxy)
  }  
  
  if(is.na(covariants[1])==F){ 
    func_cov <- paste(covariants,collapse = "+")
    func1 <- paste0(outcome,"~",exposure,"+",Eproxy,"+",func_cov)  
    func2 <- paste0(Oproxy,"~",exposure,"+",Eproxy,"+",func_cov)
  }
  
  model1_x <- summary(lm(func1,data = datause))$coefficients[2,1]
  model1_e <- summary(lm(func1,data = datause))$coefficients[3,1]
  model2_x <- summary(lm(func2,data = datause))$coefficients[2,1]
  model2_e <- summary(lm(func2,data = datause))$coefficients[3,1]
  
  estimate <- model1_x - (model1_e/model2_e)*model2_x
  
  varmodel1_x <- (summary(lm(func1,data=datause))$coefficients[2,2])^2
  varmodel1_e <- (summary(lm(func1,data=datause))$coefficients[3,2])^2
  varmodel2_x <- (summary(lm(func2,data=datause))$coefficients[2,2])^2
  varmodel2_e <- (summary(lm(func2,data=datause))$coefficients[3,2])^2
  
  a_ratio <- model1_e/model2_e
  var_ratio <- ((model1_e/model2_e)^2)*((varmodel1_e/(model1_e)^2)+(varmodel2_e/(model2_e)^2))
  
  varestimate <- varmodel1_x+((a_ratio^2) * varmodel2_x + (model2_x^2) * var_ratio + varmodel2_x*var_ratio)
  
  results <- data.frame(matrix(c(estimate,varestimate),1,2))
  colnames(results) <- c("est","var_est")
  return(results)
}
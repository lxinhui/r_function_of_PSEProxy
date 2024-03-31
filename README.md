# r_function_of_PSEProxy

PSE-Proxy is a multiple mediation analysis method based on proxy variable by controlling unmeasured confounders. There are differences in the specific effect identification process for different causal graph models. Therefore, we provide the basis functions for the estimation of causal pathes between variables and present the estimation functions for direct and two indirect effects for ordered two mediations' case. 
 
The basis function of causal edge estimation between variables is based on three estimation methods: identification formula method, P2SLS method and GMM method. Corresponding function respectively: PSE_Proxy_IF_est (), PSE_Proxy_P2SLS_est () and PSE_Proxy_GMM_est (). It is worth mentioning that PSE_Proxy_GMM_est() is derived from the GMM estimation function proposed by Miao et al（2018）. with a slight modification. Further, taking the case of ordered two mediators as an example, we provide the R function for mediation analysis twocausallynonordered(): 
 
twocausallynonordered(nameuse,data), 
 
The parameter nameuse is a vector of column names for exposure, mediator 1, mediator 2, outcome, and proxy variables, e.g. nameuse <- c ("BMI","TC","PLT","lbapwv","TB"); 
 
The parameter data is the input data set dataframe. The first column should be the individual identification number, and the column name should be named "id_num". 




Reference:

Miao W, Shi X, Tchetgen E T. A confounding bridge approach for double negative control inference on causal effects[J]. arXiv preprint arXiv:1808.04945, 2018.

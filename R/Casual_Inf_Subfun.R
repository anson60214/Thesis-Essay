
#2022 new measure and ultra-high dim
#with D considerable

choose_p_D <- function(test_data,cov_e){
  data = test_data
  xi = c()
  res = c()
  d = 0
  C = 0
  n = nrow(data)
  p = ncol(data)-2
  
  index0 <- which(data$D == 0)
  index1 <- which(data$D == 1)
  alpha = length(index0)/n
  for( i in 1:p){
    w_0 = data[index0,i]
    var_e = cov_e[i,i]
    var_w_0 = var(w_0)
    x_0 = mean(w_0) + as.numeric( (var_w_0-var_e)*solve(var_w_0) ) *(w_0-mean(w_0))
    
    w_1 = data[index1,i]
    var_e = cov_e[i,i]
    var_w_1 = var(w_1)
    x_1 = mean(w_1) + as.numeric( (var_w_1-var_e)*solve(var_w_1) ) *(w_1-mean(w_1))
    
    C = alpha*calculateXI(x_0 , data$Y[index0] ,seed = sample(1:9999,1)) + 
      (1-alpha)*calculateXI(x_1 , data$Y[index1],
                            seed =sample(1:9999,1))
    xi = c(xi , C ) 
  }
  
  var.list = colnames(data)[1:p]
  names(xi) <- c( var.list )
  
  d = round(n/log(n),0)
  selected_Var = sort(xi,decreasing = TRUE)[1:d]
  return(selected_Var)
}

#2022 new measure and ultra-high dim

choose_p <- function(test_data,cov_e){
  data = test_data
  xi = c()
  res = c()
  d = 0
  C = 0
  n = nrow(data)
  p = ncol(data)-2
  
  for( i in 1:p){
    w = data[,i]
    var_e = cov_e[i,i]
    var_w = var(w)
    x = mean(w) + as.numeric( (var_w-var_e)*solve(var_w) )*(w-mean(w))
    C = calculateXI(x,  data$Y ,seed = sample(1:9999,1))
    xi = c(xi , C) 
  }
  
  var.list = colnames(data)[1:p]
  names(xi) <- c( var.list )
  
  d = round(n/log(n),0)
  selected_Var = sort(xi,decreasing = TRUE)[1:d]
  return(selected_Var)
}

## define some functions for generating data, ATE estimates, and the wAMD, 

expit = function(x){ 
  pr = ( exp(x) / (1+exp(x)) ) 
  return(pr)
}
ATE_est = function(fY,fw,fA,fp){
  t_ATE = fY*fw
  tt_ATE = ( ( sum(t_ATE[fA==1])  / sum(fw[fA==1]) ) - ( sum(t_ATE[fA==0]) /  sum(fw[fA==0]) ) )
  return(tt_ATE) 
}


create_weights = function(fp,fA){
  fw = (fp)^(-1)
  fw[fA==0] = (1 - fp[fA==0])^(-1)
  return(fw)
}
wAMD_function = function(DataM,varlist,trt.var,wgt,beta){
  trt = untrt = diff_vec = rep(NA,length(beta)) 
  names(trt) = names(untrt) = names(diff_vec) = varlist
  for(jj in 1:length(varlist)){ 
    this.var = paste("w",varlist[jj],sep="") 
    DataM[,this.var] = DataM[,varlist[jj]] * DataM[,wgt] 
    trt[jj] = sum( DataM[DataM[,trt.var]==1, this.var ]) / sum(DataM[DataM[,trt.var]==1, wgt]) 
    untrt[jj] = sum(DataM[DataM[,trt.var]==0, this.var]) / sum(DataM[DataM[,trt.var]==0, wgt]) 
    diff_vec[jj] = abs( trt[jj] - untrt[jj] ) 
  } 
  wdiff_vec = diff_vec * abs(beta) 
  wAMD = c( sum(wdiff_vec))
  ret = list( diff_vec = diff_vec, wdiff_vec = wdiff_vec, wAMD = wAMD )
  return(ret) 
}

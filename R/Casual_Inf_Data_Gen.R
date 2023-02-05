#Data Generate

Data_Gen <- function(n,p,mean_x,sig_x,rho,alpha,beta,theta,a,var_e,e_distr="normal",num_pi,delta,linearY,typeY){
  Data = NULL
  n=n 
  p=p 
  mean_x=mean_x  
  sig_x=sig_x  
  rho=rho
  delta=delta 
  theta=theta
  a=a
  var_e=var_e
  var.list2=c()
  
  scale_list = as.data.frame(rbind(alpha,beta))
  
  colnames(scale_list)[ which(scale_list[1,]*scale_list[2,]!=0 ) ] = "Xc"
  colnames(scale_list)[ which(scale_list[1,]==0 ) ] = "Xp"
  colnames(scale_list)[ which(scale_list[2,]==0 ) ] = "Xi"
  
  pC = sum(colnames(scale_list)=="Xc" )  # pC: associate with both
  pP = sum(colnames(scale_list)=="Xp" )  # pP: associate with outcome
  pI = sum(colnames(scale_list)=="Xi" )  # pI: associate with treatment
  
  scale_1 = c( unlist( scale_list[1,colnames(scale_list)=="Xc"] ) )
  scale_2 = c( unlist( scale_list[1,colnames(scale_list)=="Xi"] ) )
  scale_3 = c( unlist( scale_list[2,colnames(scale_list)=="Xc"] ) )
  scale_4 = c( unlist( scale_list[2,colnames(scale_list)=="Xp"] ) )
  
  pS = p - (pC+pI+pP) # pS: associate with neither
  
  var.list2 = c( paste( "Xc", 1:pC, sep ="") , paste( "Xp", 1:pP, sep ="") , 
                 paste( "Xi", 1:pI, sep =""), paste( "XS", 1:pS, sep =""))
  
  Beta0 = rgamma(n,5,5)-delta
  Gamma0 = rgamma(n,5,5)-delta #delta=1 has equal num 1 0  ,default 0.5 since has small pi10 
  pi_10 = 1-exp(Beta0)/ (1+exp(Beta0))
  pi_01 = 1-exp(Gamma0)/ (1+exp(Gamma0))
  
  if(num_pi==1){pi_10=pi_01}
  
  #generate data 
  Sigma_x = matrix(rho*sig_x^2,nrow=length(var.list2),ncol=length(var.list2)) 
  diag(Sigma_x) = sig_x^2
  
  Mean_x = rep(mean_x,length(var.list2))
  X = as.matrix(mvrnorm(n = n,mu=Mean_x,Sigma = Sigma_x,empirical = FALSE))
  colnames(X) = var.list2
  
  #generate measurement error for X
  cov_e = matrix(0,p,p)
  diag(cov_e) = var_e
  colnames(cov_e) = var.list2
  mean_e = rep(0,p)
  if(e_distr == "normal"){
    e = mvrnorm(n = n, mu = mean_e, Sigma = cov_e, empirical = FALSE)
    }
  else{
    e = matrix(0,n,p)
    for(i in 1:p){
      e[,i] = rt(n,e_distr)
      } 
    }
  W = var_e*a + X + e
  
  # set associate with treatment 
  prob =  rowSums( cbind( X[,grepl("Xc", colnames(X))]*scale_1 , X[,grepl("Xi", colnames(X))]*scale_2 ) )
  prob = exp(prob)/(1+exp(prob))
  D = (prob >0.5)*1
  #D = rbinom(n,1,prob)
  D_star = rep(0,n)
  for( i in 1:n){
    pi_matrix = matrix( c(1-pi_10[i],pi_10[i],pi_01[i],1-pi_01[i]) ,nrow=2,ncol=2)
    prob_star = pi_matrix %*% rbind(prob[i],1-prob[i])
    D_star[i] = (prob_star[1,]>prob_star[2,])*1
  }
  
  if(linearY == TRUE){
    CP_vec = rowSums( cbind( X[,grepl("Xc", colnames(X))]*scale_3 , X[,grepl("Xp", colnames(X))]*scale_4) )
  }
  else{
    CP_vec =  rowSums( cbind(sin(X[,grepl("Xc", colnames(X))])*scale_3 , exp(X[,grepl("Xp", colnames(X))])*scale_4))
  }
  
  #set associate with outcome
  if(typeY == "binary") {
    prob = exp(D*theta + CP_vec ) / 
      (1+exp(D*theta + CP_vec ))
    y = rbinom(n, size = 1,prob ) 
  }
  else if(typeY == "pois"){
    prob = abs(D*theta + CP_vec )
    y = rpois(n,prob)
  }
  else {
    y = D*theta + CP_vec +rnorm(n,0,1)
  }
  
  Data = as.data.frame(cbind(X,D,y) )
  colnames(Data) <- c(var.list2,"D","Y")
  Error_Data = as.data.frame(cbind(W,D_star,y) )
  colnames(Error_Data) <-  c(var.list2,"D","Y")
  return( list(Data=Data, Error_Data = Error_Data, Pi = cbind(pi_10,pi_01 ), cov_e = cov_e) )
}


# for user to input X as cavariate

Data_GenX <- function(X,alpha,beta,theta,a,var_e,e_distr="normal",num_pi,delta,linearY,typeY){
  Data = NULL
  X=X
  n=nrow(X)
  p=ncol(X)
  delta=delta 
  theta=theta
  a=a
  var_e=var_e
  var.list2=c()
  
  scale_list = as.data.frame(rbind(alpha,beta))
  
  colnames(scale_list)[ which(scale_list[1,]*scale_list[2,]!=0 ) ] = "Xc"
  colnames(scale_list)[ which(scale_list[1,]==0 ) ] = "Xp"
  colnames(scale_list)[ which(scale_list[2,]==0 ) ] = "Xi"
  
  pC = sum(colnames(scale_list)=="Xc" )  # pC: associate with both
  pP = sum(colnames(scale_list)=="Xp" )  # pP: associate with outcome
  pI = sum(colnames(scale_list)=="Xi" )  # pI: associate with treatment
  
  scale_1 = c( unlist( scale_list[1,colnames(scale_list)=="Xc"] ) )
  scale_2 = c( unlist( scale_list[1,colnames(scale_list)=="Xi"] ) )
  scale_3 = c( unlist( scale_list[2,colnames(scale_list)=="Xc"] ) )
  scale_4 = c( unlist( scale_list[2,colnames(scale_list)=="Xp"] ) )
  
  pS = p - (pC+pI+pP) # pS: associate with neither
  
  var.list2 = c( paste( "Xc", 1:pC, sep ="") , paste( "Xp", 1:pP, sep ="") , 
                 paste( "Xi", 1:pI, sep =""), paste( "XS", 1:pS, sep =""))
  
  Beta0 = rgamma(n,5,5)-delta
  Gamma0 = rgamma(n,5,5)-delta #delta=1 has equal num 1 0  ,default 0.5 since has small pi10 
  pi_10 = 1-exp(Beta0)/ (1+exp(Beta0))
  pi_01 = 1-exp(Gamma0)/ (1+exp(Gamma0))
  
  if(num_pi==1){pi_10=pi_01}
  
  #generate data 
  colnames(X) = var.list2
  
  #generate measurement error for X
  cov_e = matrix(0,p,p)
  diag(cov_e) = var_e
  colnames(cov_e) = var.list2
  mean_e = rep(0,p)
  if(e_distr == "normal"){
    e = mvrnorm(n = n, mu = mean_e, Sigma = cov_e, empirical = FALSE)
  }
  else{
    e = matrix(0,n,p)
    for(i in 1:p){
      e[,i] = rt(n,e_distr,lower)
    } 
  }
  W = var_e*a + X + e
  
  
  # set associate with treatment 
  prob =  rowSums( cbind( X[,grepl("Xc", colnames(X))]*scale_1 , X[,grepl("Xi", colnames(X))]*scale_2 ) )
  prob = exp(prob)/(1+exp(prob))
  D = (prob >0.5)*1
  D_star = rep(0,n)
  for( i in 1:n){
    pi_matrix = matrix( c(1-pi_10[i],pi_10[i],pi_01[i],1-pi_01[i]) ,nrow=2,ncol=2)
    prob_star = pi_matrix %*% rbind(prob[i],1-prob[i])
    D_star[i] = (prob_star[1,]>prob_star[2,])*1
  }
  
  if(linearY == TRUE){
    CP_vec = rowSums( cbind( X[,grepl("Xc", colnames(X))]*scale_3 , X[,grepl("Xp", colnames(X))]*scale_4) )
  }
  else{
    CP_vec =  rowSums( cbind(sin(X[,grepl("Xc", colnames(X))])*scale_3 , exp(X[,grepl("Xp", colnames(X))])*scale_4))
  }
  
  #set associate with outcome
  if(typeY == "binary") {
    prob = exp(D*theta + CP_vec ) / 
      (1+exp(D*theta + CP_vec ))
    y = rbinom(n, size = 1,prob ) 
  }
  else if(typeY == "pois"){
    prob = abs(D*theta + CP_vec )
    y = rpois(n,prob)
  }
  else {
    y = D*theta + CP_vec +rnorm(n,0,1)
  }
  
  Data = as.data.frame(cbind(X,D,y) )
  colnames(Data) <- c(var.list2,"D","Y")
  Error_Data = as.data.frame(cbind(W,D_star,y) )
  colnames(Error_Data) <-  c(var.list2,"D","Y")
  return( list(Data=Data, Error_Data = Error_Data, Pi = cbind(pi_10,pi_01), cov_e = cov_e) )
}

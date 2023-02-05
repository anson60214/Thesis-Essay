source("Casual_Inf_Subfun.R")
source("Casual_Inf_Data_Gen.R")

##install and load the "MASS" package
#install.packages("MASS")
library(MASS) # version 3.3.1
library(dplyr)
library(grpreg)
library(Ball)
## procedure to install and load the "lqa" package from GitHub in two steps

# first install and load the "devtools" package 
#install.packages("devtools") 
library(devtools)

# now install and load "lqa" package from GitHub
#install_github("cran/lqa")  
library(lqa)

Causal.cor <- function(x, y, z, distance = FALSE) {
  x <- as.matrix(x)
  y <- as.matrix(y)
  index0 <- which(z == 0)
  index1 <- which(z == 1)
  alpha = length(index0)/length(z)
  if (distance == TRUE) {
    x0 <- x[index0, index0]
    y0 <- y[index0, index0]
    x1 <- x[index1, index1]
    y1 <- y[index1, index1]
    #see definition 4
    Causal.cor <- alpha*bcov(x0, x0, distance = TRUE)^2 + (1-alpha)*bcov(x1, y1, distance = TRUE)^2
  } else {
    x0 <- x[index0, ]
    y0 <- y[index0, ]
    x1 <- x[index1, ]
    y1 <- y[index1, ]
    #see definition 4
    Causal.cor <- alpha*bcov(x0, y0)^2 + (1-alpha)*bcov(x1, y1)^2
  }
  return(Causal.cor)
}


Ball_GOAL <- function(Data){
  X = Data[,1: (ncol(Data)-2)]
  D = Data[,(ncol(Data)-1)]
  Y = Data[,(ncol(Data))]
  p = ncol(X)
  n = nrow(X)
  threshold = min(30,p) 
  #set threshold for screening.
  ballcor = rep(NA, p)
  for (j in 1:p){
    # calculate conditional ball covariance for each variable.
    ballcor[j] <- Causal.cor(X[,j],Y,D)
  }
  
  # screening procedure 
  names(ballcor) = names(X)
  ballorder = order(ballcor, decreasing=TRUE)
  ballcor = sort(ballcor,decreasing = TRUE)[1:threshold]
  kersballcorye <- ballcor[sort(names(ballcor))]
  weight = ballcor
  var.list = names(ballcor)
  
  # weight for each variable for refined selection
  betaXY = weight
  betaXY = weight/max(weight)
  
  # the data matrix after screening
  Data = X[,var.list]
  Data = as.data.frame(Data)
  Data$A = D
  Data$Y = Y
  p = ncol(Data)-2
  n.p=n+p
  
  ## set vector of possible lambda's to try
  # set lambda values for grid search.
  
  lambda_vec =  c(-10, -5, -2, -1, -0.75, -0.5, -0.25, 0.25, 0.49)
  names(lambda_vec) = as.character(lambda_vec)
  
  ## lambda_n (n)^(gamma/2 - 1) = n^(gamma_convergence_factor)
  gamma_convergence_factor = 2
  
  ## get the gamma value for each value in the lambda vector that corresponds to convergence factor
  gamma_vals = 2*(gamma_convergence_factor - lambda_vec + 1)
  names(gamma_vals) = names(lambda_vec)
  
  ## normalize covariates to have mean 0 and standard deviation 1
  temp.mean = colMeans(Data[,var.list])
  Temp.mean = matrix(temp.mean,ncol=length(var.list),nrow=nrow(Data),byrow=TRUE)
  Data[,var.list] = Data[,var.list] - Temp.mean
  temp.sd = apply(Data[var.list],FUN=sd,MARGIN=2)
  Temp.sd = matrix(temp.sd,ncol=length(var.list),nrow=nrow(Data),byrow=TRUE)
  Data[var.list] = Data[,var.list] / Temp.sd
  rm(list=c("temp.mean","Temp.mean","temp.sd","Temp.sd"))
  

  ## want to save ATE, wAMD and propensity score coefficients for each lambda value
  ATE = wAMD_vec =coeff_XA=NULL
  ATE = wAMD_vec = rep(NA, length(lambda_vec))
  names(ATE) = names(wAMD_vec) = names(lambda_vec)
  coeff_XA = as.data.frame(matrix(NA,nrow=length(var.list),ncol=length(lambda_vec)))
  names(coeff_XA) = names(lambda_vec)
  rownames(coeff_XA) = var.list
  
  ## set the possible lambda2 value (taken from Zou and Hastie (2005))
  S_lam=c(0,10^c(-2,-1.5,-1,-0.75,-0.5,-0.25,0,0.25,0.5,1))
  
  ## want to save ATE, wAMD and propensity score coefficients for each lambda2 value
  WM_N=M_N=S_wamd=rep(NA,length(S_lam))
  M_mat=matrix(NA,length(S_lam),p)
  colnames(M_mat)= var.list
  
  
  for (m in 1:length(S_lam)) {
    ## create augmented A and X 
    lambda2=S_lam[m]
    I=diag(1,p,p)
    Ip=sqrt(lambda2)*I
    Anp=c(Data$A,rep(0,p))
    Xnp=matrix(0,n+p,p)
    X=Data[,var.list]
    for (j in 1:p){
      Xnp[,j]=c(X[,j],Ip[,j])
    }
    newData=as.data.frame(Xnp)
    names(newData)=var.list
    newData$A=Anp
    
    ## want to save ATE, wAMD and propensity score coefficients for each lambda value
    ATE_try=ATE = wAMD_vec =coeff_XA=NULL;
    ATE_try=ATE = wAMD_vec = rep(NA, length(lambda_vec))
    names(ATE) = names(wAMD_vec) = names(lambda_vec)
    coeff_XA = as.data.frame(matrix(NA,nrow=length(var.list),ncol=length(lambda_vec)))
    names(coeff_XA) = names(lambda_vec)
    rownames(coeff_XA) = var.list
    
    ## run GOAL with lqa using augmented data
    # weight model with all possible covariates included, this is passed into lasso function
    
    for( lil in names(lambda_vec) ){
      il = lambda_vec[lil]
      ig = gamma_vals[lil]
      
      # create the outcome adaptive lasso penalty with coefficient specific weights determined by outcome model
      oal_pen = adaptive.lasso(lambda=n.p^(il),al.weights = abs(betaXY)^(-ig) )
      # run outcome-adaptive lasso model with appropriate penalty
      X=as.matrix(newData[var.list]);y=as.vector(newData$A);
      logit_oal = lqa.default ( X,y, penalty=oal_pen, family=binomial(logit) )
      # save propensity score coefficients    !!!!adaptive elastic net
      coeff_XA[var.list,lil] = (1+lambda2)*coef(logit_oal)[var.list]
      # generate propensity score 
      Data[,paste("f.pA",lil,sep="")]=
        expit(as.matrix(cbind(rep(1,n),Data[var.list]))%*%as.matrix((1+lambda2)*coef(logit_oal)))
      
      # create inverse probability of treatment weights for ATE
      Data[,paste("w",lil,sep="")] = create_weights(fp=Data[,paste("f.pA",lil,sep="")],fA=Data$A)
      # estimate weighted absolute mean different over all covaraites using this lambda to generate weights
      wAMD_vec[lil] = wAMD_function(DataM=Data,varlist=var.list,trt.var="A",
                                    wgt=paste("w",lil,sep=""),beta=betaXY)$wAMD
      # save ATE estimate for this lambda value
      ATE[lil] = ATE_est(fY=Data$Y,fw=Data[,paste("w",lil,sep="")],fA=Data$A)
      
      
    } # close loop through lambda values
    
    # print out wAMD for all the lambda values evaluated
    wAMD_vec
    # find the lambda value that creates the smallest wAMD
    tt = which.min(wAMD_vec)
    # print out ATE corresponding to smallest wAMD value
    ATE[tt][[1]]
    # save the coefficients for the propensity score that corresponds to smallest wAMD value 
    GOAL.lqa.c=coeff_XA[,tt]
    names(GOAL.lqa.c) = var.list
    
    Coefficients_propensity_score_table = as.data.frame( matrix(GOAL.lqa.c,length(GOAL.lqa.c),3) )
    colnames(Coefficients_propensity_score_table) = c("Covariate","Coefficients Propensity Score","Significantness")
    Coefficients_propensity_score_table[,1] = var.list  
    Coefficients_propensity_score_table[,3] = ifelse(abs(Coefficients_propensity_score_table[,3])> 10^(-8),"***"," ")
    
    # check which covariates are selected
    M_mat[m,]=ifelse(abs(coeff_XA[,tt])> 10^(-8),1,0)
    # save the ATE corresponding to smallest wAMD value
    GOAL_est1=ATE[tt][[1]]
    GOAL.lqa.ate=GOAL_est1
    
    M_N[m]=GOAL.lqa.ate
    names(M_N) <- S_lam
    # save the smallest wAMD value
    WM_N[m]=wAMD_vec[tt][[1]]
  }
  return( list( ATE = GOAL.lqa.ate , wAMD = wAMD_vec, 
                Coefficients_propensity_score = Coefficients_propensity_score_table,
                Ball = ballcor ) )
} 
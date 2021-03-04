#------------------------------------------------------------------------
#Generating design matrix x, kind=1：Random Gaussian matrix (big p) 
get.x1=function(n,p,cor){
  X1=matrix(rnorm(n*p),nr=n,nc=p)
  X2=scale(X1, center = TRUE, scale = TRUE)
  X3=sqrt((n-1)/n)*X2
  data1=matrix(0,nr=n,nc=p)
  data1[,1]=X3[,1]
  data1[,p]=X3[,p]
  for(j in 2:(p-1)){
    data1[,j]=X3[,j]+cor*(X3[,j+1]+X3[,j-1])
  }
  return(data1)
}

#------------------------------------------------------------------------
#Generating design matrix x,kind=2：Classical Gaussian matrix (small p)
get.x2=function(n,p,cor){
mean=rep(0,p)
V=matrix(0,p,p)
for(i in 1:p){
for(j in 1:p){
V[i,j]=cor^abs(i-j)
}
}
library(MASS)
x=mvrnorm(n,mean,V)
return(x)
}

#------------------------------------------------------------------------
getdata=function(n,p,T1,cor,R,seednum,kind){
#========================================================================
# Generate data for simulation
#------------------------------------------------------------------------
# INPUT:
# 	n                sample size
# 	p                dimension 
# 	T1               sparsity level of true beta
#	cor              correlation parameter of X 
# 	R                ratio of the largest non-zero absolute element to the smallest non-zero absolute element in beta
# 	seednum          the seed number for repeatability and reproducibility
# 	kind    	     design matrix type:  1=Random Gaussian matrix (big p); 2=Classical Gaussian matrix (small p)                                       
# OUTPUT:
# 	betature           the true value of beta 
#     suppture           the true support of betature 
# 	X                  design matrix 
# 	y                  response vector 
#------------------------------------------------------------------------
# Written by:
# 	Yuling Jiao (yulingjiaomath@whu.edu.cn)
# 	Lican Kang (kanglican@whu.edu.cn)                              
# This version: May 2, 2020   
#======================================================================== 
# fix seed 
set.seed(seednum)
#create beta
m1=5*sqrt(2*log(p)/n)#1
m2=R*m1
set.seed(seednum)
b1=runif(T1,m1,m2)
supptrue=sample(p,T1)###真实支撑集
betatrue=rep(0,p)
betatrue[supptrue]=b1
prob.hat=numeric(n)
y=numeric(n)
if(kind==1){
x=get.x1(n,p,cor)
}else{
x=get.x2(n,p,cor)
}
for(i in 1:n){
prob.hat[i]=1/(1+exp(-sum(x[i,]*betatrue)))
y[i]=rbinom(1,1,prob.hat[i])
}
return(list(betatrue,supptrue,x,y))
}






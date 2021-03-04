#===================================================================
# Demo: solving L0 regularizatio problem via GSDAR for given sparsity level
#-------------------------------------------------------------------
# Written by:
#              Yuling Jiao (yulingjiaomath@whu.edu.cn)                         
#              Lican Kang (kanglican@whu.edu.cn)                                
# This version: May 2, 2020
#===================================================================
gc(rm(list=ls())) 
# parameter of generating data 
T1=10
n=300
p=5000
cor=0.2
R=100
seednum=2
kind=1
# parameter of GSDAR Algorithm
beta0=rep(0,p)
w=0.9
itermax=30
#------------------------------------------------------------------------
# Generate data 
source("gdataSimu.R")
res=getdata(n,p,T1,cor,R,seednum,kind)
betatrue=res[[1]]
supptrue=res[[2]]
x=res[[3]]
y=res[[4]]
#------------------------------------------------------------------------
# GSDAR computing
source("GSDAR.R")
ptm=proc.time()
result=GSDAR(beta0,T1,x,y,w,itermax)
Time=proc.time()-ptm
betahat=result[[1]]
sethat=result[[2]]
#
set.seed(seednum)
#Train set
train.num=sample(n,round(n*0.8))
x.train=x[train.num,]
y.train=y[train.num]
#Test set
test.num=(1:n)[-train.num]
x.test=x[test.num,]
y.test=y[test.num]
rate=apply(x.test,1,function(z){exp(sum(z*betahat))/(1+exp(sum(z*betahat)))})
#rate[which(is.na(rate))]=rep(1,length(which(is.na(rate))))
y.hat=as.numeric(rate>0.5)
#
CRP=sum(y.test==y.hat)/nrow(x.test)
ReErr=sqrt(sum((betahat-betatrue)^2))/(sqrt(sum(betatrue^2)))









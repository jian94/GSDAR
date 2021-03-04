#-----------------------------------------------------------------
# Negative log likelihood function
F.function=function(belta,x,y){
  f=sum(log(1+exp(x%*%belta)))-t(y)%*%x%*%belta
  return(f)
}

# Derivate of the negative log likelihood function about beta
DF.function=function(belta,x,y){
  pw=numeric(length(y))
  for(i in 1:length(y)){
  pw[i]=1/(1+exp(-sum(x[i,]*belta))) 
  }
  fd=t(x)%*%(pw-y)
  return(fd)
}

#-----------------------------------------------------------------
SQN1=function(x,beta0,y){
maxoutiter = 20
#maxiniter = 20 
stoprule = 0.01
# maxoutiter = 300
# maxiniter = 30
# stoprule = 0.001
yita=0.4;rho=0.55; 
Hk=diag(ncol(x))  
outiter=0
repeat{
 gk=DF.function(beta0,x,y)    
 if(sqrt(sum(gk^2))< stoprule | outiter > maxoutiter) break   ##belta=beltao
 dk=-Hk%*%gk   
 initer=0; mk=0
 repeat{
  beltanew=beta0 + rho^initer*dk
  newPNLL=F.function(beltanew,x,y)
  oldPNLL=F.function(beta0,x,y)
  temp=yita*(rho^initer)*sum(gk*dk)
  mk=initer
  initer=initer+1
  if(newPNLL<(oldPNLL+temp)) break
  }
##======== DFP==========##
belta=beta0+(rho^mk)*dk
##--------- update Hk----------##
sk=belta-beta0     ##displacement shift
yk=DF.function(belta,x,y)-gk   ##gradient shift
if(sum(sk*yk)>0) Hk=Hk-(Hk%*%(yk%*%t(yk))%*%Hk)/as.numeric(t(yk)%*%Hk%*%yk)+(sk%*%t(sk))/as.numeric(t(sk)%*%yk)
outiter=outiter+1
beta0=belta
}
return(beta0)
}

#-----------------------------------------------------------------
GSDAR=function(beta0,T1,x,y,w,itermax){
#=================================================================
#computating L0 regularization problem via GSDAR for given sparsity level
#------------------------------------------------------------------
# INPUT:
# 	x	design matrix      
# 	y   	response vector 
# 	algorithm parameters: 
# 		              beta0   	     intial value of beta (default: 0)
# 		              T1              sparsity level of true beta
#                          w              weight between primal and dual part 
# 		              itermax         maximum number of iterations  
# OUTPUT:
# 	betahat       estimator of beta
# 	sethat        estimated support set of betahat 
#-----------------------------------------------------------------
# Written by:
#           Yuling Jiao (yulingjiaomath@whu.edu.cn)                         
#           Lican Kang (kanglican@whu.edu.cn)                            
# This version: May 2, 2020                                             
#=================================================================
  A=list()
  I=list()
  d=list()
  belta1=list()
  belta1[[1]]=beta0
  d[[1]]=-DF.function(beta0,x,y)
  d.T0=sort(abs(w*belta1[[1]]+(1-w)*d[[1]]),decreasing =T)[T1]
  A[[1]]=which(abs(w*belta1[[1]]+(1-w)*d[[1]])>=d.T0)
  I[[1]]=which(abs(w*belta1[[1]]+(1-w)*d[[1]])<d.T0)
  b=A[[1]]
  bc=I[[1]]
  belta1[[2]]=numeric(ncol(x))
  belta1[[2]][b]=SQN1(x[,b],beta0[b],y)
  belta1[[2]][bc]=rep(0,length(bc))
  d[[2]]=numeric(ncol(x))
  d[[2]][b]=rep(0,length(b))
  X.I=x[,bc]
  pw=numeric(nrow(x))
  for(i in 1:nrow(x)){
  pw[i]=exp(sum((x[i,][b])*(belta1[[2]][b])))/(1+exp(sum((x[i,][b])*(belta1[[2]][b]))))
  }
  d[[2]][bc]=t(X.I)%*%(y-pw)
  d.T1=sort(abs(w*belta1[[2]]+(1-w)*d[[2]]),decreasing =T)[T1]
  A[[2]]=which(abs(w*belta1[[2]]+(1-w)*d[[2]])>=d.T1)
  I[[2]]=which(abs(w*belta1[[2]]+(1-w)*d[[2]])<d.T1)
  k=1
  while(setequal(A[[k]],A[[k+1]])==0 && k<itermax){
    k=k+1
    d.T0=sort(abs(w*belta1[[k]]+(1-w)*d[[k]]),decreasing =T)[T1]
    A[[k]]=which(abs(w*belta1[[k]]+(1-w)*d[[k]])>=d.T0)
    I[[k]]=which(abs(w*belta1[[k]]+(1-w)*d[[k]])<d.T0)
    b=A[[k]]
    bc=I[[k]]
    belta1[[k+1]]=numeric(ncol(x))
    belta1[[k+1]][b]=SQN1(x[,b],beta0[b],y)
    belta1[[k+1]][bc]=rep(0,length(bc))
    d[[k+1]]=numeric(ncol(x))
    d[[k+1]][b]=rep(0,length(b))
    X.I=x[,bc]
    pw=numeric(nrow(x))
    for(i in 1:nrow(x)){
    pw[i]=exp(sum((x[i,][b])*(belta1[[k+1]][b])))/(1+exp(sum((x[i,][b])*(belta1[[k+1]][b]))))
    }
    d[[k+1]][bc]=t(X.I)%*%(y-pw)
    d.T1=sort(abs(w*belta1[[k+1]]+(1-w)*d[[k+1]]),decreasing =T)[T1]
    A[[k+1]]=which(abs(w*belta1[[k+1]]+(1-w)*d[[k+1]])>=d.T1)
    I[[k+1]]=which(abs(w*belta1[[k+1]]+(1-w)*d[[k+1]])<d.T1)
  }
  betahat=belta1[[length(belta1)]]
  sethat=which(betahat!=0)
  return(list(betahat,sethat))
}







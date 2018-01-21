##############ornstein-uhlenbeck simulation######################
ornstein_uhlenbeck <- function(mu=.1,theta=1,sigma=1,from=0,to=1,n=len,...){
  
  dt  <- (to-from)/(n-1)
  dw  <- rnorm(n+20, 0, sqrt(dt))
  x <- from
  for (i in 2:(n+20)) {
    x[i]  <-  x[i-1] + theta*(mu-x[i-1])*dt + sigma*dw[i-1]*sqrt(dt)
  }
  return(x[-c(1:20)])
}

##Global false discovery rate
## set the tolerance for integration
rel.tol=1e-8
##set the tolerance for uniroot
tol=1e-6
options(digits = 7)

##calculate T-statistics 
Tstat<-function(x,y){
  x=x[is.finite(x)];y=y[is.finite(y)]
  n1=length(x);n2=length(y)
  tvalue=(mean(x)-mean(y))/sqrt(((n1-1)*var(x)+(n2-1)*var(y))/(n1+n2-2)*(1/n1+1/n2)) 
  return(list(qt=tvalue,pt=pt(tvalue,df=n1+n2-2)))
}

##Local false discovery rate
# h<-function(x,tau,parm) (new.tau[1]*new.f0(x,tau,parm))/new.fw(x,tau,parm)
h<-function(x){
  # A11=integrate(f11,0,1,parm)$value  
  # new.tau=c(1-tau[2]*A11,tau[2]*A11)
  return((new.tau[1]*new.f0(x))/new.fw(x))
}

##use untilted distribution to find the symetric piont t to x
froot1 <- function(t,x) Fw(t)+Fw(x)-1
##use tilted distribution to find the symetric piont t to x
froot2 <- function(t,x) {
  return(ifelse(t!=0,integrate(kernel_fw,0,t,rel.tol = rel.tol,stop.on.error = FALSE)$value,0)-
           ifelse(x!=1,integrate(kernel_fw,x,1,rel.tol = rel.tol,stop.on.error = FALSE)$value,0))
}

##EM algorithm for two-point mixture empirical null
fw <- function(x) tau[1]+tau[2]*dbeta(x,parm[1],parm[2])
Fw <- function(x) tau[1]*x+tau[2]*pbeta(x,parm[1],parm[2])

#============================================#
# Method of Moment to estimate the mixture of
# uniform and beta distributions     
#============================================#
MOM<-function(x,start=c(0.9,1,2),...){
  parms=c(mean(x),mean(x^2),mean(x^3))
  myfun<-function(x,parms){
    F1=1/2*x[1]+x[2]/sum(x[2:3])*(1-x[1])-parms[1]
    F2=1/3*x[1]+(x[2]*(x[2]+1))/(sum(x[2:3])*(sum(x[2:3])+1))*(1-x[1])-parms[2]
    F3=1/4*x[1]+((x[2]+2)*(x[2]+1)*x[2])/((sum(x[2:3])+2)*(sum(x[2:3])+1)*sum(x[2:3]))*(1-x[1])-parms[3]
    c(F1,F2,F3)
  }
  res=nleqslv(x=start,myfun,parms=mu,method="Newton",
          control = list(btol=1e-8,allowSingular=TRUE))$x
  return(list(p=c(res[1],1-res[1]),parm=res[2:3]))
}

##########################################################################
##########  Using EM algorithm to fit the original model  ################  
empiri.null<-function(y,parm,tau,tol=1e-8,iterations=1000,...){
  # ldiff=NULL
  loglike=NULL
  old.tau=tau;old.parm=parm
  itr=1
  while(itr<iterations){
    fy=rbind(dunif(y),dbeta(y,old.parm[1],old.parm[2]))
    ##calculate new tau
    p=old.tau*fy;p.scale=apply(p,2,function(t) t/sum(t))
    lold=sum(log(colSums(p))) ##the total observed loglike
    new.tau=apply(p.scale,1,mean)  ##the updated weights
    
    ##calculate new parm
    logQ<-function(parm) { ##only account for the beta distribution
      sum(p.scale[2,]*dbeta(y,parm[1],parm[2],log = TRUE)) 
    }
    new.parm=tryCatch(maxNR(logQ,start=old.parm)$estimate,
                      error=function(e) old.parm) ## updated parm for beta distribution
    
    ##iteration diagnose 
    lnew=sum(log(colSums(new.tau*rbind(dunif(y),dbeta(y,new.parm[1],new.parm[2])))))
    loglike=c(loglike,lold)
    # ldiff=c(ldiff,lold-lnew)
    if(abs(lold-lnew)<tol){
      break
    }else{
      old.parm=new.parm
      old.tau=new.tau
      itr=itr+1  
    }
  }
  out=list(tau=new.tau,parm=new.parm,itr=itr,loglike=loglike)
  return(out)
}

##EM algorithm for three-point mixture empirircal null
empiri.null1<-function(y,parm1,parm2,tau){
  # y=simdata;parm1=c(5,1);parm2=c(1,5);tau=c(.7,.15,.15)
  ldiff=NULL
  old.tau=tau;old.parm1=parm1;old.parm2=parm2
  itr=1
  while(itr<6000){
    # print(paste("the",itr," iteration is running!"))
    fy=rbind(dunif(y),dbeta(y,old.parm1[1],old.parm1[2]),dbeta(y,old.parm2[1],old.parm2[2]))
    ##calculate new tau
    p=old.tau*fy;p.scale=apply(p,2,function(t) t/sum(t))
    lold=sum(log(colSums(p))) ##the total observed loglike
    new.tau=apply(p.scale,1,mean)  ##the updated weights
    
    ##optimize new parm1
    logQ1<-function(parm) { ##only account for the beta distribution
      sum(p.scale[2,]*log(dbeta(y,parm[1],parm[2]))) 
    }
    
    ##optimize new parm2
    logQ2<-function(parm) { ##only account for the beta distribution
      sum(p.scale[3,]*log(dbeta(y,parm[1],parm[2]))) 
    }
    
    new.parm1=maxNR(logQ1,start=old.parm1)$estimate ## updated parm for beta distribution
    new.parm2=maxNR(logQ2,start=old.parm2)$estimate ## updated parm for beta distribution
    
    ##iteration diagnose 
    lnew=sum(log(colSums(new.tau*rbind(dunif(y),dbeta(y,new.parm1[1],new.parm1[2]),
                                       dbeta(y,new.parm2[1],new.parm2[2])))))
    ldiff=c(ldiff,lold-lnew)
    if(abs(lold-lnew)<1e-8){
      break
    }else{
      old.parm1=new.parm1
      old.parm2=new.parm2
      old.tau=new.tau
      itr=itr+1  
    }
  }
  out=list(tau=new.tau,parm1=new.parm1,parm2=new.parm2,itr=itr,logdiff=ldiff)
  return(out)
}

# simdata=sapply(factor(sample(0:2,3000,replace = TRUE,prob=c(0.8,.1,.1)),
#                       labels = c("f0","f1","f2")),
#                function(x) switch(as.character(x),f0=runif(1),f1=rbeta(1,10,0.5),f2=rbeta(1,0.5,10)))
# emnull1=empiri.null1(simdata,c(10,1),c(1,10),c(.9,.05,.05))

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##
## ## ## ##   adjusted two-point Uniform-Beta mixture model    ## ## ## ##
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##

##new empiri.null model
f11<-function(x) ifelse(dbeta(x,parm[1],parm[2])<1,0,dbeta(x,parm[1],parm[2])-1)
f10<-function(x) ifelse(dbeta(x,parm[1],parm[2])<1,dbeta(x,parm[1],parm[2]),1)

new.f1<-function(x){
  # A11=integrate(f11,0,1,parm,rel.tol = rel.tol)$value   ##provide as a global constant in each CV
  return(f11(x)/A11) 
}

new.f0<-function(x) {
  # A11=integrate(f11,0,1,parm,rel.tol = rel.tol)$value
  # new.tau=c(1-tau[2]*A11,tau[2]*A11)
  return(tau[1]/new.tau[1]*dunif(x)+tau[2]/new.tau[1]*f10(x))
}

new.fw<-function(x) {
  # A11=integrate(f11,0,1,parm,rel.tol = rel.tol)$value
  # new.tau=c(1-tau[2]*A11,tau[2]*A11)
  return(new.tau[1]*new.f0(x)+new.tau[2]*new.f1(x))
} 


## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##
## ## ##   Tilted two-point Uniform-Beta mixture   ## ## ##
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##

##calculate the optimal theta such that E[h(x)]=T where T is
##average of global fdr using the verification data
##m(theta)=numerator(the)/denomenator(the)
TiltGmean<-function(theta){
  numerator <- function(x,the) h(x)*exp(the*h(x))*fw(x)
  denomenator <- function(x,the) exp(the*h(x))*fw(x)
  
  return(sapply(theta,function(t) 
    log(integrate(numerator,0,1,t,rel.tol = 1e-5,stop.on.error = FALSE)$value)-
      log(integrate(denomenator,0,1,t,rel.tol=1e-5,stop.on.error=FALSE)$value)-log(Tverf)))
}

##tilted empirical null distribution(may not need them)
##calculate the tilted tau 
##calculate E[exp(theta*h(x))] with respect to the empirical null
ctheta<-function(theta){
  kernel_fw <- function(x,the) exp(the*h(x))*fw(x)
  return(sapply(theta,function(the)
    integrate(kernel_fw,0,1,the,rel.tol = 1e-6,stop.on.error = FALSE)$value))
}

##not set the rel.tol too low; otherwise the error of non-finit function values might arise
tilt.tau<-function(theta){
  temp <- function(x) exp(theta*h(x))*f11(x)
  kernel_fw <- function(x) exp(theta*h(x))*fw(x)
  new.tau2=tau[2]*integrate(temp,0,1,rel.tol = 1e-5,stop.on.error = FALSE)$value/
    integrate(kernel_fw,0,1,rel.tol = 1e-5,stop.on.error = FALSE)$value
  return(c(1-new.tau2,new.tau2))
}

##global environment functions/values
q=seq(0.01,0.99,0.01)
foonames1=c("h","f11","f10","new.f0","new.f1","new.fw","fw","Fw","froot1","froot2",
            "A11","new.tau","tol","rel.tol","tau","parm","q")
foonames=c("h","f11","f10","new.f0","new.f1","new.fw","fw","Fw","froot1","froot2",
           "A11","new.tau","tol","rel.tol","tau","parm","tilted.tau")

##ROC1 is the untilted one
# x=quantv[1:4];theta=theta.opt;tau=em.tau;parm=em.parm
ROC1<-function(x){
  ##false discovery rate
  # cl <- makeCluster(4)
  # registerDoParallel(cl)
  res=foreach(i=1:length(x),.combine=rbind,.multicombine=T,.export=foonames1) %dopar%{
    tryCatch({
      ##use the above froot1 function to get the reflex point to the observed ptvalue
      sym.point=tryCatch(uniroot(froot1,x[i],interval=c(0,1),tol=rel.tol)$root,
                         error=function(e) 1-x[i])
      
      left=min(x[i],sym.point)  ##denote the left point wrt "pvalue"
      right=max(x[i],sym.point) ##denote the right point wrt "pvalue"
      
      ##probability of rejection based on untilt.f0 and untilt.fw
      if(left==0 && right==1){
        F0_x=new.f0(x[i])
        Fw_x=new.fw(x[i])
      }else if(left==0 && right!=1){
        F0_x=2*integrate(new.f0,right,1,stop.on.error = FALSE)$value
        Fw_x=2*integrate(new.fw,right,1,stop.on.error = FALSE)$value
      }else if(left!=0 && right==1){
        F0_x=2*integrate(new.f0,0,left,stop.on.error = FALSE)$value
        Fw_x=2*integrate(new.fw,0,left,stop.on.error = FALSE)$value
      }else {
        F0_x=1-integrate(new.f0,left,right,stop.on.error = FALSE)$value
        Fw_x=1-integrate(new.fw,left,right,stop.on.error = FALSE)$value
      }
      
      ##probability of non-rejection based on new.f0 and new.fw(fw)
      F0c_x = 1 - F0_x 
      Fwc_x = 1 - Fw_x
      
      ##calculate type I and type II errors
      if(left!=right){
        Fa_x=integrate(new.f1,left,right,stop.on.error = FALSE)$value
      }else{
        Fa_x=0
      }
      
      c(FDR=new.tau[1]*F0_x/Fw_x,NPV=new.tau[1]*F0c_x/Fwc_x,TypeI=F0_x,TypeII=Fa_x)
    },
    error=function(e) NA
    )
  }
  # stopCluster(cl)
  rownames(res)=NULL
  res
}

##ROC2 is for the tilted empirical null
##better not set rel.tol in ROC2

ROC2<-function(x,...){
  kernel_fw <- function(x) exp(theta*h(x))*fw(x)
  kernel_f0 <- function(x) exp(theta*h(x))*new.f0(x)
  kernel_f1 <- function(x) exp(theta*h(x))*new.f1(x)

  const0=integrate(kernel_f0,0,1,stop.on.error = FALSE)$value
  const1=integrate(kernel_f1,0,1,stop.on.error = FALSE)$value
  constw=integrate(kernel_fw,0,1,stop.on.error = FALSE)$value
  
  ##false discovery rate
  res=foreach(i=1:length(x),.combine=rbind,.multicombine=T,.export=foonames1)%dopar%{
    ##use the above froot2 function to get the reflex point to the observed ptvalue
    tryCatch({
      sym.point=tryCatch(uniroot(froot2,x[i],interval=c(0,1),tol=rel.tol)$root,
                         error=function(e) 1-x[i])
      left=min(x[i],sym.point)  ##denote the left point wrt "pvalue"
      right=max(x[i],sym.point) ##denote the right point wrt "pvalue"
      
      ##probability of rejection based on tilt.f0 and tilt.fw
      if(left==0 && right==1){
        F0_x=kernel_f0(x[i])/const0
        Fw_x=kernel_fw(x[i])/constw
      }else if(left==0 && right!=1){
        F0_x=2*integrate(kernel_f0,right,1,stop.on.error = FALSE)$value/const0
        Fw_x=2*integrate(kernel_fw,right,1,stop.on.error = FALSE)$value/constw
      }else if(left!=0 && right==1){
        F0_x=2*integrate(kernel_f0,0,left,stop.on.error = FALSE)$value/const0
        Fw_x=2*integrate(kernel_fw,0,left,stop.on.error = FALSE)$value/constw
      }else {
        F0_x=1-integrate(kernel_f0,left,right,stop.on.error = FALSE)$value/const0
        Fw_x=1-integrate(kernel_fw,left,right,stop.on.error = FALSE)$value/constw
      }
      
      ##probability of non-rejection based on tilt.f0 and tilt.fw
      F0c_x = 1 - F0_x
      Fwc_x = 1 - Fw_x
      
      ##calculate type I and type II errors
      if(left!=right){
        Fa_x=integrate(kernel_f1,left,right,stop.on.error = FALSE)$value/const1
      }else{
        Fa_x=0
      }
      
      c(FDR=tilted.tau[1]*F0_x/Fw_x,NPV=tilted.tau[1]*F0c_x/Fwc_x,
        TypeI=F0_x,TypeII=Fa_x)
    },
    error=function(e) NA
    )
  }
  rownames(res)=NULL
  res
}

##ROC4 is for the tilted empirical null
##use local false discovery rate 
ROC4<-function(x){
  # x=quantv[1:4];theta=theta.opt;tau=em.tau;parm=em.parm
  kernel_fw <- function(x) exp(theta*h(x))*fw(x)
  kernel_f0 <- function(x) exp(theta*h(x))*new.f0(x)
  kernel_f1 <- function(x) exp(theta*h(x))*new.f1(x)
  
  const0=integrate(kernel_f0,0,1,stop.on.error = FALSE)$value
  const1=integrate(kernel_f1,0,1,stop.on.error = FALSE)$value
  constw=integrate(kernel_fw,0,1,stop.on.error = FALSE)$value
  
  ##false discovery rate
  # cl <- makeCluster(4)
  # registerDoParallel(cl)
  res<-foreach(i=1:length(x),.combine=rbind,.multicombine=T,.export=foonames) %dopar%{
    ##use the above froot2 function to get the reflex point to the observed ptvalue
    tryCatch({
      c(pvalue=x[i],
        fdr1=new.tau[1]*new.f0(x[i])/fw(x[i]),
        fdr2=new.tau[1]*kernel_f0(x[i])/kernel_fw(x[i]))
    },
    error=function(e) NA
    )
  }
  # stopCluster(cl)
  res=data.frame(res)
  colnames(res)=c("pvalue","fdr","adj.fdr")
  rownames(res)=1:length(x)
  res
}


#===================================================================#
#                      power diagnosis                              #
#===================================================================#

##=============================================##
##new power function... -_- 
##Date: June 14, 2017   -_-
##=============================================##
##define a combine function
cfun<-function(a,b) mapply(cbind,a,b,SIMPLIFY = FALSE)
cfun1<-function(a,b) mapply(rbind,a,b,SIMPLIFY = FALSE)

##================================================##
##        use untilted f0,f1,and fw               ##
##================================================##
mypower1<-function(q){
  ##use fdr threshold
  fdr<-function(x,target=0){
    sapply(x,function(t){
      temp=(new.tau[1]*new.f0(t))/new.fw(t)
      ifelse(temp>1,1,ifelse(temp<0,0,temp))-target
    })
  }
  
  ##use Fdr threshold
  Fdr=function(x,target=0){
    sapply(x,function(t){
      tryCatch({
        sym.point=tryCatch(uniroot(froot1,t,interval=c(0,1),tol=tol)$root,
                           error=function(e) 1-t)
        left=min(t,sym.point)  ##denote the left point wrt "pvalue"
        right=max(t,sym.point) ##denote the right point wrt "pvalue"
        
        ##probability of rejection based on untilt.f0 and untilt.fw
        if(left==0 && right==1){
          F0_x=new.f0(t)
          Fw_x=new.fw(t)
        }else if(left==0 && right!=1){
          F0_x=2*integrate(new.f0,right,1,stop.on.error = FALSE)$value
          Fw_x=2*integrate(new.fw,right,1,stop.on.error = FALSE)$value
        }else if(left!=0 && right==1){
          F0_x=2*integrate(new.f0,0,left,stop.on.error = FALSE)$value
          Fw_x=2*integrate(new.fw,0,left,stop.on.error = FALSE)$value
        }else {
          F0_x=integrate(new.f0,0,left,stop.on.error = FALSE)$value+
            integrate(new.f0,right,1,stop.on.error = FALSE)$value
          Fw_x=integrate(new.fw,0,left,stop.on.error = FALSE)$value+
            integrate(new.fw,right,1,stop.on.error = FALSE)$value
        }
        
        temp=new.tau[1]*F0_x/Fw_x
        
        ifelse(temp>1,1,ifelse(temp<0,0,temp))-target
      },
      error=function(e) NA
      )
    })
  }
  
  ##find endpoints
  res=foreach(i=1:length(q),.combine = cfun1,.packages = c("rootSolve"),
              .export = foonames1) %dopar%{
                
                ##use fdr as threshold
                fdr_endpoints=tryCatch({uniroot.all(fdr,target=q[i],
                                                    interval = c(1e-7,1-1e-7),
                                                    tol = 1e-6)},
                                       error=function(e) NA)
                
                if(length(fdr_endpoints)==0){
                  fdr_power.h1=fdr_power.h0=NA
                  
                }else if(length(fdr_endpoints)==1){
                  ##power related to new.f1
                  fdr_power.h1=ifelse(fdr_endpoints<0.5,
                                integrate(new.f1,0,fdr_endpoints,rel.tol=1e-5,
                                          stop.on.error = FALSE)$value,
                                integrate(new.f1,fdr_endpoints,1,rel.tol=1e-5,
                                          stop.on.error = FALSE)$value)
                  ##power related to new.f0
                  fdr_power.h0=ifelse(fdr_endpoints<0.5,
                                  integrate(new.f0,0,fdr_endpoints,rel.tol=1e-5,
                                            stop.on.error = FALSE)$value,
                                  integrate(new.f0,fdr_endpoints,1,rel.tol=1e-5,
                                            stop.on.error = FALSE)$value)
                  
                }else{
                  ##power related to new.f1
                  fdr_power.h1=integrate(new.f1,0,min(fdr_endpoints),rel.tol=1e-6,
                                         stop.on.error = FALSE)$value+
                    integrate(new.f1,max(fdr_endpoints),1,rel.tol=1e-6,
                              stop.on.error = FALSE)$value
                    
                  ##power related to new.f0
                  fdr_power.h0=integrate(new.f0,0,min(fdr_endpoints),rel.tol=1e-6,
                                         stop.on.error = FALSE)$value+
                    integrate(new.f0,max(fdr_endpoints),1,rel.tol=1e-6,
                              stop.on.error = FALSE)$value
                }
                
                EPower_fdr=fdr_power.h1
                CPower_fdr=new.tau[2]*fdr_power.h1/
                  (new.tau[1]*fdr_power.h0+new.tau[2]*fdr_power.h1)
                
                ##unconditional error I/II
                ETypeI_fdr=fdr_power.h0; ETypeII_fdr=1-fdr_power.h1
                
                ##conditional error I/II
                CTypeI_fdr=1-CPower_fdr
                CTypeII_fdr=new.tau[2]*ETypeII_fdr/
                  (new.tau[1]*(1-ETypeI_fdr)+new.tau[2]*ETypeII_fdr)
                
                
                ##use Fdr as threshold
                if(q[i]>=new.tau[1]){
                  Fdr_power.h0=Fdr_power.h1=1
                }else{
                  Fdr_endpoints=tryCatch({uniroot.all(Fdr,target=q[i],
                                                      interval = c(1e-7,1-1e-7),
                                                      tol = 1e-6)}, 
                                         error=function(e) NA)
                  Fdr_endpoints
                  
                  if(length(Fdr_endpoints)==0){
                    Fdr_power.h1=Fdr_power.h0=NA
                  }else if(length(Fdr_endpoints)==1){
                    Fdr_power.h1=ifelse(Fdr_endpoints<0.5,
                                  integrate(new.f1,0,Fdr_endpoints,rel.tol=1e-5,
                                            stop.on.error = FALSE)$value,
                                  integrate(new.f1,Fdr_endpoints,1,rel.tol=1e-5,
                                            stop.on.error = FALSE)$value)
                    
                    Fdr_power.h0=ifelse(Fdr_endpoints<0.5,
                                        integrate(new.f0,0,Fdr_endpoints,rel.tol=1e-5,
                                                  stop.on.error = FALSE)$value,
                                        integrate(new.f0,Fdr_endpoints,1,rel.tol=1e-5,
                                                  stop.on.error = FALSE)$value)
                  }else{
                    Fdr_power.h1=integrate(new.f1,0,min(Fdr_endpoints),rel.tol=1e-6,
                                           stop.on.error = FALSE)$value+
                      integrate(new.f1,max(Fdr_endpoints),1,rel.tol=1e-6,
                                stop.on.error = FALSE)$value
                    
                    Fdr_power.h0=integrate(new.f0,0,min(Fdr_endpoints),rel.tol=1e-6,
                                           stop.on.error = FALSE)$value+
                      integrate(new.f0,max(Fdr_endpoints),1,rel.tol=1e-6,
                                stop.on.error = FALSE)$value
                  }
                }
                EPower_Fdr=Fdr_power.h1
                CPower_Fdr=new.tau[2]*Fdr_power.h1/
                  (new.tau[1]*Fdr_power.h0+new.tau[2]*Fdr_power.h1)
                
                ##unconditional error I/II
                ETypeI_Fdr=Fdr_power.h0; ETypeII_Fdr=1-Fdr_power.h1
                
                ##conditional error I/II
                CTypeI_Fdr=1-CPower_Fdr
                CTypeII_Fdr=new.tau[2]*ETypeII_Fdr/
                  (new.tau[1]*(1-ETypeI_Fdr)+new.tau[2]*ETypeII_Fdr)
                
                list(c(EPower_fdr,CPower_fdr,ETypeI_fdr,ETypeII_fdr,CTypeI_fdr,CTypeII_fdr),
                     c(EPower_Fdr,CPower_Fdr,ETypeI_Fdr,ETypeII_Fdr,CTypeI_Fdr,CTypeII_Fdr))
              }
  # rownames(res)=NULL
  res
}


##================================================##
##          use tilted f0,f1 and fw               ##
##================================================##
mypower2<-function(q){
  kernel_fw <- function(x) exp(theta*h(x))*fw(x)
  kernel_f0 <- function(x) exp(theta*h(x))*new.f0(x)
  kernel_f1 <- function(x) exp(theta*h(x))*new.f1(x)
  
  const0=integrate(kernel_f0,0,1,stop.on.error = FALSE)$value
  const1=integrate(kernel_f1,0,1,stop.on.error = FALSE)$value
  constw=integrate(kernel_fw,0,1,stop.on.error = FALSE)$value
  
  ##use fdr threshold
  fdr<-function(x,target=0) {
    sapply(x,function(t){
      temp=(tilted.tau[1]*kernel_f0(t)/const0)/(kernel_fw(t)/constw)
      ifelse(temp>1,1,ifelse(temp<0,0,temp))-target
    })
  }
  
  ##use Fdr threshold
  Fdr=function(x,target=0){
    sapply(x,function(t){
      tryCatch({
        sym.point=tryCatch(uniroot(froot2,t,interval=c(0,1),tol=tol)$root,
                           error=function(e) 1-t)
        left=min(t,sym.point)  ##denote the left point wrt "pvalue"
        right=max(t,sym.point) ##denote the right point wrt "pvalue"
        
        ##probability of rejection based on tilt.f0 and tilt.fw
        if(left==0 && right==1){
          F0_x=kernel_f0(t)
          Fw_x=kernel_fw(t)
        }else if(left==0 && right!=1){
          F0_x=2*integrate(kernel_f0,right,1,stop.on.error = FALSE)$value
          Fw_x=2*integrate(kernel_fw,right,1,stop.on.error = FALSE)$value
        }else if(left!=0 && right==1){
          F0_x=2*integrate(kernel_f0,0,left,stop.on.error = FALSE)$value
          Fw_x=2*integrate(kernel_fw,0,left,stop.on.error = FALSE)$value
        }else {
          F0_x=integrate(kernel_f0,0,left,stop.on.error = FALSE)$value+
            integrate(kernel_f0,right,1,stop.on.error = FALSE)$value
          Fw_x=integrate(kernel_fw,0,left,stop.on.error = FALSE)$value+
            integrate(kernel_fw,right,1,stop.on.error = FALSE)$value
        }
        
        temp=new.tau[1]*F0_x/Fw_x
        
        ifelse(temp>1,1,ifelse(temp<0,0,temp))-target
      },
      error=function(e) NULL
      )
    })
  }
  
  ##find endpoints
  res=foreach(i=1:length(q),.combine = cfun1,.packages = c("rootSolve"),
              .export = foonames) %dopar%{
                
                ##use fdr as threshold
                fdr_endpoints=tryCatch({uniroot.all(fdr,target=q[i],
                                                    interval = c(1e-7,1-1e-7),
                                                    tol = 1e-6)}, 
                                       error=function(e) NA)
                if(length(fdr_endpoints)==0){
                  fdr_power.h1=fdr_power.h0=NA
                }else if(length(fdr_endpoints)==1){
                  ##power related to new.f1
                  fdr_power.h1=ifelse(fdr_endpoints<0.5,
                                      integrate(kernel_f1,0,fdr_endpoints,rel.tol=1e-5,
                                                stop.on.error = FALSE)$value/const1,
                                      integrate(kernel_f1,fdr_endpoints,1,rel.tol=1e-5,
                                                stop.on.error = FALSE)$value/const1)
                  ##power related to new.f0
                  fdr_power.h0=ifelse(fdr_endpoints<0.5,
                                      integrate(kernel_f0,0,fdr_endpoints,rel.tol=1e-5,
                                                stop.on.error = FALSE)$value/const0,
                                      integrate(kernel_f0,fdr_endpoints,1,rel.tol=1e-5,
                                                stop.on.error = FALSE)$value/const0)
                  
                }else{
                  
                  ##power related to new.f1
                  fdr_power.h1=integrate(kernel_f1,0,min(fdr_endpoints),
                                           rel.tol=1e-6,stop.on.error = FALSE)$value/const1+
                    integrate(kernel_f1,max(fdr_endpoints),1,
                              rel.tol=1e-6,stop.on.error = FALSE)$value/const1
                  
                  ##power related to new.f0
                  fdr_power.h0=integrate(kernel_f0,0,min(fdr_endpoints),
                                         rel.tol=1e-6,stop.on.error = FALSE)$value/const0+
                    integrate(kernel_f0,max(fdr_endpoints),1,
                              rel.tol=1e-6,stop.on.error = FALSE)$value/const0
                }
                
                EPower_fdr=fdr_power.h1
                CPower_fdr=tilted.tau[2]*fdr_power.h1/
                  (tilted.tau[1]*fdr_power.h0+tilted.tau[2]*fdr_power.h1)
                
                ##unconditional error I/II
                ETypeI_fdr=fdr_power.h0; ETypeII_fdr=1-fdr_power.h1
                
                ##conditional error I/II
                CTypeI_fdr=1-CPower_fdr
                CTypeII_fdr=tilted.tau[2]*ETypeII_fdr/
                  (tilted.tau[1]*(1-ETypeI_fdr)+tilted.tau[2]*ETypeII_fdr)
                
                ##use Fdr as threshold
                if(q[i]>=tilted.tau[1]){
                  Fdr_power.h0=Fdr_power.h1=1
                }else{
                  Fdr_endpoints=tryCatch({uniroot.all(Fdr,target=q[i],
                                                      interval = c(1e-7,1-1e-7),
                                                      tol = 1e-6)}, 
                                         error=function(e) NA)
                  
                  if(length(Fdr_endpoints)==0){
                    Fdr_power.h1=Fdr_power.h0=NA
                  }else if(length(Fdr_endpoints)==1){
                    Fdr_power.h1=ifelse(Fdr_endpoints<0.5,
                                        integrate(kernel_f1,0,Fdr_endpoints,rel.tol=1e-5,
                                                  stop.on.error = FALSE)$value/const1,
                                        integrate(kernel_f1,Fdr_endpoints,1,rel.tol=1e-5,
                                                  stop.on.error = FALSE)$value)/const1
                    
                    Fdr_power.h0=ifelse(Fdr_endpoints<0.5,
                                        integrate(kernel_f0,0,Fdr_endpoints,rel.tol=1e-5,
                                                  stop.on.error = FALSE)$value/const0,
                                        integrate(kernel_f0,Fdr_endpoints,1,rel.tol=1e-5,
                                                  stop.on.error = FALSE)$value/const0)
                  }else{
                    ##power related to new.f1
                    Fdr_power.h1=integrate(kernel_f1,0,min(Fdr_endpoints),
                                           rel.tol=1e-6,stop.on.error = FALSE)$value/const1+
                      integrate(kernel_f1,max(Fdr_endpoints),1,
                                rel.tol=1e-6,stop.on.error = FALSE)$value/const1
                    
                    ##power related to new.f0
                    Fdr_power.h0=integrate(kernel_f0,0,min(Fdr_endpoints),
                                           rel.tol=1e-6,stop.on.error = FALSE)$value/const0+
                      integrate(kernel_f0,max(Fdr_endpoints),1,
                                rel.tol=1e-6,stop.on.error = FALSE)$value/const0
                  }
                }
                EPower_Fdr=Fdr_power.h1
                CPower_Fdr=tilted.tau[2]*Fdr_power.h1/
                  (tilted.tau[1]*Fdr_power.h0+tilted.tau[2]*Fdr_power.h1)
                
                ##unconditional error I/II
                ETypeI_Fdr=Fdr_power.h0; ETypeII_Fdr=1-Fdr_power.h1
                
                ##conditional error I/II
                CTypeI_Fdr=1-CPower_Fdr
                CTypeII_Fdr=tilted.tau[2]*ETypeII_Fdr/
                  (tilted.tau[1]*(1-ETypeI_Fdr)+tilted.tau[2]*ETypeII_Fdr)
                
                list(c(EPower_fdr,CPower_fdr,ETypeI_fdr,ETypeII_fdr,CTypeI_fdr,CTypeII_fdr),
                     c(EPower_Fdr,CPower_Fdr,ETypeI_Fdr,ETypeII_Fdr,CTypeI_Fdr,CTypeII_Fdr))
              }
  # rownames(res)=NULL
  res
}


##=================================================##
##get combined endpoints use untilted f0,f1 and fw ##
##=================================================##
my_endpoints1<-function(q){
  
  ##use fdr threshold
  fdr<-function(x,target=0){
    sapply(x,function(t){
      temp=(new.tau[1]*new.f0(t))/new.fw(t)
      ifelse(temp>1,1,ifelse(temp<0,0,temp))-target
    })
  }
  
  ##use Fdr threshold
  Fdr=function(x,target=0){
    sapply(x,function(t){
      tryCatch({
        sym.point=tryCatch(uniroot(froot1,t,interval=c(0,1),tol=tol)$root,
                           error=function(e) 1-t)
        left=min(t,sym.point)  ##denote the left point wrt "pvalue"
        right=max(t,sym.point) ##denote the right point wrt "pvalue"
        
        ##probability of rejection based on untilt.f0 and untilt.fw
        if(left==0 && right==1){
          F0_x=new.f0(t)
          Fw_x=new.fw(t)
        }else if(left==0 && right!=1){
          F0_x=2*integrate(new.f0,right,1,stop.on.error = FALSE)$value
          Fw_x=2*integrate(new.fw,right,1,stop.on.error = FALSE)$value
        }else if(left!=0 && right==1){
          F0_x=2*integrate(new.f0,0,left,stop.on.error = FALSE)$value
          Fw_x=2*integrate(new.fw,0,left,stop.on.error = FALSE)$value
        }else {
          F0_x=integrate(new.f0,0,left,stop.on.error = FALSE)$value+
            integrate(new.f0,right,1,stop.on.error = FALSE)$value
          Fw_x=integrate(new.fw,0,left,stop.on.error = FALSE)$value+
            integrate(new.fw,right,1,stop.on.error = FALSE)$value
        }
        
        temp=new.tau[1]*F0_x/Fw_x
        
        ifelse(temp>1,1,ifelse(temp<0,0,temp))-target
      },
      error=function(e) NA
      )
    })
  }
  
  ##find endpoints
  res=foreach(i=1:length(q),.combine = cfun1,.packages = c("rootSolve"),
              .export = foonames1) %dopar%{
                
                fdr_endpoints=Fdr_endpoints=NULL
                ##use fdr as threshold
                ##tol=1e-6;tol1=1e-7;tol2=1e-8;tol3=1e-9
                tol1=1e-7
                while((length(fdr_endpoints)!=2 || any(is.na(fdr_endpoints))) && 
                      tol1 >= 1e-10){
                  fdr_endpoints=tryCatch({uniroot.all(fdr,target=q[i],
                                                      interval = c(tol1,1-tol1),
                                                      tol = 1e-6)},
                                         error=function(err) c(NA,NA)
                                         )
                  tol1=tol1*0.1
                }
                
                ##use Fdr as threshold
                tol1=1e-6
                while((length(Fdr_endpoints)!=2 || any(is.na(Fdr_endpoints))) &&
                      tol1 >= 1e-10){
                  Fdr_endpoints=tryCatch({uniroot.all(Fdr,target=q[i],
                                                      interval = c(tol1,1-tol1),
                                                      tol = 1e-6)},
                                         error=function(e) return(c(NA,NA)))
                  tol1=tol1*0.1
                }
                
                ##error handling 
                if(length(fdr_endpoints)==2 && all(!is.na(fdr_endpoints))){
                  fdr_endpoints=sort(fdr_endpoints)
                }else if(length(fdr_endpoints)==1){
                  fdr_endpoints=sort(c(fdr_endpoints,1-fdr_endpoints))
                }else {
                  fdr_endpoints=c(NA,NA)
                }
                
                if(length(Fdr_endpoints)==2 && all(!is.na(Fdr_endpoints))){
                  Fdr_endpoints=sort(Fdr_endpoints)
                }else if(length(Fdr_endpoints)==1){
                  Fdr_endpoints=sort(c(Fdr_endpoints,1-Fdr_endpoints))
                }else {
                  Fdr_endpoints=c(NA,NA)
                }
                list(fdr_endpoints,Fdr_endpoints)
              }
  # rownames(res)=NULL
  res
}

##new power function use endpoints; use untilted f0,f1,fw
newPower1<-function(EndPoints,ncore=4){
  ## endpoints is a list of two matrix;
  ## one is for using fdr and the other for Fdr
  cl<-makeCluster(ncore)
  registerDoParallel(cl)
  ##find endpoints
  res=foreach(i=1:nrow(EndPoints[[1]]),.combine = cfun1,
              .packages = c("rootSolve"),
              .export = foonames1) %dopar%{
                
                ##use fdr as threshold
                fdr_endpoints=EndPoints[[1]][i,]
                if(sum(is.na(fdr_endpoints))==2){
                  fdr_power.h1=fdr_power.h0=NA
                }else if(sum(is.na(fdr_endpoints))==1){
                  fdr_endpoints=na.omit(fdr_endpoints)
                  ##power related to new.f1
                  fdr_power.h1=ifelse(fdr_endpoints<0.5,
                                      integrate(new.f1,0,fdr_endpoints,rel.tol=1e-5,
                                                stop.on.error = FALSE)$value,
                                      integrate(new.f1,fdr_endpoints,1,rel.tol=1e-5,
                                                stop.on.error = FALSE)$value)
                  ##power related to new.f0
                  fdr_power.h0=ifelse(fdr_endpoints<0.5,
                                      integrate(new.f0,0,fdr_endpoints,rel.tol=1e-5,
                                                stop.on.error = FALSE)$value,
                                      integrate(new.f0,fdr_endpoints,1,rel.tol=1e-5,
                                                stop.on.error = FALSE)$value)
                  
                }else{
                  ##power related to new.f1
                  fdr_power.h1=integrate(new.f1,0,min(fdr_endpoints),rel.tol=1e-6,
                                         stop.on.error = FALSE)$value+
                    integrate(new.f1,max(fdr_endpoints),1,rel.tol=1e-6,
                              stop.on.error = FALSE)$value
                  
                  ##power related to new.f0
                  fdr_power.h0=integrate(new.f0,0,min(fdr_endpoints),rel.tol=1e-6,
                                         stop.on.error = FALSE)$value+
                    integrate(new.f0,max(fdr_endpoints),1,rel.tol=1e-6,
                              stop.on.error = FALSE)$value
                }
                
                EPower_fdr=fdr_power.h1
                CPower_fdr=new.tau[2]*fdr_power.h1/
                  (new.tau[1]*fdr_power.h0+new.tau[2]*fdr_power.h1)
                
                ##unconditional error I/II
                ETypeI_fdr=fdr_power.h0; ETypeII_fdr=1-fdr_power.h1
                
                ##conditional error I/II
                CTypeI_fdr=1-CPower_fdr
                CTypeII_fdr=new.tau[2]*ETypeII_fdr/
                  (new.tau[1]*(1-ETypeI_fdr)+new.tau[2]*ETypeII_fdr)
                
                
                ##use Fdr as threshold
                if(q[i]>=new.tau[1]){
                  Fdr_power.h0=Fdr_power.h1=1
                }else{
                  Fdr_endpoints=EndPoints[[2]][i,]
                  if(sum(is.na(Fdr_endpoints))==2){
                    Fdr_power.h1=Fdr_power.h0=NA
                  }else if(sum(is.na(Fdr_endpoints))==1){
                    Fdr_endpoints=na.omit(Fdr_endpoints)
                    
                    Fdr_power.h1=ifelse(Fdr_endpoints<0.5,
                                        integrate(new.f1,0,Fdr_endpoints,rel.tol=1e-5,
                                                  stop.on.error = FALSE)$value,
                                        integrate(new.f1,Fdr_endpoints,1,rel.tol=1e-5,
                                                  stop.on.error = FALSE)$value)
                    
                    Fdr_power.h0=ifelse(Fdr_endpoints<0.5,
                                        integrate(new.f0,0,Fdr_endpoints,rel.tol=1e-5,
                                                  stop.on.error = FALSE)$value,
                                        integrate(new.f0,Fdr_endpoints,1,rel.tol=1e-5,
                                                  stop.on.error = FALSE)$value)
                  }else{
                    Fdr_power.h1=integrate(new.f1,0,min(Fdr_endpoints),rel.tol=1e-6,
                                           stop.on.error = FALSE)$value+
                      integrate(new.f1,max(Fdr_endpoints),1,rel.tol=1e-6,
                                stop.on.error = FALSE)$value
                    
                    Fdr_power.h0=integrate(new.f0,0,min(Fdr_endpoints),rel.tol=1e-6,
                                           stop.on.error = FALSE)$value+
                      integrate(new.f0,max(Fdr_endpoints),1,rel.tol=1e-6,
                                stop.on.error = FALSE)$value
                  }
                }
                EPower_Fdr=Fdr_power.h1
                CPower_Fdr=new.tau[2]*Fdr_power.h1/
                  (new.tau[1]*Fdr_power.h0+new.tau[2]*Fdr_power.h1)
                
                ##unconditional error I/II
                ETypeI_Fdr=Fdr_power.h0; ETypeII_Fdr=1-Fdr_power.h1
                
                ##conditional error I/II
                CTypeI_Fdr=1-CPower_Fdr
                CTypeII_Fdr=new.tau[2]*ETypeII_Fdr/
                  (new.tau[1]*(1-ETypeI_Fdr)+new.tau[2]*ETypeII_Fdr)
                
                list(c(EPower_fdr,CPower_fdr,ETypeI_fdr,ETypeII_fdr,CTypeI_fdr,CTypeII_fdr),
                     c(EPower_Fdr,CPower_Fdr,ETypeI_Fdr,ETypeII_Fdr,CTypeI_Fdr,CTypeII_Fdr))
              }
  stopCluster(cl)
  # rownames(res)=NULL
  res
}


##=================================================##
##  get combined endpoints use tilted f0,f1 and fw ##
##=================================================##
my_endpoints2<-function(q){
  kernel_fw <- function(x) exp(theta*h(x))*fw(x)
  kernel_f0 <- function(x) exp(theta*h(x))*new.f0(x)
  kernel_f1 <- function(x) exp(theta*h(x))*new.f1(x)
  
  const0=integrate(kernel_f0,0,1,stop.on.error = FALSE)$value
  const1=integrate(kernel_f1,0,1,stop.on.error = FALSE)$value
  constw=integrate(kernel_fw,0,1,stop.on.error = FALSE)$value
  
  ##use fdr threshold
  fdr<-function(x,target=0) {
    sapply(x,function(t){
      temp=(tilted.tau[1]*kernel_f0(t)/const0)/(kernel_fw(t)/constw)
      ifelse(temp>1,1,ifelse(temp<0,0,temp))-target
    })
  }
  
  ##use Fdr threshold
  Fdr=function(x,target=0){
    sapply(x,function(t){
      tryCatch({
        sym.point=tryCatch(uniroot(froot2,t,interval=c(0,1),tol=tol)$root,
                           error=function(e) 1-t)
        left=min(t,sym.point)  ##denote the left point wrt "pvalue"
        right=max(t,sym.point) ##denote the right point wrt "pvalue"
        
        ##probability of rejection based on tilt.f0 and tilt.fw
        if(left==0 && right==1){
          F0_x=kernel_f0(t)
          Fw_x=kernel_fw(t)
        }else if(left==0 && right!=1){
          F0_x=2*integrate(kernel_f0,right,1,stop.on.error = FALSE)$value
          Fw_x=2*integrate(kernel_fw,right,1,stop.on.error = FALSE)$value
        }else if(left!=0 && right==1){
          F0_x=2*integrate(kernel_f0,0,left,stop.on.error = FALSE)$value
          Fw_x=2*integrate(kernel_fw,0,left,stop.on.error = FALSE)$value
        }else {
          F0_x=1-integrate(kernel_f0,left,right,stop.on.error = FALSE)$value
          Fw_x=1-integrate(kernel_fw,left,right,stop.on.error = FALSE)$value
        }

        temp=new.tau[1]*F0_x/Fw_x
        
        ifelse(temp>1,1,ifelse(temp<0,0,temp))-target
      },
      error=function(e) NULL
      )
    })
  }
  
  ##find endpoints
  res=foreach(i=1:length(q),.combine = cfun1,.packages = c("rootSolve"),
              .export = foonames1) %dopar%{
                
                ##use fdr as threshold
                fdr_endpoints=tryCatch({uniroot.all(fdr,target=q[i],interval = c(1e-6,1-1e-6),
                                                    tol = 1e-6)},
                                       error=function(e) NA)
                ##use Fdr as threshold
                Fdr_endpoints=tryCatch({uniroot.all(Fdr,target=q[i],interval = c(1e-6,1-1e-6),
                                                    tol = 1e-6)},
                                       error=function(e) NA)
                list(fdr_endpoints,Fdr_endpoints)
              }
  # rownames(res)=NULL
  res
}

##new power function use endpoints; use tilted f0,f1,fw
newPower2<-function(endpoints){
  ## endpoints is a list of two matrix;
  ## one is for using fdr and the other for Fdr
  endpoints1=endpoints[[1]] ## fdr
  endpoints2=endpoints[[2]] ## Fdr
  n=dim(endpoints1)[1]
  
  ##find endpoints
  res=foreach(i=1:n,.combine = cfun1,.packages = c("rootSolve"),
              .export = foonames1) %dopar%{
                
                ##use fdr as threshold
                fdr_endpoints=endpoints1[i,]
                if(sum(is.na(fdr_endpoints))==2){
                  fdr_power.h1=fdr_power.h0=NA
                }else if(sum(is.na(fdr_endpoints))==1){
                  
                  fdr_endpoints=is.na(fdr_endpoints)
                  ##power related to new.f1
                  fdr_power.h1=ifelse(fdr_endpoints<0.5,
                                      integrate(kernel_f1,0,fdr_endpoints,rel.tol=1e-5,
                                                stop.on.error = FALSE)$value/const1,
                                      integrate(kernel_f1,fdr_endpoints,1,rel.tol=1e-5,
                                                stop.on.error = FALSE)$value/const1)
                  ##power related to new.f0
                  fdr_power.h0=ifelse(fdr_endpoints<0.5,
                                      integrate(kernel_f0,0,fdr_endpoints,rel.tol=1e-5,
                                                stop.on.error = FALSE)$value/const0,
                                      integrate(kernel_f0,fdr_endpoints,1,rel.tol=1e-5,
                                                stop.on.error = FALSE)$value/const0)
                  
                }else{
                  ##power related to new.f1
                  fdr_power.h1=1-integrate(kernel_f1,min(fdr_endpoints),max(fdr_endpoints),
                                           rel.tol=1e-6,stop.on.error = FALSE)$value/const1
                  ##power related to new.f0
                  fdr_power.h0=1-integrate(kernel_f0,min(fdr_endpoints),max(fdr_endpoints),
                                           rel.tol=1e-6,stop.on.error = FALSE)$value/const0
                }
                
                EPower_fdr=fdr_power.h1
                CPower_fdr=new.tau[2]*fdr_power.h1/
                  (new.tau[1]*fdr_power.h0+new.tau[2]*fdr_power.h1)
                
                ##unconditional error I/II
                ETypeI_fdr=fdr_power.h0; ETypeII_fdr=1-fdr_power.h1
                
                ##conditional error I/II
                CTypeI_fdr=1-CPower_fdr
                CTypeII_fdr=tilted.tau[2]*ETypeII_fdr/
                  (tilted.tau[1]*(1-ETypeI_fdr)+tilted.tau[2]*ETypeII_fdr)
                
                
                ##use Fdr as threshold
                if(q[i]>=new.tau[1]){
                  Fdr_power.h0=Fdr_power.h1=1
                }else{
                  Fdr_endpoints=endpoints2[i,]
                  
                  if(sum(is.na(Fdr_endpoints))==2){
                    Fdr_power.h1=Fdr_power.h0=NA
                  }else if(sum(is.na(Fdr_endpoints))==1){
                    Fdr_endpoints=is.na(Fdr_endpoints)
                    
                    Fdr_power.h1=ifelse(Fdr_endpoints<0.5,
                                        integrate(kernel_f1,0,Fdr_endpoints,rel.tol=1e-5,
                                                  stop.on.error = FALSE)$value/const1,
                                        integrate(kernel_f1,Fdr_endpoints,1,rel.tol=1e-5,
                                                  stop.on.error = FALSE)$value)/const1
                    
                    Fdr_power.h0=ifelse(Fdr_endpoints<0.5,
                                        integrate(kernel_f0,0,Fdr_endpoints,rel.tol=1e-5,
                                                  stop.on.error = FALSE)$value/const0,
                                        integrate(kernel_f0,Fdr_endpoints,1,rel.tol=1e-5,
                                                  stop.on.error = FALSE)$value/const0)
                  }else{
                    Fdr_power.h1=1-integrate(kernel_f1,min(Fdr_endpoints),max(Fdr_endpoints),
                                             rel.tol=1e-6,stop.on.error = FALSE)$value/const1
                    
                    Fdr_power.h0=1-integrate(kernel_f0,min(Fdr_endpoints),max(Fdr_endpoints),
                                             rel.tol=1e-6,stop.on.error = FALSE)$value/const0
                  }
                }
                EPower_Fdr=Fdr_power.h1
                CPower_Fdr=new.tau[2]*Fdr_power.h1/
                  (new.tau[1]*Fdr_power.h0+new.tau[2]*Fdr_power.h1)
                
                ##conditional error I/II
                CTypeI_Fdr=1-CPower_Fdr
                CTypeII_Fdr=tilted.tau[2]*ETypeII_Fdr/
                  (tilted.tau[1]*(1-ETypeI_Fdr)+tilted.tau[2]*ETypeII_Fdr)
                
                list(c(EPower_fdr,CPower_fdr,ETypeI_fdr,ETypeII_fdr,CTypeI_fdr,CTypeII_fdr),
                     c(EPower_Fdr,CPower_Fdr,ETypeI_Fdr,ETypeII_Fdr,CTypeI_Fdr,CTypeII_Fdr))
              }
  # rownames(res)=NULL
  res
}


##geometric mean
geomean<-function(x) 2*x[1]*x[2]/sum(x)

##approximation to theoretical power function
##not run!! results bad!!!
# mypower2<-function(q){
#   # t=c(seq(1e-4,1e-1,1e-4),seq(0.1,0.9,0.01),seq(0.9+1e-4,1-1e-4,1e-4))
#   t=seq(1e-3,1-1e-3,1e-3)
#   fdr<-function(x) new.tau[1]*new.f0(x)/new.fw(x)
#   fdr=Vectorize(fdr)
#   wh=which(fdr(t)<q)
#   ifelse(length(wh)==0,0,sum(new.f1(t)[wh])*0.001)
# }
# mypower2=Vectorize(mypower2)











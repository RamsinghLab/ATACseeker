#poisson - lognormal integral
poiscompv <- function(mu,sig,histr,histn=NULL){
  if(is.null(histn)){histn=0:length(histr)}
  t=mu+sqrt(2*sig)*gknots$nodes
  -sum(sapply(1:length(histr),function(i){
    histr[i]*log(sum(dpois(histn[i],exp(t))*gknots$weights/sqrt(pi)))
  }))
}

#poisson lognormal density
dpoiscomp <- function(d,mu,sig){
  t=mu+sqrt(2*sig)*gknots$nodes
  log(sum(dpois(d,exp(t))*gknots$weights/sqrt(pi)))
}

#negative binomial 
fbinom<-function(x,lmu,size,histn=NULL){
  if(is.null(histn)){histn=0:(length(histr)-1)}  
  -sum(x*dnbinom(histn,mu=exp(lmu),size=size,log=T))
}

#optimizer over pois lognormal
fitPoisLognorm<-function(tdd,lcset){
  mcur=log(sum(lcset*tdd)/sum(tdd))
  sigcur=optimize(poiscompv,mu=mcur,histr=tdd,interval=c(1e-5,5))
  for(i in 1:100){
    mcur=optimize(poiscompv,sig=sigcur$minimum,histr=tdd,histn=lcset,c(-10,5))
    sigcur=optimize(poiscompv,mu=mcur$minimum,histr=tdd,histn=lcset,c(0,10))
  }
  c(mcur$minimum,sigcur$minimum)
}

#####
# LCD fit

#fit log-concave density using proximal gradient 
fitPLCD<-function(lcset,lcinit,minval=-10,maxval=max(log(lcset))*2,gsize=1000,step=500,maxn=max(lcset+50,500),tol=1e-5,maxit=5000,lambda=0,debug=F){
  r=c(minval,maxval)
  t=seq(r[1],r[2],length=gsize)
  eps=step/sum(lcinit)#1e-5
  eld=rep(-log(length(t)),length(t))
  eldx=rep(-log(length(t)),length(t))
  eldxp=rep(0,length(t))
  if(debug){  par(mfrow=c(2,1))}
  err=Inf
  ctr=0
  while(err>tol & maxit > ctr){
    #print(ctr)
    sevun=sapply(1:maxn,function(i){
      evv=exp(eld+dpois(i-1,exp(t),log=T))
      evv
    })
    snn=t(t(sevun)/(colSums(sevun)+1e-200))
    cdelt=snn[,lcset+1]%*%lcinit-rowSums(sevun)/sum(sevun)*sum(lcinit)
    slc=eld+eps*(cdelt-lambda)
    crs=conreg(t,slc,maxit=c(5000,20))
    if(debug){
      plot(crs$yf,type='l')
      plot(cdelt,type='l')
    }
    eldt=crs$yf
    eldx=eldt-log(sum(exp(eldt)))
    eld=eldx+(ctr-1)/(ctr+2)*(eldx-eldxp)
    err=sum((eldx-eldxp)^2)
    eldxp=eldx
    if(debug){
      print(log(lcinit/sum(lcinit))-log(colSums(sevun)/sum(sevun))[lcset+1])
      print(c(err,ctr))
    }
    ctr=ctr+1
  }
  list(t,eld,sevun)
}

####
# continuous mapping

getflx <- function(flcd,xseq,nint=10000){
  approx(cumsum(colSums(flcd[[3]])/sum(flcd[[3]])),xseq,xout=seq(0,1,length=nint),rule=0)
}

getMapfun <- function(flcd,mm){
  xseq=seq(0,(ncol(flcd[[3]])-1),length=ncol(flcd[[3]]))
  df = dfun(log(xseq),exp(mm)); df=exp(df)/sum(exp(df))
  bfun=approx(cumsum(df),xseq,xout=seq(0,1,length=10000),rule=0)
  mapfun=approxfun(getflx(flcd,xseq,10000)$y,bfun$y,rule=0)
  mapfun
}


eval = function(x,lcinit,lcset){
  lh = log(lcinit/sum(lcinit))
  map=x[-1]
  llamb=x[1]
  lamb = exp(llamb)
  mv=cumsum(c(lcset[1],map))
  lhv= -lamb-lgamma(mv+1)+mv*log(lamb)
  sum(exp(lh)*(lh-lhv)^2)
}

grad = function(x,lcinit,lcset){
  lh = log(lcinit/sum(lcinit))
  map=x[-1]
  llamb=x[1]
  lamb = exp(llamb)
  mv=cumsum(c(lcset[1],map))
  lhv= -lamb-lgamma(mv+1)+mv*log(lamb)
  lhd= -2*(lh-lhv)*exp(lh)
  lambdiff = -exp(llamb) + mv
  lambgrad=sum(lhd*lambdiff)
  dgrad = (log(lamb) - psigamma(1+mv,0))*lhd
  dgdiff=  rev(cumsum(rev(dgrad[-1])))
  c(lambgrad,dgdiff)
}

dfun=function(x,l){exp(x)*log(l)-l-lgamma(exp(x)+1)}

fittrunc <- function(lcinit,lcset,flcd){
#find trunctaion point
  dkk=sapply(2:min(length(lcset),100),function(i){
  #print(i)
  # i = candidate truncation point
    tlcset=lcset
    ip1=min(i+1,length(lcinit))
  # set truncated counts
    tlcset[ip1:length(lcinit)]=lcset[i]
  # set truncated freqs
    tlcinit=lcinit
    tlcinit[i]=lcinit[i]+sum(lcinit[ip1:(length(lcinit))])
    tlcinit[ip1:(length(lcinit))]=0
    tmu=sum(tlcset*tlcinit)/sum(tlcinit) #mean truncated pois
    dpsub=dpois(lcset[1:i],tmu,log=T) #likelihood of truncated pois
    k1=sum(lcinit/sum(lcinit)*((dpsub[match(tlcset,lcset)]-log(lcinit/sum(lcinit)))^2))
    c(k1)
  })
}

fitdiscrete<- function(lcinit,lcset,flcd){
  maps=sapply(seq(-10,5,length=500),function(llinit){
    linit=exp(llinit)
    flsums=log(colSums(flcd[[3]])/sum(flcd[[3]])+1e-10)
    dfsums=dfun(seq(-10,5,length=200),linit)
    if(any(is.finite(dfsums))){
      os=outer(dpois(0:200,linit,log=T)-max(dfsums),(flsums-max(flsums))[lcset+1],'-')^2
      lcmap=c(1,sapply(2:ncol(os),function(i){which.min(os[2:min(i,nrow(os)),i,drop=F])})+1)
      lcmap[1]=lcset[1]+1
      lcmap[2]=lcset[2]+1
      lipost=sum((lcmap-1)*lcinit)/sum(lcinit)
      fseqerr=sum(lcinit/sum(lcinit)*(dpois(floor(lcmap-1),lambda=lipost,log=T)-log(lcinit/sum(lcinit)))^2)
      list(fseqerr,lcmap)
    }else{
    list(Inf,lcinit)
  }
  })
  minmap=maps[,which.min(maps[1,])]
  minmap
}

truncate<-function(lcset,trunc.fitted){
  trp=lcset
  trp[trp>(which.min(trunc.fitted))]=which.min(trunc.fitted)
  trp
}


fitcont.ll <- function(lcinit,lcset,flcd){
  cmaps=sapply(seq(-3,3,length=100),function(llinit){
    mapfun=getMapfun(flcd,llinit)
    fseqerr=sum(lcinit/sum(lcinit)*(dfun(log(mapfun(lcset)),sum(mapfun(lcset)*lcinit)/sum(lcinit)) - log(lcinit/sum(lcinit)))^2)
    list(fseqerr,mapfun(lcset),mapfun)
  })
  #
  elp=exp(flcd[[2]])
  mm=flcd[[1]][which(cumsum(elp)>sum(elp)/2)[1]]
  mediansel=which.min((mm-seq(-3,3,length=100))^2)
  #
  cminmap=cmaps[,which.min(cmaps[1,])]
  cminmap[[2]]=cminmap[[2]]
  cminmap.median=cmaps[,mediansel]
  cminmap.median[[2]]=cminmap.median[[2]]
  #
  par=c(sum(lcset*lcinit)/sum(lcinit),diff(cminmap[[2]]*sum(lcset*lcinit)/sum(lcinit)/(sum(cminmap[[2]]*lcinit)/sum(lcinit))))
  op=optim(par,eval,grad,control=list(trace=0,maxit=10000),lower=c(-Inf,rep(0.00001,length(par))),method='L-BFGS',lcset=lcset,lcinit=lcinit)
  #
  newmap = c(lcset[1],cumsum(op$par[-1]))
  lh = log(lcinit/sum(lcinit))
  mulfact=exp(op$par[1])/sum(newmap*exp(lh))
  newmap=newmap*mulfact
  list(cminmap,cminmap.median,newmap,op)
}

getContMap<-function(cont.map){
  list(cont.map[[1]][[2]]-cont.map[[1]][[2]][1],cont.map[[2]][[2]]-cont.map[[2]][[2]][1],cont.map[[3]]-cont.map[[3]][1])
}

#####
# helper functions

tabout<-function(x,w){
  rx=range(x)
  zxs=rep(0,max(rx)-min(rx)+1)
  for(i in 1:length(zxs)){
    zxs[i]=sum(w[x==min(rx)+i-1])
  }
  cbind(min(rx):max(rx),zxs)
}

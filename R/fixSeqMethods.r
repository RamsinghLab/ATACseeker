#!/usr/bin/Rscript

source('functions.r')

set.seed(1)

#Required packages
library("statmod")
library('cobs')
library('logcondens')

## default 
gknots = gauss.quad(40,kind="hermite")

###
# Simulations (nbin) (generates input.csv)
if(!file.exists('input.csv')){
  print('NO INPUT.CSV DETECTED: GENERATING SIMULATED INPUT')
  sizetr=0.1
  mutr=0.1
  countin=rnbinom(80000,mu=mutr,size=sizetr)
  evalpts=seq(-20,10,length=500)
  lcinit=table(as.vector(countin))
  lcset=as.double(names(lcinit))
  write.csv(data.frame(frequency=as.double(lcinit),counts=lcset),file='input.csv',row.names=F)
}

input=read.csv('input.csv')
lcsub=input[,2]<=500
lcset=sort(input[lcsub,2])
lcinit=input[lcsub,1][order(input[lcsub,2])]

###
# fit log-concave distr.
flcd = fitPLCD(lcset,lcinit,maxit=100,lambda=1e-5)

trunc.fitted = fittrunc(lcinit,lcset,flcd)
trp = truncate(lcset,trunc.fitted)

###
# rounding / count adjustment

discrete.fitted = fitdiscrete(lcinit,lcset,flcd)

continuous.fitted=fitcont.ll(lcinit,lcset,flcd)
cont.maps = getContMap(continuous.fitted)

op=continuous.fitted[[4]]
pdf('truncation-loss.pdf')
plot((2:min(length(lcset),100))-1,trunc.fitted,ylim=range(c(discrete.fitted[[1]],op$value,trunc.fitted)),log='y',xlab='truncation point (reads / bp)',ylab='expected error in log-likelihood',type='l')
abline(h=discrete.fitted[[1]],col='blue')
abline(h=op$value,col='red')
legend('topleft',lwd=1,col=c('black','blue','red','maroon'),legend=c('truncation','integer log-concave mapping','weighted log-concave mapping'))
dev.off()

write.csv(data.frame(counts=lcset,weighted=cont.maps[[1]],weighted.medianrule=cont.maps[[2]],weighted.logloss=cont.maps[[3]],integer=(discrete.fitted[[2]]-1),optimal.truncation=trp),file='output.csv')

###
# other distribution fits
tdd=lcinit
mcur=log(sum(lcinit*lcset)/sum(lcinit))
sigcur=optimize(poiscompv,mu=mcur,histr=tdd,interval=c(1e-5,5))
for(i in 1:100){
  mcur=optimize(poiscompv,sig=sigcur$minimum,histr=tdd,histn=lcset,c(-10,5))
  sigcur=optimize(poiscompv,mu=mcur$minimum,histr=tdd,histn=lcset,c(0,10))
}

muin=log(sum(lcinit*lcset)/sum(lcinit))
for(i in 1:100){
ops=optimize(fbinom,lmu=muin,x=lcinit,histn=lcset,interval=c(0,10))
op=optimize(fbinom,size=ops$minimum,x=lcinit,histn=lcset,interval=c(-10,10))
muin=op$minimum
}
prob=ops$minimum/(ops$minimum+exp(muin))
size=ops$minimum

sevun=flcd[[3]]

pdf('llhfit.pdf')
plot(lcset,log(lcinit/sum(lcinit)),type='l',xlab='count',ylab='log-density',main='goodness of fit of various distributions',log='x')
maxn=max(lcset)
points(0:maxn,dnbinom(0:maxn,prob=prob,size=size,log=T),type='l',col='red')
points(0:maxn,sapply(0:maxn,function(i){dpoiscomp(i,mcur$minimum,sigcur$minimum)}),col='blue',type='l')
points(0:(ncol(sevun)-1),log(colSums(sevun)/sum(sevun)),col='purple',type='l')
points(0:maxn,dpois(0:maxn,sum(lcinit*lcset)/sum(lcinit),log=T),type='l',col='green')
abline(v=which.min(trunc.fitted),col='maroon')
legend('topright',legend=c('true','normal','gamma','l-concave','poisson','truncation point'),col=c('black','blue','red','purple','green','maroon'),lwd=1)
#
trtrunc=tabout(trp,lcinit)
plot(trtrunc[,1],log(trtrunc[,2]/sum(trtrunc[,2])),type='l',main='goodness of fit of a poisson fit (green) optimal truncation (black)',xlab='count',ylab='log likelihood')
points(0:max(trp),dpois(0:max(trp),sum(trtrunc[,1]*trtrunc[,2])/sum(trtrunc[,2]),log=T),type='l',col='green')
#
trmap=tabout((discrete.fitted[[2]]-1),lcinit)
plot(trmap[,1],log(trmap[,2]/sum(trmap[,2])),type='l',main='goodness of fit of poisson (green) to optimal integer mapping (black)',xlab='count',ylab='log-likelihood')
points(0:max(trmap[,1]),dpois(0:max(trmap[,1]),sum(trmap[,1]*trmap[,2])/sum(trmap[,2]),log=T),type='l',col='green')
#
dev.off()

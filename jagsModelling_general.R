library(rjags)
library(coda)

rm(list=ls())

outDir<-paste(sep="/","output",Sys.Date())
if(!dir.exists(outDir)){dir.create(outDir)}

##--------------------
## set seed
##

todayDate<-as.integer(unlist(strsplit(split="-",as.character(Sys.Date()))))
yy<-todayDate[1]
mm<-todayDate[2]
dd<-todayDate[3]
randSeed<-rpois(1,lambda=sample(size=1,x=c(10,50,100,500,1000)))*dd+10*mm+yy

set.seed(randSeed) # we will reset the seed before every call to jags.model()


##--------------------
## read data
##

tt<-read.csv("data/TTRInvASensitivity20170724_corrected.csv")
tt<-tt[!is.na(tt$Sample.ID) & !is.na(tt$Visit..),]


##--------------------
## read recruitment dates and compute actual timepoints (in months)
##

dates<-read.csv("data/salexpoLIMSDataSetComplete.csv")
dates$monthYear<-gsub("^[0-9]+/","",dates$DATE) # dates run from 08/2013 to 12/2014
datesDictionary<-data.frame(monthYear=c(paste(sep="/",c("08","09","10","11","12"),"2013"),paste(sep="/",c(paste(sep="","0",1:9),paste(sep="","1",0:2)),"2014")),t=1:17)
dates$t<-datesDictionary$t[match(dates$monthYear,datesDictionary$monthYear)]

##--------------------
## manual correction for now of the visit number for one sample - DOUBLE CHECK WITH ANGE AS THIS IS NOT 100% OBVIOUS - cf email to Ange from 17/04/2018
tt$Visit..[tt$Sample.ID==59 & tt$Visit..==9 & tt$TTR=="neg"]<-8
## end manual correction

tt$t<-dates$t[match(paste(sep="::",tt$Sample.ID,tt$Visit..),paste(sep="::",dates$ID,gsub(pattern="Visit ",replacement="",dates$Visit)))]
tt$monthYear<-datesDictionary$monthYear[match(tt$t,datesDictionary$t)]


##--------------------
## manual imputation (it's always very obvious) of "missing" dates for visits
##

# despite being so obvious: confirm with Ange
tt$t[tt$Sample.ID==8 & tt$Visit..==11]<-12
tt$t[tt$Sample.ID==10 & tt$Visit..==10]<-11
tt$t[tt$Sample.ID==14 & tt$Visit..==7]<-8
tt$t[tt$Sample.ID==14 & tt$Visit..==10]<-12
tt$t[tt$Sample.ID==15 & tt$Visit..==7]<-8
tt$t[tt$Sample.ID==17 & tt$Visit..==12]<-13
tt$t[tt$Sample.ID==28 & tt$Visit..==5]<-7
tt$t[tt$Sample.ID==51 & tt$Visit..==8]<-12
tt$t[tt$Sample.ID==51 & tt$Visit..==11]<-15
tt$t[tt$Sample.ID==53 & tt$Visit..==10]<-14


##--------------------
## recode t so that time proceeds monthly from the first visit - this is to avoid duplicated t values for a given individuals (problem caused by one visit early one month and another late the same month)
##
tt$t2<-tt$t
tt$t2[tt$Sample.ID==3 & tt$Visit..==9]<-9
tt$t2[tt$Sample.ID==3 & tt$Visit..==10]<-10
tt$t2[tt$Sample.ID==11 & tt$Visit..==12]<-13
tt$t2[tt$Sample.ID==21 & tt$Visit..==7]<-8
tt$t2[tt$Sample.ID==24 & tt$Visit..==11]<-13
tt$t2[tt$Sample.ID==40 & tt$Visit..==9]<-12
tt$t2[tt$Sample.ID==40 & tt$Visit..==10]<-13    
tt$t2[tt$Sample.ID==40 & tt$Visit..==11]<-14
tt$t2[tt$Sample.ID==41 & tt$Visit..==7]<-10
tt$t2[tt$Sample.ID==49 & tt$Visit..==6]<-10
tt$t2[tt$Sample.ID==49 & tt$Visit..==7]<-11
tt$t2[tt$Sample.ID==49 & tt$Visit..==9]<-13
    
##--------------------
## recode test results in a better, cleaner way
##

tt2<-tt
tt2$TTR<-ifelse(tolower(tt2$TTR)=="pos",1,ifelse(tolower(tt2$TTR)=="neg" | tolower(tt2$TTR)=="ind",0,NA))
tt2$InvA<-ifelse(tolower(tt2$InvA)=="pos",1,ifelse(tolower(tt2$InvA)=="neg" | tolower(tt2$InvA)=="ind",0,NA))
tt2$Culture<-ifelse(tolower(tt2$Culture)=="pos",1,ifelse(tolower(tt2$Culture)=="neg",0,NA))
tt2$TAC.TTR<-ifelse(tolower(tt2$TAC.TTR)=="pos",1,ifelse(tolower(tt2$TAC.TTR)=="neg" | tolower(tt2$TAC.TTR)=="ind",0,NA))
tt2$TAC.InvA<-ifelse(tolower(tt2$TAC.InvA)=="pos",1,ifelse(tolower(tt2$TAC.InvA)=="neg" | tolower(tt2$TAC.InvA)=="ind",0,NA))


##--------------------
## reformat for the way the JAGS/BUGS model file is set-up
##

uniqTimes<-sort(unique(tt2$t2))
nT<-length(uniqTimes)
uniqIds<-sort(unique(tt2$Sample.ID))
nN<-length(uniqIds)
uniqTests<-c("Culture","TTR","InvA","TAC.TTR","TAC.InvA")
nK<-length(uniqTests)
  

dat<-list()
dat$tests<-array(dim=c(nT,nN,nK))
dimnames(dat$tests)<-list(paste(sep="","t",uniqTimes),paste(sep="","ind",uniqIds),uniqTests)
for(k in 1:nK){
  dat$tests[,,uniqTests[k]]<-matrix(nrow=nT,ncol=nN)
}

for(i_t in 1:nT){
  for(i_n in 1:nN){
    t<-uniqTimes[i_t]
    id<-uniqIds[i_n]
    for(i_k in 1:nK){
      k<-uniqTests[i_k]
      dat$tests[i_t,i_n,i_k]<-ifelse(sum(tt2$Sample.ID==id & tt2$t2==t)==1,tt2[[k]][tt2$Sample.ID==id & tt2$t2==t],NA)
    }
  }
}

dat$N<-nN
dat$T<-nT
dat$K<-nK


##--------------------
## run the models
##

set.seed(randSeed)
jagsTimeHomoModel <- jags.model('scripts/jagsModelFile_TimeHomo_general.jags', data=dat, n.chains = 4, n.adapt = 1000)
#jagsTimeHomoModel <- jags.model('GordonMelita_StatsSupportAngeziwe/jagsModelFile_TimeHomo.jags', data=dat, n.chains = 4, n.adapt = 0, inits=list(pInfInit=0.05, pUninf2Inf=0.05, pInf2Inf=0.08, crpCult=c(1-0,0.5),crpTtrStd=c(1-0.9,0.95),crpInvaStd=c(1-0.85,0.9),crpTtrTac=c(1-0.9,0.95),crpInvaTac=c(1-0.85,0.9)))
#update(jagsTimeHomoModel,n.iter=50000) # note technicaly not required as we use coda.samples() below which also updates the model, but we use this as burn-in - so in total we will get 1,000 (adaptation iterations) + 10,000 (burn-in) + 100,000 (sampling) iterations
parsTimeHomoModel<-coda.samples(model=jagsTimeHomoModel,variable.names=c('pInfInit','pUninf2Inf','pInf2Inf','crp'),n.iter=10000) # 5*2+2 = 12 parameters (each of the crp parameters is a vector of length 2); pUninf2Inf fixed to pInfInit
save(list=ls(),file=paste(sep="/",outDir,"jags_AngeData_general_TimeHomoModel_TimeHomoModel.RData"))

set.seed(randSeed)
jagsTimeHeteroModel <- jags.model('scripts/jagsModelFile_TimeHetero_general.jags', data=dat, n.chains = 4, n.adapt = 1000)
#jagsTimeHeteroModel <- jags.model('GordonMelita_StatsSupportAngeziwe/jagsModelFile_TimeHetero.jags', data=dat, n.chains = 4, n.adapt = 0, inits=list(pInfInit=0.05, pUninf2Inf=0.05, pInf2Inf=0.08, crpCult=c(1-0,0.5),crpTtrStd=c(1-0.9,0.95),crpInvaStd=c(1-0.85,0.9),crpTtrTac=c(1-0.9,0.95),crpInvaTac=c(1-0.85,0.9)))
#update(jagsTimeHeteroModel,n.iter=50000) # note technicaly not required as we use coda.samples() below which also updates the model, but we use this as burn-in - so in total we will get 1,000 (adaptation iterations) + 10,000 (burn-in) + 100,000 (sampling) iterations
parsTimeHeteroModel<-coda.samples(model=jagsTimeHeteroModel,variable.names=c('pInfInit','pUninf2Inf','pInf2Inf','crp'),n.iter=10000) # 5*2+(13-1)*2+1 = 35 parameters (each of the crp parameters is a vector of length 2); pUninf2Inf fixed to pInfInit
save(list=ls(),file=paste(sep="/",outDir,"jags_AngeData_general_TimeHeteroModel_TimeHeteroModel.RData"))

set.seed(randSeed)
jagsTimeHomoMixedModel <- jags.model('scripts/jagsModelFile_TimeHomoMixed_general.jagss', data=dat, n.chains = 4, n.adapt = 1000)
#jagsTimeHeteroMixedModel <- jags.model('GordonMelita_StatsSupportAngeziwe/jagsModelFile_TimeHeteroMixed.jags', data=dat, n.chains = 4, n.adapt = 1000)
update(jagsTimeHomoMixedModel,n.iter=10000)
parsTimeHomoMixedModel<-coda.samples(model=jagsTimeHomoMixedModel,variable.names=c('pInfInit','pUninf2Inf','pInf2Inf','crpNoRand','interceptStateCrp','slopeRandCrp'),n.iter=20000)
save(list=setdiff(ls(),c("jagsTimeHomoModel","parsTimeHomoModel")),file=paste(sep="/",outDir,"jags_AngeData_general_TimeHomoMixedModel.RData"))

#jagsTimeHeteroModel <- jags.model('GordonMelita_StatsSupportAngeziwe/jagsModelFile_TimeHetero.jags', data=dat, n.chains = 4, n.adapt = 1000)


##--------------------
## model comparison: DIC & PED
##

## A) standard DIC

dicTimeHomoModel<-dic.samples(jagsTimeHomoModel,n.iter=5000,type="pD")
dicTimeHeteroModel<-dic.samples(jagsTimeHeteroModel,n.iter=5000,type="pD")
dicTimeHomoMixedModel<-dic.samples(jagsTimeHomoMixedModel,n.iter=5000,type="pD") # the penalty term fails to evaluate for some samples, hence need to deal with this; we do this in 2 ways, see below

diffdic(dicTimeHomoModel,dicTimeHeteroModel)

# way 1: impute NaNs in penalty terms to average penalty
dicTimeHomoMixedModelNoMiss<-dicTimeHomoMixedModel
dicTimeHomoMixedModelNoMiss$penalty[is.nan(dicTimeHomoMixedModelNoMiss$penalty)]<-mean(dicTimeHomoMixedModelNoMiss$penalty[!is.nan(dicTimeHomoMixedModelNoMiss$penalty)])

diffdic(dicTimeHomoMixedModelNoMiss,dicTimeHomoModel)

# way 2: compare only non-NaN samples
idxNaN<-which(is.nan(dicTimeHomoMixedModel$penalty))
dd<-dicTimeHomoModel
ddMix<-dicTimeHomoMixedModel
if(length(idxNaN)>0){
  dd$penalty<-dd$penalty[-idxNaN]
  dd$deviance<-dd$deviance[-idxNaN]
  ddMix$penalty<-ddMix$penalty[-idxNaN]
  ddMix$deviance<-ddMix$deviance[-idxNaN]
}

diffdic(ddMix,dd) # prefer the mixed model; just about


## B) Penalized Expected Deviance (PED)

pedTimeHomoModel<-dic.samples(jagsTimeHomoModel,n.iter=5000,type="popt")
pedTimeHeteroModel<-dic.samples(jagsTimeHeteroModel,n.iter=5000,type="popt")
pedTimeHomoMixedModel<-dic.samples(jagsTimeHomoMixedModel,n.iter=5000,type="popt") # the penalty term fails to evaluate for some samples, hence need to deal with this; we do this in 2 ways, see below

# way 1: impute NaNs in penalty terms to average penalty
pedTimeHomoMixedModelNoMiss<-pedTimeHomoMixedModel
pedTimeHomoMixedModelNoMiss$penalty[is.nan(pedTimeHomoMixedModelNoMiss$penalty)]<-mean(pedTimeHomoMixedModelNoMiss$penalty[!is.nan(pedTimeHomoMixedModelNoMiss$penalty)])

diffdic(pedTimeHomoMixedModelNoMiss,pedTimeHomoModel) # prefer the mixed model; somewhat

# way 2: compare only non-NaN samples
idxNaN<-which(is.nan(pedTimeHomoMixedModel$penalty))
dd<-pedTimeHomoModel
ddMix<-pedTimeHomoMixedModel
if(length(idxNaN)>0){
  dd$penalty<-dd$penalty[-idxNaN]
  dd$deviance<-dd$deviance[-idxNaN]
  ddMix$penalty<-ddMix$penalty[-idxNaN]
  ddMix$deviance<-ddMix$deviance[-idxNaN]
}

diffdic(ddMix,dd) # prefer the standard model, somewhat



##--------------------
## MCMC diagnostics
##

print(summary(parsTimeHomoModel))
print(gelman.diag(parsTimeHomoModel))

pdf(paste(sep="/",outDir,"jags_AngeData_general_TimeHomoModel_MCMCdiagnostics.pdf"))
plot(parsTimeHomoModel)
dev.off()

print(summary(parsTimeHeteroModel))
print(gelman.diag(parsTimeHeteroModel))

pdf(paste(sep="/",outDir,"jags_AngeData_general_TimeHeteroModel_MCMCdiagnostics.pdf"))
plot(parsTimeHeteroModel)
dev.off()

print(gelman.diag(parsTimeHomoMixedModel))
print(summary(parsTimeHomoMixedModel))

pdf(paste(sep="/",outDir,"jags_AngeData_general_TimeHomoMixedModel_MCMCdiagnostics_new.pdf"))
plot(parsTimeHomoMixedModel)
dev.off()


##--------------------
## save everything in one RData object
##

save(list=ls(),file=paste(sep="/",outDir,"jags_AngeData_general_everything.RData"))
#setwd("~/../work/MLW_LSTM")
setwd("~/work")



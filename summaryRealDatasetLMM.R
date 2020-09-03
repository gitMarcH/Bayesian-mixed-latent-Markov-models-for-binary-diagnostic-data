library(MASS)
library(RColorBrewer)
library(HDInterval)
library(rjags)

rm(list=ls())

# input arguments

# args<-commandArgs(TRUE)
# rdataJagsFile<-args[1] # typically "output/ToBeUsedForPaper/Real/jags_AngeData_general_TimeHomoModel.RData"
# objName<-args[2] # typically "parsTimeHomoModel"
# if(!is.na(args[3])){plotLims<-as.numeric(unlist(strsplit(split=",",args[3])))}else{plotLims<-NA}
# if(!is.na(args[4])){testNames<-unlist(strsplit(split=",",args[4]))}else{testNames<-NA}; testNames[testNames=="NA"]<-NA
# if(!is.na(args[5])){parsNames<-unlist(strsplit(split=",",args[5]))}else{parsNames<-NA}; parsNames[parsNames=="NA"]<-NA
# if(!is.na(args[6])){trueSens<-as.numeric(unlist(strsplit(split=",",args[6])))}else{trueSens<-NA}
# if(!is.na(args[7])){trueSpec<-as.numeric(unlist(strsplit(split=",",args[7])))}else{trueSpec<-NA}
# if(!is.na(args[8])){outPrefixAdd<-args[8]}else{outPrefixAdd<-""}

rdataJagsFile<-"output/ToBeUsedForPaper/Real/jags_AngeData_general_TimeHomoModel.RData"
objName<-"parsTimeHomoModel"
plotLims<-NA
testNames<-c("Culture","TTR","InvA","TAC-TTR","TAC-InvA")
parsNames<-NA
trueSens<-NA
trueSpec<-NA
outPrefixAdd<-"_manuscriptComments"

outPrefix<-paste(sep="",gsub(rdataJagsFile,pattern=".RData",replacement=paste(sep="","_",gsub(Sys.Date(),pattern="-",replacement=""))),outPrefixAdd)


# load input
load(rdataJagsFile)
mcmcObj<-get(objName)

# define needed subroutines
addWhisker<-function(low,high,anchor,whiskWidth=0.005,vertical=T,col="black"){
      # low = lower whisker height; can be a vector
      # high = upper whisker height; can be a vector
      # anchor = position on which to anchor the whisker [x position if vertical==T, else y position]
      # whiskWidth = width of whisker bars
      # vertical = logical; if TRUE whiskers are plotted on the y axis, otherwise on the x axis

    if(vertical){
        segments(x0=anchor,x1=anchor,y0=low,y1=high,col=col)
        segments(x0=anchor-whiskWidth,x1=anchor+whiskWidth,y0=low,y1=low,col=col)
        segments(x0=anchor-whiskWidth,x1=anchor+whiskWidth,y0=high,y1=high,col=col)
    }else{
        segments(y0=anchor,y1=anchor,x0=low,x1=high,col=col)
        segments(y0=anchor-whiskWidth,y1=anchor+whiskWidth,x0=low,x1=low,col=col)
        segments(y0=anchor-whiskWidth,y1=anchor+whiskWidth,x0=high,x1=high,col=col)
    }
}

darkenCol<-function(color,factor=1.5){
    col<-col2rgb(color)
    col<-col/factor
    col<-rgb(t(col),maxColorValue=255)
    col
}

addGrid<-function(doX=TRUE,doY=TRUE,xExt=c(0,1),xBwMinor=0.05,xBwMajor=0.1,yExt=xExt,yBwMinor=xBwMinor,yBwMajor=xBwMajor){
    if(doX){
        for(xTmp in seq(xExt[1],xExt[2],by=xBwMinor)){abline(v=xTmp,col="lightgray",lty=2,lwd=0.1)}
        for(xTmp in seq(xExt[1],xExt[2],by=xBwMajor)){abline(v=xTmp,col="gray",lty=2,lwd=0.3)}
    }
    if(doY){
        for(yTmp in seq(yExt[1],yExt[2],by=yBwMinor)){abline(h=yTmp,col="lightgray",lty=2,lwd=0.1)}
        for(yTmp in seq(yExt[1],yExt[2],by=yBwMajor)){abline(h=yTmp,col="gray",lty=2,lwd=0.3)}
    }
}

addAlpha <- function(COLORS, ALPHA){ # courtesy of https://menugget.blogspot.com/2012/04/adding-transparent-image-layer-to-plot.html#more
     if(missing(ALPHA)) stop("provide a value for alpha between 0 and 1")
      RGB <- col2rgb(COLORS, alpha=TRUE)
      RGB[4,] <- round(RGB[4,]*ALPHA)
      NEW.COLORS <- rgb(RGB[1,], RGB[2,], RGB[3,], RGB[4,], maxColorValue = 255)
      return(NEW.COLORS)
}

getMaP<-function(dens,xPosterior,yPosterior,is.dens=T){
    if(!is.dens){
        dens<-kde2d(xPosterior,yPosterior,n=200)
    }
    idxMax<-which(dens$z==max(dens$z))
    idxMaxX<-idxMax %% length(dens$y)
    if(idxMaxX==0){idxMaxX<-length(dens$y)}
    idxMaxY<-floor(idxMax / length(dens$x)) + 1
    return(c(dens$x[idxMaxX],dens$y[idxMaxY]))
}


# sensitivities vs specificities
nChains<-length(mcmcObj)
nTests<-sum(grepl(pattern="crp",colnames(mcmcObj[[1]])))/2
if(is.na(testNames[1])){testNames<-paste(sep="","test",1:nTests)}
testsSens<-vector("list",nTests); names(testsSens)<-testNames
testsSpec<-vector("list",nTests); names(testsSpec)<-testNames
for(j in 1:nChains){
    colnames(mcmcObj[[j]])<-gsub(pattern="NoRand",replacement="",colnames(mcmcObj[[j]]))
    for(k in 1:nTests){
        testsSens[[k]]<-c(testsSens[[k]],mcmcObj[[j]][,paste(sep="","crp[",k,",2]")])
        testsSpec[[k]]<-c(testsSpec[[k]],1-mcmcObj[[j]][,paste(sep="","crp[",k,",1]")])
    }
}

testsDens<-vector("list",nTests); names(testsDens)<-testNames
for(k in 1:nTests){
    testsDens[[k]]<-kde2d(testsSens[[k]],testsSpec[[k]],n=200)
}

cols<-colorRampPalette(brewer.pal(7,"Spectral")[-c(4,5)])(nTests)
colsDarker<-cols
for(i in 1:length(cols)){colsDarker[i]<-darkenCol(cols[i],factor=1.5)}

if(is.na(plotLims[1])){plotLims<-c(0.35,1,0.8,1)}

png(paste(sep="",outPrefix,"_sensVSspec.png"),width=12,height=6,units="in",res=600)
plot(xlim=plotLims[1:2],ylim=plotLims[3:4],main="",xlab="sensitivity",ylab="specificity",0:1,0:1,col=cols,pch=20,type="n",asp=1)
addGrid()
par(xpd=T)

for(k in 1:nTests){
  paletteCust<-function(n){
    addAlpha(colorRampPalette(c("white",cols[k],cols[k],darkenCol(cols[k],factor=1.1)))(n),ALPHA=seq(0,1,length=n))
  }
  
  .filled.contour(testsDens[[k]]$x,testsDens[[k]]$y,testsDens[[k]]$z,col=paletteCust(50),levels=seq(0,max(testsDens[[k]]$z),length=50))
  contour(testsDens[[k]],col=addAlpha(colorRampPalette(c("white",cols[k]))(25),ALPHA=seq(0.2,0.8,length=25)),add=T,drawlabels=F,levels=seq(4,max(testsDens[[k]]$z),length=25))
}

fullMCMCtable<-mcmcObj[[1]]
for(i in 2:length(mcmcObj)){
  fullMCMCtable<-rbind(fullMCMCtable,mcmcObj[[i]])
}
fullMCMCtable<-as.data.frame(fullMCMCtable)
nc<-ncol(fullMCMCtable)

fullMCMCtable$ppvCult<-fullMCMCtable$`crp[1,2]`*fullMCMCtable$pInfInit / (fullMCMCtable$`crp[1,2]`*fullMCMCtable$pInfInit + fullMCMCtable$`crp[1,1]`*(1-fullMCMCtable$pInfInit))
fullMCMCtable$ppvTTR<-fullMCMCtable$`crp[2,2]`*fullMCMCtable$pInfInit / (fullMCMCtable$`crp[2,2]`*fullMCMCtable$pInfInit + fullMCMCtable$`crp[2,1]`*(1-fullMCMCtable$pInfInit))
fullMCMCtable$ppvInvA<-fullMCMCtable$`crp[3,2]`*fullMCMCtable$pInfInit / (fullMCMCtable$`crp[3,2]`*fullMCMCtable$pInfInit + fullMCMCtable$`crp[3,1]`*(1-fullMCMCtable$pInfInit))
fullMCMCtable$ppvTacTTR<-fullMCMCtable$`crp[4,2]`*fullMCMCtable$pInfInit / (fullMCMCtable$`crp[4,2]`*fullMCMCtable$pInfInit + fullMCMCtable$`crp[4,1]`*(1-fullMCMCtable$pInfInit))
fullMCMCtable$ppvTacInvA<-fullMCMCtable$`crp[5,2]`*fullMCMCtable$pInfInit / (fullMCMCtable$`crp[5,2]`*fullMCMCtable$pInfInit + fullMCMCtable$`crp[5,1]`*(1-fullMCMCtable$pInfInit))

fullMCMCtable$npvCult<-(1-fullMCMCtable$`crp[1,1]`)*(1-fullMCMCtable$pInfInit) / ((1-fullMCMCtable$`crp[1,2]`)*(1-fullMCMCtable$pInfInit) + (1-fullMCMCtable$`crp[1,1]`)*(1-fullMCMCtable$pInfInit))
fullMCMCtable$npvTTR<-(1-fullMCMCtable$`crp[2,1]`)*(1-fullMCMCtable$pInfInit) / ((1-fullMCMCtable$`crp[2,2]`)*(1-fullMCMCtable$pInfInit) + (1-fullMCMCtable$`crp[2,1]`)*(1-fullMCMCtable$pInfInit))
fullMCMCtable$npvInvA<-(1-fullMCMCtable$`crp[3,1]`)*(1-fullMCMCtable$pInfInit) / ((1-fullMCMCtable$`crp[3,2]`)*(1-fullMCMCtable$pInfInit) + (1-fullMCMCtable$`crp[3,1]`)*(1-fullMCMCtable$pInfInit))
fullMCMCtable$npvTacTTR<-(1-fullMCMCtable$`crp[4,1]`)*(1-fullMCMCtable$pInfInit) / ((1-fullMCMCtable$`crp[4,2]`)*(1-fullMCMCtable$pInfInit) + (1-fullMCMCtable$`crp[4,1]`)*(1-fullMCMCtable$pInfInit))
fullMCMCtable$npvTacInvA<-(1-fullMCMCtable$`crp[5,1]`)*(1-fullMCMCtable$pInfInit) / ((1-fullMCMCtable$`crp[5,2]`)*(1-fullMCMCtable$pInfInit) + (1-fullMCMCtable$`crp[5,1]`)*(1-fullMCMCtable$pInfInit))

sensEstimates<-data.frame(matrix(nrow=nTests,ncol=6))
rownames(sensEstimates)<-testNames
colnames(sensEstimates)<-c("lower95SymCI","median","upper95SymCI","lower95HdiCI","MaP","upper95HdiCI")
specEstimates<-sensEstimates
ppvEstimates<-sensEstimates
npvEstimates<-specEstimates
hdiTmp<-hdi(mcmcObj,credMass=0.95)
hdiTmp2<-hdi(fullMCMCtable[,(nc+1):(nc+2*nTests)],credMass=0.95)
for(k in 1:nTests){
    sensEstimates[k,1:3]<-summary(mcmcObj)$quantiles[paste(sep="","crp[",k,",2]"),c("2.5%","50%","97.5%")]
    specEstimates[k,1:3]<-1-summary(mcmcObj)$quantiles[paste(sep="","crp[",k,",1]"),c("97.5%","50%","2.5%")]
    map<-getMaP(dens=testsDens[[k]])
    sensEstimates$MaP[k]<-map[1]
    specEstimates$MaP[k]<-map[2]
    sensEstimates$lower95HdiCI[k]<-hdiTmp["lower",paste(sep="","crp[",k,",2]")]
    sensEstimates$upper95HdiCI[k]<-hdiTmp["upper",paste(sep="","crp[",k,",2]")]
    specEstimates$lower95HdiCI[k]<-1-hdiTmp["upper",paste(sep="","crp[",k,",1]")]
    specEstimates$upper95HdiCI[k]<-1-hdiTmp["lower",paste(sep="","crp[",k,",1]")]
    
    ppvEstimates[k,1:3]<-quantile(probs=c(0.025,0.5,0.975),fullMCMCtable[,nc+k])
    npvEstimates[k,1:3]<-quantile(probs=c(0.025,0.5,0.975),fullMCMCtable[,nc+nTests+k])
    map<-getMaP(is.dens=F,xPosterior=fullMCMCtable[,nc+k],yPosterior=fullMCMCtable[,nc+nTests+k])
    ppvEstimates$MaP[k]<-map[1]
    npvEstimates$MaP[k]<-map[2]
    ppvEstimates$lower95HdiCI[k]<-hdiTmp2["lower",k]
    ppvEstimates$upper95HdiCI[k]<-hdiTmp2["upper",k]
    npvEstimates$lower95HdiCI[k]<-hdiTmp2["lower",nTests+k]
    npvEstimates$upper95HdiCI[k]<-hdiTmp2["upper",nTests+k]
}

addWhisker(low=sensEstimates$lower95HdiCI,high=sensEstimates$upper95HdiCI,anchor=specEstimates$MaP,vertical=F,col=cols,whiskWidth=0.005/((1-0.35)/(1-0.85)))
addWhisker(low=specEstimates$lower95HdiCI,high=specEstimates$upper95HdiCI,anchor=sensEstimates$MaP,vertical=T,col=cols)

points(sensEstimates$MaP,specEstimates$MaP,col=cols,pch=19,cex=1.5)

if(!is.na(trueSens[1]) & !is.na(trueSpec[1])){points(trueSens,trueSpec,col=cols,pch=4,cex=1)}

legend(x="bottomleft",pch=c(rep(19,5)),col=cols,bty="n",legend=testNames)
dev.off()


# pUninf2Inf, pInf2Inf
pUninf2InfMCMC<-numeric(0)
pInf2InfMCMC<-numeric(0)
for(j in 1:nChains){
    pUninf2InfMCMC<-c(pUninf2InfMCMC,mcmcObj[[j]][,"pUninf2Inf"])
    pInf2InfMCMC<-c(pInf2InfMCMC,mcmcObj[[j]][,"pInf2Inf"])
}

#pUninf2InfDens<-density(pUninf2InfMCMC)
#pInf2InfDens<-density(pInf2InfMCMC)
mapEstimatesTransitionMatrix<-getMaP(is.dens=F,xPosterior=pUninf2InfMCMC,yPosterior=pInf2InfMCMC)

# caterpillar plots, posterior density estimates, prior distributions
#nPars<-ncol(mcmcObj[[1]])
#if(is.na(parsNames[1])){parsNames<-colnames(mcmcObj[[1]])}


# summary table of model parameters for paper: sensitivity, specificity, PPV, NPV (posterior median, MAP, 95% HDI CrI)
estTable<-matrix(nrow=nTests,ncol=4*3)
rownames(estTable)<-testNames
colnames(estTable)<-apply(FUN=paste,MARGIN=1,collapse="_",expand.grid(c("postMedian","MAP","HDI"),c("sens","spec","ppv","npv"))[,2:1])
estTable<-as.data.frame(estTable)
for(k in 1:nTests){
  estTable$sens_postMedian[k]<-format(nsmall=4,round(digits=4,sensEstimates$median[k]))
  estTable$sens_MAP[k]<-format(nsmall=4,round(digits=4,sensEstimates$MaP[k]))
  estTable$sens_HDI[k]<-paste(sep="","(",paste(collapse=",",format(nsmall=4,round(digits=4,as.matrix(sensEstimates[k,c("lower95HdiCI","upper95HdiCI")])))),")")
  
  estTable$spec_postMedian[k]<-format(nsmall=4,round(digits=4,specEstimates$median[k]))
  estTable$spec_MAP[k]<-format(nsmall=4,round(digits=4,specEstimates$MaP[k]))
  estTable$spec_HDI[k]<-paste(sep="","(",paste(collapse=",",format(nsmall=4,round(digits=4,as.matrix(specEstimates[k,c("lower95HdiCI","upper95HdiCI")])))),")")
  
  estTable$ppv_postMedian[k]<-format(nsmall=4,round(digits=4,ppvEstimates$median[k]))
  estTable$ppv_MAP[k]<-format(nsmall=4,round(digits=4,ppvEstimates$MaP[k]))
  estTable$ppv_HDI[k]<-paste(sep="","(",paste(collapse=",",format(nsmall=4,round(digits=4,as.matrix(ppvEstimates[k,c("lower95HdiCI","upper95HdiCI")])))),")")
  
  estTable$npv_postMedian[k]<-format(nsmall=4,round(digits=4,npvEstimates$median[k]))
  estTable$npv_MAP[k]<-format(nsmall=4,round(digits=4,npvEstimates$MaP[k]))
  estTable$npv_HDI[k]<-paste(sep="","(",paste(collapse=",",format(nsmall=4,round(digits=4,as.matrix(npvEstimates[k,c("lower95HdiCI","upper95HdiCI")])))),")")
}

write.table(estTable,sep="\t",row.names=F,col.names=T,quote=F,file=paste(sep="",outPrefix,"_sensSpecPpvNpv_Table.tab"))

# PPV vs. NPV plot
plotLims<-c(0.2,1,0.6,1)
#pdf(paste(sep="",outPrefix,"_ppvVSnpv.pdf"),width=12,height=8)
png(paste(sep="",outPrefix,"_ppvVSnpv.png"),width=12,height=8,units="in",res=600)
plot(xlim=plotLims[1:2],ylim=plotLims[3:4],main="",xlab="positive predictive value",ylab="negative predictive value",0:1,0:1,col=cols,pch=20,type="n",asp=1)
addGrid()
par(xpd=T)

testsDens<-vector("list",nTests); names(testsDens)<-testNames
for(k in 1:nTests){
  testsDens[[k]]<-kde2d(fullMCMCtable[,nc+k],fullMCMCtable[,nc+nTests+k],n=250)
}
for(k in 1:nTests){
  paletteCust<-function(n){
    addAlpha(colorRampPalette(c("white",cols[k],cols[k],darkenCol(cols[k],factor=1.1)))(n),ALPHA=seq(0,1,length=n))
  }
  
  .filled.contour(testsDens[[k]]$x,testsDens[[k]]$y,testsDens[[k]]$z,col=paletteCust(50),levels=seq(0,max(testsDens[[k]]$z),length=50))
  contour(testsDens[[k]],col=addAlpha(colorRampPalette(c("white",cols[k]))(25),ALPHA=seq(0.2,0.8,length=25)),add=T,drawlabels=F,levels=seq(8,max(testsDens[[k]]$z),length=25))
}

addWhisker(low=ppvEstimates$lower95HdiCI,high=ppvEstimates$upper95HdiCI,anchor=npvEstimates$MaP,vertical=F,col=cols,whiskWidth=0.005/((1-0.35)/(1-0.85)))
addWhisker(low=npvEstimates$lower95HdiCI,high=npvEstimates$upper95HdiCI,anchor=ppvEstimates$MaP,vertical=T,col=cols)

points(ppvEstimates$MaP,npvEstimates$MaP,col=cols,pch=19,cex=1.5)

legend(x="bottomleft",pch=c(rep(19,5)),col=cols,bty="n",legend=testNames)
dev.off()


save(list=ls(),file=paste(sep="",outPrefix,"_summaryMetrics.RData"))

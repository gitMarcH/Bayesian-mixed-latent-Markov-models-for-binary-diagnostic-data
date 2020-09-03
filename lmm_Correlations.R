library(LMest)
library(RColorBrewer)
library(tidyverse)
library(corrplot)
library(rmcorr)

rm(list=ls())
setwd("~/../work/MLW_LSTM/")


##--------------------
## read data
##

tt<-read.csv("GordonMelita_StatsSupportAngeziwe/data/TTRInvASensitivity20170724_corrected.csv")
tt<-tt[!is.na(tt$Sample.ID) & !is.na(tt$Visit..),]


##--------------------
## read recruitment dates and compute actual timepoints (in months)
##

dates<-read.csv("GordonMelita_StatsSupportAngeziwe/data/salexpoLIMSDataSetComplete.csv")
dates$monthYear<-gsub("^[0-9]+/","",dates$DATE) # dates run from 08/2013 to 12/2014
datesDictionary<-data.frame(monthYear=c(paste(sep="/",c("08","09","10","11","12"),"2013"),paste(sep="/",c(paste(sep="","0",1:9),paste(sep="","1",0:2)),"2014")),t=1:17)
dates$t<-datesDictionary$t[match(dates$monthYear,datesDictionary$monthYear)]
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
## recode test results in a better, cleaner way
##

tt2<-tt
tt2$TTR<-ifelse(tolower(tt2$TTR)=="pos",1,ifelse(tolower(tt2$TTR)=="neg" | tolower(tt2$TTR)=="ind",0,NA))
tt2$InvA<-ifelse(tolower(tt2$InvA)=="pos",1,ifelse(tolower(tt2$InvA)=="neg" | tolower(tt2$InvA)=="ind",0,NA))
tt2$Culture<-ifelse(tolower(tt2$Culture)=="pos",1,ifelse(tolower(tt2$Culture)=="neg",0,NA))
tt2$TAC.TTR<-ifelse(tolower(tt2$TAC.TTR)=="pos",1,ifelse(tolower(tt2$TAC.TTR)=="neg" | tolower(tt2$TAC.TTR)=="ind",0,NA))
tt2$TAC.InvA<-ifelse(tolower(tt2$TAC.InvA)=="pos",1,ifelse(tolower(tt2$TAC.InvA)=="neg" | tolower(tt2$TAC.InvA)=="ind",0,NA))


##--------------------
## reformat data for LMest and fit basic LMM (no covariates adjustment) model
##

ttt<-long2wide(data=tt2,nameid="Sample.ID",namet="t",coly=c("TTR","InvA","Culture","TAC.TTR","TAC.InvA"),colx=NULL,full=NA)

##--------------------
## load Ct values data
##

ctSingleAssay<-read.csv("GordonMelita_StatsSupportAngeziwe/data/TTR & InvA master file Ct for correlation.csv",na.strings=c("Undetermined"))
ctTAC<-read.csv("GordonMelita_StatsSupportAngeziwe/data/TAC Results_4Marc TAC TTR TAC InvA Ct For Correlation.csv",na.strings=c("0","Undetermined"))
ctSingleAssay[is.na(ctSingleAssay)]<-40
ctTAC[is.na(ctTAC)]<-40
# check with Ange that all NA/0/Undetermined values can be considered to have Ct>40
# check with Ange the TAC-InvA Ct value for sample p55v11: 55.037??

##--------------------
##  combine culture, Ct values for TTR, InvA, TAC-TTR, TAC-InvA, pos/neg for TTR, InvA, TAC-TTR, TAC-InvA 
##

tt2$SID<-paste(sep="","p",formatC(tt2$Sample.ID,width=2,format="d",flag="0"),"v",formatC(tt2$Visit..,width=2,format="d",flag="0"))
ctSingleAssay$SID<-paste(sep="","p",formatC(ctSingleAssay$Sample.ID,width=2,format="d",flag="0"),"v",formatC(ctSingleAssay$Visit..,width=2,format="d",flag="0"))
ctTAC$SID<-paste(sep="","p",formatC(ctTAC$Participant.ID,width=2,format="d",flag="0"),"v",formatC(ctTAC$Visit,width=2,format="d",flag="0"))

uniqueSID<-sort(unique(c(tt2$SID,ctSingleAssay$SID,ctTAC$SID)))

mergedDat<-cbind(
  NA, NA, NA,
  tt2[match(uniqueSID,tt2$SID),c("SID","Culture","TTR","InvA","TAC.TTR","TAC.InvA")],
  ctSingleAssay[match(uniqueSID,ctSingleAssay$SID),c("SID","TTR.Ct.Mean","InvA.Ct.Mean")],
  ctTAC[match(uniqueSID,ctTAC$SID),c("SID","Ct_TTR","Ct_InvA")]
)

colnames(mergedDat)<-c(
  "Participant.ID","Visit","Sample.ID",
  "Sample.ID.binData","Culture_bin","TTR_bin","InvA_bin","TAC.TTR_bin","TAC.InvaA_bin",
  "Sample.ID.ctSingle","TTR_ct","InvA_ct",
  "Sample.ID.ctTAC","TAC.TTR_ct","TAC.InvA_ct"
)

mergedDat$Sample.ID<-ifelse(!is.na(mergedDat$Sample.ID.binData),mergedDat$Sample.ID.binData,ifelse(!is.na(mergedDat$Sample.ID.ctSingle),mergedDat$Sample.ID.ctSingle,mergedDat$Sample.ID.ctTAC))
mergedDat$Participant.ID<-substr(mergedDat$Sample.ID,1,3)
mergedDat$Visit<-substr(mergedDat$Sample.ID,4,6)

mergedDat<-mergedDat[order(mergedDat$Sample.ID),]


##--------------------
#-- naive way of doing things: ignore right-censoring (use Ct=40 for all 0 detections), ignore longitudinal nature of data
##

cor(mergedDat[,c("Culture_bin","TTR_ct","InvA_ct","TAC.TTR_ct","TAC.InvA_ct")],use="pairwise.complete.obs")

##--------------------
#-- less naive way of doing things: i) convert data to ranks to deal with right-censoring at Ct=40, ii) use rmcorr
##

tests<-c("Culture_bin","TTR_ct","InvA_ct","TAC.TTR_ct","TAC.InvA_ct")
for(var in tests[-1]){
  mergedDat[,paste(sep="_",var,"rank")]<-rank(na.last="keep",ties.method="max",-mergedDat[,var]) # taking rank of negative Ct values (and with ties.method="max" rather than "min") so that we get positive correlations
}

rmCensCorMat<-matrix(nrow=5,ncol=5)
colnames(rmCensCorMat)<-tests
rownames(rmCensCorMat)<-tests
diag(rmCensCorMat)<-1
for(j in 1:length(colnames(rmCensCorMat))){
  for(i in 1:length(rownames(rmCensCorMat))){
    if(i!=1){
      var1<-paste(sep="",rownames(rmCensCorMat)[i],"_rank")
    }else{
      var1<-"Culture_bin"
    }
    if(j!=1){
      var2<-paste(sep="",colnames(rmCensCorMat)[j],"_rank")
    }else{
      var2<-"Culture_bin"
    }
    if(i!=j){rmCensCorMat[i,j]<-rmcorr(participant=factor(Participant.ID),measure1=var1,measure2=var2,dataset=mergedDat)$r}
  }
}

ppvMat<-matrix(rep(0,5*5),nrow=5,ncol=5)
colnames(ppvMat)<-gsub(pattern="\\.",replacement="-",colnames(tt2)[c(5,3,4,6,7)])
rownames(ppvMat)<-gsub(pattern="\\.",replacement="-",colnames(tt2)[c(5,3,4,6,7)])
ppvMatCount<-ppvMat
npvMat<-ppvMat
npvMatCount<-ppvMat

ppvMatList<-list()
npvMatList<-list()
for(t in unique(tt$t)){
  ppvMatList[[t]]<-matrix(nrow=5,ncol=5)
  colnames(ppvMatList[[t]])<-gsub(pattern="\\.",replacement="-",colnames(tt2)[c(5,3,4,6,7)])
  rownames(ppvMatList[[t]])<-gsub(pattern="\\.",replacement="-",colnames(tt2)[c(5,3,4,6,7)])
  npvMatList[[t]]<-ppvMatList[[t]]
  
  for(j in 1:length(colnames(ppvMatList[[t]]))){
    for(i in 1:length(rownames(ppvMatList[[t]]))){
      var1<-gsub(pattern="-",replacement="\\.",rownames(ppvMatList[[t]])[i])
      var2<-gsub(pattern="-",replacement="\\.",colnames(ppvMatList[[t]])[j])
      ppvMatList[[t]][i,j]<-sum(tt2$t==t & tt2[,var1]==1 & tt2[,var2]==1)/sum(tt2$t==t & tt2[,var1]==1)
      npvMatList[[t]][i,j]<-sum(tt2$t==t & tt2[,var1]==0 & tt2[,var2]==0)/sum(tt2$t==t & tt2[,var1]==0)
      
      if(!is.nan(ppvMatList[[t]][i,j])){
        ppvMat[i,j]<-ppvMat[i,j]+ppvMatList[[t]][i,j]*sum(tt2$t==t)
        ppvMatCount[i,j]<-ppvMatCount[i,j]+sum(tt2$t==t)
      }
      
      if(!is.nan(npvMatList[[t]][i,j])){
        npvMat[i,j]<-npvMat[i,j]+npvMatList[[t]][i,j]*sum(tt2$t==t)
        npvMatCount[i,j]<-npvMatCount[i,j]+sum(tt2$t==t)
      }
    }
  }

}
ppvMat<-ppvMat/ppvMatCount
npvMat<-npvMat/npvMatCount

colnames(rmCensCorMat)<-gsub(pattern="_bin",replacement="",gsub(pattern="_ct",replacement=" (Ct)",gsub(pattern="\\.",replacement="-",colnames(rmCensCorMat))))
rownames(rmCensCorMat)<-gsub(pattern="_bin",replacement="",gsub(pattern="_ct",replacement=" (Ct)",gsub(pattern="\\.",replacement="-",rownames(rmCensCorMat))))

png("GordonMelita_StatsSupportAngeziwe/output/lmmFits_CorrelationBetweenTestsRepeatedMeasurements_Ct.png",width=10,height=10,units="in",res=600)
corrplot(rmCensCorMat,is.corr=F,type="lower",cl.lim=c(0,1),method="number",tl.pos="lt",cl.pos="r",title="          repeated measures rank correlation",mar=c(0,0,1,0),cex.main=1.5,tl.cex=1.5,cl.cex=1.5,number.cex=1.5)
corrplot(rmCensCorMat,is.corr=F,type="upper",cl.lim=c(0,1),tl.pos="n",add=T,cl.pos="n")
dev.off()

png("GordonMelita_StatsSupportAngeziwe/output/lmmFits_PositiveDiagnosticConcordance.png",width=10,height=10,units="in",res=600)
corrplot(ppvMat,is.corr=F,type="full",cl.lim=c(0,1),method="circle",tl.pos="lt",cl.pos="r",title="          concordance of positive diagnosis",mar=c(0,0,1,0),cex.main=1.5,tl.cex=1.5,cl.cex=1.5)
for(i in 1:nrow(ppvMat)){
  for(j in 1:ncol(ppvMat)){
    text(adj=c(0.5,0.5),x=j,y=nrow(ppvMat)-i+1,font=2,col="gray40",labels=round(digits=2,ppvMat[i,j]),cex=1.25)
  }
}
dev.off()

png("GordonMelita_StatsSupportAngeziwe/output/lmmFits_NegativeDiagnosticConcordance.png",width=10,height=10,units="in",res=600)
corrplot(npvMat,is.corr=F,type="full",cl.lim=c(0,1),method="circle",tl.pos="lt",cl.pos="r",title="          concordance of negative diagnosis",mar=c(0,0,1,0),cex.main=1.5,tl.cex=1.5,cl.cex=1.5)
for(i in 1:nrow(npvMat)){
  for(j in 1:ncol(npvMat)){
    text(adj=c(0.5,0.5),x=j,y=nrow(npvMat)-i+1,font=2,col="gray40",labels=round(digits=2,npvMat[i,j]),cex=1.25)
  }
}
dev.off()

png("GordonMelita_StatsSupportAngeziwe/output/lmmFits_CorrelationAndConcordance.png",width=30,height=10,units="in",res=600)
par(mfrow=c(1,3))

corrplot(rmCensCorMat,is.corr=F,type="lower",cl.lim=c(0,1),method="number",tl.pos="lt",cl.pos="r",title="          repeated measures rank correlation",mar=c(0,0,1,0),cex.main=2,tl.cex=2,cl.cex=2,number.cex=2)
corrplot(rmCensCorMat,is.corr=F,type="upper",cl.lim=c(0,1),tl.pos="n",add=T,cl.pos="n")

corrplot(ppvMat,is.corr=F,type="full",cl.lim=c(0,1),method="circle",tl.pos="lt",cl.pos="r",title="          concordance of positive diagnosis",mar=c(0,0,1,0),cex.main=2,tl.cex=2,cl.cex=2)
for(i in 1:nrow(ppvMat)){
  for(j in 1:ncol(ppvMat)){
    text(adj=c(0.5,0.5),x=j,y=nrow(ppvMat)-i+1,font=2,col="gray40",labels=round(digits=2,ppvMat[i,j]),cex=2)
  }
}

corrplot(npvMat,is.corr=F,type="full",cl.lim=c(0,1),method="circle",tl.pos="lt",cl.pos="r",title="          concordance of negative diagnosis",mar=c(0,0,1,0),cex.main=2,tl.cex=2,cl.cex=2)
for(i in 1:nrow(npvMat)){
  for(j in 1:ncol(npvMat)){
    text(adj=c(0.5,0.5),x=j,y=nrow(npvMat)-i+1,font=2,col="gray40",labels=round(digits=2,npvMat[i,j]),cex=2)
  }
}
dev.off()




#########
## OLD ##
#########

##--------------------
#-- look at agreements between diagnostic tests (binary results only)
##
# agreeMat<-matrix(nrow=5,ncol=5)
# tests<-colnames(agreeMat)
# colnames(agreeMat)<-gsub(pattern="\\.",replacement="-",colnames(tt2)[3:7])
# rownames(agreeMat)<-gsub(pattern="\\.",replacement="-",colnames(tt2)[3:7])
# corMat<-agreeMat
# ppvMat<-agreeMat
# diag(corMat)<-1
# for(j in 1:length(colnames(agreeMat))){
#   for(i in 1:length(rownames(agreeMat))){
#     var1<-gsub(pattern="-",replacement="\\.",rownames(agreeMat)[i])
#     var2<-gsub(pattern="-",replacement="\\.",colnames(agreeMat)[j])
#     agreeMat[i,j]<-sum(tt2[,var1]==tt2[,var2])/nrow(tt2)
#     ppvMat[i,j]<-sum(tt2[,var1]==1 & tt2[,var2]==1)/sum(tt2[,var1]==1)
#     if(i!=j){corMat[i,j]<-rmcorr(participant=factor(Sample.ID),measure1=var1,measure2=var2,dataset=tt2)$r}
#   }
# }
# 
# pdf("GordonMelita_StatsSupportAngeziwe/lmmFits_overallAgreementBetweenTests.pdf",width=10,height=10)
# corrplot(agreeMat,is.corr=F,type="lower",cl.lim=c(0.85,1),method="number",tl.pos="lt",cl.pos="r")
# corrplot(agreeMat,is.corr=F,type="upper",cl.lim=c(0.85,1),tl.pos="n",add=T)
# dev.off()
# 
# pdf("GordonMelita_StatsSupportAngeziwe/lmmFits_CorrelationBetweenTestsRepeatedMeasurements.pdf",width=10,height=10)
# corrplot(corMat,is.corr=F,type="lower",cl.lim=c(0,1),method="number",tl.pos="lt",cl.pos="r")
# corrplot(corMat,is.corr=F,type="upper",cl.lim=c(0,1),tl.pos="n",add=T)
# dev.off()
# 
# png("GordonMelita_StatsSupportAngeziwe/lmmFits_overallAgreementBetweenTests.png",width=10,height=10,units = "in", res = 300)
# corrplot(agreeMat,is.corr=F,type="lower",cl.lim=c(0.85,1),method="number",tl.pos="lt",cl.pos="r")
# corrplot(agreeMat,is.corr=F,type="upper",cl.lim=c(0.85,1),tl.pos="n",add=T)
# dev.off()
# 
# png("GordonMelita_StatsSupportAngeziwe/lmmFits_CorrelationBetweenTestsRepeatedMeasurements.png",width=10,height=10,units = "in", res = 300)
# corrplot(corMat,is.corr=F,type="lower",cl.lim=c(0,1),method="number",tl.pos="lt",cl.pos="r")
# corrplot(corMat,is.corr=F,type="upper",cl.lim=c(0,1),tl.pos="n",add=T)
# dev.off()

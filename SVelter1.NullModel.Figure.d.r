#!R
Args <- commandArgs()
NullModel=Args[6]
outputname=Args[7]
boxcolor=Args[8]
linecolor=Args[9]
outputname2=Args[10]
outputFileType=strsplit(outputname,'[.]')[[1]][length(strsplit(outputname,'[.]')[[1]])]
dataNull=read.table(file=NullModel,header=F)
dataNull=dataNull[order(dataNull[,1],decreasing=F),]

if(sum(dataNull[,2])>2^31-1){
shrink=int(sum(dataNull[,2])/2^31-1)+3
dataNull[,2]=dataNull[,2]/shrink
}

MedianData=median(dataNull[,1])
for(i in nrow(dataNull):1){
    if(dataNull[i,1]>3*MedianData){
        dataNull=dataNull[-i,]
    }
}

StatType=strsplit(strsplit(NullModel,'/')[[1]][length(strsplit(NullModel,'/')[[1]])],'[.]')[[1]][1]
if (StatType=='RDNull'){
data_RD=c()
for(i in 1:nrow(dataNull)){
    data_RD=c(data_RD,rep(dataNull[i,1],dataNull[i,2]))
}
RD_Mean=mean(data_RD)
RD_Median=median(data_RD)
RD_STD=sd(data_RD)
data_RD=data_RD[data_RD<2*RD_Mean]
p=RD_Mean/RD_STD^2
r=RD_Mean*p/(1-p)
xrange=range(0,2*RD_Mean)
yrange=range(0,max(dataNull[,2]/sum(dataNull[,2])))
if(outputFileType=='pdf'){
    pdf(outputname)
}
if(outputFileType=='jpg'){
    jpeg(outputname)
}
if(outputFileType=='png'){
    png(outputname)
}
par(mfrow=c(1,2))
breakNumber=as.integer(nrow(dataNull)/2)
hist(data_RD,freq=F,col=boxcolor,breaks=breakNumber,xlab='Read Depth',ylab='Frequency',xlim=xrange,main='Fitted By Negative Binomial Distribution')
lines(dnbinom(seq(range(data_RD)[1],range(data_RD)[2],1),r,p),type='l',lwd=2,col=linecolor)
hist(data_RD,freq=F,col=boxcolor,breaks=breakNumber,xlab='Read Depth',ylab='Frequency',xlim=xrange,main='Fitted By Normal Distribution')
lines(dnorm(seq(range(data_RD)[1],range(data_RD)[2],1),RD_Mean,RD_STD),type='l',lwd=2,col=linecolor)
dev.off()
StatMatrix=data.frame(Mean=RD_Mean,Median=RD_Median,STD=RD_STD)
write.table(file=outputname2,StatMatrix,quote=F,row.name=F)
}


if (StatType=='ILNull'){
library(mixtools)
library(fitdistrplus)
data_IL=c()
for(i in 1:nrow(dataNull)){
    data_IL=c(data_IL,rep(dataNull[i,1],dataNull[i,2]))
}
IL_Mean=mean(data_IL)
IL_Median=median(data_IL)
IL_STD=sd(data_IL)
data_IL=data_IL[data_IL<2*IL_Mean]
data_IL=data_IL[data_IL>0.1*IL_Mean]
mixmdl = normalmixEM(data_IL)
normdl=fitdist(data_IL,"norm",method="mge",gof="CvM")
IL_Normal_Mean=mean(normdl$estimate[1])
IL_Normal_STD=mean(normdl$estimate[2])
IL_STD=IL_Normal_STD
if(outputFileType=='pdf'){
    pdf(outputname)
}
if(outputFileType=='jpg'){
    jpeg(outputname)
}
if(outputFileType=='png'){
    png(outputname)
}
xrange=range(0,2*IL_Mean)
yrange=range(0,max(dataNull[,2]/sum(dataNull[,2])))
par(mfrow=c(1,2))
hist(data_IL,freq=F,col=boxcolor,breaks=as.integer(nrow(dataNull)/5),xlab='Insert Length',ylab='Frequency',xlim=xrange,main='Fitted By Bimodal Distribution')
xin=c(1:dataNull[nrow(dataNull),1])
yin=mixmdl$lambda[1] * dnorm(xin,mean=mixmdl$mu[1],sd=mixmdl$sigma[1])+
mixmdl$lambda[2] * dnorm(xin,mean=mixmdl$mu[2],sd=mixmdl$sigma[2])
lines(xin,yin,pch='l',col=linecolor)
hist(data_IL,freq=F,col=boxcolor,breaks=as.integer(nrow(dataNull)/5),xlab='Insert Length',ylab='Frequency',xlim=xrange,main='Fitted By Normal Distribution')
lines(dnorm(seq(range(data_IL)[1],range(data_IL)[2],1),IL_Normal_Mean,IL_Normal_STD),type='l',lwd=2,col=linecolor)
dev.off()
StatMatrix1=data.frame(Mean=IL_Mean, Median=IL_Median, STD=IL_STD)
StatMatrix2=data.frame(Bimodal1=mixmdl$lambda[1],Mean1=mixmdl$mu[1], STD1=mixmdl$sigma[1])
StatMatrix3=data.frame(Bimodal2=mixmdl$lambda[2],Mean2=mixmdl$mu[2], STD2=mixmdl$sigma[2])
StatMatrix4=data.frame(Normal=1, Mean=IL_Normal_Mean, STD=IL_Normal_STD)
StatMatrix=rbind(colnames(StatMatrix1),c(StatMatrix1[1,]),
colnames(StatMatrix2),c(StatMatrix2[1,]),
colnames(StatMatrix3),c(StatMatrix3[1,]),
colnames(StatMatrix4),c(StatMatrix4[1,]))
write.table(file=outputname2,StatMatrix,quote=F,col.name=F,row.name=F)
}


if (StatType=='TBNull'){
library(mixtools)
library(fitdistrplus)
data_TB=c()
for(i in 1:nrow(dataNull)){
    data_TB=c(data_TB,rep(dataNull[i,1],dataNull[i,2]))
}
TB_Mean=mean(data_TB)
TB_Median=median(data_TB)
TB_STD=sd(data_TB)
data_TB=data_TB[data_TB<2*TB_Mean]
data_TB=data_TB[data_TB>TB_Mean/10]
mixmdl = normalmixEM(data_TB)
normdl=fitdist(data_TB,"norm",method="mge",gof="CvM")
TB_Normal_Mean=mean(normdl$estimate[1])
TB_Normal_STD=mean(normdl$estimate[2])
if(outputFileType=='pdf'){
    pdf(outputname)
}
if(outputFileType=='jpg'){
    jpeg(outputname)
}
if(outputFileType=='png'){
    png(outputname)
}
xrange=range(0,2*TB_Mean)
yrange=range(0,max(dataNull[,2]/sum(dataNull[,2])))
par(mfrow=c(1,2))
hist(data_TB,freq=F,col=boxcolor,breaks=as.integer(nrow(dataNull)/3),xlab='Read Pairs Through Break Point',ylab='Frequency',xlim=xrange,main='Fitted By Bimodal Distribution')
xin=c(1:dataNull[nrow(dataNull),1])
yin=mixmdl$lambda[1] * dnorm(xin,mean=mixmdl$mu[1],sd=mixmdl$sigma[1])+
mixmdl$lambda[2] * dnorm(xin,mean=mixmdl$mu[2],sd=mixmdl$sigma[2])
lines(xin,yin,pch='l',col=linecolor)
hist(data_TB,freq=F,col=boxcolor,breaks=as.integer(nrow(dataNull)/3),xlab='Read Pairs Through Break Point',ylab='Frequency',xlim=xrange,main='Fitted By Normal Distribution')
lines(dnorm(seq(range(data_TB)[1],range(data_TB)[2],1),TB_Normal_Mean,TB_Normal_STD),type='l',lwd=2,col=linecolor)
dev.off()
StatMatrix1=data.frame(Mean=TB_Mean, Median=TB_Median, STD=TB_STD)
StatMatrix2=data.frame(Bimodal1=mixmdl$lambda[1],Mean1=mixmdl$mu[1], STD1=mixmdl$sigma[1])
StatMatrix3=data.frame(Bimodal2=mixmdl$lambda[2],Mean2=mixmdl$mu[2], STD2=mixmdl$sigma[2])
StatMatrix4=data.frame(Normal=1, Mean=TB_Normal_Mean, STD=TB_Normal_STD)
StatMatrix=rbind(colnames(StatMatrix1),c(StatMatrix1[1,]),
colnames(StatMatrix2),c(StatMatrix2[1,]),
colnames(StatMatrix3),c(StatMatrix3[1,]),
colnames(StatMatrix4),c(StatMatrix4[1,]))
write.table(file=outputname2,StatMatrix,quote=F,col.name=F,row.name=F)
}



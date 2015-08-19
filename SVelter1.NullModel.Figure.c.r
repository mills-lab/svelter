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

if (StatType=='ILNull'){
data_IL=c()
for(i in 1:nrow(dataNull)){
    data_IL=c(data_IL,rep(dataNull[i,1],dataNull[i,2]))
}
if(outputFileType=='pdf'){
    pdf(outputname)
}
if(outputFileType=='jpg'){
    jpeg(outputname)
}
if(outputFileType=='png'){
    png(outputname)
}
xrange=range(0,dataNull[,1])
yrange=range(0,max(dataNull[,2]/sum(dataNull[,2])))
hist(data_IL,freq=F,col=boxcolor,breaks=as.integer(nrow(dataNull)/5),xlab='Insert Length',ylab='Frequency',xlim=xrange,main='Null Distribution of Insert Length')
dev.off()
}


if (StatType=='TBNull'){
data_IL=c()
for(i in 1:nrow(dataNull)){
    data_IL=c(data_IL,rep(dataNull[i,1],dataNull[i,2]))
}
if(outputFileType=='pdf'){
    pdf(outputname)
}
if(outputFileType=='jpg'){
    jpeg(outputname)
}
if(outputFileType=='png'){
    png(outputname)
}
xrange=range(0,dataNull[,1])
yrange=range(0,max(dataNull[,2]/sum(dataNull[,2])))
hist(data_IL,freq=F,col=boxcolor,breaks=as.integer(nrow(dataNull)/3),xlab='Read Pairs Through Break Point',ylab='Frequency',xlim=xrange,main='Null Distribution of Number of Read Pairs Going Through Break Poins')
dev.off()
}

if (StatType=='RDNull'){
data_RD=c()
for(i in 1:nrow(dataNull)){
    data_RD=c(data_RD,rep(dataNull[i,1],dataNull[i,2]))
}
RD_Mean=mean(data_RD)
RD_Median=median(data_RD)
RD_STD=sd(data_RD)
p=RD_Mean/RD_STD^2
r=RD_Mean*p/(1-p)
xrange=range(0,dataNull[,1])
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
hist(data_RD,freq=F,col=boxcolor,breaks=as.integer(nrow(dataNull)/3),xlab='Read Depth',ylab='Frequency',xlim=xrange, main='Null Distribution of Read Length')
lines(dnbinom(seq(range(data_RD)[1],range(data_RD)[2],1),r,p),type='l',lwd=2,col=linecolor)
dev.off()
}



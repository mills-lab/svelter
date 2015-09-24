#!/usr/bin/env R

#!R
#Usage: Rscript Code1.NullModel.Figure.r 
Args <- commandArgs()
NullModel=Args[6]
outputname=Args[7]
boxcolor=Args[8]
outputname2=Args[9]
outputFileType=strsplit(outputname,'[.]')[[1]][length(strsplit(outputname,'[.]')[[1]])]
dataNull=read.table(file=NullModel,header=F)
dataNull=dataNull[order(dataNull[,1],decreasing=F),]

datatotal=c(rep(dataNull[1,1],dataNull[1,2]))
if(nrow(dataNull)>1){
for(i in 2:nrow(dataNull)){
datatotal=c(datatotal, c(rep(dataNull[i,1],dataNull[i,2])))
}}

datatotal2=c(rep(dataNull[2,1],dataNull[2,2]))
if(nrow(dataNull)>2){
for(i in 3:nrow(dataNull)){
datatotal2=c(datatotal2, c(rep(dataNull[i,1],dataNull[i,2])))
}}

if(strsplit(strsplit(NullModel,'/')[[1]][length(strsplit(NullModel,'/')[[1]])],'[.]')[[1]][1]=="InsertLenNull"){
xlabn='Number of Read Pairs with Aberrant Insert Length'
ylabn='frequency'
}

if(strsplit(strsplit(NullModel,'/')[[1]][length(strsplit(NullModel,'/')[[1]])],'[.]')[[1]][1]=="DirectionNull"){
xlabn='Number of Read Pairs with Aberrant Direction'
ylabn='frequency'
}

if(strsplit(strsplit(NullModel,'/')[[1]][length(strsplit(NullModel,'/')[[1]])],'[.]')[[1]][1]=="SplitNull"){
xlabn='Number of Split Reads'
ylabn='frequency'
}

if(outputFileType=='pdf'){
pdf(file=outputname)
hist(datatotal,col=boxcolor,xlim=range(dataNull[,1]),xlab=xlabn,ylab=ylabn,freq=T, main='histogram')
dev.off()
pdf(file=outputname2)
hist(datatotal2,col=boxcolor,xlim=range(dataNull[,1]),xlab=xlabn,ylab=ylabn,freq=T, main='histogram')
dev.off()
}

if(outputFileType=='jpg'){
jpeg(file=outputname)
hist(datatotal,col=boxcolor,xlim=range(dataNull[,1]),xlab=xlabn,ylab=ylabn,freq=T, main='histogram')
dev.off()
jpeg(file=outputname2)
hist(datatotal2,col=boxcolor,xlim=range(dataNull[,1]),xlab=xlabn,ylab=ylabn,freq=T, main='histogram')
dev.off()
}

if(outputFileType=='png'){
png(file=outputname)
hist(datatotal,col=boxcolor,xlim=range(dataNull[,1]),xlab=xlabn,ylab=ylabn,freq=T, main='histogram')
dev.off()
png(file=outputname2)
hist(datatotal2,col=boxcolor,xlim=range(dataNull[,1]),xlab=xlabn,ylab=ylabn,freq=T, main='histogram')
dev.off()
}


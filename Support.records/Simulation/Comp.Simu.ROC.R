setwd('~/Box Sync/git-repo/svelter.xuefang.github/svelter/Support.records/Simulated/')
data1=read.table('./comp.homo/Integrated.CSV.stat')
par(mar=c(3,3,3,3))
par(mfrow=c(2,2))
SV_type_names=c('INV_DUP' ,'DEL_INV' ,'DEL_DUP', 'DEL_INV_DUP' )
SV_name=sort(unique(data1[,3]))
for(k0 in unique(data1[,3])){
  pdf(paste('./Comp.homo.CSV.V2.',k0,'.Stat.pdf',sep=''))
  #par(mfrow=c(3,2))
  par(fig=c(0.1,0.7,0.4,1))
  plot(c(0,1),c(0,1),mgp=c(1.5,0.5,0),xlab='False Discovery Rate',ylab='Sensitivity',cex.lab=0.8,type='n')
  title(main=SV_type_names[match(k0,SV_name)],cex.main=0.9)
  cols=c('red','blue','seagreen','purple','orange','black')
  algorithms=c('SVelter','Delly','Lumpy','Pindel','erds')
  for(k1 in unique(data1[,1])){
    color=cols[match(k1,algorithms)]
    color2=paste(color,c('','1','2','3','4'),sep='')
    temp=data1[data1[,1]==k1 & data1[,3]==k0,]
    temp=temp[order(temp[,2]),]
    points(1-temp[,5],temp[,4],col=color2,bg=color2,pch=21,cex=seq(0.6,1.4,by=0.2))
  }
  legend(c(0.8,1),c(.75,1),algorithms,col=cols,cex=0.8,pch=20)
  fig_num=c(c(.7,1,.7,1),c(.7,1,.4,.7),c(.7,1,.1,.4),c(.4,.7,.1,.4),c(.1,.4,.1,.4))
  rec=-1
  for(k1 in unique(data1[,1])){
    rec=rec+1
    par(mar=rep(1.5,4))
    par(fig=fig_num[c((rec*4+1):(rec*4+4))],new=TRUE)
    color=cols[match(k1,algorithms)]
    color2=paste(color,c('','1','2','3','4'),sep='')
    temp=data1[data1[,1]==k1 & data1[,3]==k0,]
    temp=temp[order(temp[,2]),]
    plot(range(1-temp[,5]),range(temp[,4]),mgp=c(1.5,0.5,0),xlab='False Discovery Rate',ylab='Sensitivity',cex.lab=0.8,type='n')
    points(1-temp[,5],temp[,4],col=color2,bg=color2,pch=21,cex=seq(0.6,1.4,by=0.2))
  }
  dev.off()}

data1=read.table('./comp.het/Integrated.CSV.stat')
par(mar=c(3,3,3,3))
par(mfrow=c(2,2))
SV_type_names=c('INV_DUP' ,'DEL_INV' ,'DEL_DUP', 'DEL_INV_DUP' )
SV_name=sort(unique(data1[,3]))
for(k0 in unique(data1[,3])){
  pdf(paste('./Comp.het.CSV.V2.',k0,'.Stat.pdf',sep=''))
  #par(mfrow=c(3,2))
  par(fig=c(0.1,0.7,0.4,1))
  plot(c(0,1),c(0,1),mgp=c(1.5,0.5,0),xlab='False Discovery Rate',ylab='Sensitivity',cex.lab=0.8,type='n')
  title(main=SV_type_names[match(k0,SV_name)],cex.main=0.9)
  cols=c('red','blue','seagreen','purple','orange','black')
  algorithms=c('SVelter','Delly','Lumpy','Pindel','erds')
  for(k1 in unique(data1[,1])){
    color=cols[match(k1,algorithms)]
    color2=paste(color,c('','1','2','3','4'),sep='')
    temp=data1[data1[,1]==k1 & data1[,3]==k0,]
    temp=temp[order(temp[,2]),]
    points(1-temp[,5],temp[,4],col=color2,bg=color2,pch=21,cex=seq(0.6,1.4,by=0.2))
  }
  legend(c(0.8,1),c(.75,1),algorithms,col=cols,cex=0.8,pch=20)
  fig_num=c(c(.7,1,.7,1),c(.7,1,.4,.7),c(.7,1,.1,.4),c(.4,.7,.1,.4),c(.1,.4,.1,.4))
  rec=-1
  for(k1 in unique(data1[,1])){
    rec=rec+1
    par(mar=rep(1.5,4))
    par(fig=fig_num[c((rec*4+1):(rec*4+4))],new=TRUE)
    color=cols[match(k1,algorithms)]
    color2=paste(color,c('','1','2','3','4'),sep='')
    temp=data1[data1[,1]==k1 & data1[,3]==k0,]
    temp=temp[order(temp[,2]),]
    plot(range(1-temp[,5]),range(temp[,4]),mgp=c(1.5,0.5,0),xlab='False Discovery Rate',ylab='Sensitivity',cex.lab=0.8,type='n')
    points(1-temp[,5],temp[,4],col=color2,bg=color2,pch=21,cex=seq(0.6,1.4,by=0.2))
  }
  dev.off()}

data1=read.table('./comp.homo/Integrated.CSV.stat')
data2=read.table('./comp.het/Integrated.CSV.stat')
data1[,6]=data1[,6]+data2[,6]
data1[,7]=data1[,7]+data2[,7]
data1[,8]=data1[,8]+data2[,8]
data1[,4]=data1[,6]/data1[,7]
data1[,5]=data1[,6]/data1[,8]
par(mar=c(3,3,3,3))
par(mfrow=c(2,2))
SV_type_names=c('INV_DUP' ,'DEL_INV' ,'DEL_DUP', 'DEL_INV_DUP' )
SV_name=sort(unique(data1[,3]))
for(k0 in unique(data1[,3])){
  pdf(paste('./Comp.Combine.CSV.V2.',k0,'.Stat.pdf',sep=''))
  #par(mfrow=c(3,2))
  par(fig=c(0.1,0.7,0.4,1))
  plot(c(0,1),c(0,1),mgp=c(1.5,0.5,0),xlab='False Discovery Rate',ylab='Sensitivity',cex.lab=0.8,type='n')
  title(main=SV_type_names[match(k0,SV_name)],cex.main=0.9)
  cols=c('red','blue','seagreen','purple','orange','black')
  algorithms=c('SVelter','Delly','Lumpy','Pindel','erds')
  for(k1 in unique(data1[,1])){
    color=cols[match(k1,algorithms)]
    color2=paste(color,c('','1','2','3','4'),sep='')
    temp=data1[data1[,1]==k1 & data1[,3]==k0,]
    temp=temp[order(temp[,2]),]
    points(1-temp[,5],temp[,4],col=color2,bg=color2,pch=21,cex=seq(0.6,1.4,by=0.2))
  }
  legend(c(0.8,1),c(.75,1),algorithms,col=cols,cex=0.8,pch=20)
  fig_num=c(c(.7,1,.7,1),c(.7,1,.4,.7),c(.7,1,.1,.4),c(.4,.7,.1,.4),c(.1,.4,.1,.4))
  rec=-1
  for(k1 in unique(data1[,1])){
    rec=rec+1
    par(mar=rep(1.5,4))
    par(fig=fig_num[c((rec*4+1):(rec*4+4))],new=TRUE)
    color=cols[match(k1,algorithms)]
    color2=paste(color,c('','1','2','3','4'),sep='')
    temp=data1[data1[,1]==k1 & data1[,3]==k0,]
    temp=temp[order(temp[,2]),]
    plot(range(1-temp[,5]),range(temp[,4]),mgp=c(1.5,0.5,0),xlab='False Discovery Rate',ylab='Sensitivity',cex.lab=0.8,type='n')
    points(1-temp[,5],temp[,4],col=color2,bg=color2,pch=21,cex=seq(0.6,1.4,by=0.2))
  }
  dev.off()}



setwd('~/Box Sync/git-repo/svelter.xuefang.github/svelter/Support.records/Simulated/')

pdf('./Comp.homo.CSV.V1.Stat.pdf', useDingbats=FALSE)
data1=read.table('./comp.homo/Integrated.CSV.stat')
par(mar=c(3,3,3,3))
par(mfrow=c(2,2))
SV_type_names=c('INV_DUP' ,'DEL_INV' ,'DEL_DUP', 'DEL_INV_DUP' )
SV_name=sort(unique(data1[,3]))
for(k0 in unique(data1[,3])){
  plot(c(0,1),c(0,1),mgp=c(1.5,0.5,0),xlab='False Discovery Rate',ylab='Sensitivity',cex.lab=0.8,type='n')
  title(main=SV_type_names[match(k0,SV_name)],cex.main=0.9)
  cols=c('red','blue','seagreen','purple','orange','black')
  pch=22
  algorithms=c('SVelter','Delly','Lumpy','Pindel','erds')
  for(k1 in unique(data1[,1])){
    color=cols[match(k1,algorithms)]
    color2=paste(color,c('','1','2','3','4'),sep='')
      temp=data1[data1[,1]==k1 & data1[,3]==k0,]
    temp=temp[order(temp[,2]),]
    points(1-temp[,5],temp[,4],col=color,pch=c(21:25),bg=color)
  }}
dev.off()


pdf('./Comp.het.CSV.V1.Stat.pdf', useDingbats=FALSE)
data1=read.table('./comp.het/Integrated.CSV.stat')
par(mar=c(3,3,3,3))
par(mfrow=c(2,2))
SV_type_names=c('INV_DUP' ,'DEL_INV' ,'DEL_DUP', 'DEL_INV_DUP' )
SV_name=sort(unique(data1[,3]))
for(k0 in unique(data1[,3])){
  plot(c(0,1),c(0,1),mgp=c(1.5,0.5,0),xlab='False Discovery Rate',ylab='Sensitivity',cex.lab=0.8,type='n')
  title(main=SV_type_names[match(k0,SV_name)],cex.main=0.9)
  cols=c('red','blue','seagreen','purple','orange','black')
  pch=22
  algorithms=c('SVelter','Delly','Lumpy','Pindel','erds')
  for(k1 in unique(data1[,1])){
    color=cols[match(k1,algorithms)]
    color2=paste(color,c('','1','2','3','4'),sep='')
    temp=data1[data1[,1]==k1 & data1[,3]==k0,]
    temp=temp[order(temp[,2]),]
    points(1-temp[,5],temp[,4],col=color,pch=c(21:25),bg=color)
  }}

dev.off()

data1=read.table('./comp.homo/Integrated.CSV.stat')
data2=read.table('./comp.het/Integrated.CSV.stat')
data1[,6]=data1[,6]+data2[,6]
data1[,7]=data1[,7]+data2[,7]
data1[,8]=data1[,8]+data2[,8]
data1[,4]=data1[,6]/data1[,7]
data1[,5]=data1[,6]/data1[,8]
pdf('./Comp.Combine.CSV.V1.Stat.pdf', useDingbats=FALSE)
par(mar=c(3,3,3,3))
par(mfrow=c(2,2))
SV_type_names=c('INV_DUP' ,'DEL_INV' ,'DEL_DUP', 'DEL_INV_DUP' )
SV_name=sort(unique(data1[,3]))
for(k0 in unique(data1[,3])){
  plot(c(0,1),c(0,1),mgp=c(1.5,0.5,0),xlab='False Discovery Rate',ylab='Sensitivity',cex.lab=0.8,type='n')
  title(main=SV_type_names[match(k0,SV_name)],cex.main=0.9)
  cols=c('red','blue','seagreen','purple','orange','black')
  pch=22
  algorithms=c('SVelter','Delly','Lumpy','Pindel','erds')
  for(k1 in unique(data1[,1])){
    color=cols[match(k1,algorithms)]
    color2=paste(color,c('','1','2','3','4'),sep='')
    temp=data1[data1[,1]==k1 & data1[,3]==k0,]
    temp=temp[order(temp[,2]),]
    points(1-temp[,5],temp[,4],col=color,pch=c(21:25),bg=color)
  }}

dev.off()













setwd('~/Box Sync/git-repo/svelter.xuefang.github/svelter/Support.records/Simulated/')

pdf('./Comp.homo.CSV.V2.Stat.pdf', useDingbats=FALSE)
data1=read.table('./comp.homo/Integrated.CSV.stat')
par(mar=c(3,3,3,3))
par(mfrow=c(2,2))
SV_type_names=c('INV_DUP' ,'DEL_INV' ,'DEL_DUP', 'DEL_INV_DUP' )
SV_name=sort(unique(data1[,3]))
for(k0 in unique(data1[,3])){
  plot(c(0,1),c(0,1),mgp=c(1.5,0.5,0),xlab='False Discovery Rate',ylab='Sensitivity',cex.lab=0.8,type='n')
  title(main=SV_type_names[match(k0,SV_name)],cex.main=0.9)
  cols=c('red','blue','seagreen','purple','orange','black')
  pch=22
  algorithms=c('SVelter','Delly','Lumpy','Pindel','erds')
  for(k1 in unique(data1[,1])){
    color=cols[match(k1,algorithms)]
    color2=paste(color,c('','1','2','3','4'),sep='')
    temp=data1[data1[,1]==k1 & data1[,3]==k0,]
    temp=temp[order(temp[,2]),]
    points(1-temp[,5],temp[,4],col=color2,pch=c(21:25),bg=color2)
  }}
dev.off()


pdf('./Comp.het.CSV.V2.Stat.pdf', useDingbats=FALSE)
data1=read.table('./comp.het/Integrated.CSV.stat')
par(mar=c(3,3,3,3))
par(mfrow=c(2,2))
SV_type_names=c('INV_DUP' ,'DEL_INV' ,'DEL_DUP', 'DEL_INV_DUP' )
SV_name=sort(unique(data1[,3]))
for(k0 in unique(data1[,3])){
  plot(c(0,1),c(0,1),mgp=c(1.5,0.5,0),xlab='False Discovery Rate',ylab='Sensitivity',cex.lab=0.8,type='n')
  title(main=SV_type_names[match(k0,SV_name)],cex.main=0.9)
  cols=c('red','blue','seagreen','purple','orange','black')
  pch=22
  algorithms=c('SVelter','Delly','Lumpy','Pindel','erds')
  for(k1 in unique(data1[,1])){
    color=cols[match(k1,algorithms)]
    color2=paste(color,c('','1','2','3','4'),sep='')
    temp=data1[data1[,1]==k1 & data1[,3]==k0,]
    temp=temp[order(temp[,2]),]
    points(1-temp[,5],temp[,4],col=color2,pch=c(21:25),bg=color2)
  }}

dev.off()

data1=read.table('./comp.homo/Integrated.CSV.stat')
data2=read.table('./comp.het/Integrated.CSV.stat')
data1[,6]=data1[,6]+data2[,6]
data1[,7]=data1[,7]+data2[,7]
data1[,8]=data1[,8]+data2[,8]
data1[,4]=data1[,6]/data1[,7]
data1[,5]=data1[,6]/data1[,8]
pdf('./Comp.Combine.CSV.V2.Stat.pdf', useDingbats=FALSE)
par(mar=c(3,3,3,3))
par(mfrow=c(2,2))
SV_type_names=c('INV_DUP' ,'DEL_INV' ,'DEL_DUP', 'DEL_INV_DUP' )
SV_name=sort(unique(data1[,3]))
for(k0 in unique(data1[,3])){
  plot(c(0,1),c(0,1),mgp=c(1.5,0.5,0),xlab='False Discovery Rate',ylab='Sensitivity',cex.lab=0.8,type='n')
  title(main=SV_type_names[match(k0,SV_name)],cex.main=0.9)
  cols=c('red','blue','seagreen','purple','orange','black')
  pch=22
  algorithms=c('SVelter','Delly','Lumpy','Pindel','erds')
  for(k1 in unique(data1[,1])){
    color=cols[match(k1,algorithms)]
    color2=paste(color,c('','1','2','3','4'),sep='')
    temp=data1[data1[,1]==k1 & data1[,3]==k0,]
    temp=temp[order(temp[,2]),]
    points(1-temp[,5],temp[,4],col=color2,pch=c(21:25),bg=color2)
  }}

dev.off()
















setwd('~/Box Sync/git-repo/svelter.xuefang.github/svelter/Support.records/Simulated/')

pdf('./Comp.homo.CSV.V0.Stat.pdf', useDingbats=FALSE)
data1=read.table('./comp.homo/Integrated.CSV.stat')
par(mar=c(3,3,3,3))
par(mfrow=c(2,2))
SV_type_names=c('INV_DUP' ,'DEL_INV' ,'DEL_DUP', 'DEL_INV_DUP' )
SV_name=sort(unique(data1[,3]))
for(k0 in unique(data1[,3])){
  plot(c(0,1),c(0,1),mgp=c(1.5,0.5,0),xlab='False Discovery Rate',ylab='Sensitivity',cex.lab=0.8,type='n')
  title(main=SV_type_names[match(k0,SV_name)],cex.main=0.9)
  cols=c('red','blue','seagreen','purple','orange','black')
  pch=22
  algorithms=c('SVelter','Delly','Lumpy','Pindel','erds')
  for(k1 in unique(data1[,1])){
    color=cols[match(k1,algorithms)]
    color2=paste(color,c('','1','2','3','4'),sep='')
    temp=data1[data1[,1]==k1 & data1[,3]==k0,]
    temp=temp[order(temp[,2]),]
    points(1-temp[,5],temp[,4],pch=21,bg=color)
  }}
dev.off()


pdf('./Comp.het.CSV.V0.Stat.pdf', useDingbats=FALSE)
data1=read.table('./comp.het/Integrated.CSV.stat')
par(mar=c(3,3,3,3))
par(mfrow=c(2,2))
SV_type_names=c('INV_DUP' ,'DEL_INV' ,'DEL_DUP', 'DEL_INV_DUP' )
SV_name=sort(unique(data1[,3]))
for(k0 in unique(data1[,3])){
  plot(c(0,1),c(0,1),mgp=c(1.5,0.5,0),xlab='False Discovery Rate',ylab='Sensitivity',cex.lab=0.8,type='n')
  title(main=SV_type_names[match(k0,SV_name)],cex.main=0.9)
  cols=c('red','blue','seagreen','purple','orange','black')
  pch=22
  algorithms=c('SVelter','Delly','Lumpy','Pindel','erds')
  for(k1 in unique(data1[,1])){
    color=cols[match(k1,algorithms)]
    color2=paste(color,c('','1','2','3','4'),sep='')
    temp=data1[data1[,1]==k1 & data1[,3]==k0,]
    temp=temp[order(temp[,2]),]
    points(1-temp[,5],temp[,4],pch=21,bg=color)
  }}

dev.off()

data1=read.table('./comp.homo/Integrated.CSV.stat')
data2=read.table('./comp.het/Integrated.CSV.stat')
data1[,6]=data1[,6]+data2[,6]
data1[,7]=data1[,7]+data2[,7]
data1[,8]=data1[,8]+data2[,8]
data1[,4]=data1[,6]/data1[,7]
data1[,5]=data1[,6]/data1[,8]
pdf('./Comp.Combine.CSV.V0.Stat.pdf', useDingbats=FALSE)
par(mar=c(3,3,3,3))
par(mfrow=c(2,2))
SV_type_names=c('INV_DUP' ,'DEL_INV' ,'DEL_DUP', 'DEL_INV_DUP' )
SV_name=sort(unique(data1[,3]))
for(k0 in unique(data1[,3])){
  plot(c(0,1),c(0,1),mgp=c(1.5,0.5,0),xlab='False Discovery Rate',ylab='Sensitivity',cex.lab=0.8,type='n')
  title(main=SV_type_names[match(k0,SV_name)],cex.main=0.9)
  cols=c('red','blue','seagreen','purple','orange','black')
  pch=22
  algorithms=c('SVelter','Delly','Lumpy','Pindel','erds')
  for(k1 in unique(data1[,1])){
    color=cols[match(k1,algorithms)]
    color2=paste(color,c('','1','2','3','4'),sep='')
    temp=data1[data1[,1]==k1 & data1[,3]==k0,]
    temp=temp[order(temp[,2]),]
    points(1-temp[,5],temp[,4],pch=21,bg=color)
  }}

dev.off()



















setwd('~/Box Sync/git-repo/svelter.xuefang.github/svelter/Support.records/Simulated/')

pdf('./Comp.homo.CSV.Stat.pdf', useDingbats=FALSE)
data1=read.table('./comp.homo/Integrated.CSV.stat')
par(mar=c(3,3,3,3))
par(mfrow=c(2,2))
SV_type_names=c('INV_DUP' ,'DEL_INV' ,'DEL_DUP', 'DEL_INV_DUP' )
SV_name=sort(unique(data1[,3]))
for(k0 in unique(data1[,3])){
  plot(c(0,1),c(0,1),mgp=c(1.5,0.5,0),xlab='False Discovery Rate',ylab='Sensitivity',cex.lab=0.8,type='n')
  title(main=SV_type_names[match(k0,SV_name)],cex.main=0.9)
  cols=c('red','blue','seagreen','purple','orange','black')
  pch=22
  algorithms=c('SVelter','Delly','Lumpy','Pindel','erds')
  for(k1 in unique(data1[,1])){
    color=cols[match(k1,algorithms)]
    color2=paste(color,c('','1','2','3','4'),sep='')
    temp=data1[data1[,1]==k1 & data1[,3]==k0,]
    temp=temp[order(temp[,2]),]
    points(1-temp[,5],temp[,4],pch=21,col=color,bg=color)
  }}
dev.off()


pdf('./Comp.het.CSV.Stat.pdf', useDingbats=FALSE)
data1=read.table('./comp.het/Integrated.CSV.stat')
par(mar=c(3,3,3,3))
par(mfrow=c(2,2))
SV_type_names=c('INV_DUP' ,'DEL_INV' ,'DEL_DUP', 'DEL_INV_DUP' )
SV_name=sort(unique(data1[,3]))
for(k0 in unique(data1[,3])){
  plot(c(0,1),c(0,1),mgp=c(1.5,0.5,0),xlab='False Discovery Rate',ylab='Sensitivity',cex.lab=0.8,type='n')
  title(main=SV_type_names[match(k0,SV_name)],cex.main=0.9)
  cols=c('red','blue','seagreen','purple','orange','black')
  pch=22
  algorithms=c('SVelter','Delly','Lumpy','Pindel','erds')
  for(k1 in unique(data1[,1])){
    color=cols[match(k1,algorithms)]
    color2=paste(color,c('','1','2','3','4'),sep='')
    temp=data1[data1[,1]==k1 & data1[,3]==k0,]
    temp=temp[order(temp[,2]),]
    points(1-temp[,5],temp[,4],pch=21,col=color,bg=color)
  }}

dev.off()

data1=read.table('./comp.homo/Integrated.CSV.stat')
data2=read.table('./comp.het/Integrated.CSV.stat')
data1[,6]=data1[,6]+data2[,6]
data1[,7]=data1[,7]+data2[,7]
data1[,8]=data1[,8]+data2[,8]
data1[,4]=data1[,6]/data1[,7]
data1[,5]=data1[,6]/data1[,8]
pdf('./Comp.Combine.CSV.Stat.pdf', useDingbats=FALSE)
par(mar=c(3,3,3,3))
par(mfrow=c(2,2))
SV_type_names=c('INV_DUP' ,'DEL_INV' ,'DEL_DUP', 'DEL_INV_DUP' )
SV_name=sort(unique(data1[,3]))
for(k0 in unique(data1[,3])){
  plot(c(0,1),c(0,1),mgp=c(1.5,0.5,0),xlab='False Discovery Rate',ylab='Sensitivity',cex.lab=0.8,type='n')
  title(main=SV_type_names[match(k0,SV_name)],cex.main=0.9)
  cols=c('red','blue','seagreen','purple','orange','black')
  pch=22
  algorithms=c('SVelter','Delly','Lumpy','Pindel','erds')
  for(k1 in unique(data1[,1])){
    color=cols[match(k1,algorithms)]
    color2=paste(color,c('','1','2','3','4'),sep='')
    temp=data1[data1[,1]==k1 & data1[,3]==k0,]
    temp=temp[order(temp[,2]),]
    points(1-temp[,5],temp[,4],pch=21,col=color,bg=color)
  }}

dev.off()






















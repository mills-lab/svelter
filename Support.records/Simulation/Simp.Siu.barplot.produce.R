setwd('~/Box Sync/git-repo/svelter.xuefang.github/svelter/Support.records/Simulated/')
data_reorganize<-function(data){
  for(i in 1:nrow(data)){
    data[i,8]=strsplit(as.character(data[i,1]),'_')[[1]][1]
    data[i,9]=strsplit(as.character(data[i,1]),'_')[[1]][3]
    }
  data=data[,c(1,2,9,3,4,5,6,7,8)]
  return (data)
}

Input_file='Simp.Homo.Pseudo.ROC.Mappable.TRAFree.min100.max1000000000.Barplot.Stats'
data_het=read.table(Input_file,header=T)
data_new=data_reorganize(data_het)
Het_Delly=data_new[data_new[,9]=='Delly',]
for(k1 in unique(Het_Delly[,3])){
  temp1=Het_Delly[Het_Delly[,3]==k1,]
  all_dups=temp1[temp1[,2]=='DUP',4]+temp1[temp1[,2]=='DUP_TANDEM',4]
  total_dup_calls=temp1[temp1[,2]=='DUP',6]
  Het_Delly[Het_Delly[,3]==k1 & Het_Delly[,2]=='DUP',8]=0
  Het_Delly[Het_Delly[,3]==k1 & Het_Delly[,2]=='DUP_TANDEM',8]=all_dups/total_dup_calls
}

Het_Lumpy=data_new[data_new[,9]=='Lumpy',]
for(k1 in unique(Het_Lumpy[,3])){
  temp1=Het_Lumpy[Het_Lumpy[,3]==k1,]
  all_dups=temp1[temp1[,2]=='DUP',4]+temp1[temp1[,2]=='DUP_TANDEM',4]
  total_dup_calls=temp1[temp1[,2]=='DUP',6]
  Het_Lumpy[Het_Lumpy[,3]==k1 & Het_Lumpy[,2]=='DUP',8]=0
  Het_Lumpy[Het_Lumpy[,3]==k1 & Het_Lumpy[,2]=='DUP_TANDEM',8]=all_dups/total_dup_calls
}

Het_SVelter=data_new[data_new[,9]=='SVelter',]

Het_Pindel=data_new[data_new[,9]=='Pindel',]
for(k1 in unique(Het_Pindel[,3])){
  temp1=Het_Pindel[Het_Pindel[,3]==k1,]
  all_dups=temp1[temp1[,2]=='DUP',4]+temp1[temp1[,2]=='DUP_TANDEM',4]
  total_dup_calls=temp1[temp1[,2]=='DUP',6]
  Het_Pindel[Het_Pindel[,3]==k1 & Het_Pindel[,2]=='DUP',8]=0
  Het_Pindel[Het_Pindel[,3]==k1 & Het_Pindel[,2]=='DUP_TANDEM',8]=all_dups/total_dup_calls
}

Het_erds=data_new[data_new[,9]=='erds',]
for(k1 in unique(Het_erds[,3])){
  temp1=Het_erds[Het_erds[,3]==k1,]
  all_dups=temp1[temp1[,2]=='DUP',4]+temp1[temp1[,2]=='DUP_TANDEM',4]
  total_dup_calls=temp1[temp1[,2]=='DUP',6]
  Het_erds[Het_erds[,3]==k1 & Het_erds[,2]=='DUP',8]=0
  Het_erds[Het_erds[,3]==k1 & Het_erds[,2]=='DUP_TANDEM',8]=all_dups/total_dup_calls
}

pdf(file=gsub('.Stats','.Barplot.pdf',Input_file), useDingbats = FALSE )
cols_Delly=rgb(red=0, green=0, blue=1.0, alpha=0.5)
cols_Lumpy=rgb(red=0.2, green=1.0, blue=0.2, alpha=0.2)
cols=c('red','blue','darkgreen','purple','orange','black')
cols_DUP=c('red',cols_Delly,cols_Lumpy,'purple','black')

par(mfrow=c(1,1))
par(fig=c(0,0.2,0.5,0.95))
plot(c(0,1),c(0,1),type='n',xlab='',ylab='',axes=FALSE)
text(0.5,0.5,'Sensitivity',cex=1,srt=90)
par(fig=c(0,0.2,0.1,0.55),new=TRUE)
plot(c(0,1),c(0,1),type='n',xlab='',ylab='',axes=FALSE)
text(0.5,0.5,'Specificity',cex=1,srt=90)
dis_per=0.3
dis_rec=0
for(SV in unique(Het_Delly[,2])){
  dis_rec=dis_rec+1
  Methods=c('Delly','Lumpy','SVelter','Pindel','erds')
  xrange=c(1,20)
  yrange=c(-1,1)
  #plot(xrange,yrange,type='n',xlab='1-Specificity',ylab='Specificity/Sensitivity',main=paste('Het',SV,sep='_'))
  #abline(h=0,v=0)
  rec=0
  data_all1=c()
  data_all2=c()
  for(x in sort(unique(Het_Delly[,3]))){
    rec=1
    t_delly=Het_Delly[Het_Delly[,2]==SV&Het_Delly[,3]==x,]
    t_lumpy=Het_Lumpy[Het_Lumpy[,2]==SV&Het_Lumpy[,3]==x,]
    t_svelter=Het_SVelter[Het_SVelter[,2]==SV&Het_SVelter[,3]==x,]
    t_pindel=Het_Pindel[Het_Pindel[,2]==SV&Het_Pindel[,3]==x,]
    t_erds=Het_erds[Het_erds[,2]==SV&Het_erds[,3]==x,]
    data1=c(t_svelter[1,7],t_delly[1,7],t_lumpy[1,7],t_pindel[1,7],t_erds[1,7],0)
    data2=c(t_svelter[1,8],t_delly[1,8],t_lumpy[1,8],t_pindel[1,8],t_erds[1,8],0)
    data_all1=c(data_all1,data1)
    data_all2=c(data_all2,data2) 
    data_all1[is.na(data_all1)]=0
    data_all2[is.na(data_all2)]=0
   }
  if(SV=='DUP'){SVType='Dispersed Duplications'}
  if(SV=='DUP_TANDEM'){SVType='Tandem Duplications'}
  if(SV=='DEL'){SVType='Deletions'}
  if(SV=='INV'){SVType='Inversions'}
  par(fig=c(-0.1+dis_rec*0.2,0.2+dis_rec*0.2,0.5,0.95), new=TRUE)
  par(mar = rep(2.4, 4))
  barplot(data_all1,col=cols,cex.lab=0.6,ylim=c(0,1),main='',cex.names=0.8,cex.axis=0.6,mgp=c(0.05,0.25,0.01))
  axis(1,at=c(1,7,13,19,25)+3,cex.axis=0.6,label=c('10X','20X','30X','40X','50X'),mgp=c(0.05,0.25,0.01))
  #axis(2,at=c(0.2,0.4,0.6,0.8,1.0),cex.axis=0.6,label=c(0.2,0.4,0.6,0.8,1.0))
}

dis_per=0.3
dis_rec=0
for(SV in unique(Het_Delly[,2])){
  dis_rec=dis_rec+1
  xrange=c(1,20)
  yrange=c(-1,1)
  #plot(xrange,yrange,type='n',xlab='1-Specificity',ylab='Specificity/Sensitivity',main=paste('Het',SV,sep='_'))
  #abline(h=0,v=0)
  rec=0
  data_all1=c()
  data_all2=c()
  for(x in sort(unique(Het_Delly[,3]))){
    rec=1
    t_delly=Het_Delly[Het_Delly[,2]==SV&Het_Delly[,3]==x,]
    t_lumpy=Het_Lumpy[Het_Lumpy[,2]==SV&Het_Lumpy[,3]==x,]
    t_svelter=Het_SVelter[Het_SVelter[,2]==SV&Het_SVelter[,3]==x,]
    t_pindel=Het_Pindel[Het_Pindel[,2]==SV&Het_Pindel[,3]==x,]
    t_erds=Het_erds[Het_erds[,2]==SV&Het_erds[,3]==x,]
    data1=c(t_svelter[1,7],t_delly[1,7],t_lumpy[1,7],t_pindel[1,7],t_erds[1,7],0)
    data2=c(t_svelter[1,8],t_delly[1,8],t_lumpy[1,8],t_pindel[1,8],t_erds[1,8],0)
    data_all1=c(data_all1,data1)
    data_all2=c(data_all2,data2) 
    data_all1[is.na(data_all1)]=0
    data_all2[is.na(data_all2)]=0
  }
    if(SV=='DUP'){SVType='Dispersed Duplications'}
  if(SV=='DUP_TANDEM'){SVType='Tandem Duplications'}
  if(SV=='DEL'){SVType='Deletions'}
  if(SV=='INV'){SVType='Inversions'}
  par(fig=c(-0.1+dis_rec*0.2,0.2+dis_rec*0.2,0.1,0.55), new=TRUE)
  par(mar = rep(2.4, 4))
  barplot(data_all2,col=cols,cex.lab=0.6,ylim=c(0,1),main='',cex.names=0.8,cex.axis=0.6,mgp=c(0.05,0.25,0.01))
  axis(1,at=c(1,7,13,19,25)+3,cex.axis=0.6,label=c('10X','20X','30X','40X','50X'),mgp=c(0.05,0.25,0.01))
}

par(mar = rep(0.4, 4))
par(fig=c(0.1,0.4,0.85,1),new=TRUE)
plot(c(0,1),c(0,1),type='n',xlab='',ylab='',axes=FALSE)
text(0.5,0.5,'Deletion',cex=0.8)

par(fig=c(0.3,0.6,0.85,1),new=TRUE)
plot(c(0,1),c(0,1),type='n',xlab='',ylab='',axes=FALSE)
text(0.5,0.5,'Disperse Duplication',cex=0.8)

par(fig=c(0.5,0.8,0.85,1),new=TRUE)
plot(c(0,1),c(0,1),type='n',xlab='',ylab='',axes=FALSE)
text(0.5,0.5,'Tandem Duplication',cex=0.8)

par(fig=c(0.7,1,0.85,1),new=TRUE)
plot(c(0,1),c(0,1),type='n',xlab='',ylab='',axes=FALSE)
text(0.5,0.5,'Inversion',cex=0.8)

dev.off()




data_reorganize<-function(data){
  for(i in 1:nrow(data)){
    data[i,8]=strsplit(as.character(data[i,1]),'_')[[1]][1]
    data[i,9]=strsplit(as.character(data[i,1]),'_')[[1]][3]
  }
  data=data[,c(1,2,9,3,4,5,6,7,8)]
  return (data)
}
Input_file='Simp.Het.Pseudo.ROC.Mappable.TRAFree.min100.max1000000000.Barplot.Stats'
data_het=read.table(Input_file,header=T)
data_new=data_reorganize(data_het)
Het_Delly=data_new[data_new[,9]=='Delly',]
for(k1 in unique(Het_Delly[,3])){
  temp1=Het_Delly[Het_Delly[,3]==k1,]
  all_dups=temp1[temp1[,2]=='DUP',4]+temp1[temp1[,2]=='DUP_TANDEM',4]
  total_dup_calls=temp1[temp1[,2]=='DUP',6]
  Het_Delly[Het_Delly[,3]==k1 & Het_Delly[,2]=='DUP',8]=0
  Het_Delly[Het_Delly[,3]==k1 & Het_Delly[,2]=='DUP_TANDEM',8]=all_dups/total_dup_calls
}

Het_Lumpy=data_new[data_new[,9]=='Lumpy',]
for(k1 in unique(Het_Lumpy[,3])){
  temp1=Het_Lumpy[Het_Lumpy[,3]==k1,]
  all_dups=temp1[temp1[,2]=='DUP',4]+temp1[temp1[,2]=='DUP_TANDEM',4]
  total_dup_calls=temp1[temp1[,2]=='DUP',6]
  Het_Lumpy[Het_Lumpy[,3]==k1 & Het_Lumpy[,2]=='DUP',8]=0
  Het_Lumpy[Het_Lumpy[,3]==k1 & Het_Lumpy[,2]=='DUP_TANDEM',8]=all_dups/total_dup_calls
}

Het_SVelter=data_new[data_new[,9]=='SVelter',]

Het_Pindel=data_new[data_new[,9]=='Pindel',]
for(k1 in unique(Het_Pindel[,3])){
  temp1=Het_Pindel[Het_Pindel[,3]==k1,]
  all_dups=temp1[temp1[,2]=='DUP',4]+temp1[temp1[,2]=='DUP_TANDEM',4]
  total_dup_calls=temp1[temp1[,2]=='DUP',6]
  Het_Pindel[Het_Pindel[,3]==k1 & Het_Pindel[,2]=='DUP',8]=0
  Het_Pindel[Het_Pindel[,3]==k1 & Het_Pindel[,2]=='DUP_TANDEM',8]=all_dups/total_dup_calls
}

Het_erds=data_new[data_new[,9]=='erds',]
for(k1 in unique(Het_erds[,3])){
  temp1=Het_erds[Het_erds[,3]==k1,]
  all_dups=temp1[temp1[,2]=='DUP',4]+temp1[temp1[,2]=='DUP_TANDEM',4]
  total_dup_calls=temp1[temp1[,2]=='DUP',6]
  Het_erds[Het_erds[,3]==k1 & Het_erds[,2]=='DUP',8]=0
  Het_erds[Het_erds[,3]==k1 & Het_erds[,2]=='DUP_TANDEM',8]=all_dups/total_dup_calls
}


pdf(file=gsub('.Stats','.Barplot.pdf',Input_file) ,useDingbats = FALSE)
cols_Delly=rgb(red=0, green=0, blue=1.0, alpha=0.5)
cols_Lumpy=rgb(red=0.2, green=1.0, blue=0.2, alpha=0.2)
cols=c('red','blue','darkgreen','purple','orange','black')
cols_DUP=c('red',cols_Delly,cols_Lumpy,'purple','black')

par(mfrow=c(1,1))
par(fig=c(0,0.2,0.5,0.95))
plot(c(0,1),c(0,1),type='n',xlab='',ylab='',axes=FALSE)
text(0.5,0.5,'Sensitivity',cex=1,srt=90)
par(fig=c(0,0.2,0.1,0.55),new=TRUE)
plot(c(0,1),c(0,1),type='n',xlab='',ylab='',axes=FALSE)
text(0.5,0.5,'Specificity',cex=1,srt=90)
dis_per=0.3
dis_rec=0
for(SV in unique(Het_Delly[,2])){
  dis_rec=dis_rec+1
  Methods=c('Delly','Lumpy','SVelter','Pindel','erds')
  xrange=c(1,20)
  yrange=c(-1,1)
  #plot(xrange,yrange,type='n',xlab='1-Specificity',ylab='Specificity/Sensitivity',main=paste('Het',SV,sep='_'))
  #abline(h=0,v=0)
  rec=0
  data_all1=c()
  data_all2=c()
  for(x in sort(unique(Het_Delly[,3]))){
    rec=1
    t_delly=Het_Delly[Het_Delly[,2]==SV&Het_Delly[,3]==x,]
    t_lumpy=Het_Lumpy[Het_Lumpy[,2]==SV&Het_Lumpy[,3]==x,]
    t_svelter=Het_SVelter[Het_SVelter[,2]==SV&Het_SVelter[,3]==x,]
    t_pindel=Het_Pindel[Het_Pindel[,2]==SV&Het_Pindel[,3]==x,]
    t_erds=Het_erds[Het_erds[,2]==SV&Het_erds[,3]==x,]
    data1=c(t_svelter[1,7],t_delly[1,7],t_lumpy[1,7],t_pindel[1,7],t_erds[1,7],0)
    data2=c(t_svelter[1,8],t_delly[1,8],t_lumpy[1,8],t_pindel[1,8],t_erds[1,8],0)
    data_all1=c(data_all1,data1)
    data_all2=c(data_all2,data2) 
    data_all1[is.na(data_all1)]=0
    data_all2[is.na(data_all2)]=0
  }
  if(SV=='DUP'){SVType='Dispersed Duplications'}
  if(SV=='DUP_TANDEM'){SVType='Tandem Duplications'}
  if(SV=='DEL'){SVType='Deletions'}
  if(SV=='INV'){SVType='Inversions'}
  par(fig=c(-0.1+dis_rec*0.2,0.2+dis_rec*0.2,0.5,0.95), new=TRUE)
  par(mar = rep(2.4, 4))
  barplot(data_all1,col=cols,cex.lab=0.6,ylim=c(0,1),main='',cex.names=0.8,cex.axis=0.6,mgp=c(0.05,0.25,0.01))
  axis(1,at=c(1,7,13,19,25)+3,cex.axis=0.6,label=c('10X','20X','30X','40X','50X'),mgp=c(0.05,0.25,0.01))
  #axis(2,at=c(0.2,0.4,0.6,0.8,1.0),cex.axis=0.6,label=c(0.2,0.4,0.6,0.8,1.0))
}

dis_per=0.3
dis_rec=0
for(SV in unique(Het_Delly[,2])){
  dis_rec=dis_rec+1
  xrange=c(1,20)
  yrange=c(-1,1)
  #plot(xrange,yrange,type='n',xlab='1-Specificity',ylab='Specificity/Sensitivity',main=paste('Het',SV,sep='_'))
  #abline(h=0,v=0)
  rec=0
  data_all1=c()
  data_all2=c()
  for(x in sort(unique(Het_Delly[,3]))){
    rec=1
    t_delly=Het_Delly[Het_Delly[,2]==SV&Het_Delly[,3]==x,]
    t_lumpy=Het_Lumpy[Het_Lumpy[,2]==SV&Het_Lumpy[,3]==x,]
    t_svelter=Het_SVelter[Het_SVelter[,2]==SV&Het_SVelter[,3]==x,]
    t_pindel=Het_Pindel[Het_Pindel[,2]==SV&Het_Pindel[,3]==x,]
    t_erds=Het_erds[Het_erds[,2]==SV&Het_erds[,3]==x,]
    data1=c(t_svelter[1,7],t_delly[1,7],t_lumpy[1,7],t_pindel[1,7],t_erds[1,7],0)
    data2=c(t_svelter[1,8],t_delly[1,8],t_lumpy[1,8],t_pindel[1,8],t_erds[1,8],0)
    data_all1=c(data_all1,data1)
    data_all2=c(data_all2,data2) 
    data_all1[is.na(data_all1)]=0
    data_all2[is.na(data_all2)]=0
  }
  if(SV=='DUP'){SVType='Dispersed Duplications'}
  if(SV=='DUP_TANDEM'){SVType='Tandem Duplications'}
  if(SV=='DEL'){SVType='Deletions'}
  if(SV=='INV'){SVType='Inversions'}
  par(fig=c(-0.1+dis_rec*0.2,0.2+dis_rec*0.2,0.1,0.55), new=TRUE)
  par(mar = rep(2.4, 4))
  barplot(data_all2,col=cols,cex.lab=0.6,ylim=c(0,1),main='',cex.names=0.8,cex.axis=0.6,mgp=c(0.05,0.25,0.01))
  axis(1,at=c(1,7,13,19,25)+3,cex.axis=0.6,label=c('10X','20X','30X','40X','50X'),mgp=c(0.05,0.25,0.01))
}

par(mar = rep(0.4, 4))
par(fig=c(0.1,0.4,0.85,1),new=TRUE)
plot(c(0,1),c(0,1),type='n',xlab='',ylab='',axes=FALSE)
text(0.5,0.5,'Deletion',cex=0.8)

par(fig=c(0.3,0.6,0.85,1),new=TRUE)
plot(c(0,1),c(0,1),type='n',xlab='',ylab='',axes=FALSE)
text(0.5,0.5,'Disperse Duplication',cex=0.8)

par(fig=c(0.5,0.8,0.85,1),new=TRUE)
plot(c(0,1),c(0,1),type='n',xlab='',ylab='',axes=FALSE)
text(0.5,0.5,'Tandem Duplication',cex=0.8)

par(fig=c(0.7,1,0.85,1),new=TRUE)
plot(c(0,1),c(0,1),type='n',xlab='',ylab='',axes=FALSE)
text(0.5,0.5,'Inversion',cex=0.8)

dev.off()






Input_file_1='Simp.Het.Pseudo.ROC.Mappable.TRAFree.min100.max1000000000.Barplot.Stats'
Input_file_2='Simp.Homo.Pseudo.ROC.Mappable.TRAFree.min100.max1000000000.Barplot.Stats'
Input_file='Simp.Combine.Pseudo.ROC.Mappable.TRAFree.min100.max1000000000.Barplot.Stats'
data_het1=read.table(Input_file_1,header=T)
data_homo1=read.table(Input_file_2,header=T)
data_het=data_het1
data_het[,3]=data_het1[,3]+data_homo1[,3]
data_het[,4]=data_het1[,4]+data_homo1[,4]
data_het[,5]=data_het1[,5]+data_homo1[,5]
data_het[,6]=data_het[,3]/data_het[,4]
data_het[,7]=data_het[,3]/data_het[,5]
data_new=data_reorganize(data_het)
Het_Delly=data_new[data_new[,9]=='Delly',]
for(k1 in unique(Het_Delly[,3])){
  temp1=Het_Delly[Het_Delly[,3]==k1,]
  all_dups=temp1[temp1[,2]=='DUP',4]+temp1[temp1[,2]=='DUP_TANDEM',4]
  total_dup_calls=temp1[temp1[,2]=='DUP',6]
  Het_Delly[Het_Delly[,3]==k1 & Het_Delly[,2]=='DUP',8]=0
  Het_Delly[Het_Delly[,3]==k1 & Het_Delly[,2]=='DUP_TANDEM',8]=all_dups/total_dup_calls
}

Het_Lumpy=data_new[data_new[,9]=='Lumpy',]
for(k1 in unique(Het_Lumpy[,3])){
  temp1=Het_Lumpy[Het_Lumpy[,3]==k1,]
  all_dups=temp1[temp1[,2]=='DUP',4]+temp1[temp1[,2]=='DUP_TANDEM',4]
  total_dup_calls=temp1[temp1[,2]=='DUP',6]
  Het_Lumpy[Het_Lumpy[,3]==k1 & Het_Lumpy[,2]=='DUP',8]=0
  Het_Lumpy[Het_Lumpy[,3]==k1 & Het_Lumpy[,2]=='DUP_TANDEM',8]=all_dups/total_dup_calls
}

Het_SVelter=data_new[data_new[,9]=='SVelter',]

Het_Pindel=data_new[data_new[,9]=='Pindel',]
for(k1 in unique(Het_Pindel[,3])){
  temp1=Het_Pindel[Het_Pindel[,3]==k1,]
  all_dups=temp1[temp1[,2]=='DUP',4]+temp1[temp1[,2]=='DUP_TANDEM',4]
  total_dup_calls=temp1[temp1[,2]=='DUP',6]
  Het_Pindel[Het_Pindel[,3]==k1 & Het_Pindel[,2]=='DUP',8]=0
  Het_Pindel[Het_Pindel[,3]==k1 & Het_Pindel[,2]=='DUP_TANDEM',8]=all_dups/total_dup_calls
}

Het_erds=data_new[data_new[,9]=='erds',]
for(k1 in unique(Het_erds[,3])){
  temp1=Het_erds[Het_erds[,3]==k1,]
  all_dups=temp1[temp1[,2]=='DUP',4]+temp1[temp1[,2]=='DUP_TANDEM',4]
  total_dup_calls=temp1[temp1[,2]=='DUP',6]
  Het_erds[Het_erds[,3]==k1 & Het_erds[,2]=='DUP',8]=0
  Het_erds[Het_erds[,3]==k1 & Het_erds[,2]=='DUP_TANDEM',8]=all_dups/total_dup_calls
}


pdf(file=gsub('.Stats','.Barplot.pdf',Input_file), useDingbats = FALSE)
cols_Delly=rgb(red=0, green=0, blue=1.0, alpha=0.5)
cols_Lumpy=rgb(red=0.2, green=1.0, blue=0.2, alpha=0.2)
cols=c('red','blue','darkgreen','purple','orange','black')
cols_DUP=c('red',cols_Delly,cols_Lumpy,'purple','black')

par(mfrow=c(1,1))
par(fig=c(0,0.2,0.5,0.95))
plot(c(0,1),c(0,1),type='n',xlab='',ylab='',axes=FALSE)
text(0.5,0.5,'Sensitivity',cex=1,srt=90)
par(fig=c(0,0.2,0.1,0.55),new=TRUE)
plot(c(0,1),c(0,1),type='n',xlab='',ylab='',axes=FALSE)
text(0.5,0.5,'Specificity',cex=1,srt=90)
dis_per=0.3
dis_rec=0
for(SV in unique(Het_Delly[,2])){
  dis_rec=dis_rec+1
  Methods=c('Delly','Lumpy','SVelter','Pindel','erds')
  xrange=c(1,20)
  yrange=c(-1,1)
  #plot(xrange,yrange,type='n',xlab='1-Specificity',ylab='Specificity/Sensitivity',main=paste('Het',SV,sep='_'))
  #abline(h=0,v=0)
  rec=0
  data_all1=c()
  data_all2=c()
  for(x in sort(unique(Het_Delly[,3]))){
    rec=1
    t_delly=Het_Delly[Het_Delly[,2]==SV&Het_Delly[,3]==x,]
    t_lumpy=Het_Lumpy[Het_Lumpy[,2]==SV&Het_Lumpy[,3]==x,]
    t_svelter=Het_SVelter[Het_SVelter[,2]==SV&Het_SVelter[,3]==x,]
    t_pindel=Het_Pindel[Het_Pindel[,2]==SV&Het_Pindel[,3]==x,]
    t_erds=Het_erds[Het_erds[,2]==SV&Het_erds[,3]==x,]
    data1=c(t_svelter[1,7],t_delly[1,7],t_lumpy[1,7],t_pindel[1,7],t_erds[1,7],0)
    data2=c(t_svelter[1,8],t_delly[1,8],t_lumpy[1,8],t_pindel[1,8],t_erds[1,8],0)
    data_all1=c(data_all1,data1)
    data_all2=c(data_all2,data2) 
    data_all1[is.na(data_all1)]=0
    data_all2[is.na(data_all2)]=0
  }
  if(SV=='DUP'){SVType='Dispersed Duplications'}
  if(SV=='DUP_TANDEM'){SVType='Tandem Duplications'}
  if(SV=='DEL'){SVType='Deletions'}
  if(SV=='INV'){SVType='Inversions'}
  par(fig=c(-0.1+dis_rec*0.2,0.2+dis_rec*0.2,0.5,0.95), new=TRUE)
  par(mar = rep(2.4, 4))
  barplot(data_all1,col=cols,cex.lab=0.6,ylim=c(0,1),main='',cex.names=0.8,cex.axis=0.6,mgp=c(0.05,0.25,0.01))
  axis(1,at=c(1,7,13,19,25)+3,cex.axis=0.6,label=c('10X','20X','30X','40X','50X'),mgp=c(0.05,0.25,0.01))
  #axis(2,at=c(0.2,0.4,0.6,0.8,1.0),cex.axis=0.6,label=c(0.2,0.4,0.6,0.8,1.0))
}

dis_per=0.3
dis_rec=0
for(SV in unique(Het_Delly[,2])){
  dis_rec=dis_rec+1
  xrange=c(1,20)
  yrange=c(-1,1)
  #plot(xrange,yrange,type='n',xlab='1-Specificity',ylab='Specificity/Sensitivity',main=paste('Het',SV,sep='_'))
  #abline(h=0,v=0)
  rec=0
  data_all1=c()
  data_all2=c()
  for(x in sort(unique(Het_Delly[,3]))){
    rec=1
    t_delly=Het_Delly[Het_Delly[,2]==SV&Het_Delly[,3]==x,]
    t_lumpy=Het_Lumpy[Het_Lumpy[,2]==SV&Het_Lumpy[,3]==x,]
    t_svelter=Het_SVelter[Het_SVelter[,2]==SV&Het_SVelter[,3]==x,]
    t_pindel=Het_Pindel[Het_Pindel[,2]==SV&Het_Pindel[,3]==x,]
    t_erds=Het_erds[Het_erds[,2]==SV&Het_erds[,3]==x,]
    data1=c(t_svelter[1,7],t_delly[1,7],t_lumpy[1,7],t_pindel[1,7],t_erds[1,7],0)
    data2=c(t_svelter[1,8],t_delly[1,8],t_lumpy[1,8],t_pindel[1,8],t_erds[1,8],0)
    data_all1=c(data_all1,data1)
    data_all2=c(data_all2,data2) 
    data_all1[is.na(data_all1)]=0
    data_all2[is.na(data_all2)]=0
  }
  if(SV=='DUP'){SVType='Dispersed Duplications'}
  if(SV=='DUP_TANDEM'){SVType='Tandem Duplications'}
  if(SV=='DEL'){SVType='Deletions'}
  if(SV=='INV'){SVType='Inversions'}
  par(fig=c(-0.1+dis_rec*0.2,0.2+dis_rec*0.2,0.1,0.55), new=TRUE)
  par(mar = rep(2.4, 4))
  barplot(data_all2,col=cols,cex.lab=0.6,ylim=c(0,1),main='',cex.names=0.8,cex.axis=0.6,mgp=c(0.05,0.25,0.01))
  axis(1,at=c(1,7,13,19,25)+3,cex.axis=0.6,label=c('10X','20X','30X','40X','50X'),mgp=c(0.05,0.25,0.01))
}

par(mar = rep(0.4, 4))
par(fig=c(0.1,0.4,0.85,1),new=TRUE)
plot(c(0,1),c(0,1),type='n',xlab='',ylab='',axes=FALSE)
text(0.5,0.5,'Deletion',cex=0.8)

par(fig=c(0.3,0.6,0.85,1),new=TRUE)
plot(c(0,1),c(0,1),type='n',xlab='',ylab='',axes=FALSE)
text(0.5,0.5,'Disperse Duplication',cex=0.8)

par(fig=c(0.5,0.8,0.85,1),new=TRUE)
plot(c(0,1),c(0,1),type='n',xlab='',ylab='',axes=FALSE)
text(0.5,0.5,'Tandem Duplication',cex=0.8)

par(fig=c(0.7,1,0.85,1),new=TRUE)
plot(c(0,1),c(0,1),type='n',xlab='',ylab='',axes=FALSE)
text(0.5,0.5,'Inversion',cex=0.8)

dev.off()






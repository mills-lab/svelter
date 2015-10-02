#!/usr/bin/env python

#!R
#Usage: Rscript Produce.Barplot.For.Simple.Simu.R input-file1 input-file2
#eg: Input_file='/Users/xuefzhao/Box Sync/git-repo/SVelter.software/svelter/Support.Scripts/Simp.Het.Pseudo.ROC.Mappable.TRAFree.min100.max1000000000.Stats'
Args <- commandArgs()
Input_file1=Args[6]
Input_file2=Args[7]
reorganize <- function(data){
  out=data.frame('sample'=0,'sv'=0,'rd'=0,'overlap'=0, 'X.Ref.SVs'=0, 'X.Samp.SVs'=0)
  rec=0
  data=data[data[,3]!='Y',]
  data=data[data[,3]!='X',]
  for(x1 in unique(data[,1])){
    rd=as.integer(strsplit(strsplit(as.character(x1),'RD')[[1]][2],'_')[[1]][1])
    data1=data[data[,1]==x1,]
    for(x2 in unique(data1[,2])){
      rec=rec+1
      data2=data1[data1[,2]==x2,]
      out[rec,1]=x1
      out[rec,2]=x2
      out[rec,3]=rd
      out[rec,4]=sum(data2[,4])
      out[rec,5]=sum(data2[,5])
      out[rec,6]=sum(data2[,6])
    }}
  out=out[out[,3]!=5,]
  return(out)}
remove_TANDOMDUP<-function(Het_Delly){
  for(x in unique(Het_Delly[,3])){
    Het_Delly[Het_Delly[,3]==x & Het_Delly[,2]=='DUP',4]=Het_Delly[Het_Delly[,3]==x & Het_Delly[,2]=='DUP',4]+Het_Delly[Het_Delly[,3]==x & Het_Delly[,2]=='DUP_TANDEM',4]
    Het_Delly[Het_Delly[,3]==x & Het_Delly[,2]=='DUP',5]=Het_Delly[Het_Delly[,3]==x & Het_Delly[,2]=='DUP',5]+Het_Delly[Het_Delly[,3]==x & Het_Delly[,2]=='DUP_TANDEM',5]
    Het_Delly[Het_Delly[,3]==x & Het_Delly[,2]=='DUP',6]=Het_Delly[Het_Delly[,3]==x & Het_Delly[,2]=='DUP',6]+Het_Delly[Het_Delly[,3]==x & Het_Delly[,2]=='DUP_TANDEM',6]
  }
  Het_Delly=Het_Delly[Het_Delly[,2]!='DUP_TANDEM',]
  return(Het_Delly)
}
data_in=read.csv(Input_file1,sep=' ')
data_in2=read.csv(Input_file2,sep=' ')
data_in=data_in[order(data_in[,1]),]
data_in2=data_in2[order(data_in2[,1]),]
for(k1 in unique(data_in[,1])){
  for(k2 in unique(data_in[,2])){
    for(k3 in unique(data_in[,3])){
      k1b=gsub('het','homo',k1)
      data_in[data_in[,1]==k1&data_in[,2]==k2&data_in[,3]==k3,][1,4]=
        data_in[data_in[,1]==k1&data_in[,2]==k2&data_in[,3]==k3,][1,4]+
        data_in2[data_in2[,1]==k1b&data_in[,2]==k2&data_in2[,3]==k3,][1,4]
      data_in[data_in[,1]==k1&data_in[,2]==k2&data_in[,3]==k3,][1,5]=
        data_in[data_in[,1]==k1&data_in[,2]==k2&data_in[,3]==k3,][1,5]+
        data_in2[data_in2[,1]==k1b&data_in[,2]==k2&data_in2[,3]==k3,][1,5]
      data_in[data_in[,1]==k1&data_in[,2]==k2&data_in[,3]==k3,][1,6]=
        data_in[data_in[,1]==k1&data_in[,2]==k2&data_in[,3]==k3,][1,6]+
        data_in2[data_in2[,1]==k1b&data_in[,2]==k2&data_in2[,3]==k3,][1,6]
    }
  }}
Delly_het=data_in[1:480,]
Lumpy_het=data_in[481:960,]
Pindel_het=data_in[961:1440,]
SVelter_het=data_in[1441:1920,]

Het_Delly=reorganize(Delly_het)
#Het_Delly=remove_TANDOMDUP(Het_Delly)
Het_Delly[,7]=Het_Delly[,4]/Het_Delly[,5]
Het_Delly[,8]=Het_Delly[,4]/Het_Delly[,6]
Het_Delly=Het_Delly[order(Het_Delly[,1]),]
#Special process for Delly: merge all dels for specificity:
for(k1 in unique(Het_Delly[,3])){
  temp1=Het_Delly[Het_Delly[,3]==k1,]
  all_dups=temp1[temp1[,2]=='DUP',4]+temp1[temp1[,2]=='DUP_TANDEM',4]
  total_dup_calls=temp1[temp1[,2]=='DUP',6]
  Het_Delly[Het_Delly[,3]==k1 & Het_Delly[,2]=='DUP',8]=0
  Het_Delly[Het_Delly[,3]==k1 & Het_Delly[,2]=='DUP_TANDEM',8]=all_dups/total_dup_calls
}

Het_Lumpy=reorganize(Lumpy_het)
#Het_Lumpy=remove_TANDOMDUP(Het_Lumpy)
Het_Lumpy[,7]=Het_Lumpy[,4]/Het_Lumpy[,5]
Het_Lumpy[,8]=Het_Lumpy[,4]/Het_Lumpy[,6]
Het_Lumpy=Het_Lumpy[order(Het_Lumpy[,1]),]
#Special process for Lumpy: merge all dels for specificity:
for(k1 in unique(Het_Lumpy[,3])){
  temp1=Het_Lumpy[Het_Lumpy[,3]==k1,]
  all_dups=temp1[temp1[,2]=='DUP',4]+temp1[temp1[,2]=='DUP_TANDEM',4]
  total_dup_calls=temp1[temp1[,2]=='DUP',6]
  Het_Lumpy[Het_Lumpy[,3]==k1 & Het_Lumpy[,2]=='DUP',8]=0
  Het_Lumpy[Het_Lumpy[,3]==k1 & Het_Lumpy[,2]=='DUP_TANDEM',8]=all_dups/total_dup_calls
}


Het_SVelter=reorganize(SVelter_het)
#Het_SVelter=remove_TANDOMDUP(Het_SVelter)
Het_SVelter[,7]=Het_SVelter[,4]/Het_SVelter[,5]
Het_SVelter[,8]=Het_SVelter[,4]/Het_SVelter[,6]
Het_SVelter=Het_SVelter[order(Het_SVelter[,1]),]

Het_Pindel=reorganize(Pindel_het)
#Het_Pindel=remove_TANDOMDUP(Het_Pindel)
Het_Pindel[,7]=Het_Pindel[,4]/Het_Pindel[,5]
Het_Pindel[,8]=Het_Pindel[,4]/Het_Pindel[,6]
Het_Pindel=Het_Pindel[order(Het_Pindel[,1]),]
#Special process for Pindel: merge all dels for specificity:
for(k1 in unique(Het_Pindel[,3])){
  temp1=Het_Pindel[Het_Pindel[,3]==k1,]
  all_dups=temp1[temp1[,2]=='DUP',4]+temp1[temp1[,2]=='DUP_TANDEM',4]
  total_dup_calls=temp1[temp1[,2]=='DUP',6]
  Het_Pindel[Het_Pindel[,3]==k1 & Het_Pindel[,2]=='DUP',8]=0
  Het_Pindel[Het_Pindel[,3]==k1 & Het_Pindel[,2]=='DUP_TANDEM',8]=all_dups/total_dup_calls
}

Input_list1=strsplit(strsplit(Input_file1,'/')[[1]][length(strsplit(Input_file1,'/')[[1]])],'[.]')[[1]]
Input_list2=strsplit(strsplit(Input_file2,'/')[[1]][length(strsplit(Input_file2,'/')[[1]])],'[.]')[[1]]
Output_list=c()
for(x in 1:length(Input_list1)){
  if (Input_list1[x]==Input_list2[x]){
    Output_list=c(Output_list,Input_list1[x])
  }
  else{
    Output_list=c(Output_list,Input_list1[x],Input_list2[x])
  }
}
Output_file=paste(paste(strsplit(Input_file1,'/')[[1]][1:(length(strsplit(Input_file1,'/')[[1]])-1)],collapse='/'),paste(Output_list,collapse='.'),sep='/')

pdf(file=paste(gsub('.Stats','.Barplot.pdf',Output_file)))
cols_Delly=rgb(red=0, green=0, blue=1.0, alpha=0.5)
cols_Lumpy=rgb(red=0.2, green=1.0, blue=0.2, alpha=0.2)
cols=c('red','blue','darkgreen','purple','black')
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
  Methods=c('Delly','Lumpy','SVelter','Pindel')
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
    data1=c(t_svelter[1,7],t_delly[1,7],t_lumpy[1,7],t_pindel[1,7],0)
    data2=c(t_svelter[1,8],t_delly[1,8],t_lumpy[1,8],t_pindel[1,8],0)
    data_all1=c(data_all1,data1)
    data_all2=c(data_all2,data2) }
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
  cols=c('red','blue','darkgreen','purple','black')
  Methods=c('Delly','Lumpy','SVelter','Pindel')
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
    data1=c(t_svelter[1,7],t_delly[1,7],t_lumpy[1,7],t_pindel[1,7],0)
    data2=c(t_svelter[1,8],t_delly[1,8],t_lumpy[1,8],t_pindel[1,8],0)
    data_all1=c(data_all1,data1)
    data_all2=c(data_all2,data2) }
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




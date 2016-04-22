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
RDs=unique(data_new[,3])
SVs=unique(data_new[,2])
Methods=c('SVelter','Delly','Lumpy','Pindel','erds')
rec=0
out_simp=matrix(nrow=length(SVs)*length(RDs),ncol=2*length(Methods)+1)
row.names(out_simp)=rep(SVs,length(RDs))
colnames(out_simp)=c('RD',Methods,Methods)
for(k1 in RDs){
  for (k2 in SVs){
    for(k3 in Methods){
      if(nrow(data_new[data_new[,2]==k2 & data_new[,3]==k1 & data_new[,9]==k3,])==0){
        row_num=nrow(data_new)+1
        data_new[row_num,2]=k2
        data_new[row_num,3]=k1
        data_new[row_num,7]=0
        data_new[row_num,8]=0
        data_new[row_num,9]=k3
      }
    }
  }
}
for(k1 in RDs){
  for(k2 in SVs){
    rec=rec+1
    temp=data_new[data_new[,2]==k2&data_new[,3]==k1,]
    temp2=temp[,c(9,7,8)]
    temp2=temp2[order(temp2[,1]),]
    temp2=temp2[c(5,1,3,4,2),]
    out_simp[rec,1]=k1
    out_simp[rec,2:(1+nrow(temp2))]=temp2[,2]
    out_simp[rec,(2+nrow(temp2)):(1+2*nrow(temp2))]=temp2[,3]   
  }
}
out_simp_2=cbind(out_simp[,1],rownames(out_simp),out_simp[,2:11])
colnames(out_simp_2)[1:2]=c('RDs','SVs')
write.table(out_simp_2,'Integrated.compare.Simp.Homo.stat',quote=F,row.names=F)


Input_file='Simp.Het.Pseudo.ROC.Mappable.TRAFree.min100.max1000000000.Barplot.Stats'
data_het=read.table(Input_file,header=T)
data_new=data_reorganize(data_het)
RDs=unique(data_new[,3])
SVs=unique(data_new[,2])
Methods=c('SVelter','Delly','Lumpy','Pindel','erds')
rec=0
out_simp=matrix(nrow=length(SVs)*length(RDs),ncol=2*length(Methods)+1)
row.names(out_simp)=rep(SVs,length(RDs))
colnames(out_simp)=c('RD',Methods,Methods)
for(k1 in RDs){
  for (k2 in SVs){
    for(k3 in Methods){
      if(nrow(data_new[data_new[,2]==k2 & data_new[,3]==k1 & data_new[,9]==k3,])==0){
        row_num=nrow(data_new)+1
        data_new[row_num,2]=k2
        data_new[row_num,3]=k1
        data_new[row_num,7]=0
        data_new[row_num,8]=0
        data_new[row_num,9]=k3
      }
    }
  }
}
for(k1 in RDs){
  for(k2 in SVs){
    rec=rec+1
    temp=data_new[data_new[,2]==k2&data_new[,3]==k1,]
    temp2=temp[,c(9,7,8)]
    temp2=temp2[order(temp2[,1]),]
    temp2=temp2[c(5,1,3,4,2),]
    out_simp[rec,1]=k1
    out_simp[rec,2:(1+nrow(temp2))]=temp2[,2]
    out_simp[rec,(2+nrow(temp2)):(1+2*nrow(temp2))]=temp2[,3]   
  }
}
out_simp_2=cbind(out_simp[,1],rownames(out_simp),out_simp[,2:11])
colnames(out_simp_2)[1:2]=c('RDs','SVs')
write.table(out_simp_2,'Integrated.compare.Simp.Het.stat',quote=F,row.names=F)

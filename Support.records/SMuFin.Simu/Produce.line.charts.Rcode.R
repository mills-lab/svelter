cols=c('red','blue','darkgreen','purple','orange','black')
algorithms=c('SVelter','Delly','Lumpy','Pindel','erds','SMuFin')
linetype <- c(1:7) 
plotchar <- seq(18,18+7,1)

pdf('SMuFin.Simu.without.DellyPR.pdf')
par(mfrow=c(2,2))
#sensitivity
par(fig=c(0.05,0.55,0.45,0.95))
par(mar=rep(2,4))
plot(c(1,5),c(0,1),xlab='',main='Germline',bty="n",ylab='Sensitivity',xaxt='n',type='n',)
for(k1 in algorithms){
  temp=data_a[data_a[,10]==k1,]
  temp=temp[order(temp[,9]),]
  color=cols[match(k1,algorithms)]
  ltype=linetype[match(k1,algorithms)]
  ptype=plotchar[match(k1,algorithms)]
  lines(c(1:5),temp[,6],col=color,pch=ptype,lty=ltype,type='b')
}
axis(1,paste('RD',seq(10,50,by=10),sep=''),mgp=c(0.5,0.5,0.1),cex=0.8,at=c(1:5))

par(fig=c(0.5,1,0.45,0.95),new=TRUE)
plot(c(1,5),c(0,1),xlab='',main='Somatic',bty="n",ylab='Sensitivity',xaxt='n',type='n',)
for(k1 in algorithms){
  temp=data_b_new[data_b_new[,10]==k1,]
  temp=temp[order(temp[,9]),]
  color=cols[match(k1,algorithms)]
  ltype=linetype[match(k1,algorithms)]
  ptype=plotchar[match(k1,algorithms)]
  lines(c(1:5),temp[,6],col=color,pch=ptype,lty=ltype,type='b')
}
axis(1,paste('RD',seq(10,50,by=10),sep=''),mgp=c(0.5,0.5,0.1),cex=0.8,at=c(1:5))

#specificity
par(fig=c(0.05,0.55,0,0.5),new=TRUE)
plot(c(1,5),c(0,1),xlab='',bty="n",ylab='Specificity',xaxt='n',type='n',)
for(k1 in algorithms){
  temp=data_a[data_a[,10]==k1,]
  temp=temp[order(temp[,9]),]
  color=cols[match(k1,algorithms)]
  ltype=linetype[match(k1,algorithms)]
  ptype=plotchar[match(k1,algorithms)]
  lines(c(1:5),temp[,7],col=color,pch=ptype,lty=ltype,type='b')
}
axis(1,paste('RD',seq(10,50,by=10),sep=''),mgp=c(0.5,0.5,0.1),cex=0.8,at=c(1:5))

par(fig=c(0.5,1,0,0.5),new=TRUE)
plot(c(1,5),c(0,1),xlab='',bty="n",ylab='Specificity',xaxt='n',type='n',)
for(k1 in algorithms){
  temp=data_b_new[data_b_new[,10]==k1,]
  temp=temp[order(temp[,9]),]
  color=cols[match(k1,algorithms)]
  ltype=linetype[match(k1,algorithms)]
  ptype=plotchar[match(k1,algorithms)]
  lines(c(1:5),temp[,7],col=color,pch=ptype,lty=ltype,type='b')
}
axis(1,paste('RD',seq(10,50,by=10),sep=''),mgp=c(0.5,0.5,0.1),cex=0.8,at=c(1:5))

dev.off()

pdf('SMuFin.Simu.without.DellyPR.legend.pdf')
par(fig=c(0,1,0,1))
plot(c(1,5),c(0,1),xlab='',bty="n",ylab='Specificity',xaxt='n',type='n',)
legend(c(3,5),c(0,0.6),algorithms,col=cols,lty=linetype,pch=plotchar,cex=1)
dev.off()


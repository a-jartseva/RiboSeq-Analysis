bitmap("lll.lengthdist.jpg",type="jpeg",height=4,width=6,res=600,pointsize=10)

mRNAtotal = read.table("lll.lenhist.mRNA.total")

temp = c(0,1)

par(xaxs="i")
par(cex.axis=1.0)

y2 = max(mRNAtotal$V1) #,mRNAuniq$V1)

plot(temp,temp,type="n",ylab="counts",xlab="read length",main="ttt",xlim=c(25,50),ylim=c(0,y2))

lines(mRNAtotal$V2,mRNAtotal$V1,col="red")

lines(c(40,42),0.9*c(y2,y2),col="red")


text(42.5,0.9*y2,"RefSeq mRNA reads",adj=0.0)

dev.off()

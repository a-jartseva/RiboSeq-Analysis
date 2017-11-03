bitmap("lll.framing.jpg",type="jpeg",height=4,width=7,res=600,pointsize=10)

minlen = 25
maxlen = 35

framing = read.table("lll.framing.txt")
m1 = max(framing$V2)
mm = m1
col3 = c(colors()[136],colors()[591],colors()[88])

temp = c(0,1)

par(xaxs="i")
par(cex.axis=1.0)

par(mar=c(0,0,0,0))
plot(temp,type="n",ylab="",xlab="",main="",xlim=c(minlen-0.4,maxlen+1),ylim=c(0,2.4),axes=FALSE)

rect(framing$V1+framing$V3/4,1,framing$V1+(framing$V3+1)/4,1+0.8*framing$V2/mm,col=col3[framing$V3+1],border=FALSE)

text(c(minlen:maxlen)+0.15,2.0,labels=c(minlen:maxlen),adj=0.0,font=2)
text(c(minlen:maxlen)+0.45,2.0,labels="nt",adj=0.0,font=2)

text(0.5*(minlen+maxlen+0.75),2.35,labels="Phasing of reads mapping within host mRNA coding sequences",font=2,cex=1.2,adj=0.5)
text(0.5*(minlen+maxlen+0.75),2.2,labels="ttt",font=2,cex=1.2,adj=0.5)

text(minlen-0.3,1.1,srt=90,labels="total reads",adj=0.0,font=2)

dev.off()

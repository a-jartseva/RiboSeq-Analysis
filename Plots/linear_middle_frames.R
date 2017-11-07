
#Make sure sc is bigger than ylim[2] so that the plots don't get truncated.
#  Be wary of negative reads counts being truncated if sc > ylim[1].
# sc = 8
sc = 8
y2 = c(-2,10)
offset = 12
test_fwd = 0
test_rev = 0
par(mgp=c(0.0,0.4,0.2))
par(tcl=-0.2)
fcol0 = "red"
fcol1 = "green"
fcol2 = "blue"

reads_fwd_in.0 <- read.table("ddd/lll.fwd.0")
reads_fwd_in.1 <- read.table("ddd/lll.fwd.1")
reads_fwd_in.2 <- read.table("ddd/lll.fwd.2")

#------------------------------------------------------------------------------

plot(temp,temp,type="n",ylab="",xlab="",main="",axes=FALSE)
text(0.5,0.5,labels="ttt",adj=0.5,cex=cc1*.8,srt=90)

plot(temp,temp,type="n",ylab="",xlab="",main="",xlim=c(region_starts[1],region_ends[1]),ylim=y2,axes=FALSE)
abline(v=cdsstart,lty="dotted",col=col3)
abline(v=cdsend,lty="dotted",col=col3)
if (zzz) {
  rect(1,0,genomelength,-20,col="white",border=FALSE)
}
abline(h=0,col="grey")

reads.0 <- reads_fwd_in.0
reads.1 <- reads_fwd_in.1
reads.2 <- reads_fwd_in.2

fwdRmax = max(reads.0$V1,reads.1$V1,reads.2$V1)

lines(offset+reads.0$V2,sc*reads.0$V1/fwdRmax,col=fcol0)
lines(offset+reads.1$V2,sc*reads.1$V1/fwdRmax,col=fcol1)
lines(offset+reads.2$V2,sc*reads.2$V1/fwdRmax,col=fcol2)

label = signif(0.5*fwdRmax,digits=1)
axis(side=4,labels=c("0",label),at=sc*c(0,label)/fwdRmax,cex=cc1,las=2)

plot(temp,temp,type="n",ylab="",xlab="",main="",axes=FALSE)
text(0.7,0.35,labels="RPM",cex=cc1,srt=270)


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

if (file.info("ddd/lll.fwd")$size != 0) {
  reads_fwd_in <- read.table("ddd/lll.fwd")
  test_fwd = 1
}
if (file.info("ddd/lll.rev")$size != 0) {
  reads_rev_in <- read.table("ddd/lll.rev")
  test_rev = 1
}

#------------------------------------------------------------------------------

plot(temp,temp,type="n",ylab="",xlab="",main="",axes=FALSE)
text(0.5,0.5,labels="ttt",adj=0.5,cex=cc1*0.8,srt=90)

plot(temp,temp,type="n",ylab="",xlab="",main="",xlim=c(region_starts[1],region_ends[1]),ylim=y2,axes=FALSE)
abline(v=cdsstart,lty="dotted",col=col3)
abline(v=cdsend,lty="dotted",col=col3)
if (zzz) {
  rect(1,0,genomelength,-20,col="white",border=FALSE)
}
if (test_fwd) {
  reads <- reads_fwd_in[offset+reads_fwd_in$V2 >= region_starts[1] & offset+reads_fwd_in$V2 <= region_ends[1], ]
  fwdRmax = max(reads$V1)
  rect(offset+reads$V2,0,offset+reads$V2,sc*reads$V1/fwdRmax,col="red",border=NA)
  if (fwdRmax>0) {
  label = signif(0.5*fwdRmax,digits=1)
  axis(side=4,labels=c("0",label),at=sc*c(0,label)/fwdRmax,cex=cc1,las=2)
  }
}
if (test_rev) {
  reads <- reads_rev_in[offset+reads_rev_in$V2 >= region_starts[1] & offset+reads_rev_in$V2 <= region_ends[1], ]
  rect(offset+reads$V2,0,offset+reads$V2,-sc*reads$V1/fwdRmax,col="blue",border="blue")
}

plot(temp,temp,type="n",ylab="",xlab="",main="",axes=FALSE)
text(0.7,0.35,labels="RPM",cex=cc1,srt=270)

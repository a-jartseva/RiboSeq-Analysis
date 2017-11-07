#ccc = count (augmented for each member of the group to be plotted: 1 = first, 2 = second, etc.)
# if (ccc==1) {
#   o1 = -25 #each histogram will be plotted 1*10 below the previous
# } else {
#   o1 = -25*(ccc*.75)
# }

o1 = -10*ccc #each histogram will be plotted 1*10 below the previous

sc = 1.5

if (file.info("ddd/lll.fwd")$size != 0) {
  reads <- read.table("ddd/lll.fwd")
  rect(reads$V2,o1,reads$V2,o1+log(1+reads$V1)/sc,col="red",border="red") #plotting on the log(x+1) scale
}
if (file.info("ddd/lll.rev")$size != 0) {
  reads <- read.table("ddd/lll.rev")
  rect(reads$V2,o1,reads$V2,o1-log(1+reads$V1)/sc,col="blue",border="blue")
}
text(100,o1-2,labels="ttt",adj=0.0,cex=cc2, pos = 4)
text(1-genomelength/60,o1+0,labels="rev   fwd",srt=90,cex=cc1)
lines(c(1-genomelength/30,0),c(o1,o1))
text(genomelength+genomelength/75,o1+log(1+c(1,10,100,1000))/sc,labels=c("1","10",expression(10^2),expression(10^3)),cex=cc3,adj=0.0)
arrows(c(genomelength+genomelength/200),o1+log(1+c(1,10,100,1000))/sc,c(genomelength+genomelength/100),o1+log(1+c(1,10,100,1000))/sc,length=0)
text(genomelength+(genomelength/17),o1+log(10),labels="RPM",cex=cc1,srt=270)


# rect(1,o1,genomelength,o1-5,col="white",border=FALSE) #space below histograms

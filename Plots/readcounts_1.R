readcounts = read.table("readcounts.summary_pl")

# NOTE: the order in which libraries will be plotted is the REVERSE of the order in
# which they appear in the libraries.txt file that is used when generating the plot data

nsamples = nrow(readcounts)
temp = c(0,1)
par(xaxs="i")
par(cex.axis=1.0)

bitmap("readcounts.summary.jpg",type="jpeg",width=8,height=0.6*(nsamples+10.0),res=600,pointsize=10)

par(mar=c(0,0,0,0))
plot(temp,temp,type="n",ylab="",xlab="",main="",ylim=c(0.18,nsamples+2.0+0.18),xlim=c(-0.05,1.05),axes=FALSE)

#the below rect and text functions plot the legend for the figure, as well as the too-short, adaptor-only, and non-clipped reads for each sample

#All reads - this is the empty (white) rectangle
rect(0,c(1:nsamples)-0.8,readcounts$V2/readcounts$V2,c(1:nsamples)-0.2)

base = 0

#'too-short' reads
rect(base,c(1:nsamples)-0.8,base+readcounts$V3/readcounts$V2,c(1:nsamples)-0.2,col=colors()[19])
base = base + readcounts$V3/readcounts$V2 #number of too-short reads divided by total number of reads (i.e. PROPORTION of  reads that are too-short -- same for adaptor-only below, etc.)

#'adapter-only' reads
rect(base,c(1:nsamples)-0.8,base+readcounts$V4/readcounts$V2,c(1:nsamples)-0.2,col=colors()[15])
base = base + readcounts$V4/readcounts$V2

#'non-clipped' reads
rect(base,c(1:nsamples)-0.8,base+readcounts$V5/readcounts$V2,c(1:nsamples)-0.2,col=colors()[22])
base = base + readcounts$V5/readcounts$V2


text(0.0,c(1:nsamples)-0.07,labels=readcounts$V1,adj=0.0,cex=0.8) #timepoint/sample
text(1.0,c(1:nsamples)-0.07,labels=paste("[",readcounts$V2," reads ]"),adj=1.0,cex=0.8) #total no. of reads
text(0.5,nsamples+2.18,labels="Mapping statistics",adj=0.5,cex=1.2,font=2) #plot title

x1 = 0.1
cc1 = 0.8
y1 = 0.13

rect(0.0+x1,nsamples+1.65+y1,0.1+x1,nsamples+1.65-y1,col=colors()[19])
rect(0.0+x1,nsamples+1.35+y1,0.1+x1,nsamples+1.35-y1,col=colors()[15])
rect(0.0+x1,nsamples+1.05+y1,0.1+x1,nsamples+1.05-y1,col=colors()[22])
rect(0.0+x1,nsamples+0.75+y1,0.1+x1,nsamples+0.75-y1,col=colors()[315])
rect(0.0+x1,nsamples+0.45+y1,0.1+x1,nsamples+0.45-y1,col="white")

text(0.12+x1,nsamples+1.65,adj=0.0,labels="too-short reads",cex=cc1)
text(0.12+x1,nsamples+1.35,adj=0.0,labels="adapter-only reads",cex=cc1)
text(0.12+x1,nsamples+1.05,adj=0.0,labels="non-clipped reads",cex=cc1)
text(0.12+x1,nsamples+0.75,adj=0.0,labels="reverse sense matches",cex=cc1)
text(0.12+x1,nsamples+0.45,adj=0.0,labels="other/unassigned",cex=cc1)

x1 = 0.65
y0 = 1.65
yd = 1.5/nnn
y1 = 0.65/nnn

#------------------------------------------------------------------------------
#Add lines for each database (rRNA, mRNA etc) here:

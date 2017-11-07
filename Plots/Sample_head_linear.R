bitmap("Sample_linear_ggg_xxx.jpg",type="jpeg",width=6.9,height=1.2+nnn*1,res=600,pointsize=12)
# change "Sample" to your virus name, e.g. "HIV"

mar5 = 650
mar3 = 800

cc1 <- 1.0 #(relative) font size


# annotate the genome
genomelength = 8335
cdsstart = c(621,2238,5777) #start positions of gens
cdsend = c(2237,5837,7774) #end positions of genes
cdsy = c(-1,-1,1) # y coordinates at which gene rectangles will be plotted
cdsnames = c("gag","pol","env") # names of genes
namey = c(-1,-1,1) # y positions of gene names

#define the plotting region (by default, the whole genome)
region_starts = c(1)
region_ends = c(genomelength)

r1len = region_ends[1] - region_starts[1]

#initialize plot layout
layout(matrix(c(1:(3*(nnn+2))), nnn+2, 3, byrow = TRUE), c(mar5,r1len,mar3), c(0.5,1+0*c(1:nnn),0.01))

par(mar=c(0,0,0,0))
par(xaxs="i")
par(yaxs="i")
par(cex.axis=0.7)
par(cex=1.0)
cc1 = 0.7
cc2 = 0.7
cc3 = 0.7
col3 = colors()[96]


temp=c(0,1)
x2=c(1-genomelength/20,genomelength+genomelength/20)

par(mar=c(0,0,0,0))

plot(temp,temp,type="n",ylab="",xlab="",main="",axes=FALSE)
text(0.8,0.5,labels=expression(paste("5",symbol("\242"))),adj=1.0,cex=cc1)

plot(temp,temp,type="n",xlab="",ylab="",xaxt="n",yaxt="n",xlim=c(region_starts[1],region_ends[1]),ylim=c(-2.6,2.6),axes=FALSE)

abline(v=cdsstart,lty="dotted",col=col3)
abline(v=cdsend,lty="dotted",col=col3)


# black rect for genome
rect(1,-0.3,genomelength,0.3,col="black",border="black")

# rects for genes
rect(cdsstart,cdsy-0.9,cdsend,cdsy+0.9,col="light blue",border="black",lwd=1.4)
text(0.5*(cdsstart+cdsend),cdsnamey,labels=cdsnames,cex=cc1)

plot(temp,temp,type="n",ylab="",xlab="",main="",axes=FALSE)
text(0.2,0.5,labels=expression(paste("3",symbol("\242"))),adj=0.0,cex=cc1)

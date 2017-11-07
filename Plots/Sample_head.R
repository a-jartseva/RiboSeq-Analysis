bitmap("Sample_log_ggg_xxx.jpg",type="jpeg",width=6.9,height=2.8+nnn*1,res=600,pointsize=12)
# change "Sample" to your virus name, e.g. "HIV"

library(RColorBrewer)
brewer_Set2_pallette<-brewer.pal(8,"Set2")

par(xaxs="i")
par(yaxs="i")
par(cex.axis=0.7)
par(cex=1.0)
cc1 = 0.7
cc2 = 0.7
cc3 = 0.7
col3 = colors()[96]

# annotate the genome
genomelength = 8335
cdsstart = c(621,2238,5777) #start positions of gens
cdsend = c(2237,5837,7774) #end positions of genes
cdsy = c(-1,-1,1) # y coordinates at which gene rectangles will be plotted
cdsnames = c("gag","pol","env") # names of genes
namey = c(-1,-1,1) # y positions of gene names

x=c(1,genomelength)
temp=c(0,1)
par(xaxs="i")
x2=c(1-genomelength/20,genomelength+genomelength/20)
cc1 = 0.8

par(mar=c(0,0,0,0))
plot(temp,temp,type="n",xlab="",ylab="",xaxt="n",yaxt="n",xlim=x2,ylim=c(-2.6-10*nnn,2.6),axes=FALSE)

abline(v=cdsstart,lty="dotted",col=col3)
abline(v=cdsend,lty="dotted",col=col3)
rect(1,0,genomelength,5,col="white",border=FALSE)
rect(1,-0.3,genomelength,0.3,col="black",border="black")
rect(cdsstart,cdsy-0.9,cdsend,cdsy+0.9,col="light blue",border="black",lwd=1.4)
text((cdsstart+cdsend)/2,namey,labels=cdsnames,adj=0.5,cex=cc1)

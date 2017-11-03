#this script is appended to readcounts2.R, and plots the forward and reverse strand reads for each category of RNA

#ddd
fwdrev = readcounts$Vaaa+readcounts$Vbbb #total number of reads
rect(base,c(1:nsamples)-0.8,base+fwdrev/readcounts$V2,c(1:nsamples)-0.2,col=colors()[ccc])
revfrac = readcounts$Vbbb/fwdrev
rect(base,c(1:nsamples)-0.8,base+fwdrev/readcounts$V2,c(1:nsamples)-0.8+0.6*revfrac,col=colors()[315])
base = base + fwdrev/readcounts$V2
rect(0.0+x1,nsamples+y0-xxx*yd+y1,0.1+x1,nsamples+y0-xxx*yd-y1,col=colors()[ccc])
text(0.12+x1,nsamples+y0-xxx*yd,adj=0.0,labels="ddd",cex=cc1)

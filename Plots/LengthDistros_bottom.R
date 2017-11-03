#----------------------------------------------------------------------------

blank <- grid.rect(gp=gpar(col="white"))

nsamples <- nnn
heights_vec <-c()

for (i in 1:ceiling(nsamples/2)) {
  heights_vec<-c(heights_vec,1)
}

png(filename = "LD_virus_host_combined.png",width=1000,height=200*nsamples,res=150,pointsize=12)

# if there is an even number of samples, put the legend underneath in a row of its own
# if there is an odd number of samples, the legend can go on the bottom row, right hand side.

if (nsamples%%2==0) {
  plotslist[[nsamples+1]] <- blank
  plotslist[[nsamples+2]] <- legend
  grid.arrange(grobs=plotslist, widths=c(1,1), heights=c(heights_vec,0.2), ncol=2, nrow = (nsamples/2)+1)
} else {
  plotslist[[nsamples+1]] <- legend
  grid.arrange(grobs=plotslist, widths=c(1,1), heights=heights_vec, ncol=2, nrow = ceiling(nsamples/2))
}

dev.off()

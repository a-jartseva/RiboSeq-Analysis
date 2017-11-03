#----------------------------------------------------------------------------

blank <- grid.rect(gp=gpar(col="white"))

png(filename = "LD_virus_host_combined.png",width=1000,height=1200,res=150,pointsize=12)
grid.arrange(l2,l3,l4,l5,l1,legend, widths=c(1,1), heights=c(1,1,1), ncol=2, nrow = 3)
dev.off()

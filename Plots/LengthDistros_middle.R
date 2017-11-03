
lengths_table = read.table(paste("lll",".len.percs.combined",sep=""),header=T)

lengths_melt = melt(lengths_table,id="Length_nt")

ymax <- max(lengths_melt$value)

ppp<-ggplot(lengths_melt,  aes(x=Length_nt, y=value, group=variable, color=variable)) +
  geom_line(size=1.1) + scale_color_manual("",values=colours) +
  xlab("Length (nt)") + ylab("% of reads") + ylim(0,50) + xlim(25,50) +
  theme(text = element_text(size=cc1)) +
  theme(axis.text.x = element_text(colour = "black",size=cc1),
  axis.text.y = element_text(colour = "black",size=cc1)) +
  theme(axis.title.x = element_text(margin=margin(5,0,0,0)),
  axis.title.y=element_text(margin=margin(0,15,0,0), size=cc1)) +
  theme(panel.grid.minor = element_blank(),
  panel.grid.major = element_blank(),
  panel.background = element_blank(),
  panel.border = element_rect(fill = NA, colour = "black",size = 1))+
  theme(legend.position="none") +
  annotate(geom="text", x=42, y=45, label="lll",color="black", size=5)

plotslist[[nnn]] <- ppp
#extract the legend

leg_plot<-ggplot(lengths_melt,  aes(x=Length_nt, y=value, group=variable, color=variable)) +
  geom_line(size=1.1) + scale_color_manual("",values=colours) +
  theme(text = element_text(size=cc1)) +
  theme(legend.position="bottom")

g <- ggplotGrob(leg_plot)$grobs
legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]

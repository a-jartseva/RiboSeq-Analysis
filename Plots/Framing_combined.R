library(ggplot2)
library(reshape2)
library(gridExtra)
library(grid)

cc1=10

#stick with the same set of colours for each plot
colours = c(colors()[26], colors()[257], colors()[503])

framing_table <- read.table("framing_summarised.txt",header=T)

framing_melt = melt(framing_table,id="library")

virus_framing = subset(framing_melt,variable=)

virus_framing <- framing_melt[grep("virus", framing_melt$variable), ]
host_framing <- framing_melt[grep("host", framing_melt$variable), ]

p1_ymax <- max(virus_framing$value)
p2_ymax <- max(host_framing$value)

xlim <- length(framing_table$library)

p1<-ggplot(virus_framing,  aes(x=library, y=value,fill=variable)) +
  geom_bar(position="dodge", stat="identity") +
  scale_fill_manual("Codon position",values=colours,labels=c("1 (phase 0)","2 (phase +1)","3 (phase +2)")) +
  ylab("% of RPFs") + ylim(0,100) +
  theme(text = element_text(size=cc1)) +
  theme(axis.title.x = element_blank()) +
  theme(axis.title.y = element_blank()) +
  theme(axis.ticks.x = element_blank()) +
  theme(axis.text.x=element_blank())+
  theme(axis.text.y = element_text(colour = "black",size=cc1)) +
  # theme(axis.text.x = element_text(angle=45,hjust=1,colour = "black",size=cc1),
  # axis.text.y = element_text(colour = "black",size=cc1)) +
  theme(panel.grid.minor = element_blank(),
  panel.grid.major = element_blank(),
  panel.background = element_blank(),
  panel.border = element_rect(fill = NA, colour = "black",size = 1)) +
  theme(legend.position="none") +
  scale_x_discrete(labels=c("4 hpi rep. 1","4 hpi rep. 2",
    "6 hpi rep. 1", "6 hpi rep. 2")) +
  theme(plot.margin=unit(c(0.2,0.05,0.05,0), "cm")) +
  annotate("text", label = "vRNA", x = xlim-(xlim/50), y = 95,size=cc1*0.4)

p2<-ggplot(host_framing,  aes(x=library, y=value,fill=variable)) +
  geom_bar(position="dodge", stat="identity") +
  scale_fill_manual("Codon position",values=colours,labels=c("1 (phase 0)","2 (phase +1)","3 (phase +2)")) +
  ylab("% of RPFs") + ylim(0,100) +
  theme(text = element_text(size=cc1)) +
  theme(axis.title.x = element_blank()) +
  theme(axis.title.y = element_blank()) +
  theme(axis.text.x = element_text(angle=45,hjust=1,colour = "black",size=cc1),
  axis.text.y = element_text(colour = "black",size=cc1)) +
  theme(panel.grid.minor = element_blank(),
  panel.grid.major = element_blank(),
  panel.background = element_blank(),
  panel.border = element_rect(fill = NA, colour = "black",size = 1)) +
  theme(legend.position="none") +
  scale_x_discrete(labels=c("4 hpi rep. 1","4 hpi rep. 2",
    "6 hpi rep. 1", "6 hpi rep. 2"))  +
  theme(plot.margin=unit(c(0.05,0.05,0.2,0), "cm")) +
  annotate("text", label = "mRNA", x = xlim-(xlim/50), y = 95,size=cc1*0.4)

leg_plot<-ggplot(virus_framing,  aes(x=library, y=value,fill=variable)) +
  geom_bar(position="dodge", stat="identity") +
  scale_fill_manual("Codon position",values=colours,labels=c("1 (phase 0)","2 (phase +1)","3 (phase +2)")) +
  ylab("% of RPFs") + ylim(0,p1_ymax) +
  theme(text = element_text(size=cc1)) +
  theme(axis.title.x = element_blank()) +
  theme(axis.title.y = element_blank()) +
  theme(axis.text.x = element_text(angle=45,hjust=1,colour = "black",size=cc1),
  axis.text.y = element_text(colour = "black",size=cc1)) +
  theme(panel.grid.minor = element_blank(),
  panel.grid.major = element_blank(),
  panel.background = element_blank(),
  panel.border = element_rect(fill = NA, colour = "black",size = 1)) +
  theme(legend.position="right") +
  scale_x_discrete(labels=c("4 hpi rep. 1","4 hpi rep. 2",
    "6 hpi rep. 1", "6 hpi rep. 2"))

#extract the legend
g <- ggplotGrob(leg_plot)$grobs
legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]

blank <- grid.rect(gp=gpar(col="white"))

png(filename = "Combined.phasing.png",width=600,height=550,res=150,pointsize=12)

grid.arrange(p1,blank,p2,legend,
   left=textGrob("% of RPFs", gp=gpar(fontsize=cc1), rot=90,y=0.55),
   widths=c(1,0.7), heights=c(0.68,1), ncol=2, nrow = 2)

dev.off()

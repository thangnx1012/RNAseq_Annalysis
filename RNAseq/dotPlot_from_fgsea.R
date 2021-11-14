library(ggplot2)

dat <- read.delim("fgsea_output.txt")
dat$combo <- paste(dat$pathway,".",dat$CellLine)
ggplot(dat[,1:4],aes(x=CellLine, 
           y=reorder(pathway,NES), 
           colour=NES)) +
  geom_point(aes(colour=NES,size=padj)) +
  scale_size_continuous(range = c(1.0,6.0), name = "padj",limits = c(1.3,3.5)) +
  scale_color_gradient(low = "white", high = "firebrick",limits = c(1.4,2.65))+
  labs(x="Cell Line", y="Pathway", colour="NES", size="- log10 padj")

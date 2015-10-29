load(paste0(dump.dir,"alt.MGMT.Rdata"))
#head(alt.MGMT1)
library(reshape2)
library(ggplot2)

alt.MGMT1[,2]=c(rep(0,5),rep(c(rep("20 %",5),rep("40 %",5),rep("60 %",5),rep("90 %",5)),4))
alt.MGMT1[,1]=c(rep("No Action",5),rep("3 Year",20),rep("1 Year",20),rep("Sterilize",20),rep("Cull",20))

colnames(alt.MGMT1) <- c("Treatment","Percent","Year","P<","Pin","P>","Number")

plot1 <- alt.MGMT1[6:nrow(alt.MGMT1),1:3]
plot2 <- alt.MGMT1[6:nrow(alt.MGMT1),7]
plot3 <- cbind(plot1,plot2)
melt_dat=melt(plot3,id.var=c("Treatment","Percent","Year"))
#head(melt_dat)

if(DRAW==TRUE) pdf(paste(dump.dir,"cats_barplot.pdf",sep=""))
ggplot(melt_dat, aes(x = variable, y = value, fill = as.factor(Year))) + 
  geom_bar(position="dodge",stat = "identity") +
  facet_wrap(~Treatment+Percent,nrow=4) + 
  theme(axis.text.y = element_text(colour="black"),
        legend.title = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.y = element_line("black",size = 1),
        axis.ticks.x = element_blank()) +
  ylab("Number treated") + 
  xlab("Year") +
  scale_fill_grey(name=NULL)
if(DRAW==TRUE) dev.off()

load(paste(dump.dir,"alt.MGMT1.cull1st.Rdata",sep=""))
#head(alt.MGMT1.cull1st)
alt.MGMT1.cull1st[,2]=c(rep(c(rep("20 %",4),rep("40 %",4),rep("60 %",4),rep("90 %",4)),4))
alt.MGMT1.cull1st[,1]=c(rep("3 Year",16),rep("1 Year",16),rep("Sterilize",16),rep("Cull",16))

colnames(alt.MGMT1.cull1st) <- c("Treatment","Percent","Year","P<","Pin","P>","Number")

melt_dat1=melt(cbind(alt.MGMT1.cull1st[,1:3],alt.MGMT1.cull1st[,7]),
              id.var=c("Treatment","Percent","Year"))
#head(melt_dat1)

if(DRAW==TRUE) pdf(paste(dump.dir,"cats_barplot_cull1st.pdf",sep=""))
ggplot(melt_dat1, aes(x = variable, y = value, fill = as.factor(Year))) + 
  geom_bar(position="dodge",stat = "identity") +
  facet_wrap(~Treatment+Percent,nrow=4) + 
  theme(axis.text.y = element_text(colour="black"),
        legend.title = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.y = element_line("black",size = 1),
        axis.ticks.x = element_blank()) +ylab("Number treated") + 
  xlab("Year") +
  scale_fill_grey(name=NULL)
if(DRAW==TRUE) dev.off()

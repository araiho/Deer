library(ggplot2)

test.stat<-rnorm(1000000,0)
dd<- with(density(test.stat),data.frame(x,y))
tau.hat <- -.85
tau.hat.1a <- -.85

a<-ggplot(data = dd, mapping = aes(x=x, y=y)) +
  geom_line(color="red",size=1) + layer(data = dd,mapping = aes(x=ifelse(x<tau.hat,x,tau.hat),
                                   y=y), geom = "area", geom_params=list(fill="red",alpha=.3))+
  scale_y_continuous(limits = c(0,max(dd$y)),name="Probability density") +
  geom_vline(aes(xintercept=tau.hat.1a),color="black",linetype="solid") +
  theme(axis.text.x = element_blank())+
  xlab("")+ scale_x_continuous(limits = c(-5, 5))

test.stat1<-rnorm(1000000,-1)
dd1<- with(density(test.stat1),data.frame(x,y))

b<-ggplot(data = dd1, mapping = aes(x=x, y=y)) +
  geom_line(color="blue",size=1) + layer(data = dd1,mapping = aes(x=ifelse(x<tau.hat,x,tau.hat),
                                                                y=y), geom = "area", geom_params=list(fill="blue",alpha=.3))+
  scale_y_continuous(limits = c(0,max(dd1$y)),name="Probability density") +
  geom_vline(aes(xintercept=tau.hat.1a),color="black",linetype="solid") +
  theme(axis.text.x = element_blank())+
  xlab("")+ scale_x_continuous(limits = c(-5, 5))

dd2 <- rbind(dd,dd1)
group<-rep(c("A","B"),each=512)
dd2 <- as.data.frame(cbind(dd2,group))

c<-ggplot(data = dd2, mapping = aes(x=x, y=y, group = group,colour = group)) +
  geom_line(size = 1) + scale_color_manual(values=c("red", "blue")) + scale_x_continuous(limits = c(-5, 5)) +
  scale_y_continuous(limits = c(0,max(dd2$y)),name="Probability density") +
  theme(axis.text.x = element_blank(),legend.position = "none")+
  xlab("") +
  layer(data = dd2,mapping = aes(x = ifelse(x < tau.hat, x, tau.hat),
                                 y = y), geom = "area", geom_params=list(fill="blue",alpha=.3)) +
  
  layer(data = dd2[dd2$group=="A",],mapping = aes(x = ifelse(x < tau.hat, x, tau.hat),
                                                  y = y), geom = "area", geom_params=list(fill="red",alpha=.5)) +
  geom_vline(aes(xintercept=tau.hat.1a),color="black",linetype="solid")


multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

quartz()

multiplot(a,b,c,cols=1)

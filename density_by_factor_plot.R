for(i in levels(set3_2$SetName)) {
  print(i)
  lapply(set3[set3_2$SetName==i,], denisty) -> dd
  xr=quantile(unlist(lapply(dd, function(x) x$x)), probs=c(.10,.90))
  yr=quantile(unlist(lapply(dd, function(x) x$y)), probs=c(.10,.90))
  plot(0, xlim = xr, ylim =yr, type="n",main = "", ylab="density")
  lapply(dd, function(x) (lines(x, xlim=xr, ylim=yr)))
}

fil= names(which(apply(set3, 1, mad)>0))
xx<-apply(set3[fil,], 1, function(x) (x-median(x, na.rm=TRUE))/(mad(x, na.rm=TRUE)+0.1))

xx2<-cbind(set3_2[fil,1:2], t(xx))
setM<-(melt((xx2)))

p <- ggplot(data=setM ,aes(x=value))+
      geom_density(aes(x = value,color=SetName)) +
      facet_grid(Source~., scales="free")+
      theme(legend.position="none")
p
library(plotly)
p <- ggplotly(p)


# Create a shareable link to your chart
# Set up API credentials: https://plot.ly/r/getting-started
chart_link = plotly_POST(p, filename="geom_density/estimate")
chart_link

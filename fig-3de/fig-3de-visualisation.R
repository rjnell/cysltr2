

##########################################

snp_data = read.csv("fig-3de/cysltr2.tsv", stringsAsFactors = F, check.names = F, sep="\t")
colors = c("#F44336","#2196F3","#E67E22","#4CAF50")
names(colors) = c("T","C","G","A")

samples = unique(snp_data$sample)

for (sample in samples) {
  case = substr(sample,6,12)
  snp_data_sample = snp_data[which(snp_data$sample==sample),]
  
  result_file = paste0("fig-3de/res/plots/", case, ".png")
  png(result_file, res=600, width=3000, height=2000)  
  #par(mar=c(5.1, 4.5, 4.1, 2.1))
  par(mar=c(5,5,5,5))
  ylim = c(0,1)
  xlim = c(0,10)
  plot(type="n",
       x = xlim,
       y = ylim,
       axes = F,
       ylab="",
       xlab="",
       main = "",
       bty = "l", 
       xaxs = "i", 
       yaxs = "i")
  
  yat = seq(from=ylim[1],to=ylim[2],by=0.5)
  segments(xlim[1],yat,nrow(snp_data_sample)+1,col="#eeeeee",lwd=1.4,xpd=T)
  
  #segments(xat,ylim[1],xat,ylim[2],col="#eeeeee",lwd=1.4,xpd=T)
  segments(xlim[1],ylim[1],xlim[1],ylim[2],xpd=T,col="#B1B1B1",lwd=1.4)
  
  axis(side = 2,at=yat,las=2,labels=rep("",length(yat)),col.ticks = "#b1b1b1",col = "#b1b1b1", tck = -0.053,lwd=1.4)
  axis(side = 2,at=yat,las=2,labels=c("0%","50%","100%"), lwd=0, col.axis="#333333",line=-0.23)
  mtext(text = "% alleles (RNA)", side = 2, line=3.8,col="#333333")  
  
  for (i in 1:nrow(snp_data_sample)) {
    ro = snp_data_sample$RO[i]
    ao = snp_data_sample$AO[i]
    p1 = ro/(ro+ao)
    x=i-0.5
    rect(x+0.1,0,x+0.9,p1,border="white",col="#333333", lwd=1.4, xpd=T) #colors[snp_data_sample$REF[i]]
    rect(x+0.1,p1,x+0.9,1,border="white",col="#b1b1b1", lwd=1.4, xpd=T) #colors[snp_data_sample$ALT[i]]
    text(x+0.5, 0.925, col="#333333", labels=snp_data_sample$ALT[i], xpd=T,cex=0.95)  
    text(x+0.5, 0.075, col="#ffffff", labels=snp_data_sample$REF[i], xpd=T,cex=0.95)  
    
    ci = Hmisc::binconf(ro, ro+ao, 0.05)[2:3]
    arrows(x+0.5, ci[1], x+0.5, ci[2], length=0.05, angle=90, code=3, col="white", xpd=T, lwd=1.4)
    text(x+0.5, -0.1, srt=60, adj=1, col="#333333", labels=snp_data_sample$ID[i], xpd=T)  
  }
  
  segments(0.5,c(0,1),nrow(snp_data_sample)+1,col="#eeeeee",lwd=1.4,xpd=T)
  segments(xlim[1],0,nrow(snp_data_sample)+1,col="#b1b1b1",lwd=1.4,xpd=T)
  text(nrow(snp_data_sample)/2+0.5, 1.2, col="#333333", labels=case, xpd=T, font=2)  
  
  dev.off()
}








result_file = paste0("fig-3abc/fig-3c.png")
png(result_file, res=600, width=4400, height=2000)  
par(mar=c(5,5,5,5))
ylim = c(0,1)
xlim = c(0,18)
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
segments(xlim[1],yat,2+1,col="#eeeeee",lwd=1.4,xpd=T)

#segments(xat,ylim[1],xat,ylim[2],col="#eeeeee",lwd=1.4,xpd=T)
segments(xlim[1],ylim[1],xlim[1],ylim[2],xpd=T,col="#B1B1B1",lwd=1.4)

axis(side = 2,at=yat,las=2,labels=rep("",length(yat)),col.ticks = "#b1b1b1",col = "#b1b1b1", tck = -0.053,lwd=1.4)
axis(side = 2,at=yat,las=2,labels=c("0%","50%","100%"), lwd=0, col.axis="#333333",line=-0.23)
mtext(text = "% alleles", side = 2, line=3.8,col="#333333")  

#####
i = 1
data = c(.409836,.36125,.459)

x=i-0.5
rect(x+0.1,0,x+0.9,data[1],border="white",col="#333333", lwd=1.4) 
rect(x+0.1,data[1],x+0.9,1,border="white",col="#b1b1b1", lwd=1.4) 
arrows(x+0.5, data[2], x+0.5, data[3], length=0.05, angle=90, code=3, col="white", xpd=T, lwd=1.4)
text(x+0.5, -0.1, srt=60, adj=1, col="#333333", labels="DNA", xpd=T)  

i = 2
mo = 166
ro = 3
p1 = mo/(mo+ro)
x=i-0.5
rect(x+0.1,0,x+0.9,p1,border="white",col="#333333", lwd=1.4) 
rect(x+0.1,p1,x+0.9,1,border="white",col="#b1b1b1", lwd=1.4) 
ci = Hmisc::binconf(mo, mo+ro, 0.05)[2:3]
arrows(x+0.5, ci[1], x+0.5, ci[2], length=0.05, angle=90, code=3, col="white", xpd=T, lwd=1.4)
text(x+0.5, -0.1, srt=60, adj=1, col="#333333", labels="RNA", xpd=T)  
text(1.5, 1.2, col="#333333", labels="MUM-20035", xpd=T, font=2)  

segments(xlim[1],1,2+1,col="#eeeeee",lwd=1.4,xpd=T)
segments(xlim[1],0,2+1,col="#b1b1b1",lwd=1.4,xpd=T)
segments(c(0),0,c(0),1, col="#b1b1b1", lwd=1.4)

y = .5625+0.03125
rect(x+2.6,y-0.0625,x+3.4,y+0.0625,col="#333333", border=NA,xpd=T)
text(x+3.4,y,col="#333333", labels = substitute(paste(italic('CYSLTR2')^"L129Q")), pos=4, xpd=T)

y=0.4375-0.03125
rect(x+2.6,y-0.0625,x+3.4,y+0.0625,col="#b1b1b1", border=NA,xpd=T)
text(x+3.4,y,col="#333333", labels = substitute(paste(italic('CYSLTR2')^"WT")), pos=4, xpd=T)


dev.off()



# Arm-level events
arm_data_file = "fig-s2/transformed.cor.cli.txt"
arm_data = read.csv(arm_data_file, sep="\t", check.names = F, stringsAsFactors = F, row.names = 1)

exclude = c("TCGA-V4-A9ED","TCGA-VD-AA8O","TCGA-YZ-A982")
`%notin%` <- Negate(`%in%`)

arm_data_selected = arm_data[,which(colnames(arm_data) %notin% exclude)]
n = ncol(arm_data_selected)
chrs = c(1:22,"X")

result_file <<- paste0("fig-s2/fig-s2.png")
png(result_file,res=600,width=8000,height=2000)  
par(mfrow=c(1,1))
par(mar=c(5.1, 4.5, 4.1, 2.1))

par(mfrow=c(1,1),mar=c(5,5,5,5))
par(mar=c(5.1, 5, 4.1, 2.1))
ylim= c(0,1)
xlim = c(0.5,75)
plot(x = xlim,
     y = ylim,
     pch = 16,
     xlab = "",
     ylab = "", 
     ylim = ylim,
     xlim = xlim, 
     # asp=1,
     bty = "l", type="n",
     xaxs = "i", yaxs = "i",
     axes = F)
yat = seq(from=ylim[1],to=ylim[2],by=0.25)
segments(0,yat,64.5,col="#eeeeee",lwd=1.4,xpd=T)
segments(0+0.5,-0.05,0+0.5,1,xpd=T,col="#B1B1B1",lwd=1.4)

axis(side = 2,at=yat,las=2,labels=rep("",length(yat)),col.ticks = "#b1b1b1",col = "#b1b1b1", tck = -0.047,lwd=1.4)
axis(side = 2,at=yat,las=2,labels=paste0(yat*100,"%"), lwd=0, col.axis="#333333",line=-0.23)
mtext(text = "% cases", side = 2, line=3.8,col="#333333")  

colors = c("#c0392b","#3d98d3","#ffffff",rgb(0,150/255,0))

x = 1
y = -0.15
for (chr in chrs) {
  
  if (chr %notin% c(13,14,15,21,22)) {
    ploss = length(which(arm_data_selected[paste0(chr,"P loss mutation analysis"),] == paste0(chr,"P loss mutated")))/n
    pgain = length(which(arm_data_selected[paste0(chr,"P gain mutation analysis"),] == paste0(chr,"P gain mutated")))/n
    rect(x,0,x+1,pgain,col=colors[1],border="white",lwd=0.7)
    rect(x,pgain,x+1,pgain+ploss,col=colors[2],border="white",lwd=0.7) 
    text(x+.5, y, labels=paste0(chr,"p"),srt=60,adj=1,xpd=T,col="#333333",cex=0.75)
    text(x+1.5, y, labels=paste0(chr,"q"),srt=60,adj=1,xpd=T,col="#333333",cex=0.75)
    x = x+1
  }
  else{
    #text(x+.5, y, labels=paste0(chr),xpd=T,col="#333333")
    text(x+.5, y, labels=paste0(chr,"q"),srt=60,adj=1,xpd=T,col="#333333",cex=0.75)
  }
  
  qloss = length(which(arm_data_selected[paste0(chr,"Q loss mutation analysis"),] == paste0(chr,"Q loss mutated")))/n
  qgain = length(which(arm_data_selected[paste0(chr,"Q gain mutation analysis"),] == paste0(chr,"Q gain mutated")))/n
  rect(x,0,x+1,qgain,col=colors[1],border="white",lwd=0.7)
  rect(x,qgain,x+1,qgain+qloss,col=colors[2],border="white",lwd=0.7)
  x = x+1
  
  segments(x+0.5,-0.05,x+0.5,1,xpd=T,col="#B1B1B1",lwd=1.4)
  x=x+1
}
segments(0,0,64.5,0,xpd=T,col="#B1B1B1",lwd=1.4)
x = x+1
y=0.575
rect(x, y+0.05, x+1, y-0.05, col=colors[1], border="white", xpd=T,lwd=0.7)
text(x+0.8, y, pos=4, labels="CNA gain", xpd=T)

y=0.425
rect(x, y+0.05, x+1, y-0.05, col=colors[2], border="white", xpd=T,lwd=0.7)
text(x+0.8, y, pos=4, labels="CNA loss", xpd=T)

dev.off()
system(paste("open",result_file))

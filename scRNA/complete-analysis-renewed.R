# Initialize Seurat library
library(Seurat)

###
### STEP 3 - Open scaled data and determine % CYSLTR2 expressing cells per sample (Figure 4B) ###
###

# Retrieve scaled data
um_data_cysltr2_scaled = readRDS("scRNA/um_data_scaled.RDS")

# Separate per case
case_data = list()
for (case in c("UMM059","UMM061","UMM062","UMM065","UMM064","UMM066","UMM069","UMM063")) {
  case_data[[case]] = subset(um_data_cysltr2_scaled, cells = which(substr(colnames(um_data_cysltr2_scaled@assays$RNA@counts),1,6) %in% case))
}

# Calculate statistics: numbers of cells (melanocytes or other) expressing CYSLTR2
perc_cysltr2_cells = NULL
for (case in c("UMM059","UMM061","UMM062","UMM065","UMM064","UMM066","UMM069","UMM063")) {
  total = ncol(case_data[[case]]@assays$RNA@data)
  cysltr2_mel = length(which(case_data[[case]]@assays$RNA@data["CYSLTR2",] > 0 & case_data[[case]]@assays$RNA@data["MLANA",] > 0))
  cysltr2_oth = length(which(case_data[[case]]@assays$RNA@data["CYSLTR2",] > 0 & case_data[[case]]@assays$RNA@data["MLANA",] == 0))
  perc_cysltr2_cells = rbind(perc_cysltr2_cells, c(case, total, cysltr2_mel, cysltr2_oth, cysltr2_mel/total, cysltr2_oth/total))
}
colnames(perc_cysltr2_cells) = c("Case", 
                                 "Total number of cells", 
                                 "Number of CYSLTR2 expressing MLANA+ cells", 
                                 "Number of CYSLTR2 expression MLANA- cells", 
                                 "Perc CYSLTR2+/MLANA+ vs total",
                                 "Perc CYSLTR2+/MLANA- vs total")
perc_cysltr2_cells = perc_cysltr2_cells[order(perc_cysltr2_cells[,5]),]

# Generate Figure 4B
png("scRNA/figs/Figure-4B.png", res=600, width=4400, height=2000)  
{
  par(mar=c(5,5,5,5))
  ylim = c(0,.4)
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
  
  yat = seq(from=ylim[1],to=ylim[2],by=0.1)
  segments(xlim[1],yat,nrow(perc_cysltr2_cells)+1,col="#eeeeee",lwd=1.4,xpd=T)
  segments(xlim[1],ylim[1],xlim[1],ylim[2],xpd=T,col="#B1B1B1",lwd=1.4)
  
  axis(side = 2,at=yat,las=2,labels=rep("",length(yat)),col.ticks = "#b1b1b1",col = "#b1b1b1", tck = -0.053,lwd=1.4)
  axis(side = 2,at=yat,las=2,labels=c("0%","10%","20%","30%","40%"), lwd=0, col.axis="#333333",line=-0.23)
  mtext(text = "% cells", side = 2, line=3.8,col="#333333")  
  
  for (i in 1:nrow(perc_cysltr2_cells)) {
    p1 = as.numeric(perc_cysltr2_cells[i,5])
    p2 = as.numeric(perc_cysltr2_cells[i,6])
    x=i-0.5
    rect(x+0.1,0,x+0.9,p1,border="white",col=RColorBrewer::brewer.pal(9,"Set1")[1], lwd=1.4)
    rect(x+0.1,p1,x+0.9,p1+p2,border="white",col=RColorBrewer::brewer.pal(9,"Set1")[2], lwd=1.4)
    
    text(x+0.5, -0.1/4, srt=60, adj=1, col="#333333", labels=perc_cysltr2_cells[i,1], xpd=T)  
  }
  
  segments(0.5,c(0,1),nrow(perc_cysltr2_cells)+1,col="#eeeeee",lwd=1.4,xpd=T)
  segments(xlim[1],ylim[1],nrow(perc_cysltr2_cells)+1,ylim[1],xpd=T,col="#B1B1B1",lwd=1.4)
  text(nrow(perc_cysltr2_cells)/2+1, 1.2, col="#333333", labels=case, xpd=T, font=2)  
  
  y = .225+0.0125
  rect(x+2.6,y-.025,x+3.4,y+.025,col=RColorBrewer::brewer.pal(9,"Set1")[1], border=NA)
  text(x+3.4,y,col="#333333", labels = substitute(paste(italic('CYSLTR2')^"+", italic('MLANA')^"+")), pos=4)
  
  y=0.175-0.0125
  rect(x+2.6,y-.025,x+3.4,y+.025,col=RColorBrewer::brewer.pal(9,"Set1")[2], border=NA)
  text(x+3.4,y,col="#333333", labels = substitute(paste(italic('CYSLTR2')^"+", italic('MLANA')^"−")), pos=4)
}
dev.off()



###
### STEP 4 - Generate TSNE plot with expression of relevant genes per tumour ###
###

# Initialize Seurat library
library(Seurat)

# Iterate through cases, run PCA, find neighbours and clusters and run TSNE
for (case in c("UMM059","UMM061","UMM062","UMM065","UMM064","UMM066","UMM069","UMM063")) {
  case_data[[case]] = RunPCA(case_data[[case]], features = VariableFeatures(object = case_data[[case]]))
  case_data[[case]] = FindNeighbors(case_data[[case]], dims = 1:10)
  case_data[[case]] = FindClusters(case_data[[case]], resolution = 0.5)
  case_data[[case]] = RunTSNE(case_data[[case]], dims = 1:10)
}

# Gene list to plot
genes = c("CYSLTR2","MLANA","CD14","CD3D","ALOX5AP","MITF","TYR","ETV5","DUSP4","ALOX5") 
 
# Palette 'Blues' 
palette = colorRampPalette(c("#DDDDDD",RColorBrewer::brewer.pal(9,"Blues")[4:9]))(101)

# Iterate through cases and generate multiplot TSNE with gene expression
for (case in c("UMM059","UMM061","UMM062","UMM065","UMM064","UMM066","UMM069","UMM063")) {
  graphics.off()
  
  # Multiplot
  png(paste0("scRNA/figs/Figure-S3-",case,"-overview.png"),res=600,width=4000,height=2000)
  {
    par(mfrow=c(2,5), mar=c(0.5,0.5,3.5,0.5))
    for(gene in genes) {
      tsne = case_data[[case]]@reductions$tsne@cell.embeddings
      plot(tsne, type="n", axes=F, main="", xlab="", ylab="")
      title(main = gene, font.main=3)
      ge = case_data[[case]]@assays$RNA@data[gene,]
      points(tsne, pch=16, cex=0.4, col=palette[floor(ge/max(ge)*100)+1])
    }
  }
  dev.off()
}



###
### STEP 5 - Generate TSNE plot with expression of selected genes for highlighted tumour UMM066 ###
###

# Select highlighted case
case = "UMM066"

# Select highlighted genes
genes = c("CYSLTR2","MITF","ETV5","ALOX5AP")

# Iterate through genes to plot expression
for (gene in genes) {
  png(paste0("scRNA/figs/Figure-4C-",case,"-",gene,".png"),res=600,width=2750,height=3000)
  {
    tsne = case_data[[case]]@reductions$tsne@cell.embeddings
    plot(tsne, type="n", axes=F, main="", xlab="", ylab="")
    ge = case_data[[case]]@assays$RNA@data[gene,]
    points(tsne, pch=16, cex=0.4, col=palette[floor(ge/max(ge)*100)+1])
  }
  dev.off()
}
graphics.off()

# DimPlot of TSNE to annotate clusters
png(paste0("scRNA/figs/Supporting-Figure-",case,"-dimplot.png"),res=600,width=2750,height=3000)
DimPlot(case_data[[case]], reduction = "tsne", label=T)
dev.off()

# Generate tsnepalette
tsnepalette = rep("#FFFFFF",11)
names(tsnepalette) = 1:11

# Clusters 5, 2, 4, 10 and 9 represent T-cells (green)
tsnepalette[c(5,2,4,10,9)+1] = rep(RColorBrewer::brewer.pal(10,"Paired")[4],5)

# Clusters 8, 0, 3 and 1 represent melanoma cells (red)
tsnepalette[c(8,0,3,1)+1] = rep(RColorBrewer::brewer.pal(10,"Paired")[6],4)

# Clusters 6 and 7 represent monocytes/macrophages (purple)
tsnepalette[c(6,7)+1] = rep(RColorBrewer::brewer.pal(10,"Paired")[10],2)

# Generate TSNE plot (1.5x)
png(paste0("scRNA/figs/Figure-4C-",case,"-TSNE.png"),res=600,width=4125,height=4500)
{
  tsne = case_data[[case]]@reductions$tsne@cell.embeddings
  plot(tsne, type="n", axes=F, xlab="", ylab="")
  cluster = as.numeric(case_data[[case]]@meta.data$seurat_clusters)
  points(tsne, pch=16, cex=0.6, col=tsnepalette[cluster])
}
dev.off()


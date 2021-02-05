### R script to visualize CYSLTR2 expression in TCGA cohort. ###

'%!in%' <- function(x,y)!('%in%'(x,y))

# Load gene_expression
RSEM_h = drop(read.csv("fig-4a/RSEM_genes_normalized.txt", sep="\t", stringsAsFactors = F, check.names = F, header = F, nrows = 1))
RSEM = read.csv("fig-4a/RSEM_genes_normalized.txt", sep="\t", stringsAsFactors = F, check.names = F, skip = 2, header = F)
RSEM_genes = RSEM[,1]

# Link unique gene names
RSEM_genes_names = sapply(strsplit(RSEM_genes, "\\|"), "[", 1)
RSEM_genes_names[which(duplicated(RSEM_genes_names))] = RSEM_genes[which(duplicated(RSEM_genes_names))]

# Initialize RSEM data matrix
RSEM = as.matrix(RSEM[,2:ncol(RSEM)])
rownames(RSEM) = RSEM_genes_names
colnames(RSEM) = substr(as.character(RSEM_h[1,2:81]),1,12)
head(RSEM)

# Separate CYSLTR2 mutant and wild-type cases and order by CYSLTR2 expression [log2(RSEM+1)]
mutant = c("TCGA-YZ-A982","TCGA-V4-A9ED","TCGA-VD-AA8O")
o = names(sort(log2(RSEM[c("CYSLTR2"),which(colnames(RSEM) %!in% mutant)]+1)))

# Collect immune signature data
data = log2(RSEM[c("CYSLTR2","CD14","CD163","CD3D","CD8A"),o]+1)

# Initialize palette
palette = rev(colorRampPalette(RColorBrewer::brewer.pal(9,"RdYlBu"))(101))

# Create heatmap (figure 4A)
png("fig-4a/fig-4a.png", res=600, height=2500, width=6000)
{
  plot(c(-10,100),
       c(0,8),
       type="n",ylab="",xlab="",
       axes=F)
  
  y = 6
  for (i in 1:nrow(data)) {
    d = round((data[i,]-min(data[i,]))/(max(data[i,])-min(data[i,]))*100)
    for (x in 1:ncol(data)) {
      rect(x, y, x+1, y+1, col=palette[d[x]+1], border=NA)
    }
    y = y-1
    if (y == 5) { y = y-0.5}
  }
  rect(1,6,78,7, lwd=1.4, border="#b1b1b1")
  rect(1,5.5,78,1.5, lwd=1.4, border="#b1b1b1")
  
  for (i in 1:101) {
    rect(80,7-5.5/101*i,82,7-5.5/101*(i-1), lwd=1.5,col = (rev(palette)[i]),border=NA)
  }
}
dev.off()

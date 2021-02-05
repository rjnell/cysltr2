# Function to obtain relevant data from VCF file
#install.packages("vcfR")
read_vcf = function(file, case, region) {
  vcf_file = vcfR::read.vcfR(file)  
  cols = c("GT","DP","AD","RO","QR","AO","QA","GL","POS_ID")
  vcf_table = cbind(rep(case, nrow(vcf_file@fix)), rep(region, nrow(vcf_file@fix)), vcf_file@fix, matrix(data=NA, nrow=nrow(vcf_file@fix), ncol=length(cols)))
  colnames(vcf_table) = c("sample", "region_target", colnames(vcf_file@fix), cols)
  for (index in row(vcf_file@fix)[,1]) {
    formats = strsplit(vcf_file@gt[index, 1], split=":")[[1]]  
    values = strsplit(vcf_file@gt[index, 2], split=":")[[1]]  
    vcf_table[index,formats] = values
    vcf_table[index,"POS_ID"] = paste(vcf_table[index,c("CHROM","POS","REF","ALT")],collapse ="-")
  }  
  return(vcf_table)
}

# Open SNP array copy number data
snp_array_data = read.csv("fig-s3/SNP_array.txt", sep="\t", check.names = F, stringsAsFactors = F)

# Only analyse CYSLTR2 mutant cases
cases = c("TCGA-YZ-A982","TCGA-V4-A9ED","TCGA-VD-AA8O")

# Iterate through cases
for (case in cases) {

  # Initialize sample-specific SNP array data
  sample_snp_array_data = snp_array_data[which(substr(snp_array_data$Sample,1,15) == paste0(case,"-01")),]
  
  # Initialize data savers
  data = list()
  total = 0
  
  # Iterate through chromosomes
  chrs = paste0("chr",c(1:22,"X"))
  
  for (chr in chrs) {
    
    # Read normal vcf and select heterozygous SNPs
    vcf_normal = read_vcf(paste0("fig-s3/",case,"/",case,"-",chr,"-normal.vcf"), paste0(case, "-", "normal"), "chr1")
    AO = as.numeric(vcf_normal[,"AO"])
    RO = as.numeric(vcf_normal[,"RO"])
    w = which(vcf_normal[,"GT"] == "0/1" & abs(.5-AO/(AO+RO)) < .2)
    vcf_normal_selection = vcf_normal[w,]
    
    # Read tumor vcf and analyse selected SNPS and save data
    vcf_tumor = read_vcf(paste0("fig-s3/",case,"/",case,"-",chr,"-tumor.vcf"), paste0(case, "-", "normal"), "chr1")
    results = NULL
    for (vcf_i in 1:nrow(vcf_tumor)) {
      vcf_item = vcf_tumor[vcf_i,]
      w = which(vcf_normal_selection[,"POS_ID"] == vcf_item["POS_ID"])
      if (length(w) == 1) {
        res = c(vcf_item[c("CHROM","POS","REF","ALT","POS_ID","RO","AO")])
        results = rbind(results, res)
      }
    }
    MAF = as.numeric(results[,"AO"]) / (as.numeric(results[,"RO"]) + as.numeric(results[,"AO"]))
    results = cbind(results, MAF)
    w = which((as.numeric(results[,"RO"]) + as.numeric(results[,"AO"])) > 50)
    res = results[w,]
    data[[chr]] = res
    total = total+nrow(res)
    
  }
  
  # Generate plot for file
  file = paste0("fig-s3/chr-plot-",case,".png")
  png(file, res=600, height=3000,width=6000)
  {
    plot(c(1,total+5),c(-1.5,1.1),type="n",axes=F,ylab='',yaxs="i",xaxs="i", main = case, xlab="")
    abline(h=c(2/3, 1/3, 1/2,1,0)-1.5, lty=2, col="#aaaaaaaa")
    
    axis(side = 2,at=c(0,0.3333,0.5,0.66667,1)-1.5, las=2,labels=rep("",length(yat)),col.ticks = "#b1b1b1",col = "#b1b1b1", tck = -0.023,lwd=1.4)
    axis(side = 2,at=c(0,0.3333,0.5,0.66667,1)-1.5,las=2,labels=paste0(c(0,33,50,67,100),"%"), lwd=0, col.axis="#333333",line=-0.33, cex.axis=0.9)
    text(-(total+5)/12.55,-1,srt=90,col="#333333",labels="B-allele fraction",xpd=T, cex=0.9)
    
    # axis(2, las = 2, cex.axis = 1, at=c(0,0.3333,0.5,0.66667,1)-1.5, labels=paste0(c(0,33,50,67,100),"%"), cex.axis=0.9)
    vals = c(5/3, 4/3, 3/3, 2/3)
    values = c("5/3","4/3","3/3", "2/3")
    segments(total+5,0,total+5,1,xpd=T,lwd=1.4,col="#b1b1b1")
    #axis(4, las = 2, cex.axis = 1, at=log2(vals)/2+0.5, labels=values, cex.axis=0.9)
    axis(side = 4,at=log2(vals)/2+0.5, las=2,labels=rep("",length(vals)),col.ticks = "#b1b1b1",col = "#b1b1b1", tck = -0.023,lwd=1.4)
    axis(side = 4,at=log2(vals)/2+0.5, las=2,labels=values, lwd=0, col.axis="#333333",line=-0.33, cex.axis=0.9)
    abline(h=log2(vals)/2+0.5, lty=2, col="#aaaaaaaa")
    
    vals = c(4/2, 3/2, 2/2, 1/2)
    values = c("4/2","3/2","2/2", "1/2")
    # axis(2, las = 2, cex.axis = 1, at=log2(vals)/2+0.5, labels=values, cex.axis=0.9)
    
    axis(side = 2,at=log2(vals)/2+0.5, las=2,labels=rep("",length(vals)),col.ticks = "#b1b1b1",col = "#b1b1b1", tck = -0.023,lwd=1.4)
    axis(side = 2,at=log2(vals)/2+0.5, las=2,labels=values, lwd=0, col.axis="#333333",line=-0.33, cex.axis=0.9)
    text(-(total+5)/12.55,0.5,srt=90,col="#333333",labels="Copy number ratio",xpd=T, cex=0.9)
    #text(-(total+5)/12.55,-0.25,srt=0,col="#333333",labels="Chr.",xpd=T, cex=0.9,pos=4)
    
    
    abline(h=log2(vals)/2+0.5, lty=2, col="#aaaaaaaa")
    
    x=1
    for (chr in chrs) {
      points(x:(x+nrow(data[[chr]])-1),as.numeric(data[[chr]][,"MAF"])-1.5, pch=16, col="#33333333",xpd=T,cex=0.25)
      
      chr2 = substr(chr,4,nchar(chr))
      if (chr == "chrX") {
        chr2 = 23
      }
      chr_snp_array_data = sample_snp_array_data[which(sample_snp_array_data$Chromosome == chr2),]
      i = 0
      for (pos in as.numeric(data[[chr]][,"POS"])) {
        w = which(chr_snp_array_data$Start < pos & chr_snp_array_data$End > pos)    
        if (length(w) == 1) {
          segments(x+i, chr_snp_array_data$Segment_Mean[w]/2+0.5, x+i+1,lwd=1.6,col="darkred")
        }
        i = i +1
      }
      text(x+nrow(data[[chr]])/2,-.25,labels=substr(chr,4,nchar(chr)),xpd=T,cex=0.9,col="#333333")
      x =x+nrow(data[[chr]])
      segments(x,0,x,1,lwd=1.4,col="#b1b1b1")
      segments(x,-1.5,x,-0.5,lwd=1.4,col="#b1b1b1")
    }
  }
  dev.off()
  system(paste("open",file))
}

# Load vcfR library to read VCF files
library(vcfR)

# Define position to search (CYSLTR2 locus)
pos = "chr13:48649049-48718052"

# Set results directory
res_dir = "fig-3cd/res"

# Open samplesheet with RNA-seq info from TCGA
samplesheet = read.csv("K:/TCGA/RNA-Seq/gdc_sample_sheet.2019-08-03.tsv", sep="\t", check.names=F, stringsAsFactors = F)

# Initiate complete table
table_complete = NULL

for (i in 1:80) {
  
  # Collect RNA-seq bam file location for each TCGA case
  input_bam = paste0("/mnt/k/TCGA/RNA-Seq/", samplesheet$`File ID`[i], "/", samplesheet$`File Name`[i])
  TCGA_ID = samplesheet$`Case ID`[i]
  
  # Create a BAM file with only reads in CYSLTR2 locus
  output_bam = paste0("/mnt/f/cysltr2/", res_dir, "/", TCGA_ID, "-CYSLTR2.bam")
  cmd = paste0('wsl samtools view -b ', input_bam, ' "', pos, '" > ', output_bam)
  system(cmd)
  
  # Find variants in CYSLTR2 BAM file using freebayes
  output_vcf = paste0("/mnt/f/cysltr2/", res_dir, "/", TCGA_ID, "-CYSLTR2.vcf")
  cmd = paste0("wsl freebayes -f /mnt/f/hg38/38.fa ", output_bam, " > ", output_vcf)
  system(cmd)
  
  # Open CYSLTR2 VCF file
  vcf_win = paste0("F:/cysltr2/", res_dir, "/", TCGA_ID, "-CYSLTR2.vcf")
  vcf_file = read.vcfR(vcf_win) 
  
  # Obtain relevant information from VCF file and save as 'vcf_table'
  cols = c("GT","DP","AD","RO","QR","AO","QA","GL")
  vcf_table = cbind(rep(TCGA_ID, nrow(vcf_file@fix)), rep("CYSLTR2", nrow(vcf_file@fix)), vcf_file@fix, matrix(data=NA, nrow=nrow(vcf_file@fix), ncol=length(cols)))
  colnames(vcf_table) = c("sample", "gene_target", colnames(vcf_file@fix), cols)
  for (index in row(vcf_file@fix)[,1]) {
    formats = strsplit(vcf_file@gt[index, 1], split=":")[[1]]  
    values = strsplit(vcf_file@gt[index, 2], split=":")[[1]]  
    vcf_table[index,formats] = values
  }
  
  # Add 'vcf_table' from this sample to the complete table
  table_complete = rbind(table_complete, vcf_table) 
}

snps = read.csv("fig-3/snps.txt", stringsAsFactors = F, check.names = F, sep="\t")

select = NULL
for (i in 1:nrow(table_complete)) {
  if (table_complete[i,"POS"] %in% snps$pos & (as.numeric(table_complete[i,"RO"]) > 10 | as.numeric(table_complete[i,"AO"]) > 10)) {
    table_complete[i,"ID"] = snps$id[which(snps$pos == table_complete[i,"POS"])]
    select = c(select,i)  
  }
}

View(table_complete[select,])

# Save complete table
write.table(x = table_complete[select,], file="fig-3/cysltr2.tsv", sep="\t")


#Identification of differentially expressed genes
#Requires a input text file containing raw read counts


#Md Shamimuzzaman
# Postdoctoral Research Computational Biologist at USDA-ARS
# Date: 10 Nov, 2016

source("https://bioconductor.org/biocLite.R")
biocLite("edgeR")
library("edgeR")
#Reading read counts data from a txt file
y = read.table("Input.txt", header=TRUE, row.names=1)
#Grouping read count columns into control(ctrl) and treatment as there are 3 replicates  
group=c("ctrl", "ctrl", "ctrl", "treatment","treatment","treatment")
dge = DGEList(counts = y, group = group)
#Estimate common dispersion
dge = estimateCommonDisp(dge)
# Estimate Tagwise dispersion
dge = estimateTagwiseDisp(dge)
#Perform pair-wise tests for differential expression between two groups.
et = exactTest(dge)
#topTags returns the top differentially expressed genes
resEdgeR = topTags(et)
#Wanted to see the results for all of my genes (1100)
resEdgeR = topTags(et, n=1100)
resEdgeR=data.frame(resEdgeR)
# Subsetting to get DE genes with corrected/adjusted p-value <= 0.05, here it is referred as FDR.
resEdgeR=subset(resEdgeR, resEdgeR$FDR<=0.05)
head(resEdgeR)
write.csv(resEdgeR,"DE_genes_edgeR.csv")



Without replicate
y = read.table("Input.txt", header=TRUE, row.names=1)
group=c("ctrl", "treatment")
dge = DGEList(counts = y, group = group)
dge = estimateGLMCommonDisp(dge,method="deviance",robust=TRUE,subset = NULL)
dge = estimateGLMTagwiseDisp(dge)
et = exactTest(dge)
resEdgeR = topTags(et, n=1100)
write.csv(resEdgeR,"DE_genes_edgeR.csv")

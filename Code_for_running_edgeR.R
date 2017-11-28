#Identification of differentially expressed genes
#Requires a input text file containing raw read counts


#Md Shamimuzzaman
# Postdoctoral Research Computational Biologist at USDA-ARS
# Date: 10 Nov, 2016

source("https://bioconductor.org/biocLite.R")
biocLite("edgeR")
library("edgeR")
y = read.table("dexample.txt", header=TRUE, row.names=1)  
group=c("ctrl", "ctrl", "ctrl", "treatment","treatment","treatment")
dge = DGEList(counts = y, group = group)
dge = estimateCommonDisp(dge)
dge = estimateTagwiseDisp(dge)
et = exactTest(dge)

## Extract results from edgeR analysis
resEdgeR = topTags(et)
resEdgeR
resEdgeR = topTags(et, n=201)
write.csv(resEdgeR,"Avg_24H_edgeR_output.csv")

With Replicate
y = read.csv("72H_Reps.csv", header=TRUE, row.names=1)  
group=c("ctrl", "ctrl", "ctrl", "treatment","treatment","treatment")
dge = DGEList(counts = y, group = group)
dge = estimateCommonDisp(dge)
dge = estimateTagwiseDisp(dge)
et = exactTest(dge)
resEdgeR = topTags(et, n=201)
write.csv(resEdgeR,"Reps_72H_edgeR_output.csv")


Without replicate
y = read.csv("24H.csv", header=TRUE, row.names=1) 
group=c("ctrl", "treatment")
dge = DGEList(counts = y, group = group)
dge = estimateGLMCommonDisp(dge,method="deviance",robust=TRUE,subset = NULL)
dge = estimateGLMTagwiseDisp(dge)
et = exactTest(dge)
resEdgeR = topTags(et, n=201)
write.csv(resEdgeR,"24H_edgeR_output.csv")

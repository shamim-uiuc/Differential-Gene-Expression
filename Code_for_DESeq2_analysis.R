DESeq2 for Differential expression:
source("http://bioconductor.org/biocLite.R")
biocLite("DESeq2")
library("DESeq2")

#Method1:
raw.data = read.table(file= "rexample.txt",header = TRUE)
counts = raw.data[, -1]
rownames(counts)=raw.data[,1]
colData <- DataFrame(condition=factor(c("ctrl", "ctrl", "ctrl", "treatment","treatment","treatment")))
dds <- DESeqDataSetFromMatrix(counts, colData, formula(~ condition))
dds = estimateSizeFactors(dds)
dds = estimateDispersions(dds)
dds = nbinomWaldTest(dds)
resDESeq2 = results(dds)
head(resDESeq2)
write.csv(resDESeq2, file = "all_genes_expression.csv")


#Mthod-2 : Slightly different than method 1

CountData = read.table(file= "rexample.txt",header = TRUE,row.names=1)
# Filter data where you only have 0 or 1 read count across all samples.
CountData = CountData[rowSums(CountData)>1, ]
head(CountData)

# Import metadata
colData = read.table("meta.txt", row.names=1)
colData
##               condition
## SRR493366 control_sirna
## SRR493367 control_sirna
## SRR493368 control_sirna
## SRR493369      hoxa1_kd
## SRR493370      hoxa1_kd
## SRR493371      hoxa1_kd
# Set up the DESeqDataSet Object and run the DESeq pipeline
dds = DESeqDataSetFromMatrix(countData=CountData,
                              colData=colData,
                              design=~condition)
dds = DESeq(dds)
res = results(dds, contrast=c("condition", "ctrl", "treatment"))
write.csv(res, file = "all_genes_expression.csv")



#Method-3: Pairwise comparision included
> Library("DESeq2")                
 > CountTable  = read.table("dexample.txt", header=TRUE, row.names=1)                
 > head(CountTable)            
              
 gene    tissue1    tissue2    tissue3    tissue4
 gene1    233    91    17    593
 gene2    1011    0    7    1
 gene3    963    2    3    66
 gene4    908    41    1    74
 gene5    596    50    26    328
 gene6    1    0    0    0
                 
colData = data.frame(row.names= colnames(CountTable),condition = c("ctrl", "ctrl", "ctrl", "treatment","treatment","treatment"))
dds <- DESeqDataSetFromMatrix( countData = CountTable, colData = colData,design = ~ condition)                
dds <- DESeq(dds)   
res <- results(dds)   



For pairwise comparision, you have to define

res.A.B <- results(dds1, contrast=c("condition","A","B"))
res.A.C <- results(dds1, contrast=c("condition","A","C"))
res.B.C <- results(dds1, contrast=c("condition","B","C"))




#####Method for DESeq2 without replicates:
library(DESeq2)
raw.data = read.csv(file= "Avg_24H.csv",header = TRUE)
counts = raw.data[, -1]
rownames(counts)=raw.data[,1]
colData <- DataFrame(condition=factor(c("ctrl","treatment")))
dds <- DESeqDataSetFromMatrix(counts, colData, formula(~ condition))
colData(dds)$condition<-factor(colData(dds)$condition,levels=c("ctrl","treatment"))
dds <- DESeq(dds)
res<-results(dds)
write.csv(res, file = "all_genes_expression.csv")



##### DESeq code for data without replicte
library(DESeq)
y = read.csv("24H_Reps.csv", header=TRUE, row.names=1) 
conds=c("ctrl", "treatment")
cds <- newCountDataSet( y, conds )
cds <- estimateSizeFactors( cds )
cds<-estimateDispersions(cds,method="blind",sharingMode="fit-only",fitType="local") 
res<-nbinomTest(cds,"ctrl","treatment")
write.csv(res,"DESeq_24H_output.csv")


##### DESeq code for data with replicate
library(DESeq)
y = read.csv("72H_Reps.csv", header=TRUE, row.names=1) 
conds=c("ctrl", "ctrl","ctrl","treatment","treatment","treatment")
cds <- newCountDataSet( y, conds )
cds <- estimateSizeFactors( cds )
cds<-estimateDispersions(cds) 
res<-nbinomTest(cds,"ctrl","treatment")
write.csv(res,"Reps_DESeq_72H_output.csv")


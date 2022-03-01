suppressPackageStartupMessages({
   library(tximport)
   library(SummarizedExperiment)
})
t2g <- read.delim("/reference/Homo_sapiens/Ensembl/GRCh38.p13/Annotation/Release_105-2021/tx2gene", header=FALSE)
lf <- list.files(pattern="quant\\.sf", recursive=TRUE, full=TRUE)
names(lf) <- basename(dirname(lf))
ti <- tximport(lf, type="salmon", tx2gene=t2g)
se <- SummarizedExperiment(list(counts=ti$counts, TPM=ti$abundance))
rowData(se)$medianLength <- matrixStats::rowMedians(ti$length)
se$sample_id <- sapply(strsplit(colnames(se),"_"), FUN=function(x) x[2])
se$batch <- paste0("set",1+(as.numeric(se$sample_id)>13)+
  (as.numeric(se$sample_id)>26) +
  (as.numeric(se$sample_id)>39))

conds <- lapply(list(
"18h untreated"=c(1,14,27,40),
"18h DMSO"=c(2,15,28,41),
"18h KH-103"=c(3,16,29,42),
"18h MIF"=c(4,17,30,43),
"18h Cort113"=c(5,18,31,44),
"16h DMSO > 2h DEX+DMSO"=c(6,19,32,45),
"16h KH-103 > 2h DEX+KH-103"=c(7,20,33,46),
"16h MIF  > 2h DEX+MIF"=c(8,21,34,47),
"16h Cort113 > 2h DEX+Cort113"=c(9,22,35,48),
"2h DEX > 16h DEX+DMSO"=c(10,23,36,49),
"2h DEX > 16h DEX+KH-103"=c(11,24,37,50),
"2h DEX > 16h DEX+MIF"=c(12,25,38,51),
"2h DEX > 16h DEX+Cort113"=c(13,26,39,52)
), as.character)

s2cond <- setNames(rep(names(conds),lengths(conds)),unlist(conds))
se$condition <- s2cond[se$sample_id]

saveRDS(se, file="geneLevel.SE.rds")


#!/usr/bin/Rscript
args = commandArgs(trailingOnly=TRUE)


basename = args[1]
count_file = args[2]
sample_info = args[3]
diff_type = args[4]
reps = args[5]



library(DESeq2)

merge <- read.delim(count_file,sep=",")
#Drop the tenth column
#merge <- merge[,-14]
head(merge,1)
data_file <- read.delim(sample_info,sep=",")
coldata <- as.data.frame(apply(data_file,2,as.factor))

if (diff_type == "TE") {
	#TE analysis
	mat <- merge[,-1]
	rownames(mat) <- merge[,1]
	ddsMat <- DESeqDataSetFromMatrix(countData = mat,colData = coldata, design = ~ Condition + SeqType + Condition:SeqType)
	ddsMat$SeqType = relevel(ddsMat$SeqType,"rnaseq")
	ddsMat <- DESeq(ddsMat)
	print(resultsNames(ddsMat))
	res <- results(ddsMat, contrast=list("Conditiontreatment.SeqTyperiboseq"))
	write.table(res, paste(basename,"DESeq2_TE.txt",sep="_"),sep=",",quote=FALSE,col.names=NA)

	if (reps == "2") {
		#RNA-Seq only analysis
		mat <- merge[c(4,5,8,9)]
		rnacoldata = coldata[c(3,4,7,8),]
		rownames(mat) <- merge[,1]
		ddsMat <- DESeqDataSetFromMatrix(countData = mat,colData = rnacoldata, design = ~ Condition)
		ddsMat <- DESeq(ddsMat)
		print(resultsNames(ddsMat))
		res <- results(ddsMat, contrast=list("Condition_treatment_vs_control"))
		write.table(res, paste(basename,"DESeq2_RNASEQ.txt",sep="_"),sep=",",quote=FALSE,col.names=NA)
		
		#Riboseq only analysis
		mat <- merge[c(2,3,6,7)]
		ribocoldata = coldata[c(1,2,5,6),]
		rownames(mat) <- merge[,1]
		ddsMat <- DESeqDataSetFromMatrix(countData = mat,colData = ribocoldata, design = ~ Condition)
		ddsMat <- DESeq(ddsMat)
		print(resultsNames(ddsMat))
		res <- results(ddsMat, contrast=list("Condition_treatment_vs_control"))
		write.table(res, paste(basename,"DESeq2_RIBOSEQ.txt",sep="_"),sep=",",quote=FALSE,col.names=NA)
		}

	if (reps == "3") {
		#RNA-Seq only analysis
		mat <- merge[c(4,5,8,9,12,13)]
		rnacoldata = coldata[c(3,4,7,8,11,12),]
		rownames(mat) <- merge[,1]
		ddsMat <- DESeqDataSetFromMatrix(countData = mat,colData = rnacoldata, design = ~ Condition)
		ddsMat <- DESeq(ddsMat)
		print(resultsNames(ddsMat))
		res <- results(ddsMat, contrast=list("Condition_treatment_vs_control"))
		write.table(res, paste(basename,"DESeq2_RNASEQ.txt",sep="_"),sep=",",quote=FALSE,col.names=NA)

		#Riboseq only analysis
		mat <- merge[c(2,3,6,7,10,11)]
		ribocoldata = coldata[c(1,2,5,6,9,10),]
		rownames(mat) <- merge[,1]
		ddsMat <- DESeqDataSetFromMatrix(countData = mat,colData = ribocoldata, design = ~ Condition)
		ddsMat <- DESeq(ddsMat)
		print(resultsNames(ddsMat))
		res <- results(ddsMat, contrast=list("Condition_treatment_vs_control"))
		write.table(res, paste(basename,"DESeq2_RIBOSEQ.txt",sep="_"),sep=",",quote=FALSE,col.names=NA)
		}
	
	if (reps == "4") {
		#RNA-Seq only analysis
		mat <- merge[c(4,5,8,9,12,13,16,17)]
		rnacoldata = coldata[c(3,4,7,8,11,12,15,16),]
		rownames(mat) <- merge[,1]
		ddsMat <- DESeqDataSetFromMatrix(countData = mat,colData = rnacoldata, design = ~ Condition)
		ddsMat <- DESeq(ddsMat)
		print(resultsNames(ddsMat))
		res <- results(ddsMat, contrast=list("Condition_treatment_vs_control"))
		write.table(res, paste(basename,"DESeq2_RNASEQ.txt",sep="_"),sep=",",quote=FALSE,col.names=NA)
		
		#Riboseq only analysis
		mat <- merge[c(2,3,6,7,10,11,14,15)]
		ribocoldata = coldata[c(1,2,5,6,9,10,13,14),]
		rownames(mat) <- merge[,1]
		ddsMat <- DESeqDataSetFromMatrix(countData = mat,colData = ribocoldata, design = ~ Condition)
		ddsMat <- DESeq(ddsMat)
		print(resultsNames(ddsMat))
		res <- results(ddsMat, contrast=list("Condition_treatment_vs_control"))
		write.table(res, paste(basename,"DESeq2_RIBOSEQ.txt",sep="_"),sep=",",quote=FALSE,col.names=NA)
		}	
	
	
	
	
}




if (diff_type == "Riboseq") {
	#Riboseq only analysis
	mat = merge[,-1]
	rownames(mat) <- merge[,1]
	ddsMat <- DESeqDataSetFromMatrix(countData = mat,colData = coldata, design = ~ Condition)
	ddsMat <- DESeq(ddsMat)
	print(resultsNames(ddsMat))
	res <- results(ddsMat, contrast=list("Condition_treatment_vs_control"))
	write.table(res, paste(basename,"DESeq2_RIBOSEQ.txt",sep="_"),sep=",",quote=FALSE,col.names=NA)
	}


if (diff_type == "Rnaseq") {
	#Rnaseq only analysis
	mat = merge[,-1]
	rownames(mat) <- merge[,1]
	head(mat,10)
	ddsMat <- DESeqDataSetFromMatrix(countData = mat,colData = coldata, design = ~ Condition)
	ddsMat <- DESeq(ddsMat)
	print(resultsNames(ddsMat))
	res <- results(ddsMat, contrast=list("Condition_treatment_vs_control"))
	write.table(res, paste(basename,"DESeq2_RNASEQ.txt",sep="_"),sep=",",quote=FALSE,col.names=NA)
	foo <- counts(ddsMat, normalized = FALSE)
	write.table(foo, paste(basename,"DESeq2_non_norm_counts.txt",sep="_"),sep=",",quote=FALSE,col.names=NA)
	foo <- counts(ddsMat, normalized = TRUE)
	write.table(foo, paste(basename,"DESeq2_norm_counts.txt",sep="_"),sep=",",quote=FALSE,col.names=NA)
	}

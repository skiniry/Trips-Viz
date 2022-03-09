#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)



basename = args[1]
count_file = args[2]
sample_info = args[3]
diff_type = args[4]
reps = args[5]
maxp = as.numeric(args[6])


library(anota2seq)
infile <- read.csv(file=count_file, row.names = 1, header=TRUE, sep=",")
if (diff_type == "TE") {
	anota2seq_data_P <- infile[,c(1,2,5,6,9,10)]
	anota2seq_data_T <- infile[,c(3,4,7,8,11,12)]
	head(anota2seq_data_P)
	head(anota2seq_data_T)
	anota2seq_pheno_vec <- c("control","treatment","control","treatment","control","treatment")
	ads <- anota2seqDataSetFromMatrix(dataP = anota2seq_data_P,dataT = anota2seq_data_T,phenoVec = anota2seq_pheno_vec,dataType = "RNAseq",normalize = TRUE)
	ads <- anota2seqRun(ads)

	translation_res <- anota2seqGetOutput(ads, analysis = "translation",output = "full",selContrast = 1,getRVM = TRUE)
	ribo_fc_res <- anota2seqGetDeltaData(ads, output = "full", analysis="translated mRNA", selContrast = 1)
	ribo_fc_res <- as.data.frame(ribo_fc_res)
	translation_fc_res <-  merge(translation_res,ribo_fc_res["deltaP"],by="row.names",all.x=TRUE)
	write.table(translation_fc_res, paste(basename,"anota2seq_RIBOSEQ.csv",sep="_"),sep=",",quote=FALSE,col.names=NA)

	mRNA_abundance_res <- anota2seqGetOutput(ads, analysis = "total mRNA", output = "full",selContrast = 1,getRVM = TRUE)
	mrna_fc_res <- anota2seqGetDeltaData(ads, output = "full", analysis="total mRNA", selContrast = 1)
	mrna_fc_res <- as.data.frame(mrna_fc_res)
	mrna_res <- merge(mRNA_abundance_res,mrna_fc_res["deltaT"],by="row.names",all.x=TRUE)
	head(mRNA_abundance_res)
	head(mrna_fc_res)
	head(mrna_res)
	write.table(mrna_res, paste(basename,"anota2seq_mRNA.csv",sep="_"),sep=",",quote=FALSE,col.names=NA)


	buffered_res <- anota2seqGetOutput(ads, analysis = "buffering",output = "full",selContrast = 1,getRVM = TRUE)
	write.table(buffered_res, paste(basename,"anota2seq_BUFFERED.txt",sep="_"),sep=",",quote=FALSE,col.names=NA)

	print ("Before sig genes")
	ads <- anota2seqSelSigGenes(ads, useRVM = TRUE,selContrast = 1,maxPAdj = maxp)
	ads <- anota2seqRegModes(ads)
	sig_translated_df <- anota2seqGetOutput(ads, analysis = "translation", output = "selected", getRVM = TRUE,selContrast = 1)
	sig_translated_genes <- row.names(sig_translated_df)
	write.table(sig_translated_genes, paste(basename,"anota2seq_sig_translated_genes.csv",sep="_"),sep=",",quote=FALSE,col.names=NA)
	print ("translated genes written")
	sig_buffered_df <- anota2seqGetOutput(ads, analysis = "buffering", output = "selected", getRVM = TRUE,selContrast = 1)
	sig_buffered_genes <- row.names(sig_buffered_df)
	write.table(sig_buffered_genes, paste(basename,"anota2seq_sig_buffered_genes.csv",sep="_"),sep=",",quote=FALSE,col.names=NA)
	print ("buffered genes written")
	sig_rna_df <- anota2seqGetOutput(ads, analysis = "mRNA abundance", output = "selected", getRVM = TRUE,selContrast = 1)
	print ("sig rna dif:")
	sig_rna_genes <- row.names(sig_rna_df$"abundance translated mRNA")
	write.table(sig_rna_genes, paste(basename,"anota2seq_sig_rna_genes.csv",sep="_"),sep=",",quote=FALSE,col.names=NA)
	print ("abundance genes written")

	#ribo_fc_res <- anota2seqGetDeltaData(ads, output = "full", analysis="translated mRNA", selContrast = 1)
	#ribo_fc_res <- as.data.frame(ribo_fc_res)
	#write.table(ribo_fc_res, paste(basename,"anota2seq_ribo_FC.csv",sep="_"),sep=",",quote=FALSE,col.names=NA)
	#fc_res <- merge(mrna_fc_res,ribo_fc_res["deltaP"],by="row.names",all.x=TRUE)
	#write.table(fc_res, paste(basename,"anota2seq_FC.csv",sep="_"),sep=",",quote=FALSE,col.names=NA)
	}

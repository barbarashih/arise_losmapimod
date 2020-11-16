#### STEP 0: import libraries ####
if(!require("DESeq2",character.only = TRUE)){
        source("https://bioconductor.org/biocLite.R")
        biocLite("DESeq2")
}
if(!require("edgeR",character.only = TRUE)){
        source("https://bioconductor.org/biocLite.R")
        biocLite("edgeR")
}
if(!require("biomaRt",character.only = TRUE)){
        source("https://bioconductor.org/biocLite.R")
        biocLite("biomaRt")
}
if(!require("sva",character.only = TRUE)){
        source("https://bioconductor.org/biocLite.R")
        biocLite("sva")
}
#library(Rsubread)
library(edgeR)
require(plyr)
library(biomaRt)
library(sva)
library(devtools)
#library(rtracklayer)
#library(ballgown)
#library(limma)
#library("DESeq2")
library(ggplot2)

#### STEP 1: setup ####
working_dir=getwd()
phenotypic_data_fp ="sample_annotation/annotation_ARISE_all.txt"
gene_count_fp="count_matrix/gene_count_matrix.csv"
sample_gtf_dir="gtf_round2"
out_exprs_dir="exprs"
out_sig_dir="stats_sig_2019"
sig_gene_out_dir <- "stats_sig_2019/exprs_annotation"
imsig_fp <-"gene_annotation/ImSig.csv"

# Set up working enviroment
setwd(working_dir)
source("script/fun_stats_limma.r")
source("script/r_goterm_genes.r")
dir.create(file.path(working_dir, out_exprs_dir), showWarnings = FALSE) # create output directory for expression files
dir.create(file.path(working_dir, out_sig_dir), showWarnings = FALSE) # create output directory for expression files

# Import count data
countData <- as.matrix(read.csv(gene_count_fp, row.names="gene_id"))
colData <- read.table(phenotypic_data_fp, sep="\t", header=TRUE, row.names="sample.ID", stringsAsFactors=FALSE)
colData$sampleid <- row.names(colData)
countData <- countData[, rownames(colData)]

# Get gene annotation data
geneIDs = data.frame(geneIDs=rownames(countData), stringsAsFactors=FALSE)
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
G_list <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","external_gene_name", "description", "transcript_length"),values=geneIDs$geneIDs,mart= mart)
G_list <- G_list[order(G_list$ensembl_gene_id,G_list$transcript_length,decreasing = TRUE),]
G_list <- G_list[!duplicated(G_list$ensembl_gene_id),]
G_list$description <- sapply(strsplit(G_list$description," \\["), FUN=function(x){return(x[1])})
gene_annotation <- merge(geneIDs, G_list, by.x="geneIDs", by.y="ensembl_gene_id", all.x=TRUE)
colnames(gene_annotation) <- c("geneIDs", "geneNames", "description", "transcript_length")

small_RNA <- gene_annotation$geneIDs[gene_annotation$transcript_length < 200] # remove small RNA 
gene_annotation$transcript_length = NULL
# organise dge
dge <- DGEList(counts=countData, genes=rownames(countData))
countsPerMillion <- cpm(dge)
countCheck <- countsPerMillion > 1
passed_ids <- rowSums(countCheck) >= ncol(countData)*0.2
passed_ids[names(passed_ids) %in% small_RNA] <- FALSE
#passed_ids[!(names(passed_ids) %in% gene_annotation$geneIDs)] <- FALSE
keep <- which(passed_ids)
dge <- dge[keep,]
dge <- calcNormFactors(dge, method="TMM")
nrow(dge)


# Imsig geneset enrichement..etc
imsig <- read.csv(imsig_fp, stringsAsFactors=FALSE)
gene_lists <- split(imsig$Gene.name, imsig$ImSig)

#####################
#### compare pre-post losmapimod, treat vs n
# S6 - pre-los all
current_output_dir <- "contrast_fit_count_mx_treat"
dir.create(file.path(working_dir, out_sig_dir, current_output_dir), showWarnings = FALSE) # create output directory for expression files
current_out_dir=file.path(out_sig_dir,"contrast_fit_count_mx_treat")
current_comparison="subsetted_LosS6_all"
sample_annotation <- colData[colData$treatment %in% c("S6", "N"),]
sample_annotation$treatment <- factor(sample_annotation$treatment, levels=c("N","S6"))
levels(sample_annotation$treatment) <- c("N", "T")
design <- fun.makeDesign(current_samp_annotation= sample_annotation, design_formula= "~0+drug_status:treatment")
contrast.matrix<- fun.contrastTreat(design)
out_fp=file.path(working_dir, out_exprs_dir, paste(current_comparison, ".csv", sep=""))
top_genes_LosS6_all <- fun.topGeneAnalysis(dge_full = dge, gene_annotation=gene_annotation, current_samp_annotation=sample_annotation, current_design =design, current_contrast= contrast.matrix ,current_block=sample_annotation$subject,  return_exprs=TRUE)

###### START OF MROAST #######
###### mroast simplified comparison ####
#### mroast S6 ####
# pre losmapimod
current_comparison="Los_S6vsN_allRes_pre"
sample_annotation <- colData[colData$treatment %in% c("S6", "N"),]
sample_annotation <- sample_annotation[sample_annotation$drug_status == "pre",]
design <- fun.makeDesign(current_samp_annotation= sample_annotation, design_formula= "~treatment")
Los_S6vsN_allRes_pre <- fun.topGeneAnalysis(dge_full = dge, gene_annotation=gene_annotation, current_samp_annotation=sample_annotation, current_design =design,current_block=sample_annotation$subject,  return_exprs=TRUE)
Los_S6vsN_allRes_pre_mroast <- fun.mroast(current_object=Los_S6vsN_allRes_pre, gene_symbol_list=GO_genes_S6)
# post losmapimod
fc_filter <- 1.2
current_comparison="Los_S6vsN_allRes_post"
sample_annotation <- colData[colData$treatment %in% c("S6", "N"),]
sample_annotation <- sample_annotation[sample_annotation$drug_status == "post",]
design <- fun.makeDesign(current_samp_annotation= sample_annotation, design_formula= "~treatment")
Los_S6vsN_allRes_post <- fun.topGeneAnalysis(dge_full = dge, gene_annotation=gene_annotation, current_samp_annotation=sample_annotation, current_design =design,current_block=sample_annotation$subject,  return_exprs=TRUE)
Los_S6vsN_allRes_post_mroast <- fun.mroast(current_object=Los_S6vsN_allRes_post, gene_symbol_list=GO_genes_S6)
sig_path_pre <- Los_S6vsN_allRes_pre_mroast$Row.names[Los_S6vsN_allRes_pre_mroast$FDR<0.05 & Los_S6vsN_allRes_pre_mroast$average_logfc > log2(fc_filter) ]
#sig_path_pre <- Los_S6vsN_allRes_pre_mroast$Row.names[Los_S6vsN_allRes_pre_mroast$FDR<0.05 ]
Los_S6vsN_allRes_pre_mroast2 <- Los_S6vsN_allRes_pre_mroast[,c("Row.names", "FDR", "average_logfc")]
Los_S6vsN_allRes_post_mroast2 <- Los_S6vsN_allRes_post_mroast[,c("Row.names", "FDR", "average_logfc")]
Los_S6vsN_allRes_pre_mroast2$losmapimod <- "pre"
Los_S6vsN_allRes_post_mroast2$losmapimod <- "post"
S6_mroast_compare <- rbind(Los_S6vsN_allRes_pre_mroast2, Los_S6vsN_allRes_post_mroast2)
colnames(S6_mroast_compare)[1] <- "pathway"
S6_mroast_compare <- S6_mroast_compare[S6_mroast_compare$pathway %in% sig_path_pre, ]
S6_mroast_compare <- S6_mroast_compare[S6_mroast_compare$pathway %in% sig_path_pre, ]
S6_mroast_compare$FDR_binned <- "<0.001"
S6_mroast_compare$FDR_binned <- ifelse(S6_mroast_compare$FDR>= 0.001, "<0.05", S6_mroast_compare$FDR_binned) 
S6_mroast_compare$FDR_binned <- ifelse(S6_mroast_compare$FDR>= 0.01, "<0.01", S6_mroast_compare$FDR_binned) 
S6_mroast_compare$FDR_binned <- ifelse(S6_mroast_compare$FDR>= 0.05, ">=0.05", S6_mroast_compare$FDR_binned) 
S6_mroast_compare$FDR_binned <- factor(S6_mroast_compare$FDR_binned, levels= c("<0.001", "<0.01", "<0.05", ">=0.05"))
S6_mroast_compare$pathway <- sapply(strsplit(S6_mroast_compare$pathway, "_"), FUN=function(x)paste(x, collapse=" "))
S6_mroast_compare$pathway <-  sapply(S6_mroast_compare$pathway, FUN=function(x)paste(strwrap(x, width=40), collapse="\n"))
S6_mroast_compare$pathway <- factor(S6_mroast_compare$pathway)
S6_mroast_compare$losmapimod <- factor(S6_mroast_compare$losmapimod, levels = c("pre", "post"))
S6_mroast_compare$FC <- 2^S6_mroast_compare$average_logfc 
out_plot <- ggplot(data=S6_mroast_compare, aes(x=losmapimod, y=pathway)) + 
	geom_point(aes(size= average_logfc, colour=FDR_binned)) +
	scale_colour_manual(values = c("#cb181d", "#fb6a4a", "#fcae91", "#fee5d9")) +
	theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) +
	theme(plot.title = element_text(size=18, face="bold"), legend.key = element_rect(colour = "transparent", fill = "white")) +
	theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave(paste0("plots/gene_set_enrichment/S6_bothResp_logFC_prefilt",fc_filter, ".pdf"), plot = out_plot, scale = 1, width = 4.5, height = 6, units = "in", useDingbats=FALSE)
out_plot <- ggplot(data=S6_mroast_compare, aes(x=losmapimod, y=pathway)) + 
	geom_point(aes(size= FC, colour=FDR_binned)) +
	scale_colour_manual(values = c("#cb181d", "#fb6a4a", "#fcae91", "#fee5d9")) +
	theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) +
	theme(plot.title = element_text(size=18, face="bold"), legend.key = element_rect(colour = "transparent", fill = "white")) +
	theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave(paste0("plots/gene_set_enrichment/S6_bothResp_FC_prefilt",fc_filter, ".pdf"), plot = out_plot, scale = 1, width = 4.5, height = 6, units = "in", useDingbats=FALSE)
ggsave(paste0("plots/gene_set_enrichment/S6_bothResp_FC_prefilt",fc_filter, ".png"), plot = out_plot, scale = 1, width = 4.5, height = 6, units = "in")



##### STATISTICAL COMPARISONS #####
#####################
#### compare pre-post losmapimod, treat vs n
# S6 - pre-los all
current_output_dir <- "contrast_fit_count_mx_treat"
dir.create(file.path(working_dir, out_sig_dir, current_output_dir), showWarnings = FALSE) # create output directory for expression files
current_out_dir=file.path(out_sig_dir,"contrast_fit_count_mx_treat")
current_comparison="subsetted_LosS6_all"
sample_annotation <- colData[colData$treatment %in% c("S6", "N"),]
sample_annotation$treatment <- factor(sample_annotation$treatment, levels=c("N","S6"))
levels(sample_annotation$treatment) <- c("N", "T")
design <- fun.makeDesign(current_samp_annotation= sample_annotation, design_formula= "~0+drug_status:treatment")
contrast.matrix<- fun.contrastTreat(design)
out_fp=file.path(working_dir, out_exprs_dir, paste(current_comparison, ".csv", sep=""))
top_genes_LosS6_all <- fun.topGeneAnalysis(dge_full = dge, gene_annotation=gene_annotation, current_samp_annotation=sample_annotation, current_design =design, current_contrast= contrast.matrix ,current_block=sample_annotation$subject,  return_exprs=TRUE)

######## Treating improver and minimal improvers in a group and carry out the stats separately #####
# S6 improver/non-improver
current_out_dir=file.path(out_sig_dir,"contrast_fit_count_mx_treat")
current_comparison="subsetted_LosS6"
sample_annotation <- colData[colData$treatment %in% c("S6", "N"),]
sample_annotation$treatment <- factor(sample_annotation$treatment, levels=c("N","S6"))
levels(sample_annotation$treatment) <- c("N", "T")
design <- fun.makeDesign(current_samp_annotation= sample_annotation, design_formula= "~0+drug_status:response:treatment")
contrast.matrix<- fun.contrastResp_Treat(design)
out_fp=file.path(working_dir, out_exprs_dir, paste(current_comparison, ".csv", sep=""))
top_genes_LosS6 <- fun.topGeneAnalysis(dge_full = dge, gene_annotation=gene_annotation, current_samp_annotation=sample_annotation, current_design =design, current_contrast= contrast.matrix ,current_block=sample_annotation$subject, return_exprs=TRUE)
top_genes_temp <- fun.printTopTableRows(top_genes_LosS6[["top_genes"]],logfc_threshold=2, prefix=current_comparison)
fun.write_TopTables(top_genes_LosS6, out_dir=current_out_dir, output_prefix= current_comparison)



####################
#### Compare pre vs post losmapimod for each treatment for either improver or non-improver ####
current_out_dir=file.path(out_sig_dir,"contrast_fit_count_mx_treat")
current_comparison="S6_los"
sample_annotation <- colData[colData$treatment %in% c("S6", "N"),]
sample_annotation$treatment <- factor(sample_annotation$treatment, levels=c("N","S6"))
levels(sample_annotation$treatment) <- c("N", "T")
design <- fun.makeDesign(current_samp_annotation= sample_annotation, design_formula= "~0+drug_status:response:treatment")
contrast.matrix<- fun.contrastResp_Treat(design)
out_fp=file.path(working_dir, out_exprs_dir, paste(current_comparison, ".csv", sep=""))
top_genes_LosS6 <- fun.topGeneAnalysis(dge_full = dge, gene_annotation=gene_annotation, current_samp_annotation=sample_annotation, current_design =design, current_contrast= contrast.matrix ,current_block=sample_annotation$subject, return_exprs=TRUE)
top_genes_temp <- fun.printTopTableRows(top_genes_LosS6[["top_genes"]],logfc_threshold=2, prefix=current_comparison)
fun.write_TopTables(top_genes_LosS6, out_dir=current_out_dir, output_prefix= current_comparison)


current_comparison="pre_vs_post_response"
sample_annotation <- colData
sample_annotation$response <- ifelse(sample_annotation$response %in% "non-responder", "non-responder", "responder")
design <- fun.makeDesign(current_samp_annotation= sample_annotation, design_formula= "~0+drug_status : treatment : response")
contrast.matrix <- fun.contrastLos_response(design)
out_fp=file.path(working_dir, out_exprs_dir, paste(current_comparison, ".csv", sep=""))
top_genes_LosTreatResponder <- fun.topGeneAnalysis(dge_full = dge, gene_annotation=gene_annotation, current_samp_annotation=sample_annotation, current_design =design, current_contrast= contrast.matrix , current_block=sample_annotation$subject , write_exprs=out_fp)
out_geneset_fp=file.path(working_dir, out_exprs_dir, paste(current_comparison, ".csv", sep=""))
# fun.geneset(dge_full = dge, gene_annotation=gene_annotation, current_samp_annotation=sample_annotation, current_design =design, current_contrast= contrast.matrix, current_block=sample_annotation$subject, out_name=out_geneset_fp, gene_set=inflammIndex_genes)
print(current_comparison)
top_genes_temp <- fun.printTopTableRows(top_genes_LosTreatResponder[["top_genes"]],logfc_threshold=2, prefix=current_comparison)
fun.write_TopTables(stats_result=top_genes_LosTreatResponder, out_dir=out_sig_dir, output_prefix=current_comparison)

genes_losv6_resp_minresp <- top_genes_LosTreatResponder$top_genes$RespV6
genes_losv6_resp_minresp <- genes_losv6_resp_minresp[genes_losv6_resp_minresp$adj.P.Val < 0.05,]
#genes_losv6_resp_minresp <- genes_losv6_resp_minresp[abs(genes_losv6_resp_minresp$logFC) > log2(1.5),]
write.csv(genes_losv6_resp_minresp, file.path(sig_gene_out_dir, "losv6_resp_minresp.csv"),row.names=FALSE )

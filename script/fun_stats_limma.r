# this contains functions for r

# This function takes in the formula for the design, attach a ~0+ before it and returns a design matrix
fun.makeDesign <- function(current_samp_annotation, design_formula){
	design_formula_plain <- gsub(" ", "", design_formula)
	design_formula_plain <- gsub("~0+", "", design_formula_plain)
	design_formula_plain <- gsub("~", "", design_formula_plain)
	design_var <- strsplit(design_formula_plain, "\\+|:")[[1]]
	design_var <- design_var[design_var != ""]
	out_design <- model.matrix(eval(parse(text=design_formula)) , data = current_samp_annotation)
	for(i in 1:length(design_var)){
		colnames(out_design) <- gsub(design_var[i], "", colnames(out_design))
	}
	colnames(out_design) <- gsub(":", "", colnames(out_design))
	colnames(out_design) <- gsub("grouping", "", colnames(out_design))
	colnames(out_design) <- gsub("-", "", colnames(out_design))
	return(out_design)
}

#### STEP 2: differential gene expression analysis ####
# make a list of contrasts
fun.contrastLos <- function(current_design){
	out_contrasts <- makeContrasts(
					losN = postN-preN, 
					losS6 = postS6-preS6, 
					losV6 = postV6-preV6, 
					losD2 = postD2-preD2, 
					levels=current_design)
	return(out_contrasts)
}

fun.contrastLosOnly <- function(current_design){
	out_contrasts <- makeContrasts(
					losN = postN-preN, 
					losS6 = postS6-preS6, 
					losV6 = postV6-preV6, 
					losD2 = postD2-preD2, 
					levels=current_design)
	return(out_contrasts)
}
fun.contrastLosTreat <- function(current_design){
	out_contrasts <- makeContrasts(
					LosS6N = postS6 - postN - preS6 + preN, 
					LosV6N = postV6 - postN - preV6 + preN, 
					LosD2N = postD2 - postN - preD2 + preN, 
					levels=current_design)
	return(out_contrasts)
}
fun.contrastResp <- function(current_design){
	out_contrasts <- makeContrasts(
					NonResp = postnonresponder - prenonresponder, 
					Resp = postresponder - preresponder, 
					levels=current_design)
	return(out_contrasts)
}
fun.contrastTreat <- function(current_design){
	out_contrasts <- makeContrasts(
					pre = preT - preN, 
					post = postT - postN, 	
					N = postN - preN,
					T = postT - preT,
					levels=current_design)
	return(out_contrasts)
}
fun.contrastTreat2 <- function(current_design){
	out_contrasts <- makeContrasts(
					T = postT - preT,
					levels=current_design)
	return(out_contrasts)
}
fun.contrastResp_Treat <- function(current_design){
	out_contrasts <- makeContrasts(
					NonRespPre = prenonresponderT - prenonresponderN, 
					RespPre = preresponderT - preresponderN, 
					NonRespPost = postnonresponderT - postnonresponderN, 
					RespPost = postresponderT - postresponderN, 
					NonRespLosDiffN = postnonresponderN - prenonresponderN,
					RespLosDiffN = postresponderN - preresponderN,
					NonRespLosDiffT = postnonresponderT - prenonresponderT,
					RespLosDiffT = postresponderT - preresponderT,
					
					levels=current_design)
	return(out_contrasts)
}
fun.contrastLosTreatgN <- function(current_design){
	out_contrasts <- makeContrasts(
					LosS6N = postS6 - N - preS6 + N, 
					LosV6N = postV6 - N - preV6 + N, 
					LosD2N = postD2 - N - preD2 + N, 
					levels=current_design)
	return(out_contrasts)
}
fun.contrastTreatINTLos <- function(current_design){
	out_contrasts <- makeContrasts(
					preS6N = preS6-preN, 
					postS6N = postS6-postN, 
					preV6N = preV6-preN, 
					postV6N = postV6-postN, 
					preD2N = preD2-preN, 
					postD2N = postD2-postN, 
					levels=current_design)
	return(out_contrasts)
}
fun.contrastTreatINTLosgN <- function(current_design){
	out_contrasts <- makeContrasts(
					preS6N = preS6-N, 
					postS6N = postS6-N, 
					preV6N = preV6-N, 
					postV6N = postV6-N, 
					preD2N = preD2-N, 
					postD2N = postD2-N, 
					levels=current_design)
	return(out_contrasts)
}
fun.contrastTreatINTlosINTresponse <- function(current_design){
	out_contrasts <- makeContrasts(
					N = (postNresponder-preNresponder) - (postNnonresponder -preNnonresponder),
					S6 = (postS6responder-preS6responder) - (postS6nonresponder -preS6nonresponder), 
					V6 = (postV6responder-preV6responder) - (postV6nonresponder -preV6nonresponder),
					D2 = (postD2responder-preD2responder) - (postD2nonresponder -preD2nonresponder),
					levels=current_design)
	return(out_contrasts)
}
fun.contrastLos_response <- function(current_design){
	out_contrasts <- makeContrasts(
					RespN = postNresponder-preNresponder,
					RespS6 = postS6responder-preS6responder, 
					RespV6 = postV6responder-preV6responder,
					RespD2 = postD2responder-preD2responder,
					NonRespN = postNnonresponder-preNnonresponder,
					NonRespS6 = postS6nonresponder-preS6nonresponder, 
					NonRespV6 = postV6nonresponder-preV6nonresponder,
					NonRespD2 = postD2nonresponder-preD2nonresponder,
					levels=current_design)
	return(out_contrasts)
}
fun.contrastLos_NonResponder <- function(current_design){
	out_contrasts <- makeContrasts(
					losN = postNnonresponder -preNnonresponder,
					losS6 = postS6nonresponder -preS6nonresponder, 
					losV6 = postV6nonresponder -preV6nonresponder,
					losD2 = postD2nonresponder -preD2nonresponder,
					levels=current_design)
	return(out_contrasts)
}
fun.contrastVSNormal <- function(current_design){
	out_contrasts <- makeContrasts(
					PreS6Resp = preS6responder-preNresponder,
					PreS6NonResp = preS6nonresponder - preNnonresponder,
					PreV6Resp =  preV6responder-preNresponder,
					PreV6NonResp = preV6nonresponder - preNnonresponder,
					PreD2Resp = preD2responder-preNresponder,
					PreD2NonResp = preD2nonresponder - preNnonresponder,
					PostS6Resp = postS6responder-postNresponder,
					PostS6NonResp = postS6nonresponder - postNnonresponder,
					PostV6Resp =  postV6responder-postNresponder,
					PostV6NonResp = postV6nonresponder - postNnonresponder,
					PostD2Resp = postD2responder-postNresponder,
					PostD2NonResp = postD2nonresponder - postNnonresponder,
					levels=current_design)
	return(out_contrasts)
}

fun.topGeneAnalysis <- function(dge_full, gene_annotation, current_samp_annotation, current_design, current_contrast=NULL, current_block=NULL, write_exprs = NULL, current_corfit=NULL, return_exprs= FALSE){
	# organise expression data 
	current_dge <- DGEList(counts=dge_full$counts[,row.names(current_samp_annotation)], genes=rownames(dge_full))
	current_dge <- calcNormFactors(current_dge, method="TMM")
	#my_data = cpm(current_dge, log=TRUE, prior.count=2)
	current_v <- voom(current_dge, current_design, plot=TRUE, normalize.method="quantile")# no combat
	#current_v <- current_v$E
	#my_batch = current_samp_annotation$RunID
	#combat_design = model.matrix(~as.factor(treatment) + as.factor(drug_status) + as.factor(response), current_samp_annotation)
	#current_v <- ComBat(dat=my_data, batch=my_batch, mod=combat_design) # with covariate
	#current_v <- ComBat(dat=current_v, batch=my_batch) # no covariate
	# check if contrast fit is required
	if(!is.null(current_block)){
		current_corfit <- duplicateCorrelation(current_v, current_design, block=current_block)
		current_fit <- lmFit(current_v, current_design, correlation=current_corfit$consensus, block=current_block)
	} else {
		current_fit <- lmFit(current_v, current_design)
	}
	# check if it is contrast current_fit, insert the contrasts
	if(!is.null(current_contrast)){ current_fit <- contrasts.fit(current_fit, current_contrast)} 
	# fit eBayes and find toptable
		current_fit <- eBayes(current_fit) 
	if(!is.null(write_exprs)){ 
		# First write the table to the file (this is so that I don't need to worry about numeric/data.matrix format/row names)
		write.table(as.data.frame(t(current_samp_annotation)),  sep = ",", write_exprs, col.names=TRUE, row.names=TRUE)
		write.table(current_v, write_exprs, col.names=FALSE, sep = ",", append=TRUE, row.names=TRUE)
		# Add gene annotation
		temp_exprs <- read.csv(write_exprs, stringsAsFactors=FALSE) 
		colnames(temp_exprs)[1] <- "geneIDs"
		first_col <- temp_exprs$geneIDs
		temp_exprs <- join(gene_annotation,temp_exprs, type="right")
		temp_exprs$geneIDs[1:nrow(as.data.frame(t(current_samp_annotation)))] <- first_col[1:nrow(as.data.frame(t(current_samp_annotation)))]
		temp_exprs$geneIDs <- ifelse(!is.na(temp_exprs$geneNames), paste(temp_exprs$geneNames, temp_exprs$geneIDs, sep="_"), temp_exprs$geneIDs)
		write.csv(temp_exprs, write_exprs, row.names=FALSE)
	}
	# If there are more than one coefficient, print out a list of annotated topTables
	# This function takes in a fit and print out a list of annotated top tables, with each coefficient in each member of the list
	if(!is.null(current_contrast)){
		out_top_genes = list()
		for (i in 1:ncol(current_contrast)){
			current_toptable <- topTable(current_fit, adjust.method="fdr", sort.by="none", number=Inf, coef=i)
			current_toptable <- merge(gene_annotation, current_toptable, by.x="geneIDs", by.y=0)
			# annotate genes
			out_top_genes[[colnames(current_contrast)[i]]] = current_toptable
			rm(current_toptable)
		}
		current_toptable = topTable(current_fit, adjust.method="fdr", sort.by="none", number=Inf)
		current_toptable <- merge(gene_annotation, current_toptable, by.x="geneIDs", by.y=0)
		out_top_genes[["overall"]] = current_toptable
		rm(current_toptable)
	} else {
		out_top_genes <- topTable(current_fit, adjust.method="fdr", sort.by="none", number=Inf)
		out_top_genes <- merge(gene_annotation, out_top_genes, by.x="geneIDs", by.y=0)
	}
	if(return_exprs){
		return(list(top_genes = out_top_genes, analysed_geneIDs = row.names(current_v), ebayes=current_fit, exprs = current_v, dge= current_dge, design = current_design, contrast=current_contrast , block = current_block, corfit = current_corfit))
	} else {
		return(list(top_genes = out_top_genes, analysed_geneIDs = row.names(current_v), ebayes=current_fit))
	}
}

# this function takes in a list of fold changes
# loop through each element
# if the element name is not "overall"
# print out number of rows with logFC over a specified value
# If the element name is "overall"
# Find the maximum logFC and print out the number of rows with at least 1 comparison reaching that FC
fun.printTopTableRows <- function(x, logfc_threshold = 2, prefix=""){
	for (i in 1:length(x)){
		if(names(x)[i] %in% "overall"){
			if(nrow(x[["overall"]])>0){
				current_logFCs <- x[["overall"]][,names(x)[1:(length(x)-1)]]
				current_logFCs <- abs(data.matrix(current_logFCs))
				current_logFCs <- apply(current_logFCs, 1, max)
				print(length(current_logFCs[current_logFCs>log2(logfc_threshold*1)]))
				return(x[["overall"]][current_logFCs>log2(logfc_threshold*1),])
			} else {
				print(0)
			}
		} else {
			if(nrow(x[[i]])>0){
				outval <- nrow(x[[i]][((abs(as.numeric(x[[i]]$logFC))>log2(logfc_threshold)) &(x[[i]]$adj.P.Val < 0.05)),])
			} else {
				outval <- 0
			}
			print(cat(prefix, " ", names(x)[i], "\t", outval, "\t"))
		}
	}
}

fun.write_TopTables <- function(stats_result, out_dir, output_prefix){
	top_gene_result = stats_result[["top_genes"]]
	for (i in 1:length(top_gene_result)){
		if(nrow(top_gene_result[[i]]) > 0){
			output_fp=file.path(out_dir, paste(output_prefix, "_", names(top_gene_result)[i], ".csv", sep=""))
			write.csv(top_gene_result[[i]], output_fp, row.names=FALSE)
		}
	}
}


fun.mroast <- function(current_object, gene_symbol_list){
	current_v <- current_object$exprs
	matched_gene_ids <- lapply(gene_symbol_list, FUN=function(x)gene_annotation$geneIDs[gene_annotation$geneNames %in% x])
	matched_idx <- lapply(matched_gene_ids, FUN=function(x)match(x, rownames(current_v)))
	matched_idx <- lapply(matched_idx, FUN=function(x)x[!is.na(x)])
	matched_idx[sapply(matched_idx, FUN=function(x)length(x)<5)] <- NULL
	average_logfc <- lapply(matched_gene_ids, FUN=function(x)current_object$top_genes$logFC[current_object$top_genes$geneIDs %in% x])
	average_logfc[sapply(average_logfc, FUN=function(x)length(x)<5)] <- NULL
	average_logfc <- sapply(average_logfc, mean)
	average_fc_df <- data.frame(average_logfc = average_logfc)
	row.names(average_fc_df) <- names(average_logfc)
	mroast_out <- mroast(current_v, index=matched_idx, design = current_object$design, correlation = current_object$corfit$consensus, block = current_object$block, nrot=9999)
	mroast_out <- merge(mroast_out, average_fc_df, by=0)
	return(mroast_out)
}
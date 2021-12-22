###### change gene format ######
###### ID~~gene symbol #########


Clean_of_many_tags <- function(x){
	###########
	pattern = " (1 of many)."
	patterns = paste(pattern,1:100,sep='')
	patterns = c(" (1 of many)",patterns)
	class(x)
	for(p in patterns){
		x = gsub(p,'',x,fixed=T)
	}
	return(x)
}


Process_gene_features <- function(gene_features){
	ID = gene_features$V1
	symbol = gene_features$V2
	#### remove some tags #######
	symbol_cl = Clean_of_many_tags(symbol)
	####
	out = paste(ID,symbol,sep='~~')
	return(out)
}


Convert_to_seurat <- function(Mat,cell_id,gene_features){
	library(Matrix)
	########
	rownames(Mat) = gene_features
	colnames(Mat) = cell_id
	########
	library(Seurat)
	########
	Seurat_obj = CreateSeuratObject(Mat)
	######## return Seurat_obj ######## 
	return(Seurat_obj)
}

############
############


SCRNA_process_S1 <- function(matrix_folder,output_folder,output_tags,MT_tags,useProtein_coding=T){
	library(Seurat)
	library(Matrix)
	#####
	setwd(matrix_folder)
	#####
	mat = readMM('matrix.mtx.gz')
	features = read.table('features.tsv.gz',sep='\t')
	barcodes = read.table('barcodes.tsv.gz',sep='\t')
	##### keep the gene id and gene name #####
	gene_index = paste(features$V1,features$V2,sep='~~')
	cell_index = paste(barcodes$V1)
	colnames(mat) = cell_index
	rownames(mat) = gene_index
	##### extract RNA seq data #####
	##### k = which(features$V3 == 'Gene Expression')
	#####
	mat_cl = mat
	if(useProtein_coding == T){
		print('Zebrafish')
		load('/zp1/data/plyu3/NAR_paper_database/test_Zebrafish/protein_coding_genes')
		k = which(rownames(mat_cl) %in% protein_coding_genes == T)
		print(dim(mat_cl))
		mat_cl = mat_cl[k,]
		print(dim(mat_cl))
	}
	##### load into Seurat ####
	RNA_Seurat <- CreateSeuratObject(counts = mat_cl,min.cells=0,min.features=0)
	RNA_Seurat[["percent.mt"]] <- PercentageFeatureSet(RNA_Seurat, pattern = MT_tags)
	##### see the QC for the RNA ####
	setwd(output_folder)
	#####
	png_file = paste(output_tags,'rawQC.png',sep='_')
	png(png_file,height=10000,width=30000,res=72*12)
	print(VlnPlot(RNA_Seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3))
	dev.off()
	######
	###### save the raw data seurat ######
	raw_file = paste(output_tags,'Seurat_RNA_raw',sep='_')
	saveRDS(RNA_Seurat,file=raw_file)
}




SCRNA_process_S2 <- function(output_folder,output_tags){
	library(Seurat)
	#######
	setwd(output_folder)
	clean_file = paste(output_tags,'Seurat_RNA_clean_S',sep='_')
	print(clean_file)
	#######
	x = readRDS(clean_file)
	######
	library(future)
	plan("multiprocess", workers = 30)
	options(future.globals.maxSize = 10000 * 1024^2)
	######
	library(dplyr)
	x <- SCTransform(x, verbose = FALSE) %>% RunPCA() %>% RunUMAP(dims = 1:50, reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_')
	x <- FindNeighbors(x, dims = 1:50)
	###### names(x@graphs) #####
	x <- FindClusters(x,algorithm = 3,verbose = FALSE)
	######
	######
	setwd(output_folder)
	png_file = paste(output_tags,'_cluster.png',sep='')
	library(ggplot2)
	png(png_file,height=4000,width=6000,res=72*12)
	print(DimPlot(x, reduction = "umap.rna", group.by = "seurat_clusters", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("RNA_umap"))
	dev.off()
	#######
	setwd(output_folder)
	png_file = paste(output_tags,'_scrublet.png',sep='')
	png(png_file,height=4000,width=6000,res=72*12)
	print(FeaturePlot(x, reduction = "umap.rna", features=c('scrublet'), label.size = 2.5, repel = TRUE) + ggtitle("RNA_doublet"))
	dev.off()
	#######
	clean_file = paste(output_tags,'Seurat_RNA_clean_SSS',sep='_')
    setwd(output_folder)
	saveRDS(x,file=clean_file)
    print('Done!')
}





SCRNA_process_S3 <- function(output_folder,output_tags,seurat_query){
	########
	clean_file = paste(output_tags,'Seurat_RNA_clean_SSS',sep='_')
    setwd(output_folder)
	x = readRDS(file=clean_file)
	########
	DefaultAssay(x) = 'RNA'
	######## Run CCA piplines #######
	########
	library(future)
	plan("multicore", workers = 30)
	options(future.globals.maxSize = 10000 * 1024^2)
	########
	x <- NormalizeData(x, verbose = TRUE)
	########
	x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000,
        verbose = FALSE)
	x <- ScaleData(x)
	x <- RunPCA(x)
	######
	DefaultAssay(seurat_query) = 'RNA'
	seurat_query <- NormalizeData(seurat_query, verbose = TRUE)
	seurat_query <- FindVariableFeatures(seurat_query, selection.method = "vst", nfeatures = 2000,
        verbose = FALSE)
	seurat_query <- ScaleData(seurat_query)
	seurat_query <- RunPCA(seurat_query)
	######
	Xanchors <- FindTransferAnchors(reference = seurat_query, query = x,dims = 1:30, reference.reduction = "pca")
	########
	predictions <- TransferData(anchorset = Xanchors, refdata = seurat_query$celltypes,dims = 1:30)
	x <- AddMetaData(x, metadata = predictions)
	########
	### plot the x annotations ########
	########
	######## x <- RunTSNE(x, nn.name = "weighted.nn", reduction.name = "wnn.tsne", reduction.key = "wnnTSNE_",dims = 1:35)
	########
	png_file = paste(output_tags,'_celltypes_umap.png',sep='')
	library(ggplot2)
	png(png_file,height=4000,width=6000,res=72*12)
	print(DimPlot(x, reduction = "umap.rna", group.by = "predicted.id", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("umap.rna"))
	dev.off()
	########
	clean_file = paste(output_tags,'Seurat_RNA_merge_addCT',sep='_')
    setwd(output_folder)
	saveRDS(x,file=clean_file)
	########
	print('Done!!!!')
}








SCRNA_process_S4 <- function(output_folder,output_tags){
	clean_file = paste(output_tags,'Seurat_RNA_merge_addCT',sep='_')
	setwd(output_folder)
	x = readRDS(file=clean_file)	
	DefaultAssay(x) = 'SCT'
	library(future)
	plan("multicore", workers = 1)
	options(future.globals.maxSize = 10000 * 1024^2)
	x <- FindClusters(x,algorithm = 3,verbose = FALSE,resolution=5)
	#### re-cluster with high resolution ##########
	clusters = table(x$seurat_clusters)
	print(paste('trouble shooting:','Num of clusters:',length(clusters),sep=' '))
	meta_data = data.frame(x@meta.data)
	meta_data_cl = meta_data[which(meta_data$prediction.score.max > 0.5),]
	#### prediction.score.max > 0.5 ######
	clusters = table(meta_data_cl$seurat_clusters)
	x$celltypes_sm = 'Unknown'
	for (i in names(clusters)){
		meta_data_cl_sub = meta_data_cl[which(meta_data_cl$seurat_clusters == i),]
		subres = table(meta_data_cl_sub$predicted.id) / length(meta_data_cl_sub$predicted.id)
		sub_cluster = subres[which(subres==max(subres))]
		sub_celltype = names(sub_cluster)
		sub_ratio = round(as.numeric(sub_cluster),3)
		message = paste('cluster:',i,'celltype predict:',sub_celltype,'ratio',sub_ratio,sep=' ')
		print(message)
	###### if ########
		if (sub_ratio > 0.5){
			index = which(x$seurat_clusters == i)
			x$celltypes_sm[index] = sub_celltype
		}
	}
	#### plot the results: ################
	png_file = paste(output_tags,'_cluster_before.png',sep='')
	library(ggplot2)
	png(png_file,height=4000,width=6000,res=72*12)
	print(DimPlot(x, reduction = "umap.rna", group.by = "seurat_clusters", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("umap.rna"))
	dev.off()
	#### plot the results: ################
	png_file = paste(output_tags,'_sm_before.png',sep='')
	library(ggplot2)
	png(png_file,height=4000,width=6000,res=72*12)
	print(DimPlot(x, reduction = "umap.rna", group.by = "predicted.id", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("umap.rna"))
	dev.off()
	#### plot the results: ################
	png_file = paste(output_tags,'_sm_after.png',sep='')
	library(ggplot2)
	png(png_file,height=4000,width=6000,res=72*12)
	print(DimPlot(x, reduction = "umap.rna", group.by = "celltypes_sm", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("umap.rna"))
	dev.off()
	#### save the results ####
	clean_file = paste(output_tags,'Seurat_RNA_merge_CTsm',sep='_')
    setwd(output_folder)
	saveRDS(x,file=clean_file)
}





SCRNA_process_S5 <- function(output_folder,output_tags){
	clean_file = paste(output_tags,'Seurat_RNA_merge_addCT',sep='_')
	setwd(output_folder)
	x = readRDS(file=clean_file)	
	library(Seurat)
	######## percent.mt scrublet ######
	png_file = paste(output_tags,'_sm_features.png',sep='')
	library(ggplot2)
	png(png_file,height=8000,width=10000,res=72*12)
	print(FeaturePlot(x, reduction = "umap.rna", features = c('percent.mt','scrublet','nCount_RNA','nFeature_RNA'), label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("umap.rna"))
	dev.off()
	######## ################### ######
}



Get_gene_names = function(input,query){
	out = c()
	for(i in 1:length(query)){
		query_i = query[i]
		print(query_i)
		query_i = paste(query_i,'$',sep='')
		k = grep(query_i,input)
		if(length(k) > 0){
			out_cl = input[k]
			out = c(out,out_cl)
		}else{
			print('Not find')
		}
	}
	return(out)
}






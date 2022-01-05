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




SCRNA_process_S2 <- function(output_folder,output_tags,vars.to.regress=T){
	library(Seurat)
	#######
	setwd(output_folder)
	clean_file = paste(output_tags,'Seurat_RNA_clean',sep='_')
	print(clean_file)
	#######
	x = readRDS(clean_file)
	######
	library(future)
	plan("multiprocess", workers = 30)
	options(future.globals.maxSize = 10000 * 1024^2)
	######
	library(dplyr)
	######
	if(vars.to.regress){
		print('yes')
		x <- SCTransform(x, verbose = FALSE,vars.to.regress='percent.mt') %>% RunPCA() %>% RunUMAP(dims = 1:50, reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_')
	}else{
		print('no')
		x <- SCTransform(x, verbose = FALSE) %>% RunPCA() %>% RunUMAP(dims = 1:50, reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_')
	}
	
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


#### clean the scurblet ########
#### cutoff ####################


SCRNA_process_S6 <- function(output_folder,output_tags,cutoff=0.25){
	library(Seurat)
	clean_file = paste(output_tags,'Seurat_RNA_merge_addCT',sep='_')
    setwd(output_folder)
	x = readRDS(file=clean_file)
	######## ################### ######
	########
	x_cl = x[,which(x$scrublet < cutoff)]
	########
	########
	DefaultAssay(x_cl) = 'RNA'
	######### recluster #######
	library(dplyr)
	x_cl <- SCTransform(x_cl, verbose = FALSE,vars.to.regress='percent.mt') %>% RunPCA() %>% RunUMAP(dims = 1:50, reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_')
	x_cl <- FindNeighbors(x_cl, dims = 1:50)
	###### names(x@graphs) #####
	x_cl <- FindClusters(x_cl,algorithm = 3,verbose = FALSE)
	x_cl <- RunUMAP(x_cl,dims=1:50,reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_')
	setwd(output_folder)
	png_file = paste(output_tags,'_ct_filter.png',sep='')
	library(ggplot2)
	png(png_file,height=4000,width=6000,res=72*12)
	print(DimPlot(x_cl, reduction = "umap.rna", group.by = "predicted.id", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("RNA_umap"))
	dev.off()
	#####
	clean_file = paste(output_tags,'Seurat_RNA_merge_CTsm_filter',sep='_')
    setwd(output_folder)
	saveRDS(x_cl,file=clean_file)
	########
	print('Done!!!!')

}

#####

SCRNA_process_S7 <- function(output_folder,output_tags,markers_list){
	clean_file = paste(output_tags,'Seurat_RNA_merge_CTsm_filter',sep='_')
	setwd(output_folder)
	x = readRDS(file=clean_file)	
	library(Seurat)
	######## find Genes ######
	for(i in 1:length(markers_list)){
		features = Get_gene_names(rownames(x),markers_list[[i]])
		print(features)
		png_file = paste(output_tags,'_',names(markers_list)[i],'_features.png',sep='')
		print(png_file)
		library(ggplot2)
		png(png_file,height=8000,width=10000,res=72*12)
		print(FeaturePlot(x, reduction = "umap.rna", features = features, label = TRUE, label.size = 1, repel = TRUE))
		dev.off()
	}
	###################################
	png_file = paste(output_tags,'_','Marker','_seurat_clusters.png',sep='')
	print(png_file)
	library(ggplot2)
	png(png_file,height=8000,width=10000,res=72*12)
	print(DimPlot(x, reduction = "umap.rna", group.by='seurat_clusters', label = TRUE, label.size = 5, repel = TRUE))
	dev.off()
	png_file = paste(output_tags,'_','Predict','_seurat_clusters.png',sep='')
	print(png_file)
	library(ggplot2)
	png(png_file,height=8000,width=10000,res=72*12)
	print(DimPlot(x, reduction = "umap.rna", group.by='predicted.id', label = TRUE, label.size = 5, repel = TRUE))
	dev.off()
}


SCRNA_process_S9 <- function(output_folder,output_tags,markers_list){
	clean_file = paste(output_tags,'Seurat_RNA_merge_addCT_tile_filter',sep='_')
	setwd(output_folder)
	x = readRDS(file=clean_file)	
	library(Seurat)
	######## find Genes ######
	for(i in 1:length(markers_list)){
		features = Get_gene_names(rownames(x),markers_list[[i]])
		print(features)
		png_file = paste(output_tags,'_',names(markers_list)[i],'_features.png',sep='')
		print(png_file)
		library(ggplot2)
		png(png_file,height=8000,width=10000,res=72*12)
		print(FeaturePlot(x, reduction = "umap.rna", features = features, label = TRUE, label.size = 1, repel = TRUE))
		dev.off()
	}
	###################################
	png_file = paste(output_tags,'_','Marker','_seurat_clusters.png',sep='')
	print(png_file)
	library(ggplot2)
	png(png_file,height=8000,width=10000,res=72*12)
	print(DimPlot(x, reduction = "umap.rna", group.by='seurat_clusters', label = TRUE, label.size = 5, repel = TRUE))
	dev.off()
	png_file = paste(output_tags,'_','Predict','_seurat_clusters.png',sep='')
	print(png_file)
	library(ggplot2)
	png(png_file,height=8000,width=10000,res=72*12)
	print(DimPlot(x, reduction = "umap.rna", group.by='predicted.id', label = TRUE, label.size = 5, repel = TRUE))
	dev.off()
}




SCRNA_process_S8 <- function(output_folder,output_tags,celltype_list){
	clean_file = paste(output_tags,'Seurat_RNA_merge_CTsm_filter',sep='_')
	setwd(output_folder)
	x = readRDS(file=clean_file)	
	library(Seurat)
	x$new_celltypes = 'Unknown'
	######## find Genes ######
	for(i in 1:length(celltype_list)){
		print(celltype_list[[i]])
		print(names(celltype_list)[i])
		k = which(x$seurat_clusters %in% celltype_list[[i]] == T)
		x$new_celltypes[k] = names(celltype_list)[i]
	}
	###################################
	png_file = paste(output_tags,'_','Marker_new','.png',sep='')
	print(png_file)
	library(ggplot2)
	png(png_file,height=6000,width=8000,res=72*12)
	print(DimPlot(x, reduction = "umap.rna", group.by='new_celltypes', label = TRUE, label.size = 5, repel = TRUE))
	dev.off()
	###################################
	clean_file = paste(output_tags,'Seurat_RNA_merge_CTsm_filter_Markers',sep='_')
    setwd(output_folder)
	saveRDS(x,file=clean_file)
}




SCRNA_process_S3_revise <- function(output_folder,output_tags,seurat_query){
	########
	clean_file = paste(output_tags,'Seurat_RNA_clean_SSS_filter',sep='_')
    setwd(output_folder)
	x = readRDS(file=clean_file)
	########
	DefaultAssay(x) = 'SCT'
	######## Run CCA piplines #######
	########
	library(future)
	plan("multicore", workers = 30)
	options(future.globals.maxSize = 10000 * 1024^2)
	########
	DefaultAssay(seurat_query) = 'RNA'
	seurat_query = SCTransform(seurat_query, verbose = FALSE) %>% RunPCA()
	######
	#######
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




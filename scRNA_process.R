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


SCRNA_process_S1 <- function(matrix_folder,output_folder,output_tags,MT_tags){
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






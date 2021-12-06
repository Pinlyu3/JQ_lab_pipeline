#### functions for Multiome cellranger ####

#### matrix_folder: filtered_feature_bc_matrix folder output by cellranger #####
#### output_tags: sample name #####
#### MT_tags: genes which are located in Mitochondria #####

Muti_process_S1 <- function(matrix_folder,output_folder,output_tags,MT_tags){
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
	k = which(features$V3 == 'Gene Expression')
	#####
	mat_cl = mat[k,]
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


#### output_folder in step1 #####
#### output_tags in step1 #######

Muti_process_S2 <- function(output_folder,output_tags,nFeature_RNA,nCount_RNA,mt_RNA){
	library(Seurat)
	#######
	setwd(output_folder)
	raw_file = paste(output_tags,'Seurat_RNA_raw',sep='_')
	print(raw_file)
	#######
	RNA_Seurat = readRDS(raw_file)
	#######
	k_nFeature = which(RNA_Seurat$nFeature_RNA > nFeature_RNA[1] & RNA_Seurat$nFeature_RNA < nFeature_RNA[2])
	k_nCount = which(RNA_Seurat$nCount_RNA > nCount_RNA[1] & RNA_Seurat$nCount_RNA < nCount_RNA[2])
	k_mt = which(RNA_Seurat$percent.mt > mt_RNA[1] & RNA_Seurat$percent.mt < mt_RNA[2])
	#######
	k_merge = intersect(k_nFeature,k_nCount)
	k_merge = intersect(k_merge,k_mt)
	#######
	RNA_Seurat_cl1 = RNA_Seurat[,k_merge]
	#######
	print(paste('Cells remain', round(dim(RNA_Seurat_cl1)[2]/dim(RNA_Seurat)[2],2)))
	###### save the cleaned data seurat ######
	clean_file = paste(output_tags,'Seurat_RNA_clean',sep='_')
	saveRDS(RNA_Seurat_cl1,file=clean_file)
	#######
}



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


#### Seurat to scrublet input ####

Muti_process_S3 <- function(output_folder,output_tags){
	#######
	library(Seurat)
	#######
	setwd(output_folder)
	clean_file = paste(output_tags,'Seurat_RNA_clean',sep='_')
	print(clean_file)
	#######
	x = readRDS(clean_file)
	#######
	library(Matrix)
    mat = x[['RNA']]@counts
    gene = rownames(mat)
    gene = data.frame(V1=gene,V2=gene)
    barcode = colnames(mat)
    barcode = data.frame(V1=barcode,V2=barcode)
    ########
    FN1 = paste(output_tags,'scrublet_mat.mtx',sep='_')
    FN2 = paste(output_tags,'scrublet_gene.tsv',sep='_')
    FN3 = paste(output_tags,'scrublet_barcode.tsv',sep='_')
    writeMM(mat,file=FN1)
    write.table(gene,file=FN2,sep='\t',quote=F,col.names=F,row.names=F)
    write.table(barcode,file=FN3,sep='\t',quote=F,col.names=F,row.names=F)
    ###
    print('Done!')
	#######
}



#### Merge the scrublet score to seurat_obj ####

Muti_process_S4 <- function(output_folder,output_tags){
	#######
	library(Seurat)
	#######
	setwd(output_folder)
	clean_file = paste(output_tags,'Seurat_RNA_clean',sep='_')
	print(clean_file)
	#######
	x = readRDS(clean_file)
	#######
	FN4 = paste(output_tags,'scrublet_res.tsv',sep='_')
	score = read.table(FN4,sep='\t')
    ###
    m = match(colnames(x),score$V1)
    x$scrublet = score$V2[m]
    print(summary(x$scrublet))
    #######
    clean_file = paste(output_tags,'Seurat_RNA_clean_S',sep='_')
	saveRDS(x,file=clean_file)
    print('Done!')
	#######
}


####  

ArchR_zebrafish_genes <- function(gtf){
	library(rtracklayer)
	library(ArchR)
	######
	gtf <- rtracklayer::import(gtf)
	gtf_df <- as.data.frame(gtf)
	##### gene #######
	k = which(gtf_df$type == 'gene')
	gtf_df_gene = gtf_df[k,]
	chr = c(1:25)
	gtf_df_gene = gtf_df_gene[which(gtf_df_gene$seqnames %in% chr == T),]
	genes_GR = GRanges(seqnames=paste('chr',gtf_df_gene$seqnames,sep=''),IRanges(start=gtf_df_gene$start,end=gtf_df_gene$end),strand=gtf_df_gene$strand,gene_id=gtf_df_gene$gene_id,symbol=gtf_df_gene$gene_name)
	##### exons ######
	k = which(gtf_df$type == 'exon')
	gtf_df_exon = gtf_df[k,]
	chr = c(1:25)
	gtf_df_exon = gtf_df_exon[which(gtf_df_exon$seqnames %in% chr == T),]
	exons_GR = GRanges(seqnames=paste('chr',gtf_df_exon$seqnames,sep=''),IRanges(start=gtf_df_exon$start,end=gtf_df_exon$end),strand=gtf_df_exon$strand,gene_id=gtf_df_exon$gene_id,symbol=gtf_df_exon$gene_name)
	##### TSS ######
	gtf_df_trans = gtf_df[which(gtf_df$type == 'transcript'),]
	chr = c(1:25)
	gtf_df_trans = gtf_df_trans[which(gtf_df_trans$seqnames %in% chr == T),]
	trans_GR = GRanges(seqnames=paste('chr',gtf_df_trans$seqnames,sep=''),IRanges(start=gtf_df_trans$start,end=gtf_df_trans$end),strand=gtf_df_trans$strand,tx_id=gtf_df_trans$transcript_id,tx_name=gtf_df_trans$transcript_name)
	k_plus = which(as.character(strand(trans_GR)) == '+')
	k_minus = which(as.character(strand(trans_GR)) == '-')
	k_plus_TSS = start(trans_GR)[k_plus]
	k_minus_TSS = end(trans_GR)[k_minus]
	TSS_GR = trans_GR
	end(TSS_GR)[k_plus] = k_plus_TSS
	start(TSS_GR)[k_minus] = k_minus_TSS
	###### 
	geneAnnotation_GRCz11 <- ArchR::createGeneAnnotation(
  		TSS = TSS_GR,
  		exons = exons_GR,
  		genes = genes_GR
	)
	######
	return(geneAnnotation_GRCz11)
}




####


Muti_process_ATAC_S2 <- function(TSS_enrich_low,output_folder,output_tags,output_folder2,geneAnnotation,genomeAnnotation){
	#######
	library(ArchR)
	addArchRThreads(threads = 5) 
	####### on windows #####
	#######
	setwd(output_folder)
	#######
	ArrowFiles = paste(output_tags,'.arrow',sep='')
	#######
	Project_1 <- ArchRProject(
  		ArrowFiles = ArrowFiles, 
  		outputDirectory = output_folder2,
  		geneAnnotation = geneAnnotation,
  		genomeAnnotation = genomeAnnotation,
  		copyArrows = TRUE #This is recommened so that if you modify the Arrow files you have an original copy for later usage.
	)
	#######
	setwd(output_folder)
	FN = paste(output_tags,'ArchR_raw',sep='_')
	saveRDS(Project_1,file=FN)
	#######
	idxPass <- which(Project_1$TSSEnrichment >= TSS_enrich_low)
	cellsPass <- Project_1$cellNames[idxPass]
	Project_1_cl = Project_1[cellsPass, ]
	print(length(Project_1$TSSEnrichment))
	print(length(cellsPass))
	#######
	#######
	#### calculate the doublets on the filtered ArchR project #####
	#######
	doubScores <- addDoubletScores(
    	input = Project_1_cl,
    	k = 10, #Refers to how many cells near a "pseudo-doublet" to count.
    	knnMethod = "UMAP", #Refers to the embedding to use for nearest neighbor search with doublet projection.
    	LSIMethod = 1,
    	force=T
	)
	#######
	setwd(output_folder2)
	setwd('ArrowFiles')
	Project_1_cl <- ArchRProject(
  		ArrowFiles = ArrowFiles, 
  		outputDirectory = output_folder2,
  		geneAnnotation = geneAnnotation,
  		genomeAnnotation = genomeAnnotation,
  		copyArrows = FALSE 
	)
	#######
	setwd(output_folder)
	FN = paste(output_tags,'ArchR_cl',sep='_')
	saveRDS(Project_1_cl,file=FN)
	#######
	tab_Dob = data.frame(cellNames = Project_1_cl$cellNames,DoubletScore =Project_1_cl$DoubletEnrichment)
	summary(Project_1_cl$DoubletEnrichment)
	#######
	setwd(output_folder)
	FN = paste(output_tags,'ArchR_doublet.txt',sep='_')
	write.table(tab_Dob,file=FN,row.names=F,sep='\t',quote=F,col.names=F)
	####### Output the tile matrix ##########
	TileMatrix = getMatrixFromProject(
  		ArchRProj = Project_1_cl,
  		useMatrix = "TileMatrix",
  		useSeqnames = NULL,
  		verbose = TRUE,
  		binarize = T,
  		threads = getArchRThreads(),
  		logFile = createLogFile("getMatrixFromProject")
	)
	#######
	FN = paste(output_tags,'ArchR_tileMat',sep='_')
	#######
	setwd(output_folder)
	saveRDS(TileMatrix,file=FN)
	#######
	print('Done!')
	#######
}



#### merge the ArchR doublet score #####
Muti_process_S5 <- function(output_folder,output_tags){
	#######
	library(Seurat)
	#######
	setwd(output_folder)
	clean_file = paste(output_tags,'Seurat_RNA_clean_S',sep='_')
	print(clean_file)
	#######
	x = readRDS(clean_file)
	#######
	FN5 = paste(output_tags,'ArchR_doublet.txt',sep='_')
	#######
	ArchR_doub = read.table(FN5,sep='\t',header=F,comment.char = "")
	#######
    ArchR_doub$Cells = sapply(strsplit(as.character(ArchR_doub$V1),split='#',fixed=T),function(x) x[[2]])
    #######
    m = match(colnames(x),ArchR_doub$Cells)
    x$ArchR = ArchR_doub$V2[m]
    #######
    print(summary(x$ArchR))
    #######
    clean_file = paste(output_tags,'Seurat_RNA_clean_SS',sep='_')
    setwd(output_folder)
	saveRDS(x,file=clean_file)
    print('Done!')
	#######
}





# Pipelines

1. Pipelines for [Chromium Single Cell Multiome ATAC + Gene Expression](#multiome-cellranger)


# Multiome cellranger

## First processing the scRNAseq from the multiome datasets

in this example, the output of cellranger-arc is:
```r
cellranger_folder = '/zp1/data/Share/Fish/Multiome/54hrLD/'
```
we use the matrix and annotation in the filtered_feature_bc_matrix folder:
```r
matrix_folder = '/zp1/data/Share/Fish/Multiome/54hrLD/outs/filtered_feature_bc_matrix'
```

### RNA_Step1: Covert matrix to SeuratObj
Read the mat to seurat and plot the QC 

R package needed: Seurat Matrix devtools

MT_tags: mouse should be ~~mt-

MT_tags: fish should be ~~mt-

MT_tags: human should be ~~MT-

```r
#### load the functions from github ######
devtools::source_url("https://raw.githubusercontent.com/Pinlyu3/JQ_lab_pipeline/main/Muti_process.R")
#### set parameters #####
matrix_folder = '/zp1/data/Share/Fish/Multiome/54hrLD/outs/filtered_feature_bc_matrix'
output_folder = '/zp1/data/plyu3/Muti_omic/54hr_LD'
output_tags = '54hr_LD_202112'
MT_tags = '~~mt-'
#### run the Muti_process_S1 function #####
Muti_process_S1(matrix_folder,output_folder,output_tags,MT_tags)
#### see the output in the output_folder ######
#### check the QC plot: output_tags + _rawQC.png ######
```

### RNA_Step2: Clean cells by nCounts,nFeatures and %mt_RNA
```r
#### load the functions from github ######
devtools::source_url("https://raw.githubusercontent.com/Pinlyu3/JQ_lab_pipeline/main/Muti_process.R")
#### set parameters accoring to QC plot #####
nFeature_RNA = c(500,3000) #### nFeature_RNA > 500 & nFeature_RNA < 3000
nCount_RNA = c(0,10000) #### nCount_RNA > 0 & nCount_RNA < 10000
mt_RNA = c(0,5) #### mt_RNA > 0 & mt_RNA < 5
#### run the function ####
Muti_process_S2(output_folder,output_tags,nFeature_RNA = c(500,3000),nCount_RNA = c(0,10000),mt_RNA = c(0,5))
```

### RNA_Step3: Calculate doublet score using Scrublet

prepare the input to Scrublet
```r
#### load the functions from github ######
devtools::source_url("https://raw.githubusercontent.com/Pinlyu3/JQ_lab_pipeline/main/Muti_process.R")
#### set parameters #####
output_folder = '/zp1/data/plyu3/Muti_omic/54hr_LD'
output_tags = '54hr_LD_202112'
#### convert seurat to Scrublet input ####
Muti_process_S3(output_folder,output_tags)
```

then in the shell, run the script: Muti_process.py

be sure have installed these python packages: scrublet matplotlib scipy numpy os pandas

```console
#### download the python code and run ######
cd /zp1/data/plyu3/Muti_omic
curl https://raw.githubusercontent.com/Pinlyu3/JQ_lab_pipeline/main/Muti_process.py > Muti_process.py
#### run Scrublet #####
python Muti_process.py --output_folder='/zp1/data/plyu3/Muti_omic/54hr_LD' --output_tags='54hr_LD_202112'
```

then adding the doublet score to the SeuratObject
```r
#### load the functions and set the parameters ######
devtools::source_url("https://raw.githubusercontent.com/Pinlyu3/JQ_lab_pipeline/main/Muti_process.R")
output_folder = '/zp1/data/plyu3/Muti_omic/54hr_LD'
output_tags = '54hr_LD_202112'
#### convert seurat to Scrublet input ####
Muti_process_S4(output_folder,output_tags)
```

## Next processing the scATACseq from the multiome datasets

in this example, the output of cellranger-arc is:
```r
cellranger_folder = '/zp1/data/Share/Fish/Multiome/54hrLD/'
```

we use the fragments in the outs folder:
```r
fragments = "/zp1/data/Share/Fish/Multiome/54hrLD/outs/atac_fragments.tsv.gz"
```

Before running, make sure R package ArchR installed

To analysis zebrafish data, we must prepare the genome and gene for ArchR
because ArchR don't include zebrafish reference (ArchR only contain Mouse / Human reference)

if analysis zebrafish data in ArchR:
```r
devtools::source_url("https://raw.githubusercontent.com/Pinlyu3/JQ_lab_pipeline/main/Muti_process.R")
#### we need gtf file for zebrafish genes #####
gtf = '/mnt/d/ArchR_files/Danio_rerio.GRCz11.104.gtf'
#### ArchR file for gtf ############
geneAnnotation_GRCz11 = ArchR_zebrafish_genes(gtf)

#### we need genome sequence #######
#### zebrafish #####################
library("BSgenome.Drerio.UCSC.danRer11")
genomeAnnotation_GRCz11 <- createGenomeAnnotation(genome = BSgenome.Drerio.UCSC.danRer11)

#### save the reference ############
setwd('/mnt/d/ArchR_files/')
save(geneAnnotation_GRCz11,file='geneAnnotation_GRCz11')
save(genomeAnnotation_GRCz11,file='genomeAnnotation_GRCz11')
```


create arrow files from fragments:
```r
####### parameter #######
library(ArchR)
addArchRThreads(threads = 5)
output_tags = '54hr_LD_202112'
atac_fragments_file = "/mnt/d/ArchR_files/54hrLD/atac_fragments.tsv.gz"

####### if zebrafish #######
setwd('/mnt/d/ArchR_files/')
load('geneAnnotation_GRCz11')
load('genomeAnnotation_GRCz11')
geneAnnotation = geneAnnotation_GRCz11
genomeAnnotation = genomeAnnotation_GRCz11
####### if human ###########
addArchRGenome("hg38")
geneAnnotation = getGeneAnnotation()
genomeAnnotation = getGenomeAnnotation()
####### if mouse ###########
addArchRGenome("mm10")
geneAnnotation = getGeneAnnotation()
genomeAnnotation = getGenomeAnnotation()

####### arrow file output location #####
output_folder = '/home/lp123/Desktop/54hr_LD'

####### create arrow files #####
setwd(output_folder)
ArrowFiles <- createArrowFiles(
  inputFiles = atac_fragments_file,
  sampleNames = output_tags,
  minTSS = 4, #Dont set this too high because you can always increase later
  minFrags = 1000, 
  addTileMat = TRUE,
  addGeneScoreMat = FALSE,
  geneAnnotation = geneAnnotation,
  genomeAnnotation = genomeAnnotation,
  force=T
)

####### check the QC files in the output_folder ########
```

filter cells passed QC:
```r
#### filter the cell by TSS enrichment and calculate doublets ######
devtools::source_url("https://raw.githubusercontent.com/Pinlyu3/JQ_lab_pipeline/main/Muti_process.R")
TSS_enrich_low = 7
output_folder = '/home/lp123/Desktop/54hr_LD'
output_tags = '54hr_LD_202112'
#### mkdir '/home/lp123/Desktop/54hr_LD/54hr_LD_2' ######
output_folder2 = '/home/lp123/Desktop/54hr_LD/54hr_LD_tmp2'

Muti_process_ATAC_S2(TSS_enrich_low,output_folder,output_tags,output_folder2,geneAnnotation,genomeAnnotation)

#### 
```


## Last integrate the scRNAseq and scATACseq files 

copy the doublet file and tileMat file to the scRNAseq folder:
output_folder = '/zp1/data/plyu3/Muti_omic/54hr_LD'

```r
devtools::source_url("https://raw.githubusercontent.com/Pinlyu3/JQ_lab_pipeline/main/Muti_process.R")
output_folder = '/zp1/data/plyu3/Muti_omic/54hr_LD'
output_tags = '54hr_LD_202112'
Muti_process_S5(output_folder,output_tags)
```






### Jekyll Themes

Your Pages site will use the layout and styles from the Jekyll theme you have selected in your [repository settings](https://github.com/Pinlyu3/JQ_lab_pipeline/settings/pages). The name of this theme is saved in the Jekyll `_config.yml` configuration file.

### Support or Contact

Having trouble with Pages? Check out our [documentation](https://docs.github.com/categories/github-pages-basics/) or [contact support](https://support.github.com/contact) and we’ll help you sort it out.

### It works well? cheers 


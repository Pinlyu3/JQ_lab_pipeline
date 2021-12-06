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



Syntax highlighted code block

# Header 1
## Header 2
### Header 3

- Bulleted
- List

1. Numbered
2. List

**Bold** and _Italic_ and `Code` text

[Link](url) and ![Image](src)
```

For more details see [Basic writing and formatting syntax](https://docs.github.com/en/github/writing-on-github/getting-started-with-writing-and-formatting-on-github/basic-writing-and-formatting-syntax).

### Jekyll Themes

Your Pages site will use the layout and styles from the Jekyll theme you have selected in your [repository settings](https://github.com/Pinlyu3/JQ_lab_pipeline/settings/pages). The name of this theme is saved in the Jekyll `_config.yml` configuration file.

### Support or Contact

Having trouble with Pages? Check out our [documentation](https://docs.github.com/categories/github-pages-basics/) or [contact support](https://support.github.com/contact) and we’ll help you sort it out.
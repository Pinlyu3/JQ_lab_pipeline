## Pipelines

1. Pipelines for [Chromium Single Cell Multiome ATAC + Gene Expression](#multiome-cellranger)


### Multiome cellranger

first processing the scRNAseq from the multiome datasets

in this example, the output of cellranger-arc is:
```r
cellranger_folder = '/zp1/data/Share/Fish/Multiome/54hrLD/'
```
we use the matrix and annotation in the filtered_feature_bc_matrix folder:
```r
matrix_folder = '/zp1/data/Share/Fish/Multiome/54hrLD/outs/filtered_feature_bc_matrix'
```

#### RNA_Step1: covert to seurat obj
step1: read the mat to seurat and plot the QC 

package needed: Seurat Matrix devtools

MT_tags: mouse: '~~mt-'
MT_tags: fish: '~~mt-'
MT_tags: human: '~~MT-'

```r
#### load the functions from github ######


#### set parameters #####
matrix_folder = '/zp1/data/Share/Fish/Multiome/54hrLD/outs/filtered_feature_bc_matrix'
output_folder = '/zp1/data/plyu3/Muti_omic/54hr_LD'
output_tags = '54hr_LD_202112'


#### 




```





#### RNA_Step2: covert to seurat obj








```markdown
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

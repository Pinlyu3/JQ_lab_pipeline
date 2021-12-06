###### 
###### 
#### functions for the scrublets #####
import scrublet as scr
import scipy.io
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import argparse

#####  #####
p = argparse.ArgumentParser(usage = "Get parameter", description = "scrublet")

p.add_argument('--output_folder')
p.add_argument('--output_tags')

args = p.parse_args()

folder = args.output_folder
index = args.output_tags

print(folder)
print(index)


def scrublet_process(folder,index,doublet_rate=0.1):
    os.chdir(folder)
    counts_matrix = scipy.io.mmread(index + '_scrublet_mat.mtx').T.tocsc()
    genes = np.array(scr.load_genes(index + '_scrublet_gene.tsv', delimiter='\t', column=1))
    barcodes = np.array(scr.load_genes(index + '_scrublet_barcode.tsv', delimiter='\t', column=1))
    #####
    print (input_dir + index + '_scrublet_mat.mtx')
    #####
    scrub = scr.Scrublet(counts_matrix,expected_doublet_rate=doublet_rate)
    doublet_scores, predicted_doublets = scrub.scrub_doublets(min_counts=2, 
                                                          min_cells=3, 
                                                          min_gene_variability_pctl=85, 
                                                          n_prin_comps=30)
    ##### output #####
    barcodes_pd = pd.DataFrame(barcodes)
    barcodes_pd['doublet_scores'] = doublet_scores
    barcodes_pd['predicted_doublets'] = predicted_doublets
    barcodes_pd.rename(columns={0:'barcode','doublet_scores':'doublet_scores'},inplace=True)
    ###### 
    os.chdir(folder)
    Output_FN = index + '_scrublet_res.tsv'
    barcodes_pd.to_csv(Output_FN,sep='\t',index=0,header=0)

######
######

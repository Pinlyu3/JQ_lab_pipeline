#####  #####

import scvelo as scv
import os
import pandas as pd
import subprocess
import argparse
import re
#####  #####

p = argparse.ArgumentParser(usage = "Get parameter", description = "velocity")


p.add_argument('--output_folder')
p.add_argument('--metadata')
p.add_argument('--loomfile')

args = p.parse_args()

metadata = args.metadata
index = args.
output_folder = args.output_folder
loomfile =  args.loomfile
#######
#######
#######

os.chdir(output_folder)
scv.settings.verbosity = 3
scv.settings.set_figure_params('scvelo')

#### read the loom files ################ 
data = scv.read_loom(loomfile)

#### read the meta files ################
Meta_data = pd.read_csv(metadata,sep='\t')

#### filterout cells according to Meta_data #######
cell_id = Meta_data['cell_id'].tolist()

#### cell_id '-1' to x ####
cell_id_new = list()

for i in cell_id:
	i_new = i.replace('-1','x')
	cell_id_new.append(i_new)

print('trouble shooting: cell_id format in metadata is ' + cell_id_new[1])
#### 

Cell_id_all = data.obs.index
Cell_id_all = pd.Series(Cell_id_all,dtype=object)
index = Cell_id_all.isin(cell_id_new)
adata = adata[index,:]

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
p.add_argument('--output_tags')


args = p.parse_args()

metadata = args.metadata
output_folder = args.output_folder
loomfile =  args.loomfile
output_tags = args.output_tags
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

#### cell_id in data format ########

Cell_id_all = data.obs.index
Cell_id_all = pd.Series(Cell_id_all,dtype=object)
Cell_id_all_list = Cell_id_all.tolist()
Cell_id_all_new = list()
for j in Cell_id_all_list:
	j_new = j.split(':')[1]
	Cell_id_all_new.append(j_new)

print('trouble shooting: cell_id format in loom is ' + Cell_id_all_new[1])


#############################
Cell_id_all_new = pd.DataFrame(Cell_id_all_new)
index = Cell_id_all_new.isin(cell_id_new).iloc[:,0]
data_cl = data[index,:]

print('dim of cleaned loom is ' + str(data_cl.shape[0]) + ' ' + str(data_cl.shape[1]))


###### add UMAPs to data_cl ########

Cell_id_all_new_cl = Cell_id_all_new[index].iloc[:,0].tolist()
cell_id_new = list(cell_id_new)

###### 

Cell_index_FN = list()
#####for i, n in enumerate(Cell_id_all_new_cl[1:5]): 
#####	print(i)
#####	print(n) cell_id_new.index(Cell_id_all_new_cl[1])
for i, n in enumerate(Cell_id_all_new_cl): 
	Cell_index_FN.append(cell_id_new.index(n))

##### sort UMAP ##############
UMAP = Meta_data.iloc[:,0:2]
UMAP_M = UMAP.to_numpy()
UMAP_K = UMAP_M[Cell_index_FN,]
data_cl.obsm['X_umap'] = UMAP_K

########################################
#### run velocity on data_cl ###########
########################################

scv.utils.show_proportions(data_cl)
scv.pp.normalize_per_cell(data_cl)
scv.pp.filter_genes_dispersion(data_cl, n_top_genes=1000)
scv.pp.log1p(data_cl)
scv.pp.moments(data_cl, n_pcs=25, n_neighbors=50)
scv.tl.recover_dynamics(data_cl)
scv.tl.velocity(data_cl,mode='dynamical')
scv.tl.velocity_graph(data_cl)

output_file = output_tags + '_velocity.png'
print('output_file')

scv.pl.velocity_embedding_grid(data_cl,basis='umap',color="lightgrey",arrow_length=1, arrow_color='#CB2314',arrow_size=2,density=1,alpha=1,size=2.5,dpi=600,save=output_file)

####################
print('Done!!!')




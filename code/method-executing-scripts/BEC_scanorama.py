import numpy as np
import pandas as pd
import os

import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns

import scanpy as sc
import anndata as ad

import scvi
import torch.nn
from scvi.data import AnnDataManager
from scvi.data.fields import LayerField, CategoricalObsField

import time

if __name__ == '__main__':
    
    allpath = [ 
      "data"
    ]
    
    for path in allpath:    
        
        list_subfolders= [f.path for f in os.scandir(path) if f.is_dir()]
        
        for folder in list_subfolders:            

            ds_file = folder + '/' + "counts.txt"
            mt_file = folder + '/' + "cellinfo.txt"

            print(ds_file)
            print (mt_file)

            # loading data
            sdata = ad.read_csv(ds_file, delimiter='\t', first_column_names=True) # for simulations delimiter = '\t'        
            #sdata = ad.read_csv(ds_file, delimiter=' ', first_column_names=True)  # for Covid19 # for LUAD data delim=" " 

            sdata = sdata.T                
            sdata.layers['counts'] = sdata.X

            pd_obs = pd.read_csv(mt_file, delimiter="\t", header=0, index_col=0)
            #pd_obs = pd.read_csv(mt_file, delimiter=" ", header=0, index_col=0) # for Covid19, LUAD data delim=" "

            binfo = pd_obs['Batch']
            binfo = list(map(str, binfo))
            pd_obs['Batch'] = binfo
            pd_obs.index = sdata.obs_names
            sdata.obs = pd_obs

            sc.pp.filter_genes(sdata, min_cells=0)
            sc.pp.normalize_total(sdata, target_sum = 1e4)
            sc.pp.log1p(sdata)  
            
            start_time = time.time()
            # scanorama                
            scano_data = sdata.copy()

            # luad data
            #"P0018" "P0008" "P0006" "P0020" "P0030" "P0034" "P0019"
            #batch_1 = scano_data[scano_data.obs['Batch'] == 'P0018']
            #batch_2 = scano_data[scano_data.obs['Batch'] == 'P0008']
            #batch_3 = scano_data[scano_data.obs['Batch'] == 'P0006']
            #batch_4 = scano_data[scano_data.obs['Batch'] == 'P0020']
            #batch_5 = scano_data[scano_data.obs['Batch'] == 'P0030']
            #batch_6 = scano_data[scano_data.obs['Batch'] == 'P0034']
            #batch_7 = scano_data[scano_data.obs['Batch'] == 'P0019']        

            # splatter        
            batch_1 = scano_data[scano_data.obs['Batch'] == 'Batch1']
            batch_2 = scano_data[scano_data.obs['Batch'] == 'Batch2']        
            #batch_3 = scano_data[scano_data.obs['Batch'] == 'Batch3']
            #batch_4 = scano_data[scano_data.obs['Batch'] == 'Batch4']
            #batch_5 = scano_data[scano_data.obs['Batch'] == 'Batch5']
            #batch_6 = scano_data[scano_data.obs['Batch'] == 'Batch6']            
            #batch_7 = scano_data[scano_data.obs['Batch'] == 'Batch7']

            # model-free
            #batch_1 = scano_data[scano_data.obs['Batch'] == '1']
            #batch_2 = scano_data[scano_data.obs['Batch'] == '2']                

            # codvid-19            
            #batch_1 = scano_data[scano_data.obs['Batch'] == 'F']
            #batch_2 = scano_data[scano_data.obs['Batch'] == 'M']

            #list of AnnData of batches
            #list_of_ds = [batch_1, batch_2, batch_3, batch_4, batch_5, batch_6, batch_7]
            #list_of_ds = [batch_1, batch_2, batch_3, batch_4]
            list_of_ds = [batch_1, batch_2]

            # Batch correction.
            scano_corrected_exp = scanorama.correct_scanpy(list_of_ds)
            scano_corrected_data = ad.concat(scano_corrected_exp)

            print("--- %s seconds ---" % (time.time() - start_time))

            # Write normalized data to csv file
            filename = folder + '/' + "scanorama_corrected_data.csv"
            df = scano_corrected_data.to_df().T
            df.to_csv(os.path.join(filename))
            
            # recover major used memory
            del df
            del scano_corrected_exp
            del scano_corrected_data
            del batch_1, batch_2 #, batch_3, batch_4
            #del batch_6, batch_7, batch_3, batch_4, batch_5
            del list_of_ds

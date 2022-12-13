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
    #"../GBC/lowdepth/low_avg_10_LBE/data/"
    #,"../GBC/lowdepth/low_avg_10_SBE/data/"
    #"../GBC/lowdepth/avg_4_SBE/data/"
    #"../GBC/lowdepth/avg_4_LBE/data/"
    
    "Covid19"
    #"SCLung"   
    #"FPTest/Covid19"
    #"FPTest/SCLung"   
    
    #"../GBC/7batch/7b_sbe/data"
    #"../GBC/7batch/7b_lbe/data"
    ]
    
    for path in allpath:    
        
        list_subfolders= [f.path for f in os.scandir(path) if f.is_dir()]
        
        for folder in list_subfolders:            

            ds_file = folder + '/' + "counts.txt"
            mt_file = folder + '/' + "cellinfo.txt"

            print(ds_file)
            print (mt_file)

            # loading data
            #sdata = ad.read_csv(ds_file, delimiter='\t', first_column_names=True) # for simulations delimiter = '\t'        
            sdata = ad.read_csv(ds_file, delimiter=' ', first_column_names=True)  # for Covid19 # for LUAD data delim=" " 

            sdata = sdata.T                
            sdata.layers['counts'] = sdata.X

            #pd_obs = pd.read_csv(mt_file, delimiter="\t", header=0, index_col=0)
            pd_obs = pd.read_csv(mt_file, delimiter=" ", header=0, index_col=0) # for Covid19, LUAD data delim=" "

            binfo = pd_obs['Batch']
            binfo = list(map(str, binfo))
            pd_obs['Batch'] = binfo
            pd_obs.index = sdata.obs_names
            sdata.obs = pd_obs

            sc.pp.filter_genes(sdata, min_cells=0)
            sc.pp.normalize_total(sdata, target_sum = 1e4)
            sc.pp.log1p(sdata)  
            
            # scGen     
            scgen_data = sdata.copy()
            start_time = time.time()
            scgen.SCGEN.setup_anndata(scgen_data , batch_key="Batch", labels_key="Group")        
            model = scgen.SCGEN(scgen_data, n_hidden=1024, n_latent=128, n_layers=3, dropout_rate=0.2)             

            model.train(
                max_epochs=3000,
                batch_size=32,
                early_stopping=False,
                early_stopping_patience=100,
            )        
            scgen_corrected_exp = model.batch_removal()
            print("--- %s seconds ---" % (time.time() - start_time))

            filename = folder + '/' + "scgen_corrected_data.csv"
            df = pd.DataFrame(scgen_corrected_exp.X.T, columns=scgen_corrected_exp.obs_names, index=scgen_corrected_exp.var_names)
            df.to_csv(os.path.join(filename))
            
            # recover major used memory
            del df
            del model
            del scgen_corrected_exp
            del scgen_data            



## ------------------------------------ ##
## written by Hyojin Kim in Hayat group & Rafael group
## ------------------------------------ ##


import os,sys
## ---------- NOTE --------- ## 
## pip install pandas==1.3.X
## ------------------------- ##
import pandas as pd
import numpy as np
import scanpy as sc
import anndata2ri
import glob
import scanpy as sc
import scanpy.external as sce
import scvi
import datetime
import re
import anndata as ad
import matplotlib.pyplot as plt
import seaborn as sns

from natsort import natsorted

from sklearn.ensemble import RandomForestClassifier
from sklearn.neural_network import MLPClassifier
from sklearn.datasets import make_classification
from sklearn.model_selection import train_test_split
from sklearn.experimental import enable_halving_search_cv 
from sklearn.model_selection import HalvingGridSearchCV
from sklearn.model_selection import HalvingRandomSearchCV
from sklearn.metrics import confusion_matrix
from sklearn.metrics import f1_score
from lightgbm import LGBMClassifier

import subprocess
from scipy.sparse import csr_matrix
sc.settings.set_figure_params(dpi=80, facecolor='white')

from numpy import savetxt
import harmonypy as hm

from rpy2.robjects import pandas2ri
pandas2ri.activate()

import anndata2ri
from rpy2.robjects import r
import rpy2.robjects as ro
from rpy2.robjects.conversion import localconverter
anndata2ri.activate()
r('library(Seurat)')
r('library(SeuratDisk)')
r('library(symphony)')
r('library(dplyr)')
r('library(zellkonverter)')
r('library(SingleCellExperiment)')
r('library(harmony)')




## Command line
## ------------------------------------ ##
## python ./run_by_scanpy.harmony.scvi.py /home/ 401 vst_2000 zinb batch


## DEF
## ------------------------------------ ##
def read(a):
	col1= list(map(lambda x:x.split('\n')[0].strip(), open(a, "r")))
	col2= list(map(lambda x:x.split('\t'), col1))
	return (col2)


def two_list_diff(list1, list2) :
	array1 = np.array(list1)
	array2 = np.array(list2)
	subtracted_array = np.subtract(array1, array2)
	subtracted = list(subtracted_array)
	return (subtracted)


def intersection(lst1, lst2):
	lst3 = [value for value in lst1 if value in lst2]
	return (lst3)


def make_outdir(outdir, key):
	# outdir : top dir to save outputs
	# key : folder name to save outputs
	new_outdir = outdir+key
	if not os.path.exists(new_outdir):
		os.mkdir(new_outdir)
	return (new_outdir)




## INPUT 
## ------------------------------------ ##
## input_f
## ------------------------------------ ## 
## RQ --(tab)-- file_dir --------(tab)--------- cell_type ------(tab)------ batch -----(tab)------ PC --(tab)-- outer_key
## R1 --(tab)-- /home/data/data1.h5ad --(tab)-- annotation_column --(tab)-- batch_column --(tab)-- 40 --(tab)-- data1_2
## R2 --(tab)-- /home/data/data2.h5ad --(tab)-- annotation_column --(tab)-- batch_column --(tab)-- 40 --(tab)-- data1_2
## ------------------------------------ ##
input_f = sys.argv[1]  ## "/home/code/input_f.txt"  
input_meta = read(input_f) 
input_meta_pd = pd.DataFrame(input_meta[1:], columns=input_meta[0]) ## --- DO NOT TOUCH --- ##
input_meta_pd.index = input_meta_pd.RQ.values.tolist() ## --- DO NOT TOUCH --- ##

if input_f == "./input_f.txt" :
    sorter = [ "R1", "R2", "R3", "R4" ]
elif input_f == "./input_f1.txt" :
    sorter = [ "R1", "R2", "R3"]
else :
    sorter = [ "R1", "R2" ]

input_meta_pd = input_meta_pd.reindex(index=sorter) ## --- DO NOT TOUCH --- ##

if input_f == "./input_f.txt" :
    R1_f = input_meta_pd.file_dir[0] ## --- DO NOT TOUCH --- ##
    R2_f = input_meta_pd.file_dir[1] ## --- DO NOT TOUCH --- ##
    R3_f = input_meta_pd.file_dir[2] ## --- DO NOT TOUCH --- ##
    R4_f = input_meta_pd.file_dir[3] ## --- DO NOT TOUCH --- ##
elif input_f == "./input_f1.txt" :
    R1_f = input_meta_pd.file_dir[0] ## --- DO NOT TOUCH --- ##
    R2_f = input_meta_pd.file_dir[1] ## --- DO NOT TOUCH --- ##
    R3_f = input_meta_pd.file_dir[2] ## --- DO NOT TOUCH --- ##
else: 
    R1_f = input_meta_pd.file_dir[0] ## --- DO NOT TOUCH --- ##
    R2_f = input_meta_pd.file_dir[1] ## --- DO NOT TOUCH --- ##




## STEP0: SETTING : type address  
## ------------------------------------ ##
## ./current dir [1]
## 	|_ ./code ( = code_dir )  [2]
##          |_ ./this_code.py
##      |_ ./project ( = top_outdir )
## 		|_ ./outdir_key ( = project_dir ) [4]
##			|_ ./results ( = outdir ) [5]
## ------------------------------------ ##

## --- [1] --- ##
current_dir = sys.argv[2]

## --- [2] --- ## 
code_dir = current_dir + "code/"

## --- [3] --- ## 
top_outdir = current_dir + "project/" 
if not os.path.exists(top_outdir):os.mkdir(top_outdir)

## --- [4] --- ##
outdir_key = input_meta_pd.outer_key.values.tolist()[-1]
project_dir = make_outdir(top_outdir, outdir_key)

## --- [5] --- ##
outdir = make_outdir(project_dir, "/results/")





## STEP0: SETTING : objects
## ------------------------------------- ##
if input_f == "./input_f.txt" :
    R_file = [ R1_f, R2_f, R3_f, R4_f ]
elif input_f == "./input_F1.txt" :
    R_file = [ R1_f, R2_f, R3_f ]  ## --- do not touch --- ##
else:
    R_file = [ R1_f, R2_f ]  ## --- do not touch --- ##


## ------------------------------------- ##
## STEP0: SETTING : column name 
## ------------------------------------- ##
R_batch_column =    input_meta_pd.batch.values.tolist()      ## e.g., [ "batch", "Condition", "sample" ]
R_celltype_column = input_meta_pd.cell_type.values.tolist()  ## e.g., [ "cell_type_original", "Names", "cell_type" ]



## ------------------------------------- ##
## STEP0: SETTING : PC 
## ------------------------------------- ##
nPC = 50
nPC_Q = int(input_meta_pd.PC[-1])-3
data_project = outdir_key


## ------------------------------------- ##
## STEP1: read rds or h5ad and merge Ref and Query in turn  
## ------------------------------------- ##
k          = int(sys.argv[3]) 
genetype   = sys.argv[4]
model_used = sys.argv[5]
batch_used = sys.argv[6]

R=[]
for ri in range(0,len(R_batch_column)):
        D = scvi.data.read_h5ad(R_file[ri])
        del D.var
        del D.obsm
        del D.uns
        del D.layers
        D.X = D.X.tocsr() ## ---- DO NOT REMOVE ---- ##
        D.obs["RQ"] = "R"
        D.obs["RQ_N"] = "R"+str(ri) 
        D.obs["batch"] = D.obs[R_batch_column[ri]] ##.astype("string") --> it makes error !!!
        D.obs["cell_type"] = D.obs[R_celltype_column[ri]] ##.astype("string") --> it makes error !!!
        ## ---- filter out non-annotated cells : start ---- ##
        full_annot_cell_id = D.obs["cell_type"].dropna().index.tolist()
        D = D[D.obs.index.isin(full_annot_cell_id), :]
        ## ---- filter out non-annotated cells : end ---- ##
        if ". " in D.obs["cell_type"].values.tolist()[0] :
            item = [ sub.split('. ')[1].strip() for sub in D.obs["cell_type"].values.tolist() ]
            D.obs["cell_type"] = item	
        cell_type_N_combine = D.obs["cell_type"].astype("string") + "_" +D.obs["RQ_N"].astype("string")
        D.obs["cell_type_N"] = cell_type_N_combine.tolist()
        D.obs.columns = D.obs.columns.astype(str) ## ---- DO NOT REMOVE ---- ## 
        R.append(D)



RQ_L=[]
RQ_N=[]
for sub in range(0, len(R)) :
	R[sub].obs[["2000_batch"]] = R[sub].obs[["batch"]]
	obs = R[sub].obs.columns.tolist()
	obs.remove("RQ")
	obs.remove("RQ_N")
	obs.remove("cell_type")
	obs.remove("cell_type_N")
	obs.remove("batch")
	obs.remove("2000_batch")
	for sub_obs in obs :
		del R[sub].obs[sub_obs]
	RQ_L.append(R[sub])
	RQ_N.append("R"+str(sub))



total_batch=[]
for i in range(len(RQ_L)):
    total_batch.append(len(RQ_L[i].obs.2000_batch.value_counts().index.tolist()))


## ------------------------------------ ##
## STEP2 : get FEATURES
## genetype = "2000"
## ------------------------------------ ##
M = ad.AnnData.concatenate(*RQ_L)
M.obs[["cell_type"]]=M.obs[["cell_type_N"]]
M.obs[["batch"]]=M.obs[["2000_batch"]]
del M.obs["2000_batch"]
del M.obs["cell_type_N"]


adata = M
adata.layers["counts"] = M.X
adata.raw = adata
sc.pp.highly_variable_genes(adata, layer = 'counts', n_top_genes=2000, flavor='seurat_v3')
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.scale(adata, max_value=10)
sc.tl.pca(adata, svd_solver='arpack')
sc.pl.pca_variance_ratio(adata, n_pcs=50, log=True)
plt.savefig(os.path.join(outdir, "M_pca.pdf"))
plt.close()
adata = adata.raw.to_adata()
adata.write_h5ad(filename=outdir+"M"+"."+genetype+"."+str(k)+".h5ad")


## ------------------------------------ ##
## STEP3 : get REPRESENTATION
## ------------------------------------ ##
## scanpy x harmony : https://scanpy.readthedocs.io/en/stable/generated/scanpy.external.pp.harmony_integrate.html
## ------------------------------------ ##

## ---- harmony ---- ##
## NOTE : harmony installed in scanpy doesn't work for some small data, 
## but R harmony works for such a data
## ----------------- ##
h5ad = outdir+"M"+"."+genetype
data_mat = pd.DataFrame(adata.obsm["X_pca"][:,0:nPC])
meta_dat = adata.obs
data_mat.to_csv(h5ad+".pca."+str(k)+".txt", index=False)
meta_dat.to_csv(h5ad+".meta."+str(k)+".txt", index=False)
data_mat_file = h5ad+".pca."+str(k)+".txt"
meta_dat_file = h5ad+".meta."+str(k)+".txt"

r(f'my_pca_embeddings = as.matrix(read.table("{data_mat_file}", sep=",", header = TRUE))')
r(f'meta_data = read.table("{meta_dat_file}", sep=",", header = TRUE)')
my_harmony_embeddings = r('HarmonyMatrix(data_mat=my_pca_embeddings, meta_data=meta_data, vars_use="batch", do_pca=FALSE)')

## ---- remove data_mat_file, meta_dat_file --- ##
temp_file = open(outdir+"temp2.log.txt",'w')
call_R = subprocess.call([ 'rm', data_mat_file], stdout=temp_file )
call_R = subprocess.call([ 'rm', meta_dat_file], stdout=temp_file )

adata.obsm["2000_pca"] = adata.obsm["X_pca"][:,0:nPC]
adata.obsm["2000_harmony"] = my_harmony_embeddings
del adata.obsm["X_pca"]

## ---- scvi ---- ##
adata.layers["counts"] = adata.X
sc.pp.highly_variable_genes(adata, layer = 'counts', n_top_genes=2000, flavor='seurat_v3', subset=True)
scvi.model.SCVI.setup_anndata(adata, layer="counts", batch_key=batch_used)
vae = scvi.model.SCVI(adata, n_layers=2, n_latent=nPC, gene_likelihood=model_used)
vae.train()
adata.obsm["2000_scVI"] = vae.get_latent_representation(adata)


## ------------------------------------ ##
## STEP4 : save h5ad
## ------------------------------------ ##
adata.write_h5ad(filename=outdir+"M"+"."+genetype+"."+str(k)+".subset."+batch_used+"."+model_used+".h5ad")






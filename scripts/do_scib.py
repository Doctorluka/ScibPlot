import scib
import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import os
os.chdir('/home/joe/project/CYQ_ILD_PH/')

import rpy2.robjects as ro
from rpy2.robjects.packages import importr

# 手动加载utils和stats包
utils = importr('utils')
stats = importr('stats')

adata_ds = sc.read_h5ad("data/scRNA-seq/scib/11_integration.h5ad")
sc.tl.pca(adata_ds, n_comps=50)
sc.pp.neighbors(adata_ds, use_rep="X_pca")

# get unintegrated adata
adata_unint = adata_ds.raw.to_adata()
adata_unint

# CCA
adata_cca = adata_ds.copy()
# add X_emb
adata_cca.obsm['X_emb'] = adata_cca.obsm['X_CCA'].copy() 


# Harmony
adata_harmony = adata_ds.copy()
sc.pp.neighbors(adata_harmony, use_rep="X_pca_harmony", n_pcs=30)
# add X_emb
adata_harmony.obsm['X_emb'] = adata_harmony.obsm['X_pca_harmony'].copy() 


# combat
adata_combat = sc.read_h5ad("data/scRNA-seq/scib/11_combat.h5ad")
sc.pp.neighbors(adata_combat, use_rep="X_combat")
# add X_emb
adata_combat.obsm['X_emb'] = adata_combat.obsm['X_combat'].copy() 

# scanorama
adata_scanorama = adata_ds.copy()
sc.pp.neighbors(adata_scanorama, use_rep="X_scanorama")
# add X_emb
adata_scanorama.obsm['X_emb'] = adata_scanorama.obsm['X_scanorama'].copy() 

# scVI
adata_scvi = adata_ds.copy()
sc.pp.neighbors(adata_scvi, use_rep="X_scVI")
# add X_emb
adata_scvi.obsm['X_emb'] = adata_scvi.obsm['X_scVI'].copy() 


# scANVI
adata_scanvi = adata_ds.copy()
sc.pp.neighbors(adata_scanvi, use_rep="X_scANVI")
# add X_emb
adata_scanvi.obsm['X_emb'] = adata_scanvi.obsm['X_scANVI'].copy() 

# MIRA_feature
adata_mira_feature = adata_ds.copy()
sc.pp.neighbors(adata_mira_feature, use_rep="MIRA_feature")
# add X_emb
adata_mira_feature.obsm['X_emb'] = adata_mira_feature.obsm['MIRA_feature'].copy() 

# MIRA_topic
adata_mira_topic = adata_ds.copy()
sc.pp.neighbors(adata_mira_topic, use_rep="MIRA_topic")
# add X_emb
adata_mira_topic.obsm['X_emb'] = adata_mira_topic.obsm['MIRA_topic'].copy() 

# BBKNN & Cellhint
adata_bbknn = sc.read_h5ad("data/scRNA-seq/scib/11_bbknn.h5ad")
adata_cellhint = sc.read_h5ad("data/scRNA-seq/scib/11_cellhint.h5ad")

# Do benchmark
import sys
sys.path.append("./utils/")
from my_utils import run_scib_metrics


info_dict = {
    "methods": ["harmony", "cca", "scvi","scanvi","combat","scanorama","MIRA_feature","MIRA_topic", "bbknn", "cellhint"],
    "adatas": [adata_harmony, adata_cca, adata_scvi, adata_scanvi, adata_combat, adata_scanorama, adata_mira_feature, adata_mira_topic, adata_bbknn, adata_cellhint],
    "embeds": ["X_pca_harmony", "X_pca", "X_scVI", "X_scANVI", "X_combat", "X_scanorama", "MIRA_feature", "MIRA_topic", None, None],
    "method_type": ["embed", "full", "embed", "embed", "embed", "embed", "embed", "embed", "knn", "knn"],
}

metrics_list = run_scib_metrics(adata_unint=adata_unint, 
                                batch_key="sample", 
                                label_key="annotation_level2", 
                                info_dict=info_dict, 
                                organism="human")


# Concatenate metrics results
metrics = pd.concat(
    metrics_list,
    axis="columns",
)

column_names = ["Harmony", "CCA", "scVI", "scANVI", "Combat", "Scanorama", "MIRA_feature", "MIRA_topic", "BBKNN", "CellHint"]
# Set methods as column names
metrics = metrics.set_axis(
    column_names, axis="columns"
)
metrics.to_csv("data/scRNA-seq/scib/11_scib_metrics_concat.csv")




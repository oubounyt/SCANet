import time
import anndata
import numpy as np
import scanpy as sc
import pandas as pd
from anndata import AnnData
import matplotlib.pyplot as plt
from sklearn.cluster import KMeans

       
#Dimensionality Reduction
#Calculating PCA
#Calculating UMAP

class Visualization():

   #TODO CHECK IN CASE ANY PRE-PROCESSING IS MESSING AND ASK FOR IT
   #TODO CHECK IF DATA IS LOG NORMALIZED OR NOT

    @staticmethod
    def _check(adata: AnnData):
        if '__SCANclusters__' in list(adata.obs):
            return 1,1
        try:
            adata.var['highly_variable']
            return 1,1
        except:
            return -1, ValueError("Did not find highly_variable consider running `extract_highly_variable_genes` first.")
        
    
    @classmethod
    def dimensionality_reduction_parameters(cls, adata: AnnData):
        if cls._check(adata)[0] == -1:
            raise cls._check(adata)[1]

        sc.tl.pca(adata, svd_solver='arpack', n_comps=50)
        sc.pl.pca_variance_ratio(adata, log=True)

    @classmethod
    def dimensionality_reduction(cls, adata: AnnData, method: str, n_pcs: int, n_neighbors: int):
        start = time.time()
        
        if cls._check(adata)[0] == -1:
            raise cls._check(adata)[1]
        if method not in ['ALL','UMAP','t-SNE']:
            raise ValueError('Please select on of the following methods UMAP or t-SNE or ALL for both.')

            
        
        print("Calculating dimensionality reduction :"+ method)
        sc.pp.neighbors(adata, n_pcs=n_pcs, n_neighbors=n_neighbors)
        if method == 'ALL':
            sc.tl.tsne(adata)
            sc.tl.umap(adata)
        elif method == 't-SNE':
            sc.tl.tsne(adata)
        elif method == 'UMAP':
            sc.tl.umap(adata)
                
        end = time.time()
        print("it tooks :" + str((end - start)/60)+" minutes.")
        
        return adata

    @staticmethod
    def visualization(adata: AnnData, method: str, color: str):
        #print('Avialble coloring option :\n' + str(list(adata.obs)))
        print("\n Available method are: PCA, t-SNE,  UMAP, or ALL")
        #specifying a color palette to circumvent a bug that appears for older scanpy versions
        if method == 'ALL':
            sc.pl.pca_scatter(adata, color=color, palette='tab10')
            sc.pl.tsne(adata, color=color, palette='tab10')
            sc.pl.umap(adata, color=color, palette='tab10')
        elif method == 'PCA':
            sc.pl.pca_scatter(adata, color=color, palette='tab10')
        elif method == 't-SNE':
            sc.pl.tsne(adata, color=color, palette='tab10')
        elif method == 'UMAP':
            sc.pl.umap(adata, color=color, palette='tab10')
import anndata
import numpy as np
import pandas as pd
import scanpy as sc
from anndata import AnnData
import matplotlib.pyplot as plt
from sklearn.cluster import KMeans


class Pseudobulk():
    
    @staticmethod
    def representative_cells(adata: AnnData, num_rep_cells: int, cell_anno: str):
        cell_type = pd.DataFrame(adata.obs[cell_anno]) #Check cluster size
        cluster_sizes = cell_type.groupby([cell_anno]).size()
        print(f'clusters in the data {adata.obs.groupby([cell_anno]).size()}')
        
        df_tmp1 = adata.obs.groupby([cell_anno]).size()
        df_tmp1 = df_tmp1.to_frame()
        df_tmp1 = df_tmp1.rename(columns={df_tmp1.columns[0]: 'cluster-size'})
        df_tmp1 = df_tmp1.loc[df_tmp1["cluster-size"]<num_rep_cells]
        to_del = list(df_tmp1.index)
        
        if len(to_del) != 0:
            print("\nThe following cell clusters will be deleted"\
                  "as their size if less than the specified number of"\
                  "representative cells:\n" + str(to_del))
            
        for x in to_del:
            adata = adata[adata.obs[cell_anno] != x]
            
        data = pd.DataFrame.sparse.from_spmatrix(adata.X, index=list(adata.obs[cell_anno]))
        for c in list(data.columns):
            data[c] = data[c].values.to_dense().astype(np.float64)

        gene_names = list(adata.var_names)  #TODO How to get genes in a general way
        groups = data.groupby(data.index)
        
        dfs = [] #reduce cluster sizes
        for name, group in groups:
            g = group
            arr_ = g.to_numpy()
            kmeans = KMeans(n_clusters=num_rep_cells)
            kmeans.fit(arr_)
            kmeans.predict(arr_)
            centers = kmeans.cluster_centers_
            df = pd.DataFrame(data=centers, columns=gene_names)
            df.insert (0, "clusters", name)
            dfs.append(df)
        
        new_data = pd.concat(dfs, ignore_index=True)
        new_data = new_data.set_index('clusters')

        #create new anndata object
        reduced_data = anndata.AnnData(X= new_data, obs= pd.DataFrame(list(new_data.index), 
                                       index=list(new_data.index), columns=["__SCANclusters__"]), 
                                       var= pd.DataFrame(list(new_data.columns), index=list(new_data.columns), 
                                       columns=["Genes"]))
        reduced_data.var_names_make_unique()  
        reduced_data.obs_names_make_unique()  
        # NOW the data is processed dont check for HIghly variable genes anymore
        #reduced_data.obs["_processed_"] = "TRUE"
        #print(reduced_data)
        print(f'\nSelected representative cells for each cluster {reduced_data.obs.groupby(["__SCANclusters__"]).size()}')
        return reduced_data



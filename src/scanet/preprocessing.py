import pathlib
import sys
import typing
import warnings
from pathlib import Path, PurePath
from turtle import up
from typing import Optional, Tuple, Union

import matplotlib.pyplot as plt
import numpy as np
import scanpy as sc
import scipy
import seaborn as sns
from anndata import AnnData

# styles
# verbosity: errors (0), warnings (1), info (2), hints (3)

sns.set_theme(style = "white")
sc.settings.verbosity = 3 
sc.logging.print_header()
sc.settings.set_figure_params(dpi=100, facecolor='white')

# warnings
if not sys.warnoptions:
    warnings.simplefilter("ignore")  

class Preprocessing():
    
    @staticmethod
    def _intial(adata: AnnData):
        adata.var_names_make_unique()
        adata.var['mt'] = adata.var_names.str.startswith('MT-') 
        mito_genes = adata.var_names.str.startswith('MT-')
        adata.obs['percent_mito'] = np.sum(adata[:, mito_genes].X, axis=1).A1 / np.sum(adata.X, axis=1).A1  
        sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, inplace=True)
        sc.pp.filter_cells(adata, min_genes=0)
        sc.pp.filter_cells(adata, min_counts=0)
        sc.pp.filter_genes(adata, min_cells=1)
        adata.var_names_make_unique
        adata.obs_names_make_unique
        return adata
    
    @classmethod
    def read_h5ad(cls, filename: Union[Path, str], pr_process:str = "skip"):
        if pathlib.Path(filename).suffix != ".h5ad":
            raise ValueError('Please provide a .h5ad format file')
        adata = sc.read_h5ad(filename)
        if pr_process== "skip" or '__SCANclusters__' in list(adata.obs) \
        or ['mt', 'n_cells_by_counts'] in list(adata.var):
            return sc.read_h5ad(filename)
        else:
            # initial preprocessing as it is required later
            return cls._intial(adata) 
            
    @classmethod
    def plot_filter_cells(cls, adata: AnnData, max_counts: int = 0,
                          dot: Optional[int] = 50, 
                          fig_size: Optional[Tuple[float, float]]=(8,6)):
        if max_counts == 0:
            print("You did not specify a max count values ...!") 
            max_counts = np.max(adata.obs['n_counts'])
        plt.figure(figsize=(fig_size[0], fig_size[1]))
        ax = plt.axes()
        ax.set_facecolor('grey')
        max_value = min(max_counts, np.max(adata.obs['n_counts'])+10)
        plt.axhspan(0, max_value, color='white', zorder=-1)
        plt.scatter(adata.obs['n_genes'], adata.obs['n_counts'], dot)
        plt.xlabel("#N  genes")
        plt.ylabel("#N  counts")
        plt.axhline(max_value, color='red')
        plt.show()
        return print("Used max_counts is " + str(max_counts))

    @staticmethod
    def filter_cells(adata: AnnData, max_counts: int):
        print('Total number of cells: {:d}'.format(adata.n_obs))
        print('Number of cells after min genes filter: {:d}'.format(adata.n_obs))
        # filter cells with too high counts ==> doublet
        sc.pp.filter_cells(adata, max_counts=max_counts) 
        print('Number of cells after max count filter: {:d}'.format(adata.n_obs))
        return adata

    @classmethod
    def plot_cells_by_n_genes(cls, adata: AnnData, lowerbound: int, upperbound: int, 
                              nbins: Optional[int] = 50, 
                              fig_size: Optional[Tuple[float, float]]=(12,6)):
        fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=fig_size, sharey=True)
        x = adata.obs['n_genes']
        x_lowerbound = lowerbound
        x_upperbound = upperbound
        
        nbins=nbins
        sns.distplot(x, ax=ax1, norm_hist=True, bins=nbins)
        sns.distplot(x, ax=ax2, norm_hist=True, bins=nbins)
        sns.distplot(x, ax=ax3, norm_hist=True, bins=nbins)

        ax2.set_xlim(0,x_lowerbound)
        ax3.set_xlim(x_upperbound, adata.obs['n_genes'].max())

        for ax in (ax1,ax2,ax3): 
            ax.set_xlabel('')

        ax1.title.set_text('#N genes')
        ax2.title.set_text('#N genes, lower bound')
        ax3.title.set_text('#N genes, upper bound')
        fig.text(-0.01, 0.5, 'Frequency', ha='center', va='center', rotation='vertical', size='x-large')
        fig.text(0.5, 0.0, 'Genes expressed per cell', ha='center', va='center', size='x-large')

        fig.tight_layout()
        fig.show()

    @staticmethod
    def filter_cells_by_n_genes(adata: AnnData, min_n_genes: int , max_n_genes: int):
        # filter genes not present in enough cells
        print('Total number of cells: {:d}'.format(adata.n_obs))
        adata = adata[adata.obs['n_genes'] > min_n_genes, :]
        adata = adata[adata.obs['n_genes'] < max_n_genes, :]
        print('Number of cells after cell filter: {:d}'.format(adata.n_obs))
        return adata

    @classmethod    
    def plot_filter_mitochondrial(cls, adata: AnnData, thres: int,
                                  dot: Optional[int] = 50,
                                  fig_size: Optional[Tuple[float, float]]=(12,6)):
        plt.figure(figsize=(fig_size[0], fig_size[1]))
        plt.scatter(adata.obs['total_counts'], adata.obs['pct_counts_mt'], dot)
        plt.axhline(thres, color='red')
        plt.title("Filter Mitochondrial")
        plt.xlabel("Total Counts")
        plt.ylabel("% Counts Mitchondrial")

    @classmethod     
    def filter_genes(cls, adata: AnnData, min_cells: int):
        print('Total number of genes: {:d}'.format(adata.n_vars))
        sc.pp.filter_genes(adata, min_cells=min_cells)
        print('Number of genes after filtering: {:d}'.format(adata.n_vars))
        return adata

    @classmethod     
    def filter_mitochondrial(cls, adata: AnnData, max_mt_percent: int):
        print('Total number of cells: {:d}'.format(adata.n_obs))
        adata = adata[adata.obs.pct_counts_mt < max_mt_percent, :]
        print('Number of cells after MT filter: {:d}'.format(adata.n_obs))
        return adata
    
    @staticmethod
    def log_normalize(adata: AnnData, target_sum: Optional[int] = 1e6, use_log: bool = True):
        sc.pp.normalize_total(adata, target_sum=target_sum, inplace=True)
        if use_log:
            sc.pp.log1p(adata)
        return adata

    @staticmethod
    def extract_highly_variable_genes(adata: AnnData, n_top_genes: Optional[int] = 2000,  plot: bool = True):
        sc.pp.highly_variable_genes(adata, flavor='cell_ranger', n_top_genes=n_top_genes)
        print('\n','Number of highly variable genes: {:d}'.format(np.sum(adata.var['highly_variable'])))
        if plot:
            sc.pl.highly_variable_genes(adata)
        adata = adata[:, adata.var.highly_variable]
        return adata

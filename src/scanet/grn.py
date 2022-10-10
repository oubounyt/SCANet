import ast
import glob
import multiprocessing
import os
import random
import string
import subprocess
import sys
import tempfile
from typing import Optional, Tuple, Union

import loompy as lp
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
from ctxcore.rnkdb import FeatherRankingDatabase as RankingDatabase
from dask.diagnostics import ProgressBar
from pandas import DataFrame
from pyscenic.cli.utils import load_signatures
from pyvis.network import Network

letters = string.ascii_lowercase
import warnings

warnings.filterwarnings('ignore')

class GRN():
    
    global MAIN_PATH
    MAIN_PATH = os.path.dirname(__file__)

    global PATH_TF_H
    PATH_TF_H = os.path.join(MAIN_PATH, "databases/allTFs_hg38.txt")
    global PATH_TF_M
    PATH_TF_M = os.path.join(MAIN_PATH, "databases/allTFs_mm.txt")
    
    global PATH_DB 
    PATH_DB = os.path.join(MAIN_PATH, "databases")
    

    
    
    def regulators_count(modules_df: DataFrame,specie: str, figsize: Optional[Tuple[float, float]]=(8,6)):    
        
        lst = [] 
        for x in list(modules_df["Modules"]):

            Module_use = modules_df[modules_df["Modules"]==x]
            Module_genes = list(Module_use["genes"])[0]
            Module_genes = ast.literal_eval(Module_genes)

            if specie == 'human':
                f_tfs = PATH_TF_H # human
            elif specie == 'mouse':
                f_tfs = PATH_TF_M # mouse

            my_file = open(f_tfs, "r")
            content = my_file.read()
            TFs = content.split("\n")
            my_file.close()

            num_tfs_m = [x for x in Module_genes if x in TFs]
            lst.append([x,len(Module_genes),len(num_tfs_m)])

        df = pd.DataFrame(lst, columns =['Modules', 'n_genes', "n_TFs"])

        plt.figure(figsize=figsize)
        splot=sns.barplot(x="Modules",y="n_genes",data=df, color="blue")
        for p in splot.patches:
            splot.annotate(format(p.get_height(), '.0f'), 
                        (p.get_x() + p.get_width() / 2., p.get_height()), 
                        ha = 'center', va = 'center', 
                        xytext = (0, 9), 
                        textcoords = 'offset points')
            
        splot.set(ylim=(0, max(list(df["n_genes"]))+(max(list(df["n_genes"]))*0.2)))
        plt.xticks(rotation=90)
        plt.xlabel("", size=14)
        plt.ylabel("Number of genes", size=14)
        splot.spines['right'].set_visible(False)
        splot.spines['top'].set_visible(False)

        plt.figure(figsize=figsize)
        splot=sns.barplot(x="Modules",y="n_TFs",data=df, color="orange")
        for p in splot.patches:
            splot.annotate(format(p.get_height(), '.0f'), 
                        (p.get_x() + p.get_width() / 2., p.get_height()), 
                        ha = 'center', va = 'center', 
                        xytext = (0, 9), 
                        textcoords = 'offset points')

        splot.set(ylim=(0, max(list(df["n_TFs"]))++(max(list(df["n_TFs"]))*0.2)))
        plt.xticks(rotation=90)
        plt.xlabel("Modules", size=14)
        plt.ylabel("Number of TFs", size=14)
        splot.spines['right'].set_visible(False)
        splot.spines['top'].set_visible(False)
    
    
    def save_as_loom(adata, annoatation_to_use, annoatation_name, subsmp_pct, temp_dir_):

        # paramters filter 
        if annoatation_name not in list(adata.obs[annoatation_to_use]):
            raise Exception("Annotation dosent existe please chose from :" )
        
        # save data
        path_loom = os.path.join(temp_dir_,''.join(random.choice(letters) for i in range(10))+"_adata.loom")
        print(path_loom)
        adata = adata[adata.obs[annoatation_to_use]==annoatation_name]
        cells = list(adata.obs_names)

        print(f"Subsampling using {subsmp_pct}% ...")
        n = (len(cells)*subsmp_pct)/100
        cells_ = random.sample(cells, int(n))
        adata = adata[cells_,:]
        print(f"Using {len(cells_)} cells after sub-sampling ...")

        row_attrs = {
            "Gene": np.array(adata.var_names) ,
        }
        col_attrs = {
            "CellID": np.array(adata.obs_names) ,
            "nGene": np.array( np.sum(adata.X.transpose()>0 , axis=0)).flatten() ,
            "nUMI": np.array( np.sum(adata.X.transpose() , axis=0)).flatten() ,
        }
        lp.create(path_loom, adata.X.transpose(), row_attrs, col_attrs)
        return path_loom
    
    def regulon(paramters):

        loom_file = paramters[0]
        use_tf = paramters[1]
        num_workers = paramters[2]
        specie = paramters[3]

        temp_dir = tempfile.TemporaryDirectory()

        use_tf_path = os.path.join(temp_dir.name,"tf.txt")
        textfile = open(use_tf_path, "w")
        for element in use_tf:
            textfile.write(element + "\n")
        textfile.close()

        output_file_1 = os.path.join(temp_dir.name, "grn.csv") 
        output_file_2 = os.path.join(temp_dir.name, "reg_tracks.csv")    

        script_path = os.path.join(MAIN_PATH, "infer_grn.py")         
        subprocess.run(['python', script_path,
                        str(loom_file),
                        str(use_tf_path),
                        '--method', 'grnboost2',
                        '--output', str(output_file_1),
                        '--num_workers', str(num_workers),
                        '--seed', '777'])

        # ranking databases

        f_db_glob_h = os.path.join(PATH_DB, "human/hg19-*.feather")
        f_db_glob_m = os.path.join(PATH_DB, "mouse/mm9*.feather")
        db_path___ = PATH_DB
        if specie == 'human':
            f_db_glob = f_db_glob_h # human
        elif specie == 'mouse':
            f_db_glob = f_db_glob_m # mouse
        
        
        dbs = glob.glob(f_db_glob)
        print("Ranking dbs")
        print(dbs)
        # motifs
        if specie == 'human':
            MOTIF_ANNOTATIONS_FNAME = os.path.join(db_path___, "motifs-v9-nr.hgnc-m0.001-o0.0.tbl") #human
        elif specie == 'mouse':
            MOTIF_ANNOTATIONS_FNAME = os.path.join(db_path___, "motifs-v9-nr.mgi-m0.001-o0.0.tbl") #mouse

        
        cmd_ = ["pyscenic", "ctx", output_file_1,"--annotations_fname", MOTIF_ANNOTATIONS_FNAME, 
                    "--expression_mtx_fname", loom_file, "--output", output_file_2, "--num_workers", str(num_workers)]
        cmd = []
        cmd.append(cmd_[0])
        cmd.append(cmd_[1])
        cmd.append(cmd_[2])
        for i in range(len(dbs)):
            cmd.append(dbs[i])
        cmd.append(cmd_[3])
        cmd.append(cmd_[4])
        cmd.append(cmd_[5])
        cmd.append(cmd_[6])
        cmd.append(cmd_[7])
        cmd.append(cmd_[8])
        cmd.append(cmd_[9])
        cmd.append(cmd_[10])
        print("This the command")
        print(cmd)
        subprocess.run(cmd)

        reg = []
        try:
            sig = load_signatures(output_file_2)
            for i in sig:
                for j in i.genes:
                    reg.append([i.transcription_factor, j])
        except:
            pass
        
        print(f"found {len(reg)} edges")
    
        # Clean the temporary folder
        temp_dir.cleanup()

        return reg

    def grn_inference(adata_processed, modules_df, module, groupby_, anno_name, specie, 
                      subsampling_pct, n_iteration, num_workers): 
        

        Module_use = modules_df[modules_df["Modules"]== module]
        Module_genes = list(Module_use["genes"])[0]
        Module_genes = ast.literal_eval(Module_genes)

        temp_dir = tempfile.TemporaryDirectory()
        temp_dir_ = temp_dir.name

        if specie == 'human':
            f_tfs = PATH_TF_H # human
        elif specie == 'mouse':
            f_tfs = PATH_TF_M # mouse

        my_file = open(f_tfs, "r")
        content = my_file.read()
        TFs = content.split("\n")
        my_file.close()

        m_genes_ = Module_genes
        use_tf = [g for g in m_genes_ if g in TFs]

        if len(use_tf)==0:
            raise Exception("No TFs in your module")

        print(f"This Module has {len(Module_genes)} genes.")
        print(f"{len(use_tf)} Genes in this modules are regulatores (TFs).")
        print(f"TFs are: {use_tf}")
        
        # Keep only genes in the module
        adata_processed = adata_processed[:, Module_genes]
        print(f"The anndata size after gene filtring {adata_processed.shape} ...")

        adata_processed = adata_processed[adata_processed.obs[groupby_]==anno_name]
        cells = list(adata_processed.obs_names)
        print(f"Using anootation : {anno_name} found {len(cells)} cells.")
        print(f"runing : {n_iteration} iteration ...")
        args_ = []
        for i in range(n_iteration):
            loom_file = GRN.save_as_loom(adata_processed, annoatation_to_use=groupby_, annoatation_name=anno_name,
                                         subsmp_pct=subsampling_pct, temp_dir_=temp_dir_)
            args_.append([loom_file, use_tf, num_workers, specie])

        pool = multiprocessing.Pool(processes=10)
        result_list = pool.map(GRN.regulon, args_)

        temp_dir.cleanup()
        
        print(f"Analysis done ...")
        print("\n")
        for i in range(len(result_list)):
            print(f"In run {i}: {len(result_list[i])} edges where found.")
            
        l = [j for sub in result_list for j in sub]
        result_df = pd.DataFrame (l, columns = ['TF', 'TG'])
        # counting the duplicates
        df_pivot = result_df.pivot_table(index = ['TF','TG'], aggfunc ='size')
        df_pivot = df_pivot.reset_index()
        df_pivot = df_pivot.rename(columns={0:"occurrence(pct)"})
        df_pivot["occurrence(pct)"] = df_pivot["occurrence(pct)"].apply(lambda x: (x / int(n_iteration)*100))
        return df_pivot


import json
import os
import sys
from contextlib import contextmanager
from email import message
from turtle import title
from typing import List, Optional, Tuple, Union

import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd
import rpy2
import rpy2.robjects as ro
import seaborn as sns
from anndata import AnnData
from IPython.display import Image, display
from matplotlib import pyplot as plt
from pandas import DataFrame
from pyvis.network import Network
from rpy2 import robjects
from rpy2.robjects import pandas2ri
from rpy2.robjects.conversion import localconverter

# style
sns.set_theme(style="white")


class CoExpression():

    # Co-Expression analysis some partes will be done using R packages
    # creat an output folder
    
    global MAIN_PATH
    MAIN_PATH = os.path.dirname(__file__)
    
    cwd = os.getcwd() 
    path_outs = os.path.join(cwd, "outs-scanet")
    if not os.path.exists(path_outs):
          os.makedirs(path_outs)


    @contextmanager
    def suppress_stdout():
        with open(os.devnull, "w") as devnull:
            old_stdout = sys.stdout
            sys.stdout = devnull
            try:  
                yield
            finally:
                sys.stdout = old_stdout

    @classmethod
    def call_r_function(cls, func_name):
        
        path_r_script = os.path.join(MAIN_PATH,"script.R")
        with cls.suppress_stdout():
            r = ro.r
            r['source'](path_r_script)
        
        func_name_r = ro.globalenv[func_name]
        return func_name_r


    @staticmethod
    def convertor(adata_new, convertor_r):

        data = pd.DataFrame(adata_new.X, index=list(adata_new.obs_names), columns=list(adata_new.var_names))
        annotation = adata_new.obs
        print(f'Data dimension is : {data.shape[0]} Cells and {data.shape[1]} Genes')
        with localconverter(ro.default_converter + pandas2ri.converter):
            data_r = ro.conversion.py2rpy(data)
            annotation_r = ro.conversion.py2rpy(annotation)
        
        final_exp_ = convertor_r(data_r, annotation_r)
        return final_exp_


    @classmethod
    def plot_powers(cls, adata_new: AnnData, network_type: str, cor_method: str):

        convertor_r = cls.call_r_function('convertor_r')
        final_exp_ = cls.convertor(adata_new, convertor_r)

        # Loading the function we have defined in R.
        plot_power_r = cls.call_r_function('plot_power_r')
        optimal_power_ = plot_power_r(final_exp_, path=cls.path_outs,network_type=network_type,
                                      cor_method=cor_method, w=400,h=400)
        optimal_power_ = tuple(optimal_power_)
        display(Image(filename=cls.path_outs+'/power_plot.png'))   
        print(f'The optimal power to use is : {optimal_power_[0]} or visualize the plot and change the value accordingly.')
        print(f'Tip: If power values are not exponential distributed dont use the maximum value.')
        return int(optimal_power_[0])
    
    @classmethod
    def co_expression(cls, adata_new: AnnData, network_type: str, cor_method: str, 
                      power: int, module_merging_threshold: int):
        
        convertor_r = cls.call_r_function('convertor_r')
        final_exp_ = cls.convertor(adata_new, convertor_r)
         
        co_expression_r = cls.call_r_function('co_expression_r')
        net = co_expression_r(final_exp_, power, module_merging_threshold = module_merging_threshold,
                              network_type=network_type, cor_method=cor_method, w=400,h=400)
            
        return net
    
    @classmethod
    def plot_dendrogram(cls, net: robjects, fig_size: Optional[Tuple[float, float]]=(800,400)):
        print('Plotting the dendrogram')
        w = int(fig_size[0])
        h = int(fig_size[1])
        if w < 400 or h < 400:
            w=h=400
        plot_dendrogram_r = cls.call_r_function("plot_dendrogram_r")
        plot_dendrogram_r(net, cls.path_outs, w, h)
        display(Image(filename=cls.path_outs+'/dendrogram_plot.png'))

        return 
    
    @classmethod
    def plot_eigengene_network(cls, net: robjects, fig_size: Optional[Tuple[float, float]]=(400,800)):
        print('Plotting the eigengene network')
        w = fig_size[0]
        h = fig_size[1]
        w = int(fig_size[0])
        h = int(fig_size[1])
        if w < 400 or h < 400:
            w=h=400
        plot_eigengene_network_r = cls.call_r_function("plot_eigengene_network_r")
        plot_eigengene_network_r(net, cls.path_outs, w, h)
        display(Image(filename=cls.path_outs+'/eigengene_network_plot.png'))

        return 

    @staticmethod
    def mapper(df):
        mapper_modules = []
        l = df.values.tolist()
        count = 1
        for i in l:
            mapper_modules.append([i[0],i[1],i[2],f'M{count}'])
            count += 1
        return mapper_modules

    @staticmethod
    def r_2_py(df):
        data = []
        modules_ = list(set(list(df["Modules"])))
        modules_
        for x in modules_:
            genes = list(df[df["Modules"] == x]["Genes"])
            data.append([x, len(genes), json.dumps(genes)])

        out_df = pd.DataFrame (data, columns=["Module_r","n_genes","genes"])
        return out_df

    @classmethod
    def plot_modules(cls, net: robjects, figsize: Optional[Tuple[float, float]]=(6,6)):
        
        print('Plotting the modules')
        plot_modules_r = cls.call_r_function("plot_modules_r")
        genes_and_modules_df = plot_modules_r(net, cls.path_outs, w=400, h=400)
        with localconverter(ro.default_converter + pandas2ri.converter):
            genes_and_modules_df_ = ro.conversion.rpy2py(genes_and_modules_df)

        frequency_df_ = cls.r_2_py(genes_and_modules_df_)
        m = cls.mapper(frequency_df_)
        genes_frequency_df = pd.DataFrame (m, columns=["Module_r","n_genes","genes","Modules"])
        ax = genes_frequency_df.plot(kind='bar', x="Modules",y='n_genes', figsize=figsize, 
                                     legend=False, title='Number of genes per module')
        for bar in ax.patches:
            ax.annotate(format(bar.get_height(), '.0f'),
                        (bar.get_x() + bar.get_width() / 2,
                            bar.get_height()), ha='center', va='center',
                        size=12, xytext=(0,12),
                        textcoords='offset points', rotation=90)


        ax.set_yticks([])
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.spines['left'].set_visible(False)
        plt.show()
        genes_frequency_df.to_csv(cls.path_outs+"/genes_frequency_df.csv")
        
        return genes_frequency_df
    
    @classmethod
    def modules_to_annotation_cor(cls, adata_new: AnnData, net: robjects, cor_method: str,
                                  figsize: Optional[Tuple[float, float]]=(9,6)):
        
        convertor_r = cls.call_r_function('convertor_r')
        final_exp_ = cls.convertor(adata_new, convertor_r)

        module_to_annotation_cor_r = cls.call_r_function("module_to_annotation_cor_r")
        MEtrait = module_to_annotation_cor_r(final_exp_, net, cor_method=cor_method)
        with localconverter(ro.default_converter + pandas2ri.converter):
            MEtrait_ = ro.conversion.rpy2py(MEtrait)
        
        def rename_(x):
            m = pd.read_csv(cls.path_outs+"/genes_frequency_df.csv", index_col=0)
            m = m.values.tolist()
            index = [element for element in m if element[0] == x.split("ME")[1]][0]
            return index[-1]

        MEtrait_["Modules"] = MEtrait_["ME"].apply(rename_) 
        cor_heatmap = MEtrait_.pivot("Modules", "trait", "cor")
        cor_heatmap = cor_heatmap.round(3)
        f, ax = plt.subplots(figsize=figsize)
        sns.heatmap(cor_heatmap, annot=True, fmt=".3f", linewidths=.9, ax=ax, cmap = "coolwarm")
        
        MEtrait_ = MEtrait_[["Modules", "trait", "cor", "pvalue"]]
        MEtrait_ = MEtrait_.rename(columns={"trait": "annotation"})
        return MEtrait_
    
    @classmethod
    def module_to_annotation_cor(cls, cor: DataFrame, module: str, figsize: Optional[Tuple[float, float]]=(8,8)):

        df = cor[cor["Modules"]==module]
        pd.options.mode.chained_assignment = None
        df['positive'] = df['cor'] > 0
        pd.options.mode.chained_assignment = "warn"
        ax = df.plot(kind='barh', x="annotation",y='cor', figsize=figsize, legend=False, 
                     color=df.positive.map({True: 'r', False: 'b'}), width=.5,
             title='Module to annotation correlation')

        # customize the label to include the percent
        labels = [f' {v.get_width():.2f} (p={df.iloc[i, 3]:.1e})' for i, v in enumerate(ax.containers[0])]

        # set the bar label
        ax.bar_label(ax.containers[0], labels=labels, label_type='edge', size=13)
        ax.set_xticks([-1,-0,5,0,0.5,1])
        min_ = min(df["cor"])-1
        max_ = max(df["cor"])+0.5
        ax.set_xlim([min_, max_])
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.spines['left'].set_visible(False)
        plt.tight_layout()
        f = cls.path_outs+"/plot_"+module+"_trait+annotation.pdf"
        plt.savefig(f, dpi=300)
        plt.show()

    @classmethod
    def Module_Genes_Avrage_Expr(cls, module: str, adata: AnnData, figsize: Optional[Tuple[float, float]]=(8,8)):

        adata_new = adata
        df = pd.read_csv(cls.path_outs+"/genes_frequency_df.csv", index_col=0)
        genes = list(df[df["Modules"] == module]["genes"])
        genes = json.loads(genes[0])
        genes = [x for x in genes if x in list(adata_new.var_names)]
        adata_M_genes = adata_new[:, genes]   

        res = pd.DataFrame(columns=adata_M_genes.var_names, index=adata_M_genes.obs['__SCANclusters__'].cat.categories)                                                                                                 
        for clust in adata_M_genes.obs.__SCANclusters__.cat.categories: 
            res.loc[clust] = adata_M_genes[adata_M_genes.obs['__SCANclusters__'].isin([clust]),:].X.mean(0) 
        
        res = res.T

        plt.figure(figsize=figsize)
        ax = sns.violinplot(data=res)
        ax.set_xticklabels(ax.get_xticklabels(),rotation = 90)
        ax.set(xlabel='Modules', ylabel='Avrage gene expression')
        plt.show()
    
    @classmethod
    def plot_module_membership(cls, net: robjects, adata: AnnData, figsize: Optional[Tuple[float, float]]=(8,8)):
        adata_new = adata
        annotation = adata_new.obs
        def Convert_to_dict(lst):
            res_dct = {lst[i]: lst[i + 1] for i in range(0, len(lst), 2)}
            return res_dct

        print('Plotting Module Membership')
        plot_module_membership_r = cls.call_r_function("plot_module_membership_r")
        MEs = plot_module_membership_r(net)
        with localconverter(ro.default_converter + pandas2ri.converter):
            MEs_ = ro.conversion.rpy2py(MEs)

        df2 = pd.read_csv(cls.path_outs+"/genes_frequency_df.csv", index_col=0)
        df2 = df2.drop(columns="genes")
        names_ = df2.values.tolist()
        
        names__ = []
        for x in names_:
            names__.append(str("ME")+x[0])
            names__.append(x[2])

        MEs_ = MEs_.rename(columns=Convert_to_dict(names__))
        MEs_ = annotation.join(MEs_)
        MEs_.boxplot(by='__SCANclusters__', rot=90, grid=False, fontsize=15, figsize=figsize)
        title_boxplot = 'awesome title'
        plt.title( title_boxplot )
        plt.suptitle('Module Membership')
        plt.show()

        return MEs_

    @staticmethod
    def Module_Activity(module: str, df: DataFrame, plot_type: str, figsize: Optional[Tuple[float, float]]=(8,8)):
        if plot_type== "box":
            df_ = df[[module,'__SCANclusters__']]
            boxplot = df_.boxplot(by='__SCANclusters__', rot=90, grid=False, fontsize=12, figsize=figsize, return_type=None)
        elif plot_type == "violin":
            df_ = df[[module,'__SCANclusters__']]
            plt.figure(figsize=figsize)
            ax = sns.violinplot(x="__SCANclusters__", y=module, data=df_)
            ax.set_xticklabels(ax.get_xticklabels(), rotation=90)
    
    
    @classmethod
    def hub_genes(cls, adata: AnnData, net: robjects, figsize: Optional[Tuple[float, float]]=(8,8)):
        
        convertor_r = cls.call_r_function('convertor_r')
        final_exp_ = cls.convertor(adata, convertor_r)

        # Loading the function we have defined in R.
        hub_genes_r = cls.call_r_function('hub_genes_r')
        hubs_ = hub_genes_r(final_exp_, net)
        with localconverter(ro.default_converter + pandas2ri.converter):
            hubs = ro.conversion.rpy2py(hubs_)

        hub_genes = hubs.rename(columns={"Module": "Module_r"})
        df2 = pd.read_csv(cls.path_outs+"/genes_frequency_df.csv", index_col=0)
        df2 = df2.drop(columns="genes")
        df2 = df2.drop(columns="n_genes")
        hub_genes = hub_genes.merge(df2, how='inner', on='Module_r')
        hub_genes = hub_genes.drop(columns="Module_r")

        ax = hub_genes["Modules"].value_counts().plot(kind='bar', title = "Hub Gene count", figsize=figsize)
        for bar in ax.patches:
            ax.annotate(format(bar.get_height(), '.0f'),
                        (bar.get_x() + bar.get_width() / 2,
                            bar.get_height()), ha='center', va='center',
                        size=12, xytext=(0, 8),
                        textcoords='offset points')


        ax.set_yticks([])
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.spines['left'].set_visible(False)
        plt.show()


        return hub_genes

    @classmethod
    def module_to_network(cls, net: robjects, module: str, co_cutoff: float):

        module_to_network_r = cls.call_r_function('module_to_network_r')

        df2 = pd.read_csv(cls.path_outs+"/genes_frequency_df.csv", index_col=0)
        df2 = df2.drop(columns="genes")
        names_ = df2.values.tolist()
        module_ = [x[0] for x in names_ if x[2]==module]

        edges = module_to_network_r(net, module=module_[0], co_cutoff=co_cutoff)

        with localconverter(ro.default_converter + pandas2ri.converter):
            edges_filtered = ro.conversion.rpy2py(edges)
        
        print(f'Using a correlation cut-off of {co_cutoff} found {edges_filtered.shape[0]} edges')
        lis_ = list(edges_filtered["Gene1"])+list(edges_filtered["Gene2"])
        lis_ = list(set(lis_))
        print(f'There are {len(lis_)} unique genes')
        edges_filtered["Module"] = module

        return edges_filtered

    @staticmethod
    def plot_network(network, hubgenes):
        hubgenes = list(hub_genes_df[hub_genes_df['Modules'] == list(network_["Module"])[0]]['Gene'])
        print(hubgenes)
        edgelist_ = network.rename(columns={"Gene1": "source", "Gene2": "target", "Weight": "weight"})
        edgelist = edgelist_[["source","target","weight"]]
        graph = nx.from_pandas_edgelist(edgelist)
        net = Network(notebook=True)
        net.from_nx(graph)
        hubgenes = hubgenes[hubgenes["Modules"]==list(edgelist_["Module"])[0]]
        hubgenes = list(hubgenes["Gene"])
        for gene in hubgenes:
            print(gene)
            net.get_node(gene)["color"] = 'red'
        return net
    
    @classmethod
    def network_statistics(cls, network):
        network_statistics_r = cls.call_r_function('network_statistics_r')
        with localconverter(ro.default_converter + pandas2ri.converter):
            network_r = ro.conversion.py2rpy(network)
        network_statistics_r = network_statistics_r(network_r)

        return dict(zip(network_statistics_r.names, map(list,list(network_statistics_r))))
    
    @classmethod
    def module_consensus_powers(cls, adata, samples, nPerm):

        clusters_ = list(set(list(adata.obs["__clusters__"])))
        l_ = list(adata.obs["__clusters__"])
        clusters_ = sorted(set(l_), key=l_.index)
        samples_ = [clusters_.index(x)+1 for x in samples]
        print(f"Chosen samples are {samples_}")
        
        module_consensus_powers_r = cls.call_r_function('module_consensus_powers_r')
        convertor_r = cls.call_r_function('convertor_r')
        final_exp_ = cls.convertor(adata, convertor_r)

        samples__ = rpy2.robjects.IntVector(samples_)
        powers_ = module_consensus_powers_r(final_exp_, samples__, nPerm)
        
        display(Image(filename=cls.path_outs+'/power_plot_consensus_modules.png')) 
        print("Selected powers are:")
        print(powers_)

        return (int(powers_[0]), int(powers_[1]))
    
    @classmethod
    def module_consensus(cls, adata, samples, nPerm, powers):

        clusters_ = list(set(list(adata.obs["__clusters__"])))
        l_ = list(adata.obs["__clusters__"])
        clusters_ = sorted(set(l_), key=l_.index)
        samples_ = [clusters_.index(x)+1 for x in samples]
        print(f"Chosen samples are {samples_}")
        
        module_consensus_r = cls.call_r_function('module_consensus_r')
        convertor_r = cls.call_r_function('convertor_r')
        final_exp_ = cls.convertor(adata, convertor_r)

        samples__ = rpy2.robjects.IntVector(samples_)
        powers = rpy2.robjects.IntVector(powers)
        consensus_trait = module_consensus_r(final_exp_, samples__, nPerm, powers)

        with localconverter(ro.default_converter + pandas2ri.converter):
            consensus_trait_ = ro.conversion.rpy2py(consensus_trait)
        
        return consensus_trait_

    @classmethod
    def module_conservation(cls, adata, samples, nPerm, powers):

        clusters_ = list(set(list(adata.obs["__clusters__"])))
        l_ = list(adata.obs["__clusters__"])
        clusters_ = sorted(set(l_), key=l_.index)
        samples_ = [clusters_.index(x)+1 for x in samples]
        print(f"Chosen samples are {samples_}")
        module_conservation_r = cls.call_r_function('module_conservation_r')
        convertor_r = cls.call_r_function('convertor_r')
        final_exp_ = cls.convertor(adata, convertor_r)

        samples__ = rpy2.robjects.IntVector(samples_)

        display(Image(filename=cls.path_outs+'/power_plot_conservation_1.png')) 
        print("seconde one")
        display(Image(filename=cls.path_outs+'/power_plot_conservation_2.png')) 
        
        powers = rpy2.robjects.IntVector(powers)
        outputs = module_conservation_r(final_exp_, samples__, nPerm, powers)
        
        outputs_ = list(outputs)[0]
        message_ = dict(zip(outputs_.names, map(list,list(outputs_))))
        message_ = message_["message"]
        print(message_)
        print("here")
        models_ = message_[0].split(":")[1]
        models_ = models_.split("\n")[1]
        models_ = models_.split(",")
        models_ = [x.strip() for x in models_]

        # plotts etc ...
        module_conservation_plot_r = cls.call_r_function('module_conservation_plot_r')
        genes_and_modules_df = module_conservation_plot_r(outputs)
        with localconverter(ro.default_converter + pandas2ri.converter):
            genes_and_modules_df_ = ro.conversion.rpy2py(genes_and_modules_df)
    
        genes_and_modules_df_ = genes_and_modules_df_[genes_and_modules_df_['Modules'].isin(models_)]
        genes_and_modules_df_ = cls.r_2_py(genes_and_modules_df_)
        m = cls.mapper(genes_and_modules_df_)
        genes_and_modules_df_ = pd.DataFrame (m, columns=["Module_r","n_genes","genes","Modules"])
        genes_and_modules_df_.plot(kind="bar",y="n_genes", x="Modules", title="Number of genes in MODULES",
                                   figsize=(12,6))

        return genes_and_modules_df_, outputs

    @classmethod
    def consrv_module_to_network(cls, data, name, color, co_cutoff):

        consrv_module_to_network_r = cls.call_r_function('consrv_module_to_network_r')

        edges = consrv_module_to_network_r(data, module=color, co_cutoff=co_cutoff)
        display(Image(filename=cls.path_outs+'/consrv__scale_free_topology.png'))

        with localconverter(ro.default_converter + pandas2ri.converter):
            edges_filtered = ro.conversion.rpy2py(edges)
        
        print(f'Using a correlation cut-off of {co_cutoff} found {edges_filtered.shape[0]} edges')
        lis_ = list(edges_filtered["Gene1"])+list(edges_filtered["Gene2"])
        lis_ = list(set(lis_))
        print(f'There are {len(lis_)} unique genes')
        print(f'There are x hub genes genes')
        edges_filtered["Module"] = name

        return edges_filtered


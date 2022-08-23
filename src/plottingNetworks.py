from pandas import DataFrame
from turtle import shape
import pandas as pd
import networkx as nx
import requests
import json
import os

from typing import Optional, Union, Tuple
from typing import List
from pyvis.network import Network
from drug_interactions import DrugInteractions

import importlib
import drug_interactions
importlib.reload(drug_interactions)
DrugInteractions = drug_interactions.DrugInteractions

class Plot:

    cwd = os.getcwd()
    path_outs = os.path.join(cwd, "outs_scan")
    if not os.path.exists(path_outs):
        os.makedirs(path_outs)
        print(path_outs)
    
    #TODO: Talk about do we want new genes with no drug targets
    @classmethod
    def plot_gcn(cls, network_: DataFrame, hub_genes_df: DataFrame, name: str, 
                 drug_interaction: str, smooth_edges: bool = True):
        hubgenes = list(hub_genes_df[hub_genes_df['Modules'] == list(network_["Module"])[0]]['Gene'])
        print("Hub genes: " + str(hubgenes))
        graph = nx.from_pandas_edgelist(network_, "Gene1","Gene2")
        graph.remove_edges_from(nx.selfloop_edges(graph))
        network = Network(height='500px', width='800px', directed=False, notebook=True)
        network.from_nx(graph)
        
        net_genes = [x["id"] for x in network.nodes]
        for gene in hubgenes:
            if gene in net_genes:
                network.get_node(gene)["color"] = 'orangered'

        
        if drug_interaction == 'direct':
            drug_interactions, _ = DrugInteractions.get_drug_interactions(net_genes, 'degree', only_direct=True)
            for gene in drug_interactions:
                for drug in drug_interactions[gene]:
                    network.add_node(drug, shape='star', color='darkgreen')
                    network.add_edge(gene,drug)

        elif drug_interaction == 'all':
            drug_interactions, gene_edges_dir = DrugInteractions.get_drug_interactions(net_genes, 'degree', only_direct=False)
            #filter out all genes already in the network
            new_genes = [new_gene for new_gene in drug_interactions if new_gene not in net_genes]
            #new_genes = [edge[0] for edge in new_genes if new_gene not in network.get_nodes()]
            #remove directedness of gene edges
            gene_edges_undir = []
            [gene_edges_undir.append((source,target)) for (source, target) 
             in gene_edges_dir if (target,source) not in gene_edges_undir]
            network.add_nodes(new_genes, color=['green' for i in new_genes], size=[10 for i in new_genes])
            for (source, target) in gene_edges_undir:
                network.add_edge(source, target)
            for gene in drug_interactions:
                for drug in drug_interactions[gene]:
                    network.add_node(drug, shape='triangle', color='pink')
                    network.add_edge(gene,drug)


        for edge in network.edges:
            edge['physics'] = False
        if smooth_edges:
            network.set_edge_smooth('cubicBezier')
        
        #path__ = cls.path_outs + f"/{name}.html"
        #print(name)
        return  network.show(f"outs_scan/{name}.html")
    # by @Mhaned use the still and more ...

    def plot_grn(df: DataFrame, occurrence_pct: float, regulon: Optional[str]="all",
                 name: Optional[str]="GRN_network",
                 layout: Optional[str]="hierarchical", drug_interaction: Optional[str]=None):
        
        shap1 = df.shape[0]
        df = df[df["occurrence(pct)"]>=occurrence_pct]
        shap2 = df.shape[0]
        print(f"\n Out of {shap1} edges, {shap2} edges satisfying occurrence threshold {occurrence_pct}% where kept \n")
        TFs = list(set(list(df["TF"])))
        print(TFs)
        print(df.shape)

        if regulon != "all":
            if regulon in TFs:
                print(f"plottin the {regulon} network ...")
                df = df[df["TF"]==regulon]
            else:
                raise Exception(f"{regulon} is not a regulator please chose from: {TFs}") 

        graph = nx.from_pandas_edgelist(df,"TF","TG")
        graph.remove_edges_from(nx.selfloop_edges(graph))
        print(layout)
        network = Network(height='500px', width='800px', directed=True, notebook=True)
        network.from_nx(graph)

        net_genes = [x["id"] for x in network.nodes]
        if drug_interaction == 'direct':
            drug_interactions, _ = DrugInteractions.get_drug_interactions(net_genes, 'degree', only_direct=True)
            for gene in drug_interactions:
                for drug in drug_interactions[gene]:
                    network.add_node(drug, shape='star', color='darkgreen')
                    network.add_edge(drug, gene)
        elif drug_interaction == 'all':
            drug_interactions, gene_edges_dir = DrugInteractions.get_drug_interactions(net_genes, 'degree', only_direct=False)
            #filter out all genes already in the network
            new_genes = [new_gene for new_gene in drug_interactions if new_gene not in net_genes]
            #new_genes = [edge[0] for edge in new_genes if new_gene not in network.get_nodes()]
            #remove directedness of gene edges
            network.add_nodes(new_genes, color=['green' for i in new_genes], size=[10 for i in new_genes])
            for (source, target) in gene_edges_dir:
                network.add_edge(source, target)
            for gene in drug_interactions:
                for drug in drug_interactions[gene]:
                    network.add_node(drug, shape='star', color='darkgreen')
                    network.add_edge(gene,drug)

        for node in network.nodes:
            if node["label"] in TFs:
                node['shape'] = "triangleDown"
                node['color'] = "orangered"
                node['size'] = 16
            else:
                node['size'] = 14
        for edge in network.edges:
            edge['color'] = "black"
            
        
        network.set_edge_smooth('cubicBezier')
        return network.show(f"outs_scan/{name}.html")




    
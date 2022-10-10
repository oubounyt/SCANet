import importlib
import json
import os
from turtle import shape
from typing import List, Optional, Tuple, Union

import networkx as nx
import pandas as pd
import requests
from pandas import DataFrame
from pyvis.network import Network

from .drug_interactions import DrugInteractions
#DrugInteractions = drug_interactions.DrugInteractions

import warnings

warnings.filterwarnings("ignore")

class Plot:

    cwd = os.getcwd()
    path_outs = os.path.join(cwd, "outs_scan")
    if not os.path.exists(path_outs):
        os.makedirs(path_outs)
        print(path_outs)
    
    @classmethod
    def plot_gcn(cls, network_: DataFrame, hub_genes_df: DataFrame, name: str, drug_interaction: Optional[bool] = False,
                 algorithm: Optional[str] = "degree", smooth_edges: bool = True):
        
        # Make all genes uppercase
        network_ = network_.applymap(lambda s: s.upper() if type(s) == str else s)
        hub_genes_df = hub_genes_df.applymap(lambda s: s.upper() if type(s) == str else s)
        
        # print hub genes
        hubgenes = list(hub_genes_df[hub_genes_df['Modules'] == list(network_["Module"])[0]]['Gene'])
        print("Hub genes: " + str(hubgenes))
        
        # plot the gcn network 
        graph = nx.from_pandas_edgelist(network_, "Gene1","Gene2")
        graph.remove_edges_from(nx.selfloop_edges(graph))
        network = Network(height='500px', width='800px', directed=False, notebook=True)
        network.from_nx(graph)
        
        net_genes = [x["id"] for x in network.nodes]
        for gene in hubgenes:
            if gene in net_genes:
                network.get_node(gene)["color"] = 'orangered'
        
        # drug interaction
        if drug_interaction:
            drug_interactions = DrugInteractions.get_drug_interactions(net_genes, algorithm)
            print(drug_interactions)
            for inter in drug_interactions:
                network.add_node(inter[0], shape='star', color='darkgreen')
                network.add_edge(inter[0],inter[1])


        for edge in network.edges:
            edge['physics'] = False
        if smooth_edges:
            network.set_edge_smooth('cubicBezier')
        
        return  network.show(f"outs_scan/{name}.html")

    def plot_grn(df: DataFrame, occurrence_pct: float, regulon: Optional[str]="all",
                 name: Optional[str]="GRN_network", drug_interaction: Optional[bool] = False,
                 algorithm: Optional[str] = "degree", layout: Optional[str]="hierarchical"):
        
        # Make all genes uppercase
        df = df.applymap(lambda s: s.upper() if type(s) == str else s)
        
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
        # drug interaction
        if drug_interaction:
            drug_interactions = DrugInteractions.get_drug_interactions(net_genes, algorithm)
            for inter in drug_interactions:
                network.add_node(inter[0], shape='star', color='darkgreen')
                network.add_edge(inter[0],inter[1])

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




    
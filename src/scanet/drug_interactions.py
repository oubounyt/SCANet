import warnings
from typing import List, Tuple

import drugstone
from drugstone import new_task

warnings.filterwarnings("ignore")

class DrugInteractions:
    
    @classmethod
    def get_drug_interactions(cls, genes:List[str], algorithm:str):
        
        parameters = {
            "target": "drug",
            "algorithm": algorithm
        }


        task = new_task(genes, parameters)
        gene_edges_dir = []
        r = task.get_result()
        if r.get_genes() == {}:
            also_broken = cls.find_broken(genes, algorithm=algorithm)
            genes = [gene for gene in genes if gene not in also_broken]
            task = new_task(genes, parameters)
            r = task.get_result()

        drugs_ = r.get_drugs()
        genes_ = r.get_genes()

        gene_edges_list = []
        for d in drugs_.keys():
            for g in drugs_[d]['hasEdgesTo']:
                if g in genes:
                    gene_edges_list.append([d,g])

        return gene_edges_list
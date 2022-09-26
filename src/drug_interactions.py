from typing import List, Tuple
import drugstone
from drugstone import new_task

class DrugInteractions:

    @classmethod
    def get_drug_interactions(cls, genes:List[str], algorithm:str, only_direct=True) -> Tuple[dict, list]:
        parameters = {
            "target": "drug",
            "algorithm": algorithm
        }
        # TODO: Find more broken genes
        broken_genes = ['APOLD1', 'KCNK17', 'ZNF155','ZNF619', 'ZNF235', 'ZNF468', 'CENPL', 'FAM76A', 'APOBEC3H']
        genes = [gene for gene in genes if gene not in broken_genes]
        
        task = new_task(genes, parameters)
        gene_edges_dir = []
        r = task.get_result()
        if r.get_genes() == {}:
            also_broken = cls.find_broken(genes, algorithm=algorithm)
            print(also_broken)
            genes = [gene for gene in genes if gene not in also_broken]
            task = new_task(genes, parameters)
            r = task.get_result()

        drugs = r.get_drugs()
        gene_drugs = r.get_genes()
        drug_edges = {}
        gene_edges_undir = []
        if only_direct:
            gene_drugs = {gene:r.get_genes()[gene] for gene in r.get_genes().keys() if gene in genes}
            print(gene_drugs)
        else:
            
            for source in gene_drugs.keys():
                gene_edges_dir += [(source, target) for target in gene_drugs[source]['has_edges_to'] if 
                                   target in gene_drugs.keys()]

        for gene in gene_drugs.keys():
            drug_edges[gene] = {drug:drugs[drug]["score"] for drug in gene_drugs[gene]['has_edges_to'] if drug in drugs.keys()}
        

        return drug_edges, gene_edges_dir

    @staticmethod
    def find_broken(genes:List[str], algorithm:str) -> List[str]:
        parameters = {
                "target": "drug",
                "algorithm": algorithm
            }
        broken_genes= []
        for i in range(len(genes)):
            task = new_task(["CFTR", genes[i]], parameters)
            r = task.get_result()
            if r.get_genes() == {}:
                broken_genes += [genes[i]]
        return broken_genes

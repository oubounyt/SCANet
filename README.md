[![PyPI](https://img.shields.io/pypi/v/scanpy?logo=PyPI)](https://test.pypi.org/project/scanet/)




# SCANet: A workflow for single-cell co-expression based analysis 

<p align="center">
  <img src="https://github.com/oubounyt/SCAn/blob/main/logo_.png" width="400"/>
</p>

`SCANet`, a new python package that incorporates the inference of gene co-expression networks from single-cell gene expression data and a complete analysis of the identified modules through trait and cell type associations, hub genes detection, deciphering of co-regulatory signals in co-expression, and drug-gene interactions identification. This will likely accelerate network analysis pipelines and advance systems biology research.


## Installation

To install `SCANet`.
```console
pip install scanet
```
Finally, it should be possible to import drugstone to your python script.
````python
import scanet
````
You can use 
```python
import scanet as sn
```

`SCANet` officially supports Python 3.6+.

## Usage

A full detailed example of `SCANet` : Refer to this [Jupyter Notebook](https://github.com/oubounyt/SCAn/blob/main/docs/full-example.ipynb).


In this [Jupyter Notebook](https://github.com/oubounyt/SCAn/blob/main/docs/full-example.ipynb), we explored all aspects of gene coexpression networks (GCNs) using `SCANet` through a full analysis of the 3k PBMCs from 10x Genomics.

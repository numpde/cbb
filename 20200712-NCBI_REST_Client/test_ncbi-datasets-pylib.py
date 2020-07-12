# RA, 2020-07-12

"""
Experiment with the NCBI datasets client ncbi-datasets-pylib
https://github.com/ncbi/datasets/tree/master/client_docs/python

Examples:
https://github.com/ncbi/datasets/tree/master/examples/jupyter/ncbi-datasets-pylib

Installation:
    pip install ncbi-datasets-pylib
"""

import json
from ncbi import datasets

api = datasets.AssemblyDatasetDescriptorsApi()

print(
    json.dumps(
        api.assembly_descriptors_by_organism("Arabidopsis").to_dict(),
        indent=2
    )
)

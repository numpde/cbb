# RA, 2020-07-11

# Note:
# For "NCBI Datasets", use their client (cf. 20200712-NCBI_REST_Client)

"""
Assignment from
https://www.youtube.com/watch?v=s9E4Njfdib4
https://github.com/linsalrob/ComputationalGenomicsManual/tree/master/Assignments/NCBIEDirectAssignment
https://linsalrob.github.io/ComputationalGenomicsManual/Assignments/NCBIEDirectAssignment/

The file genera.txt is a plain text file that you can read on the command line using less or more.
The file has 69 different organisms, listed as genus and species, with one entry per line. Th
 assignment is to familiarize yourself with NCBIs EDirect tool and to use some advanced bash scripting.

In the first part of the assignment, you must identify which organism has the most number of genomes in the assembly database.
You should use on of the edirect scripts provided on the AWS image to complete that task as shown in the manual.

In the second part of the assignment, you must calculate the AVERAGE (mean) genome size of the genomes associated with “Prevotella buccalis”.

Hint: There is a command called countfasta.py that will likely help you with this step!
"""

import pandas as pd

from tcga.utils import download
from tcga.strings import lines
from urllib.parse import urlencode, quote

download = download.to(rel_path="cache/download")


class param:
    genera = list(lines(download(
        "https://raw.githubusercontent.com/linsalrob/ComputationalGenomicsManual/master/"
        "Assignments/NCBIEDirectAssignment/genera.txt"
    ).now.text))

    class urls:
        descriptors = "https://api.ncbi.nlm.nih.gov/datasets/v1alpha/assembly_descriptors/organism/"


for genus in param.genera:
    data = download(
        param.urls.descriptors + quote(genus) + "?"
        +
        urlencode({'returned_content': "COMPLETE", 'tax_exact_match': False})
    ).now

    if data.json:
        df = pd.DataFrame(data.json['datasets'])
        df = df.sort_values('display_name')
        print(F"{genus}, estimated genome size:", list(df.estimated_size.astype(int)))

#script to extract ENCODE accessions for a given TF

import pandas as pd
import encode_utils as eu
from encode_utils.connection import Connection

conn = Connection("prod")

#input here is a list of TFs to scrape the ENCODE portal
tf = pd.read_csv("tflist.txt", sep = "\t", header=None)
d = {} 

for i in tf[0]:

    tf_name = i.split("(")[0]
    
    query = [("assay_title","TF ChIP-seq"), ("type","Experiment"),("target.label",tf_name),('replicates.library.biosample.donor.organism.scientific_name','Homo sapiens')]

    try:
        reff = pd.json_normalize(conn.search(query))
        d[tf_name] = list(reff["accession"])
    except:
        continue

out = pd.DataFrame.from_dict({0:d},'columns')

out.to_csv("tf-encode-accessions-human.tsv", sep = "\t")
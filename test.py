import getopt
import os
import random
import sys

import pandas as pd
from py2neo import Node, Relationship
from py2neo.database import Transaction
from tqdm import tqdm

from connection import populate_db
from constants import ENCODING, DATA_DIR

# def main():
#     tx = populate_db("testReln")
#     #tx.deleteall('proxitest')
#     protacdb_df = pd.read_csv(
#         os.path.join(DATA_DIR, "protacDB.csv"),
#         usecols=[
#             "Compound ID", "Uniprot", "Target", "E3 ligase", "Name", "Smiles",
#             "DC50 (nM)", "Dmax (%)", "Assay (DC50/Dmax)", "Percent degradation (%)",
#             "Assay (Percent degradation)", "Article DOI", "Molecular Weight",
#             "Exact Mass", "logP", "logS", "Heavy Atom Count", "Ring Count",
#             "Hydrogen Bond Acceptor Count", "Hydrogen Bond Donor Count",
#             "Rotatable Bond Count", "Topological Polar Surface Area", "Molecular Formula",
#             "InChI", "InChI Key"
#         ],
#         dtype=str,
#         encoding=ENCODING
#     )
#     # pd.set_option("display.max_columns", None)
#     nodes_dict = add_nodes(tx=tx, protacdb_df=protacdb_df)
#     tx.commit()
#     #return(protacdb_df)

c = ("a","b")
print(c)
print(c[1])

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

test = pd.read_csv(
        os.path.join(DATA_DIR, "protacDB.csv"),
        usecols=[
            "Compound ID", "Uniprot", "Target", "E3 ligase", "Name", "Smiles",
            "DC50 (nM)", "Dmax (%)", "Assay (DC50/Dmax)", "Percent degradation (%)",
            "Assay (Percent degradation)", "Article DOI", "Molecular Weight",
            "Exact Mass", "logP", "logS", "Heavy Atom Count", "Ring Count",
            "Hydrogen Bond Acceptor Count", "Hydrogen Bond Donor Count",
            "Rotatable Bond Count", "Topological Polar Surface Area", "Molecular Formula",
            "InChI", "InChI Key"
        ],
        dtype=str,
        encoding=ENCODING
    )

get2nodes = test.groupby(['Target', 'E3 ligase','Uniprot']).size().reset_index().rename(columns={0: 'count'})

print(get2nodes.head(10))
print(len(get2nodes))


print(list(get2nodes.columns))


def makeNodes(tx, xyz):
    nd_dc = {
        "prtn": {},
        "e3": {}
    }

    for target,y in tqdm(get2nodes[["Target","Uniprot"]].values, total=get2nodes.shape[0]):

        if target in nd_dc["prtn"]:
            continue

        nd_dc["prtn"][target] = Node("prtn",**{"Name":target,"uniprot":f"https://www.uniprot.org/uniprot/{y}"})
        tx.create(nd_dc["prtn"][target])


    for e3 in tqdm(get2nodes["E3 ligase"], total=get2nodes.shape[0]):

        if e3 in nd_dc["e3"]:
            continue

        nd_dc["e3"][e3] = Node("e3",**{"Name":e3})
        tx.create(nd_dc["e3"][e3])
    #tx.commit()
    return(nd_dc)

t = populate_db("xyz")

xx = makeNodes(t,xyz=test)

def makeRel(tx,abc,nd):

    for row in tqdm(abc.values):
        (target,e3,y,count) = row
        rel = Relationship(nd["prtn"][target], "interacts",nd["e3"][e3],** {"Count":count,"smiles":"to be added"})
        tx.create(rel)

makeRel(t,get2nodes,xx)
t.commit()
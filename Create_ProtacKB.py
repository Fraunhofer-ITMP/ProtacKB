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


def add_nodes(tx: Transaction, protacdb_df: pd.DataFrame):
    node_dict = {
        'Protein': {},
        'E3': {}
    }

    for target, smiles, doi in tqdm(protacdb_df[["Target", "Smiles", "Article DOI"]].values, total=protacdb_df.shape[0]):

        if target in node_dict["Protein"]:
            continue

        protein_prop = {}

        if pd.notna(target):
            protein_prop["target name"] = target

        if pd.notna(smiles):
            protein_prop["smiles"] = smiles

        if pd.notna(doi):
            protein_prop["article reference"] = doi

        node_dict["Protein"][target] = Node("Protein", **protein_prop)
        tx.create(node_dict["Protein"][target])


    for ligase, mf, inchi in tqdm(protacdb_df[["E3 ligase", "Molecular Formula", "InChI"]].values, total=protacdb_df.shape[0]):

        if ligase in node_dict["E3"]:
            continue

        ligase_prop = {}

        if pd.notna(ligase):
            ligase_prop["E3 ligase"] = ligase

        if pd.notna(mf):
            ligase_prop["mol formula"] = mf

        if pd.notna(inchi):
            ligase_prop["InChi"] = inchi

        node_dict["E3"][ligase] = Node("E3", **ligase_prop)

        tx.create(node_dict["E3"][ligase])
    #print(len(node_dict['Protein']))
    #print(len(node_dict['E3']))

    return node_dict


def main():
    tx = populate_db("proxitest")
    #tx.deleteall('proxitest')
    protacdb_df = pd.read_csv(
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
    # pd.set_option("display.max_columns", None)
    nodes_dict = add_nodes(tx=tx, protacdb_df=protacdb_df)
    tx.commit()


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


#get2nodes = test.groupby(['Target', 'E3 ligase']).size().reset_index().rename(columns={0: 'count'})
#print(get2nodes.head(5))

#getComb(protacdb_df)

#getComb(protacdb_df)

def add_rel(
    tx: Transaction,
    #df: pd.DataFrame,
    get2nodes: pd.DataFrame,
    node_mapping_dict: dict
):


if __name__ == '__main__':
    #main()
    #getComb(protacdb_df)
    get2nodes = test.groupby(['Target', 'E3 ligase']).size().reset_index().rename(columns={0: 'count'})
    #print(get2nodes.head(5))
    print(list(get2nodes.columns))
    print(get2nodes['E3 ligase'][0:5])
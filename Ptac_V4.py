import pandas as pd

#import gspread
#from oauth2client.service_account import ServiceAccountCredentials

import getopt
import os
import io
from bioservices import UniProt
import uniprot
import random
import sys
from py2neo import Graph
from constants import FRAUNHOFER_ADMIN_NAME, FRAUNHOFER_ADMIN_PASS, FRAUNHOFER_URL
import urllib
import pickle
from urllib.parse import urlparse

#from constants import ADMIN_NAME, ADMIN_PASS, URL

from py2neo import Node, Relationship
#from py2neo.database import Transaction
from tqdm import tqdm

from connection import populate_db,_add_nodes
from constants import ENCODING, DATA_DIR



def create_users(
    #data_df: pd.DataFrame,
    url: str,
    name: str,
    password: str
):
    """Create automatic users."""
    graph = Graph(url, auth=(name, password), name='system')

    data_df = pd.read_csv(
        os.path.join(DATA_DIR, "userdetails.csv"),
        dtype=str,
        encoding=ENCODING
    )

    known_users = graph.run("SHOW USERS").to_series()

    for idx in tqdm(data_df['UserName'].values):
        #idx = idx.split('-')[0]

        if idx in known_users.values:  # Omit already present users
            continue

        cypher = f"CREATE USER {idx} SET PASSWORD 'abc' CHANGE REQUIRED"
        graph.run(cypher)

    # Save user details
    #data_df = data_df[['First Name', 'Last Name', 'UserName']]
    #data_df.to_csv(f'{DATA_DIR}/username_details_export.csv', sep=',', index=False)

def createNodes(tx):

    #read 3 protac databases from csv file
    # read protacdb with customized names

    ptacdb = pd.read_csv(
        os.path.join(DATA_DIR, "ptacdb_May2022.csv"),
        dtype=str,
        encoding=ENCODING
    )

    # read protacpedia with customized names

    ptacpedia = pd.read_csv(
        os.path.join(DATA_DIR, "PtacPDwithCID.csv"),
        dtype=str,
        encoding=ENCODING
    )
    # print(ptacpedia.columns)

    # read pubchemW
    pchem = pd.read_csv(
        os.path.join(DATA_DIR, "pubchem.csv"),
        dtype=str,
        encoding=ENCODING
    )
    # print(pchem.columns)

    #create protac nodes

    # add disease names as node properties for each Protein
    geneDis = pd.read_csv(
        os.path.join(DATA_DIR, "disgenet.csv"),
        usecols=["geneSymbol", "diseaseName"],
        dtype=str,
        encoding=ENCODING
    )

    # print(geneDis.columns)

    geneNames = geneDis["geneSymbol"].unique()
    geneNames = sorted(geneNames)
    # print(geneNames[1])

    gene2disease = geneDis.groupby('geneSymbol')['diseaseName'].apply(list).tolist()
    # print(gene2disease[0])
    gene2disease_list = dict(zip(geneNames, gene2disease))
    # print(gene2disease_list['APP'])
    # print(gene2disease_list['APOE'])
    # print(gene2disease_list[[0][0]])

    # read data from drugCentral
    drugTarget = pd.read_csv(
        os.path.join(DATA_DIR, "drugTarget.csv"),
        usecols=["DRUG_NAME", "GENE", "TARGET_CLASS"],
        dtype=str,
        encoding=ENCODING
    )
    # remove rows that have no protein/gene info
    drugTarget = drugTarget[~drugTarget["GENE"].isnull()]

    # create a list with drug and target class info merged together
    drugNclass = []
    for drug, target_class in drugTarget[['DRUG_NAME', 'TARGET_CLASS']].values:
        a = drug + ' (' + target_class + ')'
        drugNclass.append(a)

    drugTarget['DrugInfo'] = drugNclass

    geneNames_drugTarget = drugTarget["GENE"].unique()
    geneNames_drugTarget = sorted(geneNames_drugTarget)

    gene2drug = drugTarget.groupby('GENE')['DrugInfo'].apply(list).tolist()
    gene2drug = dict(zip(geneNames_drugTarget, gene2drug))

    # print(gene2drug['CDK2'])

    # prodb_pchem_war = pd.read_csv(
    #     os.path.join(DATA_DIR, "warhead_mapped2protac.csv"),
    #     dtype=str,
    #     encoding=ENCODING
    #     )

    warhead = pd.read_csv(
        os.path.join(DATA_DIR, "warhead_May2022.csv"),
        # usecols=["Name", "Gene", "Smiles", "IC50 (nM)", "Assay (IC50)", "Molecular Formula", "Molecular Weight",
        #          "InChI Key", "InChI", "PubChem", "ChEMBL","protac_name",'Exact Mass','logP','logS','Heavy Atom Count','Ring Count','Hydrogen Bond Acceptor Count',
        #          'Hydrogen Bond Donor Count','Rotatable Bond Count','Topological Polar Surface Area'],
        dtype=str,
        encoding=ENCODING
    )

    # remove rows with no names for warheads
    #warhead = warhead[~warhead["Name"].isnull()]

    # print(warhead["Name"])
    # print(warhead.columns)
    # print(warhead.head(20))

    #read e3 ligand from ptacdb

    e3ligand = pd.read_csv(
        os.path.join(DATA_DIR,'e3ligand_May2022.csv'),
        usecols=['Uniprot','Target','e3ligand','IC50 (nM)','Smiles','Molecular Weight','Exact Mass','logP','logS','Heavy Atom Count','Ring Count','Hydrogen Bond Acceptor Count',
                 'Hydrogen Bond Donor Count','Rotatable Bond Count','Topological Polar Surface Area','Molecular Formula','InChI','InChI Key','BindingDB','PubChem','ChEMBL','protac_name'],
        dtype=str,
        encoding=ENCODING
    )
    #print(e3lig.columns)
    #e3lig_names = set(e3lig['Name'])
    #e3lig_names = [x for x in e3lig_names if str(x) != 'nan']
    #print(len(e3lig_names))
    #print(e3lig_names)

    #breakpoint()

    # read e3 ligase file from ubinet

    ubinet_e3 = pd.read_csv(
        os.path.join(DATA_DIR, "Categorized_human_E3_ligases.csv"),
        dtype=str,
        encoding=ENCODING
    )

    # print(ubinet_e3.columns)

    # read file with e3 annotations from ubinet
    ubinet_e3_anno = pd.read_csv(
        os.path.join(DATA_DIR, "Human_E3Ligase_Annotation.csv"),
        dtype=str,
        encoding=ENCODING
    )

    # print(ubinet_e3_anno.columns)

    #read iglue file

    iglue = pd.read_csv(
        os.path.join(DATA_DIR,'iglueWithCID.csv'),
        dtype=str,
        encoding=ENCODING
    )
    # create dictionary for all nodes

    node_dict = {'Protein': {}, 'E3 ligase': {}, 'Protac': {},
                 'Warhead': {}, 'E3 binder': {},'iGLUE':{}}

    print('node_dict created')

    #1 protacdb
    #select columns and create a list of them
    cols = ["protac_name","InChI", "InChI Key", "Smiles", "Molecular Weight", "Molecular Formula",
            "Heavy Atom Count", "Ring Count", "Hydrogen Bond Acceptor Count",
            "Hydrogen Bond Donor Count", "Rotatable Bond Count","Topological Polar Surface Area", "Article DOI","CID_pchem"]
    for protac, inchi, inchikey, smiles, mw, mf, hac, rc, hbac, hbdc, rbc,tpsa, source,cid in tqdm(
            ptacdb[cols].values, total=ptacdb.shape[0],desc='ptacdb'):
        # (protac, protacsyn, inchi, inchikey, smiles, mw, mf, hac, rc, hbac, hbdc, rbc, source) = row
        if protac in node_dict["Protac"]:
            continue

        node_dict["Protac"][protac] = Node("Protac", **{"Name": protac, "InChI": inchi,
                                                        "InChI Key": inchikey,
                                                        "Smiles": smiles, "Molecular Weight": mw,
                                                        "Molecular Formula": mf, "Hydrogen Atom": hac,
                                                        "Ring Count": rc,
                                                        "Hydrogen Bond Acceptor": hbac,
                                                        "Hydrogen Bond Donor": hbdc, "Rotatable Bond": rbc,
                                                        "Source": f"https://doi.org/{source}",
                                                        #"Structure": f"https://molview.org/?q={smiles}",
                                                        "PubChem":f"https://pubchem.ncbi.nlm.nih.gov/compound/{cid}"})

        tx.create(node_dict["Protac"][protac])

    inchikeys = {val['InChI Key']: i for i, val in node_dict["Protac"].items()}
    #print(inchikeys)

    #2 Protacpedia
    cols = ["ptac_name","InChI key","PROTAC SMILES","Cells","Active/Inactive","Ligand Name","Linker Type","Hbond acceptors",
            "Hbond donors","Off Targets Reported","TPSA","Pubmed","Ligand PDB","CID_pchem"]

    for protac, inchikey, smiles, cell, status, ligname, linkertype, hba, hbd, offtar, tpsa, pubmed, ligpdb, cid in tqdm(ptacpedia[cols].values, total=ptacpedia.shape[0]):

        if inchikey in inchikeys:


            node_dict["Protac"][inchikeys[inchikey]].update({"Status": status, "Ligand": ligname, "Linker Type": linkertype, "Off targets":offtar,
                                                             "PubMed": f"https://pubmed.ncbi.nlm.nih.gov/?term={pubmed}","ProtacPedia":protac,
                                                             "Ligand PDB":f"https://www.rcsb.org/structure/{ligpdb}",
                                                             "PubChem":f"https://pubchem.ncbi.nlm.nih.gov/compound/{cid}"})


        else:
            node_dict["Protac"][protac] = Node("Protac", ** {"Name":protac,"InChI Key":inchikey,"Smiles":smiles,"Cell":cell,"Status":status,
                                                             "Ligand Name":ligname,"Linker Type":linkertype,"Hydrogen Bond Acceptor":hba,
                                                             "Hydrogen Bond Donor":hbd,"Off targets":offtar,"PubMed":f"https://pubmed.ncbi.nlm.nih.gov/?term={pubmed}",
                                                             #"Structure": f"https://molview.org/?q={smiles}",
                                                             "Ligand PDB":f"https://www.rcsb.org/structure/{ligpdb}",
                                                             "PubChem":f"https://pubchem.ncbi.nlm.nih.gov/compound/{cid}"})

            tx.create(node_dict["Protac"][protac])

    inchikeys_2 = {val['InChI Key']: i for i, val in node_dict["Protac"].items()}
    #3 pubchem

    cols_pchem = ["cmpdname","cmpdsynonym","mw","mf","polararea","hbondacc","hbonddonor","rotbonds","inchi","inchikey","isosmiles","cid"]
    for protac, ptacsyn, mw, mf, polarea, hba, hbd, rbc, inchi, inchikey, smiles, cid in tqdm(pchem[cols_pchem].values, total=pchem.shape[0]):

        if inchikey in inchikeys_2:
            node_dict["Protac"][inchikeys_2[inchikey]].update({"Protac Name":protac,"Protac Synonym":ptacsyn,"Compound":f"https://pubchem.ncbi.nlm.nih.gov/compound/{cid}","testinchi":inchikey})

        else:
            node_dict["Protac"][protac] = Node("Protac", **{"Name":protac, "Protac Synonym":ptacsyn,"InChI":inchi,"InChI Key":inchikey,
                                        "Smiles":smiles,"Molecular Weight":mw,"Molecular Formula":mf,"Ring Count":rc,
                                                        "Hydrogen Bond Acceptor Count":hba,"Hydrogen Bond Donor Count":hbd,"Rotatable Bond Count":rbc,
                                                        "PubChem":f"https://pubchem.ncbi.nlm.nih.gov/compound/{cid}","Polar Surface area":polarea,
                                                        #"Structure": f"https://molview.org/?q={smiles}"
                                                        })

    #create target nodes
    for target, uniprot in tqdm(ptacdb[["Gene","Uniprot"]].values, total=ptacdb.shape[0]):

        if target in node_dict["Protein"]:
            continue
        #add proteins in dictionary plus create a dictionary of metadata using **
        node_dict["Protein"][target] = Node("Protein", **{"Name":target,"Uniprot":uniprot,
                                                          "Uniprot link":f"https://www.uniprot.org/uniprot/{uniprot}",
                                                          "PubMed":f"https://pubmed.ncbi.nlm.nih.gov/?term={target}",
                                                          "ProteinsPlus":f"https://proteins.plus/{uniprot}",
                                                          "Alphafill":f"https://alphafill.eu/model?id={uniprot}"})

        if target in gene2disease_list:
            #node_dict["Protein"][target] = Node("Protein", **{"Diseases": geneDiseaseMapping[target]})
            node_dict["Protein"][target].update({"Diseases": gene2disease_list[target]})

        if target in gene2drug:
            node_dict["Protein"][target].update({"Drug&Class": gene2drug[target]})

    for target, uniprot in tqdm(ptacpedia[["Gene_name","uniprotid"]].values, total=ptacpedia.shape[0]):

        if target in node_dict["Protein"]:
            continue

        node_dict["Protein"][target] = Node("Protein", **{"Name": target, "Uniprot": uniprot,
                                                          "Uniprot link": f"https://www.uniprot.org/uniprot/{uniprot}",
                                                          "PubMed": f"https://pubmed.ncbi.nlm.nih.gov/?term={target}",
                                                          "ProteinsPlus":f"https://proteins.plus/{uniprot}",
                                                          "Alphafill":f"https://alphafill.eu/model?id={uniprot}"})

        if target in gene2disease_list:
            #node_dict["Protein"][target] = Node("Protein", **{"Diseases": geneDiseaseMapping[target]})
            node_dict["Protein"][target].update({"Diseases": gene2disease_list[target]})

        if target in gene2drug:
            node_dict["Protein"][target].update({"Drug&Class": gene2drug[target]})


    #create e3 nodes

    for e3 in tqdm(ptacdb["E3 ligase"].values, total=ptacdb.shape[0]):

        if e3 in node_dict["E3 ligase"]:
            continue

        node_dict["E3 ligase"][e3] = Node("E3 ligase", **{"Name": e3})

        # print(e3)

        if e3 in ubinet_e3['E3'].values:
            # print("Found")
            # print(e3)
            # print(e3 in ubinet_e3['E3'].index.values)
            pos = ubinet_e3[ubinet_e3["E3"] == e3].index
            # print(pos)

            # print(pos[0])
            # a = ubinet_e3.loc[[pos[0]]]
            a = ubinet_e3.loc[pos, "Category"].values[0]
            # print(a)

            node_dict["E3 ligase"][e3].update({"E3 Category": a})

        if e3 in ubinet_e3_anno['Gene Name'].values:
            pos = ubinet_e3_anno[ubinet_e3_anno["Gene Name"] == e3].index

            seq = ubinet_e3_anno.loc[pos, "Sequence"].values[0]
            func = ubinet_e3_anno.loc[pos, "Function [CC]"].values[0]
            gomf = ubinet_e3_anno.loc[pos, "Gene ontology (molecular function)"].values[0]
            gobp = ubinet_e3_anno.loc[pos, "Gene ontology (biological process)"].values[0]
            gocc = ubinet_e3_anno.loc[pos, "Gene ontology (cellular component)"].values[0]

            node_dict["E3 ligase"][e3].update(
                {"Sequence": seq, "Function": func, "GO molecular function": gomf, "GO biological process": gobp,
                 "GO cellular component": gocc})

    # create E3 ligase nodes from protacpedia
    for e3 in tqdm(ptacpedia["E3 Ligase"].values, total=ptacpedia.shape[0]):

        if e3 in node_dict["E3 ligase"]:
            continue

        node_dict["E3 ligase"][e3] = Node("E3 ligase", **{"Name": e3})

        # print(e3)

        if e3 in ubinet_e3['E3'].values:
            # print("Found")
            # print(e3)
            # print(e3 in ubinet_e3['E3'].index.values)
            pos = ubinet_e3[ubinet_e3["E3"] == e3].index
            # print(pos)

            # print(pos[0])
            # a = ubinet_e3.loc[[pos[0]]]
            a = ubinet_e3.loc[pos, "Category"].values[0]
            # print(a)

            node_dict["E3 ligase"][e3].update({"E3 Category": a})

        if e3 in ubinet_e3_anno['Gene Name'].values:
            pos = ubinet_e3_anno[ubinet_e3_anno["Gene Name"] == e3].index

            seq = ubinet_e3_anno.loc[pos, "Sequence"].values[0]
            func = ubinet_e3_anno.loc[pos, "Function [CC]"].values[0]
            gomf = ubinet_e3_anno.loc[pos, "Gene ontology (molecular function)"].values[0]
            gobp = ubinet_e3_anno.loc[pos, "Gene ontology (biological process)"].values[0]
            gocc = ubinet_e3_anno.loc[pos, "Gene ontology (cellular component)"].values[0]

            node_dict["E3 ligase"][e3].update(
                {"Sequence": seq, "Function": func, "GO molecular function": gomf, "GO biological process": gobp,
                 "GO cellular component": gocc})

    warhead_cols = ["warhead", "Smiles", "IC50 (nM)", "Assay (IC50)", "Molecular Formula", "Molecular Weight", "InChI Key",
                    "InChI", "PubChem", "ChEMBL"]

    #create warhead nodes
    for whead, smiles, ic50, assay, mf, mw, inchikey, inchi, pubchem, chembl in tqdm(warhead[warhead_cols].values,
                                                                                     total=warhead.shape[0]):
        if whead in node_dict["Warhead"]:
            continue

        node_dict["Warhead"][whead] = Node("Warhead",
                                           **{"Name": whead, "Smiles": smiles, "IC 50": ic50, "Assay": assay,
                                              "Molecular Formula": mf, "Molecular Weight": mw, "InChI Key": inchikey,
                                              "InChI": inchi,
                                              "PubChem": f"https://pubchem.ncbi.nlm.nih.gov/compound/{pubchem}",
                                              "ChEMBL": f"https://www.ebi.ac.uk/chembl/compound_report_card/{chembl}"})

    for prtn in tqdm(warhead['Gene'].values, total=warhead.shape[0]):

        if prtn in node_dict['Protein']:
            continue
        node_dict['Protein'][prtn] = Node('Protein',**{"Name":prtn})

    e3ligand_cols = ['Target', 'e3ligand','IC50 (nM)','Smiles', 'Molecular Weight', 'Exact Mass', 'logP', 'logS',
               'Heavy Atom Count', 'Ring Count', 'Hydrogen Bond Acceptor Count',
               'Hydrogen Bond Donor Count', 'Rotatable Bond Count', 'Topological Polar Surface Area',
               'Molecular Formula', 'InChI', 'InChI Key', 'BindingDB', 'PubChem', 'ChEMBL']

    for target, e3ligd, ic50, smiles, mw, em, logp, logs, hac, rc, hbac, hbdc, rbc, tpsa, mf, inchi, inchikey, bindingdb, pubchem, chembl in tqdm(e3ligand[e3ligand_cols].values, total=e3ligand.shape[0]):

        if e3ligd in node_dict["E3 binder"]:
            continue

        node_dict["E3 binder"][e3ligd] = Node("E3 binder",
                                                **{"Name":e3ligd, "Smiles":smiles, "IC 50": ic50, "Molecular Formula": mf, "Molecular Weight": mw,
                                                   "Exact Mass": em, "logP": logp, "logS": logs, "Hydrogen Atom count": hac, "Hydrogen Bond Acceptor Count": hbac,
                                                   "Hydrogen Bond Donor Count": hbdc, "Rotable Bond Count": rbc, "Topological Polar Surface Area": tpsa,
                                                   "InChI": inchi, "InChI Key": inchikey,
                                                   "PubChem": f"https://pubchem.ncbi.nlm.nih.gov/compound/{pubchem}",
                                                    "ChEMBL": f"https://www.ebi.ac.uk/chembl/compound_report_card/{chembl}"})


    for ligase in tqdm(e3ligand['Target'].values, total=e3ligand.shape[0],desc='read e3 ligand file'):

        if ligase in node_dict['E3 ligase']:
            continue
        node_dict['E3 ligase'][ligase] = Node('E3 ligase',**{"Name":ligase})

    ###work with iglue file
    iglue_cols = ['Catalog No.', 'Product Name', 'Clinical Research Yes/NO',
       'FDA Approved Yes/No', 'CAS Number', 'Form', 'Targets',
       'Information', 'Max Solubility in DMSO', 'URL',
       'Compound ID', 'Well', 'Rack', 'Barcode', 'CID']

    infile = open('data/cid_allProp', 'rb')
    cid_allProp = pickle.load(infile)
    infile.close()

    #print(cid_allProp.keys())

    for cat, name, clinres, fda, cas, form, target, info, dmso, url, compid, well, rack, barcode, cid in tqdm(iglue[iglue_cols].values, total=iglue.shape[0],desc= 'read iglue file'):

        if name in node_dict["iGLUE"]:
            continue

        node_dict["iGLUE"][name] = Node("iGLUE", ** {'Name':name, 'Catalog Number':cat,'Clinical Research Status':clinres,
                                                      'FDA Status':fda,'CAS':cas,'Form':form,'Targets':target,'Info':info,
                                                      'Max Solubility in DMSO':dmso,'URL':url,'Compound ID':compid,
                                                      'Well':well,'Rack':rack,'Barcode':barcode,'CID':cid})

        #print(node_dict["iGLUE"][name])

        get_cid = node_dict["iGLUE"][name]['CID']
        #print('here: '+get_cid)
        if get_cid != 'None':
            #get_cid = node_dict['iGLUE'][name]['CID']
            #schem = cid_allProp[get_cid]['SCHEMBL']
            node_dict["iGLUE"][name].update({'Canonical Smiles':cid_allProp[get_cid]['canonical_smiles'],
                                             'Exact Mass':cid_allProp[get_cid]['exact_mass'],
                                             'Hydrogen Bond Acceptor Count':cid_allProp[get_cid]['h_bond_acceptor_count'],
                                             'Hydrogen Bond Donor Count':cid_allProp[get_cid]['h_bond_donor_count'],
                                             'XlogP':cid_allProp[get_cid]['xlogp'],
                                             'Rotatable Bond Count':cid_allProp[get_cid]['rotatable_bond_count'],
                                             'Topological Polar Surface Area':cid_allProp[get_cid]['tpsa'],
                                             'InChi':cid_allProp[get_cid]['inchi'],
                                             'InChi Key':cid_allProp[get_cid]['inchikey'],
                                             'Heavy Atom Count':cid_allProp[get_cid]['heavy_atom_count'],
                                             'Isomeric Smiles':cid_allProp[get_cid]['isomeric_smiles'],
                                             'Molecular Formula':cid_allProp[get_cid]['molecular_formula'],
                                             'Molecular Weight':cid_allProp[get_cid]['molecular_weight']})

            if 'CHEMBL' in cid_allProp[get_cid]:
                chem = cid_allProp[get_cid]['CHEMBL']
                node_dict["iGLUE"][name].update({'ChEMBL':f"https://www.ebi.ac.uk/chembl/compound_report_card/{chem}"})

            if 'SCHEMBL' in cid_allProp[get_cid]:
                schem = cid_allProp[get_cid]['SCHEMBL']
                node_dict["iGLUE"][name].update({'SCHEMBL':f"https://www.surechembl.org/chemical/{schem}"})

        #print(node_dict["iGLUE"][name])

    infile = open('data/Extend_UniProtInfo', 'rb')
    uniprotDict = pickle.load(infile)
    infile.close()

    uniprot_nodes = {'Reactome': {}, 'Function': {}, 'Biological Process': {}, 'Disease': {}}

    for id in uniprotDict:

        if (uniprotDict[id]["Reactome"]):

            reactome_list = list(uniprotDict[id]["Reactome"].keys())

            for reactome in reactome_list:

                if reactome in uniprot_nodes["Reactome"]:
                    continue

                uniprot_nodes["Reactome"][reactome] = Node("Reactome", **{'Name': reactome,
                                                                          'Reactome ID': f"https://reactome.org/content/detail/{uniprotDict[id]['Reactome'][reactome]}"})

        if (uniprotDict[id]["Function"]):

            function_list = list(uniprotDict[id]['Function'].keys())

            for function in function_list:

                if function in uniprot_nodes['Function']:
                    continue

                uniprot_nodes['Function'][function] = Node('Function', **{'Name': function,
                                                                          'Gene Ontolgy Function': f"https://www.ebi.ac.uk/QuickGO/term/{uniprotDict[id]['Function'][function]}"})

        if (uniprotDict[id]['BioProcess']):

            bioprocess_list = list(uniprotDict[id]['BioProcess'].keys())

            for bioprocess in bioprocess_list:

                if bioprocess in uniprot_nodes["Biological Process"]:
                    continue

                uniprot_nodes["Biological Process"][bioprocess] = Node('Biological Process', **{'Name': bioprocess,
                                                                                                'Gene Ontolgy Biological Process': f"https://www.ebi.ac.uk/QuickGO/term/{uniprotDict[id]['BioProcess'][bioprocess]}"})

        if (uniprotDict[id]['Disease']):

            disease_list = list(uniprotDict[id]['Disease'].keys())

            for disease in disease_list:

                if disease in uniprot_nodes['Disease']:
                    continue

                uniprot_nodes['Disease'][disease] = Node('Disease', **{'Name': disease,
                                                                       'OMIM': f"https://www.omim.org/entry/{uniprotDict[id]['Disease'][disease]}"})


    node_dict.update(uniprot_nodes)



#comment following codes to be faster
    # Add also updated nodes into graph
    for node_type in node_dict:
        _add_nodes(
            node_dict=node_dict[node_type],
            tx=tx
        )

    return(node_dict)


def createReln(tx,ptacNode):

    ptacdb = pd.read_csv(
        os.path.join(DATA_DIR, "ptacdb_May2022.csv"),
        dtype=str,
        encoding=ENCODING
    )

    # print(ptacdb.columns)

    # read protacpedia with customized names

    ptacpedia = pd.read_csv(
        os.path.join(DATA_DIR, "PtacPDwithCID.csv"),
        dtype=str,
        encoding=ENCODING
    )

    warhead = pd.read_csv(
        os.path.join(DATA_DIR, "warhead_May2022.csv"),
        dtype=str,
        encoding=ENCODING
    )

    e3ligand = pd.read_csv(
        os.path.join(DATA_DIR,'e3ligand_May2022.csv'),
        dtype=str,
        encoding=ENCODING
    )

    for ligase, e3lig, protac in tqdm(e3ligand[['Target','e3ligand','protac_name']].values):
        e3ligand2ptac = Relationship(ptacNode["E3 binder"][e3lig],"isSubPartOf",ptacNode["Protac"][protac])
        e3ligand2ligase = Relationship(ptacNode["E3 binder"][e3lig],"joins",ptacNode["E3 ligase"][ligase])
        tx.create(e3ligand2ptac)
        tx.create(e3ligand2ligase)

    for target, whead, protac in tqdm(warhead[["Gene","warhead","protac_name"]].values):
        warPro = Relationship(ptacNode["Warhead"][whead], 'isApartOf', ptacNode["Protac"][protac])
        warTar = Relationship(ptacNode["Warhead"][whead], 'targets', ptacNode["Protein"][target])
        tx.create(warPro)
        tx.create(warTar)

    for target, e3, ptac in tqdm(ptacdb[["Gene","E3 ligase","protac_name"]].values):
        e3Tac = Relationship(ptacNode["E3 ligase"][e3],"binds",ptacNode["Protac"][ptac])
        targetTac = Relationship(ptacNode["Protein"][target], 'binds', ptacNode["Protac"][ptac], **{"E3 ligase":e3})
        e3Target = Relationship(ptacNode["E3 ligase"][e3], 'ubiquitinates', ptacNode["Protein"][target])
        #warPro = Relationship(ptacNode["Warhead"][headwar], 'isApartOf', ptacNode["Protac"][ptac])
        #warTar = Relationship(ptacNode["Warhead"][headwar], 'targets', ptacNode["Protein"][target])

        tx.create(e3Tac)
        tx.create(targetTac)
        tx.create(e3Target)
        #tx.create(warPro)
        #tx.create(warTar)

    for target, e3, ptac, linker in tqdm(ptacpedia[["Gene_name","E3 Ligase","ptac_name","Linker Type"]].values):
        #for protac in ptacNode["Protac"]:
            #print(ptac)
            #print(protac)
            #
            # if ptac in ptacNode["Protac"][protac].values():
            #
            #     print("yes it is")
            #     print(protac)
            #     break

        if ptac in ptacNode["Protac"]:
            e3Tac = Relationship(ptacNode["E3 ligase"][e3],"binds",ptacNode["Protac"][ptac])
            targetTac = Relationship(ptacNode["Protein"][target], 'binds', ptacNode["Protac"][ptac], **{"E3 ligase":e3})

        else:
            for protac in ptacNode["Protac"]:
                if ptac in ptacNode["Protac"][protac].values():
                    e3Tac = Relationship(ptacNode["E3 ligase"][e3], "binds", ptacNode["Protac"][protac])
                    targetTac = Relationship(ptacNode["Protein"][target], 'binds', ptacNode["Protac"][protac], **{"E3 ligase": e3})

        e3Target = Relationship(ptacNode["E3 ligase"][e3], 'ubiquitinates', ptacNode["Protein"][target])

        tx.create(e3Tac)
        tx.create(targetTac)
        tx.create(e3Target)

def createReln_uniprot(tx,ptacNode):

    infile = open('data/Extend_UniProtInfo', 'rb')
    uniprotDict = pickle.load(infile)
    infile.close()

    uprot = list(uniprotDict.keys())

    print('started to create nodes from uniprot to reactome and GOs')

    for p in ptacNode['Protein']:
        #print('start')
        #print(ptacNode['Protein'][p]['Uniprot'])
        for u in uprot:
            print('uprot from dict: ' + u)
            if u == ptacNode['Protein'][p]['Uniprot']:
                # print('check completed: '+ node_dict_new['Protein'][p]['Uniprot'])
                # print(node_dict_new['Protein'][p]['Name'])

                if uniprotDict[u]['Reactome']:
                    racts = uniprotDict[u]['Reactome'].keys()
                    for r in racts:
                        uprot2reactome = Relationship(ptacNode['Protein'][p], 'isInvolvedin', ptacNode['Reactome'][r])
                        tx.create(uprot2reactome)

                if uniprotDict[u]['Function']:
                    function = uniprotDict[u]['Function'].keys()
                    for fun in function:
                        uprot2function = Relationship(ptacNode['Protein'][p], 'hasFunction', ptacNode['Function'][fun])
                        tx.create(uprot2function)

                if uniprotDict[u]['BioProcess']:
                    bps = uniprotDict[u]['BioProcess'].keys()
                    for bp in bps:
                        uprot2bp = Relationship(ptacNode['Protein'][p], 'hasBioProcess', ptacNode['Biological Process'][bp])
                        tx.create(uprot2bp)

                if uniprotDict[u]['Disease']:
                    diseases = uniprotDict[u]['Disease'].keys()
                    for d in diseases:
                        uprot2dis = Relationship(ptacNode['Protein'][p], 'HasDisease', ptacNode['Disease'][d])
                        tx.create(uprot2dis)

    print('nodes from uniprot to reactome and GOs created')


def ExtractFromUniProt(uniprot_id):
    from bioservices import UniProt
    Uniprot_Dict = []
    # Make a link to the UniProt webservice
    service = UniProt()

    for id in uniprot_id:

        # create URL for each uniprot id
        url = 'https://www.uniprot.org/uniprot/' + id + '.txt'
        print(url)

        # Retrieve data for id in text format
        ret_uprot = urllib.request.urlopen(url)

        print(id)
        id_copy = id
        i = 0
        j = 0
        id = {}
        id['Disease'] = {}
        id['Reactome'] = {}
        id['Function'] = {}
        id['BioProcess'] = {}
        # print(id)

        # parse each line looking for info about disease, pathway, funcn, bp and so on
        for line in ret_uprot:

            line = line.decode('utf-8')

            # parse lines with disease and extract disease names and omim ids
            if '-!- DISEASE:' in line:
                if ('[MIM:' in line):
                    dis = line.split(':')
                    # dis returns list of splitted text, [1] = name of dis, [2] = OMIM ID, extra chars need cleaning
                    # print(dis[1][1:-5])
                    # print(dis[2][:-1])
                    id['Disease'].update({dis[1][1:-5]: dis[2][:-1]})

            # extract reactome ids and names
            if 'Reactome;' in line:
                ract = line.split(';')
                # ract returns list with reactome id and name, needs cleaning
                id['Reactome'].update({ract[2][1:-2]: ract[1][1:]})
                # print(ract[1][1:])
                # print(ract[2][1:-2])

            # look for functions
            if ' F:' in line:
                if j < 5:
                    # take only first 5 entries for now
                    # print(j)
                    fn = line.split(';')
                    # fn returns list with GO ids and names
                    id['Function'].update({fn[2][3:]: fn[1][1:]})
                    # print(fn[1][1:])
                    # print(fn[2][3:])
                    j += 1

            # look for biological processes
            if ' P:' in line:
                if i < 5:
                    # take only first 5 entries for now
                    # print(i)
                    bp = line.split(';')
                    # bp returns list with GO ids and names
                    id['BioProcess'].update({bp[2][3:]: bp[1][1:]})
                    # print(bp[1][1:])
                    # print(bp[2][3:])
                    i += 1

        # fetch info about gene, len, seq, mass from api directly
        # Make a query string
        query = "accession:" + str(id_copy)

        # Define a list of columns we want to retrive
        # columnlist = "id,entry name,length,mass,go(biological process),go(molecular function), pathway,feature(TOPOLOGICAL DOMAIN),comment(DISEASE),pathway"
        columnlist = "id,entry name,genes(PREFERRED),length,mass,sequence"

        # Run the remote search
        result = service.search(query, frmt="tab", columns=columnlist)

        df_result = pd.read_table(io.StringIO(result))

        geneName = {'Gene': df_result['Gene names  (primary )'][0]}
        seqLen = {'Length': df_result['Length'][0]}
        # remove comma from mass values
        mass = {'Mass': int(df_result['Mass'][0].replace(",", ""))}
        seq = {'Sequence': df_result['Sequence'][0]}

        id.update(geneName)
        id.update(seqLen)
        id.update(mass)
        id.update(seq)

        Uniprot_Dict.append(id)

    Uniprot_Dict = dict(zip(uniprot_id, Uniprot_Dict))

    return(Uniprot_Dict)

def createGraph():
    # create a new database
    #db_name, graph = populate_db("protacv4")
    graph = Graph(
        FRAUNHOFER_URL,
        auth=(FRAUNHOFER_ADMIN_NAME, FRAUNHOFER_ADMIN_PASS),
    )

    # Define the scope
    #scope = ['https://spreadsheets.google.com/feeds', 'https://www.googleapis.com/auth/drive']

    # Add credentials to the account
    #creds = ServiceAccountCredentials.from_json_keyfile_name('credentials.json', scope)

    # Authorize the clientsheet
    #client = gspread.authorize(creds)

    db_name = graph.begin()
    graph.delete_all()  # delete existing data
    getPtac = createNodes(db_name)
    getRels_uniprot = createReln_uniprot(db_name, getPtac)
    getRels = createReln(db_name,getPtac)

    graph.commit(db_name)

    #creating peronalized logins
    #create_users(url=FRAUNHOFER_URL, name=FRAUNHOFER_ADMIN_NAME, password=FRAUNHOFER_ADMIN_PASS, data_df=data_df)

    #return getPtac

createGraph()
#create_users(FRAUNHOFER_URL,FRAUNHOFER_ADMIN_NAME,FRAUNHOFER_ADMIN_PASS)


# graph = Graph(
#     FRAUNHOFER_URL,
#     auth=(FRAUNHOFER_ADMIN_NAME, FRAUNHOFER_ADMIN_PASS),
# )
# db_name = graph.begin()
# node_Dict = createNodes(db_name)
#
# outfile = open(os.path.join(DATA_DIR, "node_dict_new"),'wb')
# pickle.dump(node_Dict,outfile)
# outfile.close()



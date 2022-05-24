import pandas as pd

#import gspread
#from oauth2client.service_account import ServiceAccountCredentials

import getopt
import os
import random
import sys
from py2neo import Graph
from constants import FRAUNHOFER_ADMIN_NAME, FRAUNHOFER_ADMIN_PASS, FRAUNHOFER_URL
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

    # Form the UserName of people
    #data_df['First Name'] = data_df['First Name'].map(lambda x: x.split()[0].split('-')[0].capitalize())
    #data_df['Last Name'] = data_df['Last Name'].map(lambda x: ''.join([i.capitalize() for i in x.split()]))
    #data_df['UserName'] = data_df['First Name'] + data_df['Last Name']

    # # Replace certain characters
    # replace_char = {
    #     'ö': 'oe',
    #     'é': 'e',
    #     'í': 'i',
    #     "O'": 'o',
    #     'ø': 'o',
    #     'ä' : 'ae'
    # }
    #
    # for key, val in replace_char.items():
    #     data_df['UserName'] = data_df['UserName'].str.replace(key, val)

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


    # create dictionary for all nodes

    node_dict = {'Protein': {}, 'E3 ligase': {}, 'Protac': {}, 'Warhead': {}, 'E3 binder': {}}

    #1 protacdb
    #select columns and create a list of them
    cols = ["protac_name","InChI", "InChI Key", "Smiles", "Molecular Weight", "Molecular Formula",
            "Heavy Atom Count", "Ring Count", "Hydrogen Bond Acceptor Count",
            "Hydrogen Bond Donor Count", "Rotatable Bond Count","Topological Polar Surface Area", "Article DOI","CID_pchem"]
    for protac, inchi, inchikey, smiles, mw, mf, hac, rc, hbac, hbdc, rbc,tpsa, source,cid in tqdm(
            ptacdb[cols].values, total=ptacdb.shape[0]):
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

        #tx.create(node_dict["Protac"][protac])

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

            #tx.create(node_dict["Protac"][protac])

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


    for ligase in tqdm(e3ligand['Target'].values, total=e3ligand.shape[0]):

        if ligase in node_dict['E3 ligase']:
            continue
        node_dict['E3 ligase'][ligase] = Node('E3 ligase',**{"Name":ligase})

    # Add also updated nodes into graph
    for node_type in node_dict:
        _add_nodes(
            node_dict=node_dict[node_type],
            tx=tx
        )

    return node_dict


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
    getRels = createReln(db_name,getPtac)
    graph.commit(db_name)

    #creating peronalized logins
    #create_users(url=FRAUNHOFER_URL, name=FRAUNHOFER_ADMIN_NAME, password=FRAUNHOFER_ADMIN_PASS, data_df=data_df)

    #return getPtac

createGraph()
#create_users(FRAUNHOFER_URL,FRAUNHOFER_ADMIN_NAME,FRAUNHOFER_ADMIN_PASS)
#all over again

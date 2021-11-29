import getopt
import os
import random
import sys

import pandas as pd
from py2neo import Node, Relationship
#from py2neo.database import Transaction
from tqdm import tqdm


from connection import populate_db,_add_nodes
from constants import ENCODING, DATA_DIR

#create a new database
db_name,graph = populate_db("protacv1")

#read protacdb with customized names

ptacdb = pd.read_csv(
    os.path.join(DATA_DIR,"protacdb_withnames.csv"),
    dtype=str,
    encoding=ENCODING
)

#print(ptacdb.columns)

#read protacpedia with customized names

ptacpedia = pd.read_csv(
    os.path.join(DATA_DIR,"protacpd_withnames.csv"),
    dtype=str,
    encoding=ENCODING
)
#print(ptacpedia.columns)

#read pubchemW
pchem = pd.read_csv(
    os.path.join(DATA_DIR,"pubchem.csv"),
    dtype=str,
    encoding=ENCODING
)
#print(pchem.columns)

#add disease names as node properties for each Protein

geneDis = pd.read_csv(
    os.path.join(DATA_DIR, "disgenet.csv"),
    usecols=["geneSymbol","diseaseName"],
    dtype=str,
    encoding=ENCODING
)

#print(geneDis.columns)

geneNames = geneDis["geneSymbol"].unique()
geneNames = sorted(geneNames)
#print(geneNames[1])

gene2disease = geneDis.groupby('geneSymbol')['diseaseName'].apply(list).tolist()
#print(gene2disease[0])
gene2disease_list = dict(zip(geneNames,gene2disease))
#print(gene2disease_list['APP'])
#print(gene2disease_list['APOE'])
#print(gene2disease_list[[0][0]])

#read data from drugCentral
drugTarget = pd.read_csv(
    os.path.join(DATA_DIR,"drugTarget.csv"),
    usecols= ["DRUG_NAME","GENE","TARGET_CLASS"],
    dtype=str,
    encoding=ENCODING
)
#remove rows that have no protein/gene info
drugTarget = drugTarget[~drugTarget["GENE"].isnull()]

#create a list with drug and target class info merged together
drugNclass= []
for drug, target_class in drugTarget[['DRUG_NAME','TARGET_CLASS']].values:
    a = drug + ' ('+ target_class + ')'
    drugNclass.append(a)

drugTarget['DrugInfo'] = drugNclass

geneNames_drugTarget = drugTarget["GENE"].unique()
geneNames_drugTarget = sorted(geneNames_drugTarget)

gene2drug = drugTarget.groupby('GENE')['DrugInfo'].apply(list).tolist()
gene2drug = dict(zip(geneNames_drugTarget,gene2drug))


#print(gene2drug['CDK2'])

#read file
prodb_pchem = pd.read_csv(
    os.path.join(DATA_DIR, "protacDB_pubchem.csv"),
    dtype=str,
    encoding=ENCODING
    )

#print(prodb_pchem.columns)
#print(prodb_pchem.head(5))

prodb_pchem_war = pd.read_csv(
    os.path.join(DATA_DIR, "warhead_mapped2protac.csv"),
    dtype=str,
    encoding=ENCODING
    )

warhead = pd.read_csv(
    os.path.join(DATA_DIR, "warhead.csv"),
    usecols= ["Name","Target","Smiles","IC50 (nM)","Assay (IC50)","Molecular Formula","Molecular Weight","InChI Key","InChI","PubChem","ChEMBL"],
    dtype=str,
    encoding=ENCODING
)

#remove rows with no names for warheads
warhead = warhead[~warhead["Name"].isnull()]

#print(warhead["Name"])
#print(warhead.columns)
#print(warhead.head(20))

#read e3 ligase file from uninet

ubinet_e3 = pd.read_csv(
    os.path.join(DATA_DIR,"Categorized_human_E3_ligases.csv"),
    dtype=str,
    encoding=ENCODING
)

#print(ubinet_e3.columns)

#read file with e3 annotations from ubinet
ubinet_e3_anno = pd.read_csv(
    os.path.join(DATA_DIR,"Human_E3Ligase_Annotation.csv"),
    dtype=str,
    encoding=ENCODING
)

#print(ubinet_e3_anno.columns)

#read file sent by Andrea, PPI of computed eukaryotic proteins

# ppi_eu = pd.read_csv(
#         os.path.join(DATA_DIR,"computed_ppi.csv"),
#         dtype=str,
#         encoding=ENCODING
# )
#
# ppi_eu = ppi_eu.join(ppi_eu["pdbmono1"].str.split("_",1,expand=True).rename(columns={0:"pdb1",1:"g1"}))
# ppi_eu = ppi_eu.drop(columns = ['Unnamed: 26','Unnamed: 27'])
#
# ppi_eu = ppi_eu.join(ppi_eu["pdbmono2"].str.split("_",1,expand=True).rename(columns={0:"pdb2",1:"g2"}))
#
# print(ppi_eu.columns)
#print(ppi_eu["pdb2"].head(10))

#create protac nodes

#create dictionary for all nodes

node_dict = {'Protein': {},'E3 ligase': {},'Protac' : {},'Warhead' : {}}

def createPtac(tx):
    #create protac nodes

    #1 protacdb
    #select columns and create a list of them
    cols = ["protac_name","InChI", "InChI Key", "Smiles", "Molecular Weight", "Molecular Formula",
            "Heavy Atom Count", "Ring Count", "Hydrogen Bond Acceptor Count",
            "Hydrogen Bond Donor Count", "Rotatable Bond Count","Topological Polar Surface Area", "Article DOI"]
    for protac, inchi, inchikey, smiles, mw, mf, hac, rc, hbac, hbdc, rbc,tpsa, source in tqdm(
            ptacdb[cols].values, total=ptacdb.shape[0]):
        # (protac, protacsyn, inchi, inchikey, smiles, mw, mf, hac, rc, hbac, hbdc, rbc, source) = row
        if protac in node_dict["Protac"]:
            continue

        node_dict["Protac"][protac] = Node("Protac", **{"Protac": protac, "InChI": inchi,
                                                        "InChI Key": inchikey,
                                                        "Smiles": smiles, "Molecular Weight": mw,
                                                        "Molecular Formula": mf, "Hydrogen Atom": hac,
                                                        "Ring Count": rc,
                                                        "Hydrogen Bond Acceptor": hbac,
                                                        "Hydrogen Bond Donor": hbdc, "Rotatable Bond": rbc,
                                                        "Source": f"https://doi.org/{source}"})

        #tx.create(node_dict["Protac"][protac])

    inchikeys = {val['InChI Key']: i for i, val in node_dict["Protac"].items()}
    #print(inchikeys)

    #2 Protacpedia
    cols = ["ptac_name","InChI key","PROTAC SMILES","Cells","Active/Inactive","Ligand Name","Linker Type","Hbond acceptors",
            "Hbond donors","Off Targets Reported","TPSA","Pubmed"]

    for protac, inchikey, smiles, cell, status, ligname, linkertype, hba, hbd, offtar, tpsa, pubmed in tqdm(ptacpedia[cols].values, total=ptacpedia.shape[0]):

        if inchikey in inchikeys:


            node_dict["Protac"][inchikeys[inchikey]].update({"Cell": cell, "Status": status, "Ligand": ligname, "Linker Type": linkertype, "Off targets":offtar,"PubMed": pubmed,"ProtacPedia":protac})


        else:
            node_dict["Protac"][protac] = Node("Protac", ** {"Protac":protac,"InChI Key":inchikey,"SMILES":smiles,"Cell":cell,"Status":status,
                                                             "Ligand Name":ligname,"Linker Type":linkertype,"Hydrogen Bond Acceptor":hba,
                                                             "Hydrogen Bond Donor":hbd,"Off targets":offtar,"PubMed":pubmed})

            #tx.create(node_dict["Protac"][protac])

    inchikeys_2 = {val['InChI Key']: i for i, val in node_dict["Protac"].items()}
    #3 pubchem

    cols_pchem = ["cmpdname","cmpdsynonym","mw","mf","polararea","hbondacc","hbonddonor","rotbonds","inchi","inchikey","isosmiles","cid"]
    for protac, ptacsyn, mw, mf, polarea, hba, hbd, rbc, inchi, inchikey, smiles, cid in tqdm(pchem[cols_pchem].values, total=pchem.shape[0]):

        if inchikey in inchikeys_2:
            node_dict["Protac"][inchikeys_2[inchikey]].update({"Protac Name":protac,"Protac Synonym":ptacsyn,"Compound":f"https://pubchem.ncbi.nlm.nih.gov/compound/{cid}","testinchi":inchikey})

        else:
            node_dict["Protac"][protac] = Node("Protac", **{"Protac":protac, "Protac Synonym":ptacsyn,"InChI":inchi,"InChI Key":inchikey,
                                        "Smiles":smiles,"Molecular Weight":mw,"Molecular Formula":mf,"Ring Count":rc,
                                                        "Hydrogen Bond Acceptor Count":hba,"Hydrogen Bond Donor Count":hbd,"Rotatable Bond Count":rbc,
                                                        "Compound":f"https://pubchem.ncbi.nlm.nih.gov/compound/{cid}","Polar Surface area":polarea})


    # Add also updated nodes into graph
    for node_type in node_dict:
        _add_nodes(
            node_dict=node_dict[node_type],
            tx=tx
        )

    return node_dict

#fn to add nodes
def createNodes(tx,geneDiseaseMapping,drugInfo):

    #create dictionary of nodes
    node_dict = {
        'Protein': {},
        'E3 ligase': {},
        'Protac' : {},
        'Warhead' : {}
    }

    #create protein nodes
    for target, uniprot in tqdm(prodb_pchem[["Target","Uniprot"]].values, total=prodb_pchem.shape[0]):

        if target in node_dict["Protein"]:
            continue
        #add proteins in dictionary plus create a dictionary of metadata using **
        node_dict["Protein"][target] = Node("Protein", **{"Protein":target,"Uniprot":uniprot,"Uniprot link":f"https://www.uniprot.org/uniprot/{uniprot}",
                                                          "PubMed":f"https://pubmed.ncbi.nlm.nih.gov/?term={target}"})

        if target in geneDiseaseMapping:
            #node_dict["Protein"][target] = Node("Protein", **{"Diseases": geneDiseaseMapping[target]})
            node_dict["Protein"][target].update({"Diseases": geneDiseaseMapping[target]})

        if target in drugInfo:
            node_dict["Protein"][target].update({"Drug&Class": drugInfo[target]})

        tx.create(node_dict["Protein"][target])


    #create nodes for computed_ppi file

    #prtn_cols = ["gene name1","Uniprot1","UniFun1","KEGG1","Uniloc1","Phob1","pdbmono1","gene name2", "Uniprot2", "UniFun2", "KEGG2", "Uniloc2", "Phob2", "pdbmono2"]
    # prtn1_cols = ["gene name1", "Uniprot1", "UniFun1", "KEGG1", "Uniloc1", "Phob1", "pdb1"]
    # prtn2_cols = ["gene name2", "Uniprot2", "UniFun2", "KEGG2", "Uniloc2", "Phob2", "pdb2"]
    #
    # for p1, uniprot1, unifun1, kegg1, uniloc1, phob1, pdb1 in tqdm(ppi_eu[prtn1_cols].values, total=ppi_eu.shape[0]):
    #
    #     if p1 in node_dict["Protein"]:
    #         continue
    #
    #     node_dict["Protein"][p1] = Node("Protein", **{"Protein": p1, "Uniprot":uniprot1, "Uniprot link":f"https://www.uniprot.org/uniprot/{uniprot1}",
    #                                                   "UniFun":unifun1,"Phob":phob1,"PDB":f"https://www.rcsb.org/structure/{pdb1}"})
    #
    #     tx.create(node_dict["Protein"][p1])
    #
    # for p2, uniprot2, unifun2, kegg2, uniloc2, phob2, pdb2 in tqdm(ppi_eu[prtn2_cols].values, total=ppi_eu.shape[0]):
    #
    #     if p2 in node_dict["Protein"]:
    #         continue
    #
    #     node_dict["Protein"][p2] = Node("Protein", **{"Protein": p2, "Uniprot":uniprot2,"Uniprot link":f"https://www.uniprot.org/uniprot/{uniprot2}",
    #                                                   "UniFun":unifun2,"Phob":phob2,"PDB":f"https://www.rcsb.org/structure/{pdb2}"})
    #
    #     tx.create(node_dict["Protein"][p2])

    #create warhead nodes
    warhead_cols = ["Name","Smiles","IC50 (nM)","Assay (IC50)","Molecular Formula","Molecular Weight","InChI Key","InChI","PubChem","ChEMBL"]

    for whead, smiles, ic50, assay, mf, mw, inchikey, inchi, pubchem, chembl in tqdm(warhead[warhead_cols].values, total=warhead.shape[0]):
        if whead in node_dict["Warhead"]:
            continue

        node_dict["Warhead"][whead] = Node("Warhead", **{"Warhead":whead, "Smiles":smiles,"IC 50":ic50,"Assay":assay,"Molecular Formula":mf,"Molecular Weight": mw,"InChI Key":inchikey,"InChI":inchi,
                                    "PubChem":f"https://pubchem.ncbi.nlm.nih.gov/compound/{pubchem}","ChEMBL":f"https://www.ebi.ac.uk/chembl/compound_report_card/{chembl}"})

        tx.create(node_dict["Warhead"][whead])

    #create E3 ligase nodes
    for e3 in tqdm(prodb_pchem["E3 ligase"].values, total=prodb_pchem.shape[0]):

        if e3 in node_dict["E3 ligase"]:
            continue

        node_dict["E3 ligase"][e3] = Node("E3 ligase", **{"Name":e3})

        #print(e3)

        if e3 in ubinet_e3['E3'].values:
            #print("Found")
            #print(e3)
            #print(e3 in ubinet_e3['E3'].index.values)
            pos = ubinet_e3[ubinet_e3["E3"] == e3].index
            #print(pos)

            #print(pos[0])
            #a = ubinet_e3.loc[[pos[0]]]
            a = ubinet_e3.loc[pos,"Category"].values[0]
            print(a)

            #a=a['Category']
            #print(a)

            #print(a["Category"])
            #print(ubinet_e3.loc[[pos[0]]])
            node_dict["E3 ligase"][e3].update({"E3 Category": a})

        if e3 in ubinet_e3_anno['Gene Name'].values:

            pos = ubinet_e3_anno[ubinet_e3_anno["Gene Name"] == e3].index

            seq = ubinet_e3_anno.loc[pos,"Sequence"].values[0]
            func = ubinet_e3_anno.loc[pos,"Function [CC]"].values[0]
            gomf = ubinet_e3_anno.loc[pos,"Gene ontology (molecular function)"].values[0]
            gobp = ubinet_e3_anno.loc[pos,"Gene ontology (biological process)"].values[0]
            gocc = ubinet_e3_anno.loc[pos,"Gene ontology (cellular component)"].values[0]

            node_dict["E3 ligase"][e3].update({"Sequence": seq,"Function": func, "GO molecular function":gomf, "GO biological process": gobp,
                                               "GO cellular component":gocc})


        tx.create(node_dict["E3 ligase"][e3])

    #create protac ligase nodes
    cols = ["cmpdname","cmpdsynonym","InChI","InChI Key","Smiles","mw","mf",
            "Heavy Atom Count","Ring Count","Hydrogen Bond Acceptor Count",
            "Hydrogen Bond Donor Count","Rotatable Bond Count","Article DOI"]
    for protac, protacsyn, inchi, inchikey, smiles, mw, mf, hac, rc, hbac, hbdc, rbc, source in tqdm(prodb_pchem[cols].values, total=prodb_pchem.shape[0]):
        #(protac, protacsyn, inchi, inchikey, smiles, mw, mf, hac, rc, hbac, hbdc, rbc, source) = row
        if protac in node_dict["Protac"]:
            continue

        node_dict["Protac"][protac] = Node("Protac", **{"Protac":protac, "Protac Synonym":protacsyn,"InChI":inchi,"InChI Key":inchikey,
                                        "Smiles":smiles,"Molecular Weight":mw,"Molecular Formula":mf,"Hydrogen Atom Count":hac,"Ring Count":rc,
                                                        "Hydrogen Bond Acceptor Count":hbac,"Hydrogen Bond Donor Count":hbdc,"Rotatable Bond Count":rbc,
                                                        "Source":f"https://doi.org/{source}"})

        tx.create(node_dict["Protac"][protac])


    return(node_dict)

def createRel(tx,csvfile,passNode,warheadFile,ppiFile):

    for row in tqdm(csvfile[["cmpdname","E3 ligase","Target","Warhead_name"]].values):
        (protac, e3, target,headwar) = row
        e3Tac = Relationship(passNode["E3 ligase"][e3],'binds',passNode["Protac"][protac])
        targetTac = Relationship(passNode["Protein"][target], 'binds', passNode["Protac"][protac], **{"E3 ligase":e3})
        e3Target = Relationship(passNode["E3 ligase"][e3],'ubiquitinates',passNode["Protein"][target])
        warPro = Relationship(passNode["Warhead"][headwar],'isApartOf',passNode["Protac"][protac])
        warTar = Relationship(passNode["Warhead"][headwar], 'binds', passNode["Protein"][target])

        tx.create(e3Tac)
        tx.create(targetTac)
        tx.create(e3Target)
        tx.create(warPro)
        tx.create(warTar)

    # for name, target in tqdm(warheadFile[["Name","Target"]].values):
    #     if target in passNode["Protein"]:
    #         wheadTarget = Relationship(passNode["Warhead"][name],'binds',passNode["Protein"][target])
    #         tx.create(wheadTarget)

    # for p1, p2, ppiscore, phy_bgrid, gen_bgrid, string, pdbcomp in tqdm(ppiFile[["gene name1","gene name2","PPI Score","Phys_BioGRID",
    #                              "Gen_BioGRID","string","PDB_complex"]].values):
    #
    #     ppi = Relationship(passNode["Protein"][p1],'binds',passNode["Protein"][p2],
    #                        **{"PPI score":ppiscore,"Physical Interaction Score (BioGRID)":phy_bgrid,
    #                           "Genetic Interaction Score (BioGRID)":gen_bgrid,"String Score": string,"PDB complex":pdbcomp})
    #     tx.create(ppi)


getNodes = createPtac(db_name)

#ppi_eu not used, need to be called in the function for creating nodes and relns
#getRels = createRel(db_name,prodb_pchem_war,getNodes,warhead,ppi_eu)
#db_name.commit()
graph.commit(db_name)

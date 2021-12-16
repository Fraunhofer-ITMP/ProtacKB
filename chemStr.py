# This script is referred from http://rdkit.blogspot.jp/2015/02/new-drawing-code.html
# and http://cheminformist.itmol.com/TEST/wp-content/uploads/2015/07/rdkit_moldraw2d_2.html
# from __future__ import print_function
# from rdkit import Chem
# from rdkit.Chem.Draw import IPythonConsole
# from IPython.display import SVG
#
# from Ptac_V1 import ptacdb, ptacpedia, pchem
#
# from tqdm import tqdm
#
# import pandas as pd

# from rdkit.Chem import rdDepictor
# from rdkit.Chem.Draw import rdMolDraw2D
# def moltosvg(mol,molSize=(450,150),kekulize=True):
#     mc = Chem.Mol(mol.ToBinary())
#     if kekulize:
#         try:
#             Chem.Kekulize(mc)
#         except:
#             mc = Chem.Mol(mol.ToBinary())
#     if not mc.GetNumConformers():
#         rdDepictor.Compute2DCoords(mc)
#     drawer = rdMolDraw2D.MolDraw2DSVG(molSize[0],molSize[1])
#     drawer.DrawMolecule(mc)
#     drawer.FinishDrawing()
#     svg = drawer.GetDrawingText()
#     return svg

# def render_svg(svg):
#     # It seems that the svg renderer used doesn't quite hit the spec.
#     # Here are some fixes to make it work in the notebook, although I think
#     # the underlying issue needs to be resolved at the generation step
#     return SVG(svg.replace('svg:',''))
#
#
# ptacdb_str = []
#
# for smiles in tqdm(ptacdb["Smiles"].values, total=ptacdb.shape[0]):
#     mol = Chem.MolFromSmiles(smiles)
#     getStr = render_svg(moltosvg(mol))
#     ptacdb_str.append(getStr)

# print(len(ptacdb_str))
#
# ptacdb_test = ptacdb
# ptacdb_test["Structure"] = ptacdb_str
#
# with pd.ExcelWriter('testStr.xlsx') as writer:
#     ptacdb_test.to_excel(writer)
#
#
# def moltosvg(mol,molSize=(450,150),kekulize=True):
#     mc = Chem.Mol(mol.ToBinary())
#     if kekulize:
#         try:
#             Chem.Kekulize(mc)
#         except:
#             mc = Chem.Mol(mol.ToBinary())
#     if not mc.GetNumConformers():
#         rdDepictor.Compute2DCoords(mc)
#     drawer = rdMolDraw2D.MolDraw2DSVG(molSize[0],molSize[1])
#     drawer.DrawMolecule(mc)
#     drawer.FinishDrawing()
#     svg = drawer.GetDrawingText()
#     return svg
#
# def render_svg(svg):
#     # It seems that the svg renderer used doesn't quite hit the spec.
#     # Here are some fixes to make it work in the notebook, although I think
#     # the underlying issue needs to be resolved at the generation step
#     return SVG(svg.replace('svg:',''))


# ptacdb_str = []
# new_img = Image.new("L", (400, 150), "white")
#
# for smiles in tqdm(ptacpedia["PROTAC SMILES"].values, desc= "generating structures",total=ptacpedia.shape[0]):
#     mol = Chem.MolFromSmiles(smiles)
#     getStr = render_svg(moltosvg(mol))
#     #ptacdb_str.append(getStr)
#     new_img.putdata(new_img_list)
#     new_img.save('out.tif')
#
# #print(len(ptacdb_str))
#
# ptacpedia["Structure"] = new_img


# from __future__ import print_function
# from rdkit import Chem
# from rdkit.Chem.Draw import IPythonConsole
# from IPython.display import SVG
#
# from rdkit.Chem import rdDepictor
# from rdkit.Chem.Draw import rdMolDraw2D
#
#
# def moltosvg(mol,molSize=(450,150),kekulize=True):
#     mc = Chem.Mol(mol.ToBinary())
#     if kekulize:
#         try:
#             Chem.Kekulize(mc)
#         except:
#             mc = Chem.Mol(mol.ToBinary())
#     if not mc.GetNumConformers():
#         rdDepictor.Compute2DCoords(mc)
#     drawer = rdMolDraw2D.MolDraw2DSVG(molSize[0],molSize[1])
#     drawer.DrawMolecule(mc)
#     drawer.FinishDrawing()
#     svg = drawer.GetDrawingText()
#     return svg
#
# def render_svg(svg):
#     # It seems that the svg renderer used doesn't quite hit the spec.
#     # Here are some fixes to make it work in the notebook, although I think
#     # the underlying issue needs to be resolved at the generation step
#     return SVG(svg.replace('svg:',''))
#
#
# ptacdb_str = []
# new_img = Image.new("L", (400, 150), "white")
#
# for smiles in tqdm(ptacpedia["PROTAC SMILES"].values, desc= "generating structures",total=ptacpedia.shape[0]):
#     mol = Chem.MolFromSmiles(smiles)
#     getStr = render_svg(moltosvg(mol))
#     #ptacdb_str.append(getStr)
#     new_img.putdata(new_img_list)
#     new_img.save('out.tif')
#
# #print(len(ptacdb_str))
#
# ptacpedia["Structure"] = new_img

from Ptac_V1 import getPtac

from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import AllChem

import io
from PIL import Image

#print(getPtac.keys())
print(getPtac["Protac"]['ProtacDB 651 1419.078'].keys())
#print(getPtac["Protac"]['ProtacDB 651 1419.078']['Smiles'])


# for item in getPtac["Protac"]:
#     print(getPtac["Protac"][item]["Smiles"])
#     break

#
# for item in getPtac["Protac"]:
#     #print(getPtac["Protac"][item]["Protac"])
#     a=getPtac["Protac"][item]["Smiles"]
#
#     #mol= Chem.MolFromSmiles(a)
#     print("done")
#     break
#
# for item in getPtac["Protac"]:
#
#     path = r"C:\Users\reagon.karki\PycharmProjects\ProtacKB\data\image"
#     imageName = str(getPtac["Protac"][item]["Protac"])
#     #print(imageName)
#     fullPath = path + "\\" + imageName + ".png"
#     #print(fullPath)
#
#     a=getPtac["Protac"][item]["Smiles"]
#     mol = Chem.MolFromSmiles(a)
#     d2d = Draw.MolDraw2DCairo(500, 500)
#     d2d.DrawMolecule(mol)
#     d2d.FinishDrawing()
#     png_data = d2d.GetDrawingText()
#
#
#
#     #image.show()
#     #getPtac["Protac"][item].update({"2D Structure": image})
#
#     # save png to file
#     with open(fullPath, 'wb') as png_file:
#         png_file.write(png_data)
#

#print(imageList[0])
#print(mol)

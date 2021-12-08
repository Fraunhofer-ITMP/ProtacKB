from rdkit import Chem
from rdkit.Chem import Draw

import matplotlib.pyplot as plt
#%matplotlib inline

penicillin_g_smiles = 'CC1([C@@H](N2[C@H](S1)[C@@H](C2=O)NC(=O)Cc3ccccc3)C(=O)O)C'

penicillin_g = Chem.MolFromSmiles(penicillin_g_smiles)

#Draw.MolToMPL(penicillin_g, size=(200, 200))
AllChem.Compute2DCoords(penicillin_g)
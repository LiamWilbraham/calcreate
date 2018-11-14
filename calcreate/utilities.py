from contextlib import contextmanager
import os
from rdkit import Chem
from rdkit.Chem import AllChem
import subprocess as sp

def generate_xyz(smiles, xtb):
    mol = Chem.MolFromSmiles(smiles)
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, AllChem.ETKDG())
    print(Chem.MolToMolBlock(mol), file=open('mol','w+'))

    p = sp.Popen(['babel', 'mol', 'xyz'], stdout=sp.PIPE)
    o, e = p.communicate()

    if xtb:
        calc_params = ['xtb', 'xyz', '-opt']
        p = sp.Popen(calc_params, stdout=sp.PIPE, encoding='utf8')
        output, _ = p.communicate()

        with open('xtbopt.xyz') as f:
            return f.readlines()[2:]

    else:
        with open('xyz') as f:
            return f.readlines()[2:]

@contextmanager
def cd(newdir):
    prevdir = os.getcwd()
    os.chdir(os.path.expanduser(newdir))
    try:
        yield
    finally:
        os.chdir(prevdir)

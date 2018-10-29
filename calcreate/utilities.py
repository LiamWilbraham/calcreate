

def generate_xyz(smiles, xtb):
    #mol = Chem.MolFromSmiles(smiles)
    #mol = Chem.AddHs(mol)
    #AllChem.EmbedMolecule(mol, AllChem.ETKDG())
    #print(Chem.MolToMolBlock(mol), file=open('mol','w+'))

    #sp.call(['babel', 'mol', 'xyz'])

    if xtb:
        calc_params = ['xtb', 'xyz', '-opt']
        p = sp.Popen(calc_params, stdout=sp.PIPE, encoding='utf8')
        output, _ = p.communicate()

        with open('xtbopt.xyz') as f:
            return f.readlines()[2:]

    else:
        with open('xyz') as f:
            return f.readlines()[2:]

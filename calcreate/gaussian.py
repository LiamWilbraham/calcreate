from rdkit import Chem
from rdkit.Chem import AllChem
import subprocess as sp
import os
from .utilities import *

class Gaussian:

    '''
    A class for creating Gaussian input files from SMILES strings.
    Arbitrary keywords can be specified and initial geometries can either
    be approximated by the ETKDG method or via semi-empirical tight binding
    method, GFN-xTB (if available.)

    Parameters
    ----------
    name : :class:`string`
        Name to be given to the input and checkpoint files. The suffix '.com'
        will be automatically appended.

    smiles : :class:`string`
        SMILES string of the molecule on which calculations are to be performed.

    func : :class:`string` (default = 'pbepbe')
        Density functional to be usedself.

    basis : :class:`string` (default = '6-31G')
        Basis set to be used.

    charge : :class:`int` (default = 0)
        Molecular charge.

    mult : :class:`int` (default = 1)
        Spin multiplicity (2S + 1).

    kwds : :class:`list` (default = [])
        List of gaussian keywords to be applied (e.g. 'opt', 'td=(Nstates=10)').

    nprocs : :class:`int` (default = 8)
        Number of processors to be used.

    mem : :class:`str` (default = '1GB')
        Memory allowance.

    xtb : :class:`bool` (default = False)
        If True, xtb is used to produce a starting structure, if False, ETKDG
        is used by default.
    '''

    def __init__(self,
                 name,
                 smiles,
                 func='pbepbe',
                 basis='6-31G',
                 charge=0,
                 mult=1,
                 kwds=[],
                 nprocs=8,
                 mem='1GB',
                 xtb=False
                 ):

        self.name = name
        self.smiles = smiles
        self.func = func
        self.basis = basis
        self.charge = charge
        self.mult = mult
        self.kwds = kwds
        self.nprocs = nprocs
        self.mem = mem
        self.xtb = xtb

    def generate_input(self):

        '''
        Generates Gaussian input file.
        '''

        header = self._generate_header()
        kwds = self._generate_kwds()
        xyz = generate_xyz(self.smiles, self.xtb)
        self._write_input(header, kwds, xyz)
        self._remove_junk()


    def _generate_header(self):
        string = "%chk={}.chk\n%nprocs={}\n%mem={}\n"
        return string.format(self.name, self.nprocs, self.mem)


    def _generate_kwds(self):
        string = '{} '*len(self.kwds)+'\n'
        return string.format(*self.kwds)


    def _write_input(self, header, kwds, xyz):
        with open('{}.txt'.format(self.name), 'w') as f:
            f.write(header)
            f.write('#P {}/{} '.format(self.func, self.basis))
            f.write(kwds)
            f.write('\n{}\n\n'.format(self.name))
            f.write('{} {}\n'.format(self.charge, self.mult))
            for line in xyz:
                f.write(line)
            f.write('\n\n\n\n')


    def _remove_junk(self):
        os.system('rm mol xyz')

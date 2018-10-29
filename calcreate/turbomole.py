from .utilities import *

class Turbomole:

    def __init__(self,
                 name,
                 smiles,
                 func='pbe-pbe',
                 basis='DZP',
                 charge=0,
                 mult=1,
                 kwd_dict={'opt':300, 'rpas':5, 'exopt':300},
                 solvent=None,
                 xtb=False
                 ):

        self.name = name
        self.smiles = smiles
        self.func = func
        self.basis = basis
        self.charge = charge
        self.mult = mult
        self.kwd_dict = kwd_dict
        self.solvent = solvent
        self.xtb = xtb

    def generate_input():
        pipe = _generate_pipe()
        _generate_coord()
        _define(pipe)


    def _define():
        pass


    def _generate_pipe():
        pass

    def _generate_coord():
        pass

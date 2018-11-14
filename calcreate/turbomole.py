from .utilities import *
import subprocess as sp

class Turbomole:

    def __init__(self,
                 name,
                 smiles,
                 func='pbe-pbe',
                 basis='DZP',
                 charge=0,
                 #kwd_dict={'rpas':5},
                 solvent_epsilon=None,
                 xtb=False
                 ):

        self.name = name
        self.smiles = smiles
        self.func = func
        self.basis = basis
        self.charge = charge
        #self.kwd_dict = kwd_dict
        self.solvent_epsilon = solvent_epsilon
        self.xtb = xtb

    def generate_input(self):
        pipe = self._generate_pipe()
        generate_xyz(self.smiles, self.xtb)
        self._generate_coord()
        self._define(pipe)
        self._add_solvent()

    def _define(self, pipe):
        p = sp.Popen(['define'], stdout=sp.PIPE, stdin=sp.PIPE)
        o, e = p.communicate(pipe.encode())

    def _generate_pipe(self):
        pipe = 'a coord\nired\n*\neht\n\n'
        pipe += str(self.charge)
        pipe += '\n\n\nscf\niter\n300\n\ndft\non\nfunc '
        pipe += self.func
        pipe += '\n\n*\n\n'
        return pipe

    def _generate_coord(self):
        sp.call(['x2t', 'xyz', 'coord'])

    def _add_solvent(self):
        with open('control') as f:
            control = f.readlines()[:-1]
        control += '$cosmo\nepsilon={}\n$end'.format(self.solvent_epsilon)
        with open('control') as f:
            f.write(control)

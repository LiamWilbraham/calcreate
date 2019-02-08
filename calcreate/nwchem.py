from .utilities import *
import subprocess as sp
import errno

class NWChem:

        def __init__(self,
                     name,
                     smiles,
                     symmetry=None,
                     basis='def2-SVPD',
                     dft=False,
                     functional='b3lyp',
                     stackmem=1800,
                     heapmem=100,
                     globalmem=2100,
                     xtb=False):

            self.name = name
            self.smiles = smiles
            self.symmetry = symmetry
            self.basis = basis
            self.dft = dft
            self.functional = functional
            self.stackmem = stackmem
            self.heapmem = heapmem
            self.globalmem = globalmem
            self.xtb = xtb

        def generate_input(self):
            xyz = generate_xyz(self.smiles, self.xtb)[1:]
            header = self._header_block()
            basis = self._basisset_block()
            task = self._task()
            self._write_input(header, xyz, basis, task)

        def _geometry_block(self):
            if symmetry is not 'None':
                symmetry_block = self.symmetry

        def _header_block(self):
            string = 'START\n{}\necho\nMEMORY stack {} mb heap {} mb global {}\n'
            return string.format(self.name, self.stackmem, self.heapmem, self.globalmem)

        def _task(self):
            if self.dft is not False:
                string = 'dft\nxc {}\nend\n\n'
                return string.format(self.functional)
            return 'task scf\n'

        def _basisset_block(self):
            string = 'basis\n* library {}\nend\n\n'
            return string.format(self.basis)

        def _write_input(self, header, xyz, basis, task):
            with open('{}.nw'.format(self.name), 'w') as f:
                f.write(header)
                f.write('geometry\n')
                for line in xyz:
                    f.write(line)
                f.write('end\n\n')
                f.write(basis)
                f.write(task)
                if self.dft is not False:
                    f.write('task dft\n')

        def _remove_junk(self):
            os.system('rm mol xyz')

import numpy as np
from six import string_types
from horton import DenseLinalgFactory
from horton.gbasis.gobasis import GOBasis
from horton.gbasis.io import str_to_shell_types
from horton.meanfield.wfn import RestrictedWFN, UnrestrictedWFN

BOHR_TO_ANGSTROM = 0.529177249


class MoldenFile(object):
    def __init__(self, filename):
        if isinstance(filename, string_types):
            self.f = open(filename, 'rb')
            self.own_handle = True
        else:
            self.f = filename
            self.own_handle = False

        self.currentMode = 'NotParsing'
        self.coordFactor = 1
        self.numBasisFunctions = 0
        self.aNums = []
        self.aPos = []
        self.shellTypes = []
        self.shellNums = []
        self.shellToAtom = []
        self.a = []
        self.c = []

        self.MOcoeffs_alpha = []
        self.MOcoeffs_beta = []        
        self.MOenergy_alpha = []
        self.MOenergy_beta = []
        self.MOoccup_alpha = []
        self.MOoccup_beta = []
        
        self._read()
        self.obasis = self._buildBasis()
        
        MOcoeffs_alpha = np.array(self.MOcoeffs_alpha)
        MOcoeffs_beta = np.array(self.MOcoeffs_beta)
        MOenergy_alpha = np.array(self.MOenergy_alpha)
        MOenergy_beta = np.array(self.MOenergy_beta)
        MOoccup_alpha = np.array(self.MOoccup_alpha)
        MOoccup_beta = np.array(self.MOoccup_beta)

        if len(MOcoeffs_beta) == 0:
            raise NotImplementedError()
        else:
            nalpha = int(np.round(MOoccup_alpha.sum()))
            nbeta = int(np.round(MOoccup_beta.sum()))
            
            norb = min(len(MOcoeffs_alpha), len(MOcoeffs_beta))
            wfn = UnrestrictedWFN(DenseLinalgFactory(), self.obasis.nbasis, norb=norb)
            exp_alpha = wfn.init_exp('alpha')
            exp_alpha.coeffs[:] = MOoccup_alpha[:norb]
            exp_alpha.energies[:] = MOenergy_alpha[:norb]
            exp_alpha.occupations[:] = MOoccup_alpha[:norb]
            exp_beta = wfn.init_exp('beta')
            exp_beta.coeffs[:] = MOoccup_beta[:norb]
            exp_beta.energies[:] =  MOenergy_beta[:norb]
            exp_beta.occupations[:] = MOoccup_beta[:norb]

            self.wfn = wfn

    def _read(self):
        for line in self.f:
            try:
                self._process_line(line)
            except StopIteration:
                break

    def _process_line(self, line):
        key = line.strip()
        if len(key) == 0:
            return

        lkey = key.lower()
        words = key.split()
        
        if '[atoms]' in lkey:
            if 'angs' in lkey:
                self.coordFactor = 1.0 / BOHR_TO_ANGSTROM
            elif 'au' in lkey:
                self.coordFactor = 1.0
            self.currentMode = 'Atoms'
        elif '[gto]' in lkey:
            self.currentMode = 'GTO'
        elif '[mo]' in lkey:
            self.currentMode = 'MO'
        elif '[' in lkey:
            self.currentMode = 'NotParsing'
        else:
            if self.currentMode == 'Atoms':
                self._processAtoms(words)
            elif self.currentMode == 'GTO':
                self._processGTO(words)
            elif self.currentMode == 'MO':
                self._processMO(key)

    def _processAtoms(self, words):
        if (len(words) < 6):
            return       
        self.aNums.append(int(words[2]))
        self.aPos.append(self.coordFactor * np.array((float(words[3]), float(words[4]), float(words[5]))))
    
    def _processGTO(self, words):
        atom = int(words[0])
        key = next(self.f).strip()
        while len(key) > 0:
            words = key.split()
            shell = words[0].lower()
            shellType = None
            if 'sp' in shell:
                shellType = 'SP'
            elif 's' in shell:
                shellType = 'S'
            elif 'p' in shell:
                shellType = P;
            elif 'd' in shell:
                shellType = 'D'
            elif 'f' in shell:
                shellType = 'F'
                
            if shellType is None:
                raise ValueError('unknown shell type: %s' % shell)
            if shellType == 'SP':
                self.shellTypes.append('S')
                self.shellTypes.append('P')
                self.shellToAtom.append(atom)
                self.shellToAtom.append(atom)
            else:
                self.shellTypes.append(shellType)
                self.shellToAtom.append(atom)

            numGTOs = int(words[1])
            self.shellNums.append(numGTOs)
            if shellType == 'SP':
                self.shellNums.append(numGTOs)

            # now read all the exponents and contraction coefficients
            for i in range(numGTOs):
                words = next(self.f).strip().split()
                self.a.append(float(words[0]))
                self.c.append(float(words[1]))
                if (shellType == 'SP' and len(words) > 2):
                    self.a.append(float(words[0]))
                    self.c.append(float(words[2]))
            # start reading the next shell
            key = next(self.f).strip()

    def _processMO(self, key):
        # parse occ, spin, energy, etc.
        expansion_coeff = []
        
        energy = None
        spin = None
        occup = None
        
        while '=' in key:
            words = key.strip().split('=')
            w0 = words[0].lower()
            if w0 == 'ene':
                energy = float(words[1])
            elif w0 == 'spin':
                spin = words[1]
            elif w0 == 'sym':
                pass
            elif w0 == 'occup':
                occup = float(words[1])
            else:
                raise ValueError()
            key = next(self.f)
        while len(key) > 0 and '=' not in key:
            words = key.strip().split()
            if len(words) < 2:
                return
            expansion_coeff.append(float(words[1]))
            try:
                key = next(self.f)
            except StopIteration:
                break

        assert (energy is not None) and \
               (spin is not None) and (occup is not None)
        if spin == 'Alpha':
            self.MOcoeffs_alpha.append(np.array(expansion_coeff))
            self.MOenergy_alpha.append(energy)
            self.MOoccup_alpha.append(occup)
        elif spin == 'Beta':
            self.MOcoeffs_beta.append(np.array(expansion_coeff))
            self.MOenergy_beta.append(energy)
            self.MOoccup_beta.append(occup)

    def _buildBasis(self):
        coordinates = np.array(self.aPos)
        shell_map = np.array(self.shellToAtom) - 1
        nprims = np.array(self.shellNums)
        shell_types = np.array([str_to_shell_types(e, False)[0] for e in self.shellTypes])
        alphas = np.array(self.a)
        con_coeffs = np.array(self.c)
        return GOBasis(coordinates, shell_map, nprims, shell_types, alphas, con_coeffs)

    def close(self):
        if self.own_handle:
            self.f.close()


import numpy as np
from os.path import join
import gzip
import glob
from tempfile import NamedTemporaryFile
import molden2
from pyvdwsurface import vdwsurface

ELEMENTS = {
    1: 'H', 2: 'He', 3: 'Li', 4: 'Be', 5: 'B', 6: 'C', 7: 'N', 8: 'O', 9: 'F',
    10: 'Ne', 11: 'Na', 12: 'Mg', 13: 'Al', 14: 'Si', 15: 'P', 16: 'S',
    17: 'Cl', 18: 'Ar', 19: 'K', 20: 'Ca', 21: 'Sc', 22: 'Ti', 23: 'V',
    24: 'Cr', 25: 'Mn', 26: 'Fe', 27: 'Co', 28: 'Ni', 29: 'Cu', 30: 'Zn',
    31: 'Ga', 32: 'Ge', 33: 'As', 34: 'Se', 35: 'Br', 36: 'Kr', 37: 'Rb',
    38: 'Sr', 39: 'Y', 40: 'Zr', 41: 'Nb', 42: 'Mo', 43: 'Tc', 44: 'Ru',
    45: 'Rh', 46: 'Pd', 47: 'Ag', 48: 'Cd', 49: 'In', 50: 'Sn', 51: 'Sb',
    52: 'Te', 53: 'I', 54: 'Xe', 55: 'Cs', 56: 'Ba', 57: 'La', 58: 'Ce',
    59: 'Pr', 60: 'Nd', 61: 'Pm', 62: 'Sm', 63: 'Eu', 64: 'Gd', 65: 'Tb',
    66: 'Dy', 67: 'Ho', 68: 'Er', 69: 'Tm', 70: 'Yb', 71: 'Lu', 72: 'Hf',
    73: 'Ta', 74: 'W', 75: 'Re', 76: 'Os', 77: 'Ir', 78: 'Pt', 79: 'Au',
    80: 'Hg', 81: 'Tl', 82: 'Pb', 83: 'Bi', 84: 'Po', 85: 'At', 86: 'Rn',
    87: 'Fr', 88: 'Ra', 89: 'Ac', 90: 'Th', 91: 'Pa', 92: 'U', 93: 'Np',
    94: 'Pu', 95: 'Am', 96: 'Cm', 97: 'Bk', 98: 'Cf', 99: 'Es', 100: 'Fm',
    101: 'Md', 102: 'No', 103: 'Lr', 104: 'Rf', 105: 'Db', 106: 'Sg',
    107: 'Bh', 108: 'Hs', 109: 'Mt', 110: 'Ds', 111: 'Rg', 112: 'Cn',
    114: 'Uuq', 116: 'Uuh'}


def qc_molden_section(f):
    in_molden = False
    for line in f:
        if line.strip() == '======= MOLDEN-FORMATTED INPUT FILE FOLLOWS =======':
            in_molden = True
        elif line.strip() == '======= END OF MOLDEN-FORMATTED INPUT FILE =======':
            in_molden = False
        if in_molden:
            yield line


def compute_vectors(gzip_qcout, e_range, e_sigma):
    with gzip.open(gzip_qcout, 'rb') as f:
        mf = molden2.MoldenFile(qc_molden_section(f))

    basis = mf.obasis
    wfn = mf.wfn
    
    nuclear_coordinates = np.array(mf.aPos) * molden2.BOHR_TO_ANGSTROM
    elements = [ELEMENTS[anum] for anum in mf.aNums]
    
    # points on the vdw surface (n_surface_points x 3)
    surface = vdwsurface(nuclear_coordinates, elements)

    # indices of the orbitals to use (all of them)
    iorbs = np.arange(wfn.norb)
    # value of each of the MOs at each of the vdw surface points
    # (n_sufrace_points x n_orbitals)
    alpha_orbs = np.zeros((len(surface), len(iorbs)))
    beta_orbs = np.zeros((len(surface), len(iorbs)))
    basis.compute_grid_orbitals_exp(wfn.exp_alpha, surface, iorbs, alpha_orbs)
    basis.compute_grid_orbitals_exp(wfn.exp_beta, surface, iorbs, beta_orbs)

    weights_alpha = np.subtract.outer(e_range, wfn.exp_alpha.energies)**2 / (2*e_sigma)
    weights_alpha = (weights_alpha.sum(axis=1)**(-1))[:, np.newaxis] * weights_alpha
    weights_beta = np.subtract.outer(e_range, wfn.exp_beta.energies)**2 / (2*e_sigma)
    weights_beta = (weights_beta.sum(axis=1)**(-1))[:, np.newaxis] * weights_beta

    # n_suface_points x len(e_range)
    vectors = (np.dot(alpha_orbs**2, weights_alpha.T) +
               np.dot(beta_orbs**2, weights_beta.T))

    return vectors


if __name__ == '__main__':
    fn = join('543dd14591dd7a5ed5c13907', 'calculation.out.gz')
    e_sigma = 1
    e_range = np.linspace(-20, 2, 10)
    v = compute_vectors(fn, e_sigma=e_sigma, e_range=e_range)
    print(v.shape)

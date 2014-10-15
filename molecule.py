import os
import tempfile
import subprocess
from distutils.spawn import find_executable

from rdkit import Chem
from rdkit_utils.serial import MolReader, MolWriter
from pubchem_utils import PubChem


def cid2mol(cid):
    pc = PubChem()
    reader = MolReader()

    try:
        fd, fn = tempfile.mkstemp(suffix='.sdf.gz')
        os.close(fd)
        print('Getting record...')
        pc.get_records([cid], filename=fn, use_3d=True)
        print('Ok')
        with reader.open(fn) as f:
            return next(f.get_mols())
    finally:
        os.unlink(fn)


def protonate(mol, pH=7.4):
    babel = find_executable('obabel')
    assert babel is not None
    os.environ['BABEL_LIBDIR'] = os.path.join(os.path.dirname(babel), '..', 'lib', 'openbabel', '2.3.2')
    os.environ['BABEL_DATADIR'] = os.path.join(os.path.dirname(babel), '..', 'share', 'openbabel', '2.3.2')
    comm = subprocess.Popen([babel, '-o', 'sdf', '-i', 'sdf', '-p', str(pH)],
                            stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = comm.communicate(Chem.MolToMolBlock(mol, includeStereo=True))
    assert stderr.strip() == '1 molecule converted'
    return Chem.MolFromMolBlock(stdout)


def mol2xyz(mol):
    conf = mol.GetConformer()
    charge = sum(atom.GetFormalCharge() for atom in mol.GetAtoms())
    spin = 0.5 * (sum(atom.GetAtomicNum() for atom in mol.GetAtoms()) % 2)
    mult = int(2*spin + 1)

    lines = []
    lines.append('%s %s' % (charge, mult))
    for i, atom in enumerate(mol.GetAtoms()):
        pos = conf.GetAtomPosition(i)
        lines.append('%4s %10.5f %10.5f %10.5f' % (atom.GetSymbol(), pos.x, pos.y, pos.z))

    return lines


def qcheminpfile(mol, comment=''):
    lines = ['$comment', str(comment), '$end', '']
    lines.extend(['$molecule'] + mol2xyz(mol) + ['$end', ''])
    lines.extend([
        '$rem',
        'jobtype          sp',
        'exchange         hf',
        'basis            6-31G*',
        'print_orbitals   true',
        'molden_format    true',
        'max_scf_cycles   100',
        'thresh           7',
        'scf_convergence  4',
        '$end', ''
    ])
    return '\n'.join(lines)

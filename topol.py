
import pickle
import os
from shutil import which
from gentopol.zmat import ZMAT, get_OPT
from gentopol.BOSStools import BOSSReader, Check_H
from gentopol.BOSS2GMX import save2GMX


def TOPOL(**kwargs):
    """Configure the topology of a system.

    Keyword arguments:
    file -- input file coordinates.
    """
    file = kwargs['file']
    opt = kwargs['opt']
    charge = kwargs['charge']
    res = kwargs['resname'].lower()

    assert which('babel'), "OpenBabel is Not installed or the executable location is not accessable"

    # save file res.z
    ZMAT(file, opt, res, charge)

    # optimize geometry
    get_OPT('%s.z' % res, opt, charge)

    # Read the new geometry in a class
    mol = BOSSReader('%s.z' % res)

    assert (mol.MolData['TotalQ']['Reference-Solute'] ==
            charge), "PROPOSED CHARGE IS NOT POSSIBLE: SOLUTE MAY BE AN OPEN SHELL"

    assert(Check_H(mol.MolData['ATOMS'])
           ), "Hydrogens are not added. Please add Hydrogens"

    # saves the information in a binary
    pickle.dump(mol, open(res + ".p", "wb"))

    # Saving to GROMACS files
    save2GMX(res)
    print('Done')

    os.system(f'rm -v {res}.p')
    os.system(f'rm -v {res}.mol')
    os.system(f'rm -v {res}.z')
    mol.cleanup()

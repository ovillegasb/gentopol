"""Submodule dedicated to create a file with the Z matrix to be used by BOSS.
"""

import os
import numpy as np
from gentopol import algebra as al
import networkx as nx
from gentopol import files


def AsisZmat(file, opt, res):
    name, ext = file.split('.')

    # if file.endswith('pdb'):

    print(name, ext)


# def Gen_mol_rep(file, opt, res, charge):
#     Zmat(file, opt, res)


def load_mol(file):
    """Carga un archivo .mol."""
    mollines = open(file, 'r').readlines()
    print(mollines)

    [nats, nbonds] = map(int, (mollines[3][0:3], mollines[3][3:6]))
    print('N atoms:', nats)
    print('N bonds:', nbonds)

    coordlines = mollines[4:4 + nats]
    print(coordlines)

    coord = {}
    atsb = {}

    for i in range(nats):
        els = coordlines[i].split()
        coord[i + 1] = np.array([float(e) for e in els[0:3]])
        atsb[i + 1] = els[3]

    print(atsb)
    print(coord)

    bondlines = mollines[4 + nats:4 + nats + nbonds]
    bonds = {'BI': [], 'BJ': [], 'RIJ': [], 'UID': []}
    print(bondlines)

    for line in bondlines:
        [bi, bj] = map(int, [line[0:3], line[3:6]])
        bonds['BI'].append(bi)
        bonds['BJ'].append(bj)
        bonds['RIJ'].append(al.distance(coord[bi], coord[bj]))
        bonds['UID'].append(al.pairing_func(bi, bj))

    return coord, atsb, bonds


def ZMAT(file, opt, res, charge):
    name, ext = file.split('/')[-1].split('.')

    # create a MOL file from any file
    os.system(
        'babel -i{} {} -omol {}.mol --title {} ---errorlevel 1 -b &>LL'.format(
            ext, file, name, name
        )
    )

    # COOS, ATYPES, MolBonds = load_mol(f'{name}.mol')
    # G_mol, mol_icords = make_graphs(ATYPES, COOS, MolBonds)
    if ext.lower() == 'pdb':
        dfatoms, bonds = files.read_pdb(file)

    connect = connectivity(dfatoms, bonds)
    print(connect)

    mol_icords = connect.get_geometry()

    print(mol_icords)

    files.save_zmat(dfatoms, connect, mol_icords, res)


class connectivity(nx.DiGraph):
    """Building a class connectivity from directed graphs."""

    def __init__(self, coord, bonds):
        """Build connectivity from coordinates using nodes like atoms."""

        super().__init__()

        # Add nodes using atoms andsymbols and coordinates
        for i in coord.index:
            self.add_node(
                i,
                xyz=coord.loc[i, ['x', 'y', 'z']].values,
                atsb=coord.loc[i, 'atsb'],
                mass=coord.loc[i, 'mass'],
                atnum=coord.loc[i, 'atnum']
            )

        # Add edges like bonds
        for i, j in bonds:
            self.add_edge(i, j, distance=bonds[(i, j)][0], UID=bonds[(i, j)][1])
            self.add_edge(j, i, distance=bonds[(i, j)][0], UID=bonds[(i, j)][1])

    def nbonds(self, inode):
        """Return number of atoms in connected to iat."""
        return int(self.degree[inode] / 2)

    def get_geometry(self):
        all_ps = dict(nx.algorithms.all_pairs_shortest_path_length(self))
        print(all_ps)

        print(list(self.edges))
        print(list(self.nodes))

        all_paths = []
        for s in all_ps:
            for e in all_ps[s]:
                if all_ps[s][e] == 1:
                    all_paths += list(nx.algorithms.all_simple_paths(self, s, e, cutoff=1))
                elif all_ps[s][e] == 2:
                    all_paths += list(nx.algorithms.all_simple_paths(self, s, e, cutoff=2))
                elif all_ps[s][e] == 3:
                    all_paths += list(nx.algorithms.all_simple_paths(self, s, e, cutoff=3))

        print(all_paths)

        all_bonds = [p for p in all_paths if len(set(p)) == 2]
        new_angs = [p for p in all_paths if len(set(p)) == 3]
        new_tors = [p for p in all_paths if len(set(p)) == 4]

        print('Bonds', all_bonds)
        print('Angles', new_angs)
        print('Dihedrals', new_tors)

        dict_new_bds = {al.pairing_func(t[0], t[1]): t for t in all_bonds}
        dict_new_angs = {al.ang_id(t): t for t in new_angs}
        dict_new_tors = {al.tor_id(t): t for t in new_tors}

        print(dict_new_bds)
        print(dict_new_angs)
        print(dict_new_tors)

        imp_keys = [n for n in self.nodes() if self.nbonds(n) == 3]
        print(imp_keys)

        all_imps = {}
        for i in imp_keys:
            nei = list(self.neighbors(i))
            if self.nodes[i]['atnum'] == 6:
                all_imps[i] = [nei[0], i, nei[1], nei[2]]

        print(all_imps)

        MOL_ICOORDS = {
            'BONDS': dict_new_bds,
            'ANGLES': dict_new_angs,
            'TORSIONS': dict_new_tors,
            'IMPROPERS': all_imps
        }

        return MOL_ICOORDS


def get_OPT(zmat, opt, charge):
    """
    Do optimization

    xZCM1A - convert BOSS Z-matrix file to optimized Z-matrix.
    xSP - single-point force field calculation.
    xSPM - xSP, but prepare output plt.pdb and sum (Z-matrix).
    xOPT - optimization using the OPLS-AA force field with BFGS optimizer.

    """
    assert os.path.isfile(zmat), 'File named %10s does not exist' % zmat
    assert 'BOSSdir' in os.environ, 'Please Make sure $BOSSdir is defined \n xZCM1A and related files are in scripts directory of BOSS'
    # 1) Optimized Z-matrix
    execs = {
        # 2: os.environ['BOSSdir'] + '/scripts/xZCM1A+2 > /tmp/olog',
        # -2: os.environ['BOSSdir'] + '/scripts/xZCM1A-2 > /tmp/olog',
        0: os.environ['BOSSdir'] + '/scripts/xZCM1A > /tmp/olog',
        # 1: os.environ['BOSSdir'] + '/scripts/xZCM1A+  > /tmp/olog',
        # -1: os.environ['BOSSdir'] + '/scripts/xZCM1A-  > /tmp/olog',
    }

    print('MOLECULE HAS A CHARGE of %d' % charge)
    execfile = execs[charge]
    commad = execfile + ' ' + zmat[:-2]
    print(commad)
    os.system(commad)

    os.system('cp sum %s' % (zmat))
    # 2) Single-point, save the files plt.pdb and sum (Z-matrix)
    execfile = os.environ['BOSSdir'] + '/scripts/xSPM > /tmp/olog'
    commad = execfile + ' ' + zmat[:-2]
    print(commad)
    os.system(commad)
    os.system('cp sum %s' % (zmat))

    # 3) Optimization
    if opt > 0:
        print('Optimization level requested %d' % opt)
        for opt_lev in range(opt):
            print('Performing Stage %d of Charge Generation' % (opt_lev + 1))
            execfile = execs[charge]
            commad = execfile + ' ' + zmat[:-2]
            print(commad)
            os.system(commad)
            os.system('cp sum %s' % (zmat))

            # Optimization etape
            execfile = os.environ['BOSSdir'] + '/scripts/xOPT > /tmp/olog'
            commad = execfile + ' ' + zmat[:-2]
            os.system(commad)
            os.system('cp sum %s' % (zmat))

        # 4) final single-point
        execfile = os.environ['BOSSdir'] + '/scripts/xSPM > /tmp/olog'
        commad = execfile + ' ' + zmat[:-2]
        os.system(commad)
        os.system('cp sum %s' % (zmat))

import re
import pandas as pd
import numpy as np
import itertools as it
from numpy import linalg
from gentopol import algebra as al
import collections


Elements = {
    'H': {'mass': 1.0080, 'num': 1},
    'C': {'mass': 12.011, 'num': 6},
    'N': {'mass': 14.007, 'num': 7},
    'O': {'mass': 15.999, 'num': 8},
    'S': {'mass': 32.065, 'num': 16},
    'CL': {'mass': 35.453, 'num': 17}
}


def set_atsb(x):
    """Transform atoms with numbers like C00 to just the atomic symbol."""
    if len(x) == 1:
        return x
    else:
        return ''.join([i for i in x if not i.isdigit()])


def DBE(**kwargs):
    """Return the DBE (Double Bond Equivalent)."""
    return int((2 * kwargs['C'] + kwargs['N'] + 2 - kwargs['H']) / 2)


def Formula(**kwargs):
    """Return the molecular formule."""
    formule = ''
    for at in ['C', 'H', 'N', 'O', 'S', 'CL']:
        try:
            n = kwargs[at]
            if n == 1:
                formule += at
            else:
                formule += f'{at}{n}'
        except KeyError:
            pass

    return formule


def mu(table):
    e = 1.602176634e-19
    coord = table.loc[:, ['x', 'y', 'z']].values
    charges = table.loc[:, 'charge'].values

    q_xyz = (coord.T * charges).T * (4.8 / 1.6e-29) * e * 1e-10  # in angtroms 1e-10, in nanometer 1e-9
    q_xyz = np.sum(q_xyz, axis=0)

    return np.linalg.norm(q_xyz)


def Molecular(table):
    """Return relevant information for the structure molecular."""

    table['atsb'] = table['atsb'].apply(set_atsb)
    atoms = table.drop(['x', 'y', 'z', 'charge', 'mass'], axis=1)
    atoms = atoms.groupby('atsb').count()
    atoms = atoms.to_dict()['atnum']

    for at in ['C', 'N', 'H']:
        if at not in atoms:
            atoms[at] = 0

    info = {
        'DBE': DBE(**atoms),
        'formula': Formula(**atoms),
        'MM': table['mass'].sum(),
        'dipolar': mu(table)
    }
    return info


def read_pdb(file):
    """ Obtains the geometry of the system from pdb file."""
    coord = re.compile(r"""
            \w+\s+
            (?P<atid>\d+)\s+           # Atom id.
            (?P<atsb>\w+)\s+           # Atomic number.
            \w+\s+
            \d+\s+
            (?P<x>[+-]?\d+\.\d+)\s+    # Orthogonal coordinates for X.
            (?P<y>[+-]?\d+\.\d+)\s+    # Orthogonal coordinates for Y.
            (?P<z>[+-]?\d+\.\d+)\s+    # Orthogonal coordinates for Z.
            """, re.X)

    data = list()
    ndx_conect = dict()
    mass = list()
    atnum = list()
    symbols = list()
    with open(file, 'r') as INPUT:
        for line in INPUT:
            if coord.match(line):
                m = coord.match(line)
                data.append(m.groupdict())
                """ Adding mass """
                atsb = m.group('atsb')
                # atsb = ''.join([i for i in m.group('atsb') if not i.isdigit()])
                if 'CL' in atsb or 'Cl' in atsb:
                    atsb = 'CL'
                elif len(atsb) == 3:
                    atsb = re.sub(r'\w\w\Z', '', atsb)

                elif len(atsb) == 2:
                    atsb = re.sub(r'\w\Z', '', atsb)

                try:
                    mass.append(Elements[atsb]['mass'])
                    atnum.append(Elements[atsb]['num'])
                    symbols.append(atsb)
                except KeyError:
                    print(m.group('atsb'))
                    print(f'atom symbol {atsb} != not found')
                    exit()

            if "CONECT" in line:
                """ ndx_conect"""
                line = line.split()
                if len(line) > 2:
                    ndx_conect[int(line[1])] = [int(i) for i in line[2:]]

    dfatoms = pd.DataFrame(data)
    dfatoms['mass'] = mass
    dfatoms['atnum'] = atnum
    dfatoms['atsb'] = symbols
    dfatoms = dfatoms.astype({
        'atid': np.int64,
        'atnum': np.int64,
        'mass': np.float64,
        'x': np.float64,
        'y': np.float64,
        'z': np.float64})
    dfatoms = dfatoms.set_index('atid')
    # print(dfatoms)

    """ Listing bonds """
    # (BI, BJ) : RIJ, UID
    bonds = dict()
    for iat, jat in it.permutations([i for i in ndx_conect], 2):
        if iat < jat:
            if jat in ndx_conect[iat]:
                v = dfatoms.loc[iat, ['x', 'y', 'z']].values
                u = dfatoms.loc[jat, ['x', 'y', 'z']].values
                bonds[(iat, jat)] = list()
                bonds[(iat, jat)].append(np.linalg.norm(u - v))
                bonds[(iat, jat)].append(al.pairing_func(iat, jat))

    return dfatoms, bonds


def get_add_int(mol_icords, Z_BONDS, Z_ANGLES, Z_TORSIONS):
    all_bonds_mol, all_angles_mol, all_torsions_mol = mol_icords['BONDS'], mol_icords['ANGLES'], mol_icords['TORSIONS']
    # Bonds list
    Z_B = {al.pairing_func(i[0] - 2, i[1] - 2): [i[0] - 2, i[1] - 2] for i in Z_BONDS.values()}

    Z_A = {al.ang_id([i[0] - 2, i[1] - 2, i[2] - 2]): [i[0] - 2, i[1] - 2, i[2] - 2] for i in Z_ANGLES.values()}

    Z_T = {al.tor_id([i[0] - 2, i[1] - 2, i[2] - 2, i[3] - 2]): [i[0] - 2, i[1] - 2, i[2] - 2, i[3] - 2] for i in Z_TORSIONS.values()}

    Z_Ad_B, Z_Ad_A, Z_Ad_T = collections.OrderedDict(), collections.OrderedDict(), collections.OrderedDict()

    for b_ij in all_bonds_mol:
        # uid_b_ij = al.pairing_func(b_ij[0], b_ij[1])
        if b_ij not in list(Z_B.keys()):
            Z_Ad_B[b_ij] = [i + 2 for i in all_bonds_mol[b_ij]]

    for a_ij in all_angles_mol.keys():
        if a_ij not in list(Z_A.keys()):
            Z_Ad_A[a_ij] = [i + 2 for i in all_angles_mol[a_ij]]

    for t_ij in all_torsions_mol.keys():
        if t_ij not in list(Z_T.keys()):
            Z_Ad_T[t_ij] = [i + 2 for i in all_torsions_mol[t_ij]]

    for c in mol_icords['IMPROPERS'].values():
        Z_Ad_T["-".join(list(map(str, c)))] = [i + 2 for i in c]

    return Z_Ad_B, Z_Ad_A, Z_Ad_T


def save_zmat(df, connect, mol_icords, res):

    Z_ATOMS = {1: 'X', 2: 'X'}

    Z_NO = {1: -1, 2: -1}

    Z_BONDS = {1: (1, 0, 0.000), 2: (2, 1, 1.00), 3: (3, 2, 1.00)}

    Z_ANGLES = {1: (1, 0, 0, 0.000), 2: (2, 1, 0, 0.000),
                3: (3, 2, 1, 90.00), 4: (4, 3, 2, 90.0)}

    Z_TORSIONS = {1: (1, 0, 0, 0, 0.00), 2: (2, 1, 0, 0, 0.00), 3: (
        3, 2, 1, 0, 0.00), 4: (4, 3, 2, 1, 0.00), 5: (5, 4, 3, 2, 90.0)}

    atoms = df['atsb']

    for i in range(1, len(atoms) + 1):
        Z_ATOMS[i + 2] = atoms[i]

    for i in range(1, len(atoms) + 1):
        Z_NO[i + 2] = connect.nodes[i]['atnum']

    n_ats = 0
    B_LINK = {}
    for i in connect.nodes():
        if n_ats > 0:
            neigs = np.sort(list(connect.neighbors(i)))
            B_LINK[i] = neigs[0]
            Z_BONDS[i + 2] = (i + 2, neigs[0] + 2, connect[i][neigs[0]]['distance'])
        n_ats += 1

    n_ats = 0
    A_LINK = {}
    for i in connect.nodes():
        if n_ats > 1:
            neigs = np.sort(list(connect.neighbors(B_LINK[i])))
            A_LINK[i] = neigs[0]
            p0 = connect.nodes[i]['xyz']
            p1 = connect.nodes[B_LINK[i]]['xyz']
            p2 = connect.nodes[neigs[0]]['xyz']
            # ang = al.angle(coos[i], coos[B_LINK[i]], coos[neigs[0]])
            v0 = p0 - p1
            v1 = p2 - p1

            cos_a = np.dot(v0, v1) / linalg.norm(v0) / linalg.norm(v1)
            ang = np.arccos(round(cos_a, 3)) * 180.0 / np.pi

            Z_ANGLES[i + 2] = (i + 2, B_LINK[i] + 2, neigs[0] + 2, ang)
        n_ats += 1

    n_ats = 0
    for i in connect.nodes():
        if n_ats > 2:
            neigs = list(connect.neighbors(A_LINK[i]))
            neigs = np.array([j for j in neigs if j not in [i, B_LINK[i], A_LINK[i]]])
            neigs = np.sort(neigs)
            neigs = neigs[neigs < i]

            if len(neigs) < 1:
                neigs = [j for j in list(connect.neighbors(B_LINK[i])) if j not in [i, A_LINK[i]]]
                if (B_LINK[i] in list(mol_icords['IMPROPERS'].keys())):
                    del mol_icords['IMPROPERS'][B_LINK[i]]

            ti, tj, tk, tl = [i, B_LINK[i], A_LINK[i], neigs[0]]

            p0 = np.array(connect.nodes[i]['xyz'], dtype=float)
            p1 = np.array(connect.nodes[B_LINK[i]]['xyz'], dtype=float)
            p2 = np.array(connect.nodes[A_LINK[i]]['xyz'], dtype=float)
            p3 = np.array(connect.nodes[neigs[0]]['xyz'], dtype=float)

            v01 = p0 - p1
            v32 = p3 - p2
            v12 = p1 - p2

            v0 = np.cross(v12, v01)
            v3 = np.cross(v12, v32)

            cos_a = round(np.dot(v0, v3) / linalg.norm(v0) / linalg.norm(v3), 3)
            ang = np.arccos(cos_a) * 180.0 / np.pi

            # The cross product vectors are both normal to the axis
            # vector v12, so the angle between them is the dihedral
            # angle that we are looking for.  However, since "angle"
            # only returns values between 0 and pi, we need to make
            # sure we get the right sign relative to the rotation axis

            if np.dot(np.cross(v0, v3), v12) > 0:
                ang *= -1

            Z_TORSIONS[i + 2] = (ti + 2, tj + 2, tk + 2, tl + 2, ang)
        n_ats += 1

    Z_Ad_B, Z_Ad_A, Z_Ad_T = get_add_int(mol_icords, Z_BONDS, Z_ANGLES, Z_TORSIONS)

    # PRINTING ACTUAL Z-MATRIX
    ofile = open(f'{res}.z', 'w')
    ofile.write('BOSS Z-Matrix with LSDautozmat (written by Leela S. Dodda)\n')
    for i in range(1, len(atoms) + 3):
        ofile.write('%4d %-3s%5d%5d%5d%12.6f%4d%12.6f%4d%12.6f%4s%5d\n'
                    % (i, Z_ATOMS[i], Z_NO[i], Z_NO[i], Z_BONDS[i][1], Z_BONDS[i][-1], Z_ANGLES[i][-2], Z_ANGLES[i][-1], Z_TORSIONS[i][-2], Z_TORSIONS[i][-1], res[0:3], 1)
                    )
    ofile.write('''                    Geometry Variations follow    (2I4,F12.6)
                    Variable Bonds follow         (I4)\n''')

    for i in range(4, len(atoms) + 3):
        ofile.write('%4d\n' % i)

    ofile.write('                    Additional Bonds follow       (2I4)\n')
    if len(Z_Ad_B) > 0:
        for i in Z_Ad_B.values():
            ofile.write('%4d%4d\n' % (i[0], i[1]))

    # CREATE A FUNCTION TO DEFINE ADDITIONAL BONDS IN CASE OF RINGS
    ofile.write('''                    Harmonic Constraints follow   (2I4,4F10.4)
                    Variable Bond Angles follow   (I4)\n''')
    for i in range(5, len(atoms) + 3):
        ofile.write('%4d\n' % i)

    ofile.write('                    Additional Bond Angles follow (3I4)\n')
    if len(Z_Ad_A) > 0:
        for i in Z_Ad_A.values():
            ofile.write('%4d%4d%4d\n' % (i[0], i[1], i[2]))

    # CREATE A FUNCTION TO DEFINE ADDITIONAL BONDS IN CASE OF RINGS
    ofile.write(
        '                    Variable Dihedrals follow     (3I4,F12.6)\n')
    for i in range(6, len(atoms) + 3):
        ofile.write('%4d%4d%4d%12.6f\n' % (i, -1, -1, 0.000))
    ofile.write('                    Additional Dihedrals follow   (6I4)\n')
    if len(Z_Ad_T) > 0:
        for k in Z_Ad_T.keys():
            torsion = Z_Ad_T[k]
            ofile.write('%4d%4d%4d%4d%4d%4d\n' %
                        (torsion[0], torsion[1], torsion[2], torsion[3], -1, -1))

    ofile.write(
        '''                    Domain Definitions follow     (4I4)
                    Conformational Search (2I4,2F12.6)
                    Local Heating Residues follow (I4 or I4-I4)
                    Final blank line
''')
    ofile.close()


def save_gro(table, res, title='GRO FILE'):
    """ Save coordinate to file *.gro from dataframe with x, y, z """
    nat = len(table)
    gro = res.lower() + '.gro'
    GRO = open(gro, 'w', encoding='utf-8')
    GRO.write('%s\n' % title)
    GRO.write('%5d\n' % nat)
    for i in table.index:
        GRO.write('{:>8}{:>7}{:5d}{:8.3f}{:8.3f}{:8.3f}\n'.format(
            '1'+res.upper(),
            table.loc[i, 'atsb'],
            i,
            table.loc[i, 'x'] * 0.1,
            table.loc[i, 'y'] * 0.1,
            table.loc[i, 'z'] * 0.1))
    GRO.write('   0.00000   0.00000   0.00000\n')
    GRO.close()
    print('Saved gro file: ' + gro)

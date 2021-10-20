"""
SCRIPT TO WRITE GROMACS GRO AND ITP FILES
FROM BOSS ZMATRIX
Created on Mon Feb 15 15:40:05 2016
@author: Leela S. Dodda leela.dodda@yale.edu
@author: William L. Jorgensen Lab

Modificacionres realizadas por Orlando Villegas
2021

"""

import pickle
import pandas as pd
import numpy as np
from gentopol.BOSStools import bossPdbAtom2Element, bossElement2Mass, pairing_func, ucomb
from gentopol.files import read_pdb, save_gro


def bossData(molecule_data):
    """Extrae la data referente al tipo de atomos y los parametros de VDW."""
    ats_file = molecule_data.MolData['ATOMS']
    types = []
    for i in enumerate(ats_file):
        types.append([i[1].split()[1], 'opls_' + i[1].split()[2]])

    print(types)
    st_no = 3
    Qs = molecule_data.MolData['Q_LJ']
    print(Qs)
    assert len(Qs) == len(types), 'Please check the at_info and Q_LJ_dat files'

    num2opls = {}
    for i in range(0, len(types)):
        num2opls[i] = Qs[i][0]
    print(num2opls)

    # {i: ['AT0i', 'opls_80i', 'AT80i', 'AT', 'mass', 'atype']}
    num2typ2symb = {i: types[i] for i in range(len(Qs))}
    for i in range(len(Qs)):
        num2typ2symb[i].append(
            bossPdbAtom2Element(num2typ2symb[i][0]) + num2typ2symb[i][1][-3:]
        )
        num2typ2symb[i].append(
            bossPdbAtom2Element(num2typ2symb[i][0])
        )
        num2typ2symb[i].append(
            bossElement2Mass(num2typ2symb[i][3])
        )
        num2typ2symb[i].append(Qs[i][0])

    print(num2typ2symb)

    return types, Qs, num2opls, st_no, num2typ2symb


def boss2gmxBond(mol, st_no):
    bdat = mol.MolData['BONDS']
    bdat['cl1'] = [x - st_no if not x - st_no < 0 else 0 for x in bdat['cl1']]
    bdat['cl2'] = [x - st_no if not x - st_no < 0 else 0 for x in bdat['cl2']]
    dfbond = pd.DataFrame(bdat)
    dfbond['KIJ'] = dfbond['KIJ'] * 836.80
    dfbond['RIJ'] = dfbond['RIJ'] * 0.10
    # Revisar
    dfbond['UNQ'] = [pairing_func(i + 1, j + 1)[0]
                     for i, j in zip(dfbond.cl1, dfbond.cl2)]

    dfbond['UF'] = ((dfbond.cl1 + dfbond.cl2) *
                    (dfbond.cl1 + dfbond.cl2 + 1) * 0.5) + dfbond.cl2

    dfbond['UR'] = ((dfbond.cl1 + dfbond.cl2) *
                    (dfbond.cl1 + dfbond.cl2 + 1) * 0.5) + dfbond.cl1
    # -------

    print(dfbond)
    # connects = []
    # for ai, aj in zip(dfbond.cl1, dfbond.cl2):
    #     # ai   aj    func
    #     connects.append("{:>5}{:>5}    1".format(aj + 1, ai + 1))
    # full_bnd = bnd_df.copy()
    # bnd_df = bnd_df.drop_duplicates(['TIJ'])
    # print(dfbond.drop_duplicates(['TIJ']))
    return dfbond


def boss2gmxAngle(mol, num2opls, st_no):
    adat = mol.MolData['ANGLES']
    adat['cl1'] = [x - st_no if not x - st_no < 0 else 0 for x in adat['cl1']]
    adat['cl2'] = [x - st_no if not x - st_no < 0 else 0 for x in adat['cl2']]
    adat['cl3'] = [x - st_no if not x - st_no < 0 else 0 for x in adat['cl3']]
    dfangle = pd.DataFrame(adat)
    dfangle = dfangle[dfangle.K > 0]
    dfangle['K'] = 8.3680 * dfangle['K']
    dfangle['TY'] = np.array([num2opls[i] + '-' + num2opls[j] + '-' + num2opls[k]
                             for i, j, k in zip(dfangle.cl1, dfangle.cl2, dfangle.cl3)])
    # ang_df = ang_df.drop_duplicates(['TY'])
    # full_ang = ang_df.copy()
    print(dfangle)

    return dfangle


def boss2opmTorsion(mol, num2opls, dfbond, st_no):
    dhd = []
    for line in mol.MolData['TORSIONS']:
        dt = [float(f) for f in line]
        dhd.append(dt)

    dhd = np.array(dhd)
    dhd = dhd * 4.184  # kcal to kj conversion
    dfdih = pd.DataFrame(dhd, columns=['V1', 'V2', 'V3', 'V4'])

    ats = []
    for line in mol.MolData['ATOMS'][3:]:
        dt = [line.split()[0], line.split()[4],
              line.split()[6], line.split()[8]]
        dt = [int(d) for d in dt]
        ats.append(dt)

    for line in mol.MolData['ADD_DIHED']:
        dt = [int(l) for l in line]
        ats.append(dt)

    assert len(ats) == len(
        dhd), 'Number of Dihedral angles in Zmatrix and Out file dont match'

    ats = np.array(ats) - st_no
    for i in range(len(ats)):
        for j in range(len(ats[0])):
            if ats[i][j] < 0:
                ats[i][j] = 0

    dfat = pd.DataFrame(ats, columns=['I', 'J', 'K', 'L'])
    dfdih = pd.concat([dfdih, dfat], axis=1).reindex()

    bndlist = list(dfbond.UR) + (list(dfbond.UR))
    print(bndlist)

    # final_df['TY'] = ['Proper' if ucomb(list([final_df.I[n], final_df.J[n], final_df.K[
    #     n], final_df.L[n]]), bndlist) == 3 else 'Improper' for n in range(len(final_df.I))]

    TY = []
    for n in range(len(dfdih.I)):
        if ucomb(list([dfdih.I[n], dfdih.J[n], dfdih.K[n], dfdih.L[n]]), bndlist) == 3:
            TY.append('Proper')
        else:
            TY.append('Improper')

    dfdih['TY'] = TY
    dfdih['SumV'] = np.abs(dfdih.V1) + np.abs(dfdih.V2) + np.abs(dfdih.V3) + np.abs(dfdih.V4)

    dfdih['TI'] = [num2opls[j] for j in dfdih.I]
    dfdih['TJ'] = [num2opls[j] for j in dfdih.J]
    dfdih['TK'] = [num2opls[j] for j in dfdih.K]
    dfdih['TL'] = [num2opls[j] for j in dfdih.L]

    if len(dfdih.index) > 0:
        dfdih['NAME'] = dfdih.TI + '-' + dfdih.TJ + '-' + dfdih.TK + '-' + dfdih.TL
        dfdih = dfdih.sort_values(['NAME'])
        tor_bos = dfdih.drop(
            labels=['I', 'J', 'K', 'L', 'TI', 'TJ', 'TK', 'TL'],
            axis=1
        )
        tor_bos = tor_bos.drop_duplicates()
        df = dfdih.loc[tor_bos.index, :]
        print(dfdih)
        print(df)
        return dfdih, df
    else:
        return dfdih, dfdih


def gmxImp(df):
    odihed = []
    for pot in range(1, 5):
        if df['V' + str(pot)] != 0.00:
            odihed.append('%6d%6d%6d%6d    4     %10.3f %10.3f %5d  \n' % (
                df['I'] + 1, df['J'] + 1, df['K'] +
                1, df['L'] + 1, 180.00 * abs(pot % 2 - 1),
                float(df['V' + str(pot)]) * 0.5, pot))
    return odihed


def gmxDihed(df):
    '''
        PRINTING RYCKET-BELLMAN STYLE PARAMETERS BECAUSE FEP DOES NOT SUPPORT PROPER DIHEDRALS
    '''
    [f1, f2, f3, f4] = [df.V1, df.V2, df.V3, df.V4]
    cdat = [f2 + (f1 + f3) * 0.5, 1.5 * f3 - 0.5 * f1, 4.0 *
            f4 - f2, -2.0 * f3, -4.0 * f4, 0.0, 0.00]
    return '%5s%5s%5s%5s        3      %6.3f  %6.3f  %6.3f  %6.3f  %6.3f  %6.3f\n' % (
        df['I'] + 1, df['J'] + 1, df['K'] + 1, df['L'] + 1, cdat[0], cdat[1], cdat[2], cdat[3], cdat[4], cdat[5])


def gmxPairs(tors, ang, bond):
    atom_list = np.array(list(bond.cl1) + list(bond.cl2))
    atom_list = np.sort(atom_list)
    atom_list = np.unique(atom_list)
    dict_bond = {ano: list(bond[bond.cl1 == ano]['cl2']) +
                 list(bond[bond.cl2 == ano]['cl1']) for ano in atom_list}
    dict_ang = {}
    for ano in atom_list:
        tr_1 = []
        for bno in dict_bond[ano]:
            tr_1 += dict_bond[bno]
        tr_1 = list(set(tr_1))
        tr_1.remove(ano)
        dict_ang[ano] = tr_1
    NP_dat = []
    for a in dict_bond.keys():
        for b in dict_bond[a]:
            NP_dat.append([a, b, 1, pairing_func(a + 1, b + 1)
                           [0], pairing_func(a + 1, b + 1)[1]])
            for c in dict_ang[b]:
                NP_dat.append([a, c, 3, pairing_func(a + 1, c + 1)
                               [0], pairing_func(a + 1, c + 1)[1]])
    for a in dict_ang.keys():
        for b in dict_ang[a]:
            NP_dat.append([a, b, 2, pairing_func(a + 1, b + 1)
                           [0], pairing_func(a + 1, b + 1)[1]])
            for c in dict_bond[b]:
                NP_dat.append([a, c, 3, pairing_func(a + 1, c + 1)
                               [0], pairing_func(a + 1, c + 1)[1]])
    NP_df = pd.DataFrame(NP_dat, columns=['I', 'J', 'BSEP', 'UNQ', 'TY'])
    NP_df = NP_df[NP_df.I != NP_df.J]  # CASES WITH 3 MEMBERED RINGS
    NP_df = NP_df.sort_values(['UNQ'])
    NP_B = NP_df[NP_df.BSEP == 1]
    NP_A = NP_df[NP_df.BSEP == 2]
    NP_T = NP_df[NP_df.BSEP == 3]
    NP_T = NP_T.drop_duplicates(['UNQ'])
    NP_T = NP_T[~ NP_T.UNQ.isin(NP_B.UNQ)]
    NP_T = NP_T[~ NP_T.UNQ.isin(NP_A.UNQ)]

    return list(NP_T.TY)


def save2GMX(res):
    mol = pickle.load(open(res + ".p", "rb"))
    pdb_file = 'plt.pdb'

    types, Qs, num2opls, st_no, num2typ2symb = bossData(mol)

    # Saving ITP file.
    lines = ''
    lines += '; Writted by Orlando Villegas - orlando.villegas@univ-pau.fr\n'
    lines += '; Modifications to the LigParGen script\n'
    lines += '; Writted by Leela S. Dodda leela.dodda@yale.edu\n'
    lines += '; from Jorgensen Lab @ Yale University\n'
    lines += ';\n'

    # Atoms types.
    lines += '\n[ atomtypes ]\n'
    for i in range(len(Qs)):
        lines += '{:>10} {:>5} {:10.4f}     0.000    A    {:10.5E}\n'.format(
            num2typ2symb[i][1],
            num2typ2symb[i][2],
            num2typ2symb[i][4],
            float(Qs[i][2]) * 0.1,
            float(Qs[i][3]) * 4.184
            )

    # Molecule type.
    lines += '\n[ moleculetype ]\n'
    lines += '; Name               nrexcl\n'
    lines += '%s                   3\n' % res.upper()

    # atoms.
    lines += '\n[ atoms ]\n'
    lines += ';   nr       type  resnr residue  atom   cgnr     charge       mass  \n'
    for i in range(len(Qs)):
        lines += ' {:5d} {:>10} {:6d} {:>6} {:>5} {:6d} {:>10} {:10.4f} \n'.format(
            i + 1,
            num2typ2symb[i][1],
            1,
            res.upper(),
            num2typ2symb[i][0],
            int(((i + 1) / 33.0) + 1),  # revisar
            Qs[i][1],
            num2typ2symb[i][4]
        )

    # bonds
    dfbond = boss2gmxBond(mol, st_no)
    lines += '\n[ bonds ]\n'

    for i, bond in dfbond.iterrows():
        if bond['cl1'] < bond['cl2']:
            ai = bond['cl1'] + 1
            aj = bond['cl2'] + 1
        else:
            ai = bond['cl2'] + 1
            aj = bond['cl1'] + 1

        lines += ' {:5d} {:5d} {:5d} {:11.4f} {:10.3f} \n'.format(
            ai,
            aj,
            1,
            bond['RIJ'],
            bond['KIJ']
        )

    # angles
    dfangle = boss2gmxAngle(mol, num2opls, st_no)
    lines += '\n[ angles ]\n'
    for i, ang in dfangle.iterrows():
        if ang['cl1'] < ang['cl3']:
            ai = ang['cl1'] + 1
            aj = ang['cl2'] + 1
            ak = ang['cl3'] + 1
        else:
            ai = ang['cl3'] + 1
            aj = ang['cl2'] + 1
            ak = ang['cl1'] + 1

        lines += ' {:5d} {:5d} {:5d} {:5d} {:10.3f} {:10.3f} \n'.format(
            ai,
            aj,
            ak,
            1,
            ang['R'],
            ang['K']
        )

    # dihedrals
    dfdih, tor_df = boss2opmTorsion(mol, num2opls, dfbond, st_no)
    if len(tor_df.index) != len(dfdih.index):
        lines += '\n[ dihedrals ]\n'
        lines += '; IMPROPER DIHEDRAL ANGLES \n'
        lines += ';  ai    aj    ak    al funct            c0            c1            c2            c3            c4            c5\n'

        for i, improper in dfdih[dfdih.TY == 'Improper'].iterrows():
            for st in gmxImp(improper):
                lines += '%s' % st

        lines += '\n[ dihedrals ]\n'
        lines += '; PROPER DIHEDRAL ANGLES\n'
        lines += ';  ai    aj    ak    al funct            c0            c1            c2            c3            c4            c5\n'

        for i, dihedral in dfdih[dfdih.TY == 'Proper'].iterrows():
            for st in gmxDihed(dihedral):
                lines += '%s' % st

    # pair
    lines += '\n[ pairs ]\n'
    ppairs = gmxPairs(dfdih, dfangle, dfbond)
    for pair in ppairs:
        lines += '%s    1\n' % pair

    dfatoms, _ = read_pdb(pdb_file)
    print(dfatoms)

    # writing all
    with open("%s.itp" % res.lower(), "w") as f:
        f.write(lines)

    print("Saved itp file: %s.itp" % res.lower())

    save_gro(dfatoms, res)

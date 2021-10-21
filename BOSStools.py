import pandas as pd
import numpy as np
import os


def Refine_file(fname):
    flines = open(fname, 'r')
    lines = []
    for line in flines:
        if line.rstrip():
            line = line.rstrip()
            line = line.lstrip()
            lines.append(line)
    flines.close()
    return lines


class BOSSReader:

    def __init__(self, zmatrix):
        self.zmat = zmatrix
        self.MolData = self.get_ImpDat()

    def get_ImpDat(self):
        """Import data from the optimization output files."""
        odat = Refine_file('./out')
        sdat = Refine_file(self.zmat)

        MolData = {}
        impDat = {}
        MolData['PDB'] = Refine_file('./plt.pdb')

        for nl in range(len(odat)):
            if 'Z-Matrix for Reference Solutes' in odat[nl]:
                impDat['ATMinit'] = nl
            elif 'Net Charge' in odat[nl]:
                impDat['TotalQ'] = nl
            elif 'OPLS Force Field Parameters' in odat[nl]:
                impDat['ATMfinal'] = nl
                impDat['NBDinit'] = nl
            elif 'Fourier Coefficients' in odat[nl]:
                impDat['TORinit'] = nl
                impDat['NBDfinal'] = nl
            elif 'Bond Stretching Parameters' in odat[nl]:
                impDat['TORfinal'] = nl
                impDat['BNDinit'] = nl
            elif 'Angle Bending Parameters' in odat[nl]:
                impDat['BNDfinal'] = nl
                impDat['ANGinit'] = nl
            elif 'Non-bonded Pairs List' in odat[nl]:
                impDat['ANGfinal'] = nl
                impDat['PAIRinit'] = nl
            elif 'Solute 0:   X          Y          Z' in odat[nl]:
                impDat['XYZinit'] = nl
            elif 'Atom I      Atom J      RIJ' in odat[nl]:
                impDat['XYZfinal'] = nl
            elif 'Checking' in odat[nl]:
                impDat['PAIRfinal'] = nl

        """THIS PART IS READ FROM SUM FILE"""
        for ml in range(len(sdat)):
            if 'Additional Dihedrals follow' in sdat[ml]:
                impDat['ADDinit'] = ml
            elif 'Domain Definitions follow' in sdat[ml]:
                impDat['ADDfinal'] = ml

        MolData['ATOMS'] = self.get_atinfo(
            odat[impDat['ATMinit']:impDat['ATMfinal']]
            )
        # print(MolData['ATOMS'])
        MolData['Q_LJ'] = self.get_QLJ(
            odat[impDat['NBDinit']:impDat['NBDfinal']]
            )
        # print(MolData['Q_LJ'])
        MolData['BONDS'] = self.get_bonds(
            odat[impDat['BNDinit']:impDat['BNDfinal']]
            )
        # print(MolData['BONDS'])
        MolData['ANGLES'] = self.get_angs(
            odat[impDat['ANGinit']:impDat['ANGfinal']]
            )
        # print(MolData['ANGLES'])
        MolData['TORSIONS'] = self.get_tors(
            odat[impDat['TORinit']:impDat['TORfinal']]
            )
        # print(MolData['TORSIONS'])
        MolData['ADD_DIHED'] = self.get_addihed(
            sdat[impDat['ADDinit']:impDat['ADDfinal']]
            )
        # print(MolData['ADD_DIHED'])
        MolData['XYZ'] = self.get_XYZ(
            odat[impDat['XYZinit']:impDat['XYZfinal']]
            )
        # print(MolData['XYZ'])
        MolData['PAIRS'] = self.get_pairs(
            odat[impDat['PAIRinit']:impDat['PAIRfinal']]
            )
        # print(MolData['PAIRS'])
        MolData['TotalQ'] = self.get_charge(
            odat[impDat['TotalQ']:impDat['TotalQ'] + 4])
        # print(MolData['TotalQ'])

        return MolData

    def get_atinfo(self, data):
        ats = []
        nat = 0
        for line in data:
            if line[0].isdigit() and float(line.split()[2]) > 1:
                ats.append(line)
                nat += 1
        return (ats)

    def get_QLJ(self, data):
        qlj = []
        nqlj = 0
        for line in data:
            if 'All Solutes' in line and line[0].isalpha():
                qlj.append([line.split()[0], line.split()[2],
                            line.split()[3], line.split()[4]])
                nqlj += 1
        return (qlj)

    def get_bonds(self, data):
        bnds = {'cl1': [], 'cl2': [], 'RIJ': [], 'KIJ': [], 'TIJ': []}
        nbnd = 0
        for line in data:
            if line[0].isdigit() and float(line.split()[3]) > 0:
                word = line.split()
                bnds['cl1'].append(int(word[0]))
                bnds['cl2'].append(int(word[1]))
                bnds['RIJ'].append(float(word[2]))
                bnds['KIJ'].append(float(word[3]))
                bnds['TIJ'].append(line[-5:])
                nbnd += 1
        return (bnds)

    def get_angs(self, data):
        angs = {'cl1': [], 'cl2': [], 'cl3': [], 'R': [], 'K': []}
        nang = 0
        for line in data:
            if line[0].isdigit() and float(line.split()[4]) > 0:
                word = line.split()
                angs['cl1'].append(int(word[0]))
                angs['cl2'].append(int(word[1]))
                angs['cl3'].append(int(word[2]))
                angs['R'].append(float(word[3]))
                angs['K'].append(float(word[4]))
                nang = nang + 1
            #        print 'Total No of Non-zero Angles in BOSS is %d' % (nang)
        return (angs)

    def get_tors(self, data):
        tors = []
        ntor = 0
        for line in data:
            if 'All Solutes' in line:
                tors.append(line.split()[4:8])
                for tor in line.split()[4:8]:
                    if abs(float(tor)) > 0.0:
                        ntor = ntor + 1
        return (tors)

    def get_addihed(self, data):
        add = []
        nadd = 0
        for line in data:
            if line[0].isdigit():
                add.append(line.split()[0:4])
                nadd = nadd + 1
        return (add)

    def get_XYZ(self, data):
        XYZ = {'atnum': [], 'X': [], 'Y': [], 'Z': [], 'atsb': []}
        for line in data:
            if line[0].isdigit() and len(line.split()) == 5:
                word = line.split()
                if int(word[0]) > 0:
                    XYZ['atnum'].append(int(word[0]))
                    XYZ['X'].append(float(word[1]))
                    XYZ['Y'].append(float(word[2]))
                    XYZ['Z'].append(float(word[3]))
                    XYZ['atsb'].append(word[4])
        XYZ = pd.DataFrame(XYZ)
        return XYZ

    def get_pairs(self, data):
        data = data[1:]
        plnos = []
        for i in range(0, len(data)):
            if 'Atom' in data[i]:
                plnos.append(i)
        plnos.append(len(data))
        pair_dat = {i: ' '.join(data[plnos[i]:plnos[i + 1]])
                    for i in range(len(plnos) - 1)}
        for nu in range(len(plnos) - 1):
            pair_dat[nu] = list(pair_dat[nu][10:].split())
            pair_dat[nu] = np.array([int(a) - 2 for a in pair_dat[nu]])
        pairs = []
        for k in pair_dat.keys():
            for j in pair_dat[k]:
                pairs.append('%6d%6d%6d\n' % (k - 1, j, 1))
        return pairs

    def get_charge(self, data):
        TotQ = {}
        for line in data[1:]:
            words = line.split()
            TotQ['-'.join(words[:-1])] = round(float(words[-1]), 3)
        return TotQ

    def cleanup(self):
        os.system('rm -v sum log olog out plt.pdb')


def Check_H(atoms):
    atype = [line.split()[1][0] for line in atoms]
    ans = False
    if ('H' in atype):
        ans = True
    return ans


def bossPdbAtom2Element(attype):
    elem = ''.join([i for i in attype[:-1] if not i.isdigit()])
    return elem


def bossElement2Mass(elem):
    symb2mass = {
        'H': 1.008,
        'F': 18.998403163,
        'Cl': 35.45,
        'Br': 79.904,
        'I': 126.90447,
        'O': 15.999,
        'S': 32.06,
        'N': 14.007,
        'P': 30.973761998,
        'C': 12.011,
        'Si': 28.085,
        'Na': 22.98976928,
        'SOD': 22.98976928,
        'K': 39.0983,
        'Mg': 24.305,
        'Ca': 40.078,
        'Mn': 54.938044,
        'Fe': 55.845,
        'Co': 58.933194,
        'Ni': 58.6934,
        'Cu': 63.546,
        'Zn': 65.38, }
    try:
        mass = symb2mass[elem]
    except NameError:
        print("Mass for atom %s is not available \n add it to symb2mass dictionary" % elem)
    return mass


def pairing_func(a, b):
    ans = (a + b) * (a + b + 1) * 0.5
    if a > b:
        ans = ans + a
        pans = '%6d%6d' % (b, a)
    else:
        ans = ans + b
        pans = '%6d%6d' % (a, b)
    return (int(ans), pans)


def ucomb(vec, blist):
    res = 0
    for a in vec:
        vec.remove(a)
        for b in vec:
            ans = (a + b) * (a + b + 1) * 0.5
            if (ans + a in blist) or (ans + b in blist):
                res = res + 1

    return res

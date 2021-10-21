

import numpy as np
from numpy import linalg


def distance(u, v):
    """Return length of a vector."""

    return linalg.norm(u - v)


def angle(p0, p1, p2):
    """Return angle [0..pi] between two vectors."""
    v0 = p0 - p1
    v1 = p2 - p1

    cos_a = np.dot(v0, v1) / linalg.norm(v0) / linalg.norm(v1)

    return 180.0 * np.arccos(round(cos_a, 3)) * 7.0 / 22.0


def pairing_func(a, b):
    """The pairing function adds an identifier for each bond."""
    ans = (a + b) * (a + b + 1) * 0.5
    if a > b:
        ans = ans + a
    else:
        ans = ans + b
    return int(ans)


def tor_id(a):
    bond = pairing_func(a[1], a[2])
    ends = pairing_func(a[0], a[3])
    return '%d-%d' % (bond, ends)


def ang_id(a):
    bond_a = pairing_func(a[0], a[1])
    bond_b = pairing_func(a[1], a[2])
    return pairing_func(bond_a, bond_b)


def Mol_angle(v0, v1):
    "Return angle [0..pi] between two vectors."
    cos_a = round(np.dot(v0, v1) / linalg.norm(v0) / linalg.norm(v1), 3)
    return np.arccos(cos_a)

r"""
   _____ ______ _   _ _______ ____  _____   ____  _
  / ____|  ____| \ | |__   __/ __ \|  __ \ / __ \| |
 | |  __| |__  |  \| |  | | | |  | | |__) | |  | | |
 | | |_ |  __| | . ` |  | | | |  | |  ___/| |  | | |
 | |__| | |____| |\  |  | | | |__| | |    | |__| | |____
  \_____|______|_| \_|  |_|  \____/|_|     \____/|______|


FF formats provided :
--------------------
GROMACS      .itp & .gro

Input Files supported :
--------------------
PDB



"""

import argparse
from .topol import TOPOL


def options():
    """Generate command line interface."""

    parser = argparse.ArgumentParser(
        prog="gentopol",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        usage="%(prog)s -f file -r RES -c 0 -o 0",
        epilog="Enjoy the program!",
        description=__doc__
    )

    parser.add_argument(
        "-f", "--file",
        help="Submit file",
        type=str,
        metavar='file')

    parser.add_argument(
        "-r", "--resname",
        help="Residue name, define three letters",
        type=str,
        default='UNL',
        metavar='RES')

    parser.add_argument(
        "-o", "--opt",
        help="Optimization or Single Point Calculation",
        type=int,
        choices=[0, 1, 2, 3],
        default=0,
        metavar='N')

    parser.add_argument(
        "-c", "--charge",
        type=int,
        choices=[0, -1, 1, -2, 2],
        help="0: Neutral <0: Anion >0: Cation. Charge 1.14*CM1A ",
        default=0,
        metavar='CHARGE')

    return vars(parser.parse_args())


def main():
    args = options()
    TOPOL(**args)


# RUN
main()

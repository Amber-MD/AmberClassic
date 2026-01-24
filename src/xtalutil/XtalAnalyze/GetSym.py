#! /usr/bin/env python
import ReadAmberFiles as RAF
import sys


def NAsu(pdbfile_name):
    pdbfile = RAF.pdb(pdbfile_name)
    smtry, tr = pdbfile.Get_SMTRY()
    n_asu = len(smtry)
    return n_asu

if __name__ == "__main__":
    n_asu = NAsu(sys.argv[1])
    print(n_asu)

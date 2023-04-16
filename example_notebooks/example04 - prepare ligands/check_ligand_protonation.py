#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Check the protonation state of the molecules before docking.
In chemistry, protonation (or hydronation) is the adding of a proton (or hydron, or hydrogen cation), (H+) to an atom, molecule, or ion, forming a conjugate acid. (The complementary process, when a proton is removed from a Brønsted–Lowry acid, is deprotonation.)
Usage:
    bash check_ligand_protonation.sh <input file> <output file>
Example:
    bash check_ligand_protonation.sh 1iep_ligand_G_STI.sdf ligand_step01_protonated.sdf
"""

import os
import sys


__author__ = "Carlos Garcia-Hernandez"
__email__ = "carlos.garcia2@bsc.es"

USAGE = __doc__.format(__author__, __email__)

def check_input(args):
    """Validates user input/options.
    """
    if not os.path.isfile(args):
            emsg = 'ERROR!! File not found or not readable: \'{}\'\n'
            sys.stderr.write(emsg.format(args[0]))
            sys.stderr.write(__doc__)
            sys.exit(1)
    return args

def run(infile, outfile):
    """
    Check the protonation state of the molecules before docking.
    Parameters
    ----------
    fhandle : an iterable giving the PDB file line-by-line
    outname : str
        The base name of the output files. If None is given, tries to
        extract a name from the `.name` attribute of `fhandler`. If
        `fhandler` has no attribute name, assigns `splitchains`.
    """
    from openbabel import openbabel as ob
    from meeko import obutils

    mol = obutils.load_molecule_from_file(infile, molecule_format='SDF')

    mol.AddHydrogens()
    charge_model = ob.OBChargeModel.FindType("Gasteiger")
    charge_model.ComputeCharges(mol)

    obutils.writeMolecule(mol, fname=outfile)


def main():
    # Check Input
    args = sys.argv[1:]
    infile = check_input(args[0])
    outfile = args[1]

    # Do the job
    run(infile, outfile)

    # last line of the script
    # We can close it even if it is sys.stdin
    # infile.close()
    # outfile.close()
    sys.exit(0)


if __name__ == '__main__':
    main()

# --------------------------------------------------------------------------
# references:
# add hydrogens to molecules using openbabel
# from: https://github.com/DrrDom/rdkit-scripts/blob/master/vina_dock.py
#       https://github.com/forlilab/Meeko/blob/09610111d43fef2b11636ddd110d1423205a3d47/meeko/utils/obutils.py

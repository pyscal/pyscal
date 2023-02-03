import numpy as np
import gzip
from ase import Atom, Atoms
import gzip
import io
import os
from ase.io import write, read
import pyscal.formats.ase as ptase
import warnings

def read_snap(infile, compressed = False):
    """
    Function to read a POSCAR format.

    Parameters
    ----------
    infile : string
        name of the input file

    compressed : bool, optional
        force to read a `gz` zipped file. If the filename ends with `.gz`, use of this keyword is not
        necessary, Default False

    Returns
    -------
    atoms : list of `Atom` objects
        list of all atoms as created by user input

    box : list of list of floats
        list of the type `[[xlow, xhigh], [ylow, yhigh], [zlow, zhigh]]` where each of them are the lower
        and upper limits of the simulation box in x, y and z directions respectively.

    Examples
    --------
    >>> atoms, box = read_poscar('POSCAR')
    >>> atoms, box = read_poscar('POSCAR.gz')
    >>> atoms, box = read_poscar('POSCAR.dat', compressed=True)

    """
    aseobj = read(infile, format="vasp")
    atoms, box = ptase.read_snap(aseobj)
    return atoms, box


def write_snap(sys, outfile, comments="pyscal", species=None):
    """
    Function to read a POSCAR format.

    Parameters
    ----------
    outfile : string
        name of the input file


    """
    if species is None:
        warnings.warn("Using legacy poscar writer, to use ASE backend specify species")
        write_poscar(sys, outfile, comments=comments)
    else:
        aseobj = ptase.convert_snap(sys, species=species)
        write(outfile, aseobj, format="vasp")


def split_snaps(**kwargs):
    raise NotImplementedError("split method for mdtraj is not implemented")

def convert_snap(**kwargs):
    raise NotImplementedError("convert method for mdtraj is not implemented")

def write_poscar(sys, outfile, comments="pyscal"):
    """
    Function to read a POSCAR format.
    Parameters
    ----------
    outfile : string
        name of the input file
    """

    fout = open(outfile, 'w')

    fout.write(comments+"\n")
    fout.write("   1.00000000000000\n")

    #write box
    vecs = sys.box
    fout.write("      %1.14f %1.14f %1.14f\n"%(vecs[0][0], vecs[0][1], vecs[0][2]))
    fout.write("      %1.14f %1.14f %1.14f\n"%(vecs[1][0], vecs[1][1], vecs[1][2]))
    fout.write("      %1.14f %1.14f %1.14f\n"%(vecs[2][0], vecs[2][1], vecs[2][2]))

    atypes = sys.atoms.types
    
    tt, cc  = np.unique(atypes, return_counts=True)
    
    atomgroups = [[] for x in range(len(tt))]
    
    for count, t in enumerate(tt):
        for cc, pos in enumerate(sys.atoms.positions):
            if int(atypes[cc]) == t:
                atomgroups[count].append(pos)

    fout.write("  ")
    for c in cc:
        fout.write("%d   "%int(c))
    fout.write("\n")

    fout.write("Cartesian\n")

    for i in range(len(atomgroups)):
        for pos in atomgroups[i]:
            fout.write(" %1.14f %1.14f %1.14f\n"%(pos[0], pos[1], pos[2]))

    fout.close()

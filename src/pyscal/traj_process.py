




"""
pyscal module containing methods for processing of a trajectory. Methods for
reading of input files formats, writing of output files etc are provided in
this module.

"""
import numpy as np
import gzip
import pyscal.catom as pca
from ase import Atom, Atoms
import gzip
import io

import pyscal.trajectory.ase as ptase
import pyscal.trajectory.lammps as ptlammps
import pyscal.trajectory.mdtraj as ptmdtraj
import pyscal.trajectory.vasp as ptvasp


def read_file(filename, format="lammps-dump",
    compressed = False, customkeys=None, is_triclinic = False):
    """
    Read input file

    Parameters
    ----------
    filename : string
        name of the input file.

    format : {'lammps-dump', 'poscar', 'ase', 'mdtraj'}
        format of the input file, in case of `ase` the ASE Atoms object

    compressed : bool, optional
        If True, force to read a `gz` compressed format, default False.

    customkeys : list
        A list containing names of headers of extra data that needs to be read in from the
        input file.

    is_triclinc : bool, optional
        Only used in the case of `format='ase'`. If the read ase object is triclinic, this
        options should be set to True.

    Returns
    -------
    None
    """
    if not os.path.exists(filename):
        raise FileNotFoundError("Filename %s not found"%filename)

    if customkeys == None:
        customkeys = []
    customread = (len(customkeys) > 0)

    if format=='lammps-dump':
        atoms, box, is_triclinic = ptlammps.read_snap(filename, compressed=compressed, check_triclinic=True, customkeys=customkeys)
    elif format == 'ase':
        atoms, box = ptase.read_snap(filename)
    elif format == 'mdtraj':
        atoms, box = ptmdtraj.read_snap(filename)
    elif format == 'poscar':
        atoms, box = ptvasp.read_snap(filename, compressed=compressed)
    else:
        raise TypeError("format recieved an unknown option %s"%format)

    return atoms, box, is_triclinic        


def split_trajectory(infile, format='lammps-dump', compressed=False):
    """
    Read in a trajectory file and convert it to individual time slices.

    Parameters
    ----------

    filename : string
        name of input file

    format : format of the input file
        only `lammps-dump` is supported now.

    compressed : bool, optional
        force to read a `gz` zipped file. If the filename ends with `.gz`, use of this keyword is not
        necessary.

    Returns
    -------
    snaps : list of strings
        a list of filenames which contain individual frames from the main trajectory.

    Notes
    -----
    This is a wrapper function around `split_traj_lammps_dump` function.

    """

    snaps = []

    if format=='lammps-dump':
        snaps = ptlammps.split_snaps(infile, compressed = compressed)
    elif format == 'ase':
        snaps = ptase.split_snaps()
    elif format == 'mdtraj':
        snaps = ptmdtraj.split_snaps()
    elif format == 'poscar':
        snaps = ptvasp.split_snaps()
    else:
        raise TypeError("format recieved an unknown option %s"%format)

    return snaps



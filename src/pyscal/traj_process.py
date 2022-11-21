"""
pyscal module containing methods for processing of a trajectory. Methods for
reading of input files formats, writing of output files etc are provided in
this module.

"""
import numpy as np
import gzip
from ase import Atom, Atoms
import gzip
import io
import os
import warnings
from ase.io import write, read
import pyscal.formats.ase as ptase
import pyscal.formats.lammps as ptlammps
import pyscal.formats.mdtraj as ptmdtraj
import pyscal.formats.vasp as ptvasp
import pyscal.atoms as patoms

def read_file(filename, format="lammps-dump",
    compressed = False, customkeys=None):
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

    Returns
    -------
    None
    """

    if customkeys == None:
        customkeys = []
    customread = (len(customkeys) > 0)

    if format=='lammps-dump':
        atoms, box = ptlammps.read_snap(filename, compressed=compressed, customkeys=customkeys)
    elif format == 'ase':
        atoms, box = ptase.read_snap(filename)
    elif format == 'mdtraj':
        atoms, box = ptmdtraj.read_snap(filename)
    elif format == 'poscar':
        atoms, box = ptvasp.read_snap(filename, compressed=compressed)
    else:
        try:
            #try to use ASE backend
            aseobject = read(filename, format=format)
            atoms, box = ptase.read_snap(aseobject)
            warnings.info("Using ase backend to read file")
        except:
            raise TypeError("format recieved an unknown option %s"%format)

    pyscal_atoms = patoms.Atoms()
    pyscal_atoms.from_dict(atoms)
    return pyscal_atoms, box       

def write_file(sys, outfile, format="lammps-dump", compressed = False, 
    customkeys=None, customvals=None, timestep=0, species=None):
    """
    Write the state of the system to a trajectory file.

    Parameters
    ----------
    sys : `System` object
        the system object to be written out

    outfile : string
        name of the output file

    format : string, optional
        format of the output file

    compressed : bool, default false
        write a `.gz` format

    customkey : string or list of strings, optional
        If specified, it adds this custom column to the dump file. Default None.

    customvals : list or list of lists, optional
        If `customkey` is specified, `customvals` take an array of the same length
        as number of atoms, which contains the values to be written out.
    
    timestep: int, optional
        Specify the timestep value, default 0

    species : None, optional
        species of the atoms. Required if any format other than 'lammps-dump' is used. Required
        for convertion to ase object.


    Returns
    -------
    None

    """
    if format=='lammps-dump':
        ptlammps.write_snap(sys, outfile, compressed = compressed,
        customkeys = customkeys, customvals = customvals, timestep = timestep)
    
    elif format == 'lammps-data':
        aseobject = ptase.convert_snap(sys, species=species)
        write(outfile, aseobject, format='lammps-data')
    
    elif format == 'poscar':
        ptvasp.write_snap(sys, outfile, species=species)
        #aseobject = ptase.convert_snap(sys, species=species)
        #write(outfile, aseobject, format='vasp')

    else:
        #try a write using ase
        try:
            aseobject = ptase.convert_snap(sys, species=species)
            write(outfile, aseobject, format=format)
            warnings.info("Using ase backend to write file")
        except:
            raise TypeError("format recieved an unknown option %s"%format)


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
    if not os.path.exists(infile):
        raise FileNotFoundError("Filename %s not found"%infile)

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



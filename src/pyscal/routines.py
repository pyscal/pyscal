"""
The pyscal routines module. This will contain slightly bigger routines which needs to be used.
Most of them will not be used at the user level. User level longer routines will go to misc
module.

Imports might need specific modules!
"""

import numpy as np
import os


def get_energy_atom(outfile, species=None, pair_style=None, pair_coeff=None, mass=1.0):
    """
    Get energy per atom using a LAMMPS calculator
    """

    #file is written
    #now lammps part
    try:
        from lammps import lammps
    except ImportError:
        raise ModuleNotFoundError("energy method needs lammps compiled as a library. Either install with conda - conda install -c conda-forge lammps,\
                or see here - https://lammps.sandia.gov/doc/Howto_pylammps.html")

    #start routine
    lmp = lammps()
    lmp.command("echo log")
    lmp.command("units metal")
    lmp.command("atom_style atomic")
    lmp.command("boundary p p p")

    if (pair_style is None):
        raise ValueError("pair style has to be provided")

    if (pair_coeff is None):
        raise ValueError("pair coeff has to be provided")

    lmp.command('read_data %s'%outfile)
    lmp.command('pair_style %s'%pair_style)
    lmp.command('pair_coeff %s'%pair_coeff)
    lmp.command("mass * %f"%mass)

    lmp.command("compute 1 all pe/atom")
    lmp.command("run 0")
    eng = lmp.extract_compute("1",1,1)
    ids = lmp.extract_atom("id",0)

    natoms = lmp.get_natoms()
    eng = [eng[x] for x in range(natoms)]
    ids = [ids[x] for x in range(natoms)]

    iddict = dict(zip(np.array(ids).astype(str), eng))
    return iddict



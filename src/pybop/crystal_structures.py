"""
Small module to create structures and so on
"""

import pybop.ccore as pc
import numpy as np

def make_crystal(structure, lattice_constant = 1.00, repetitions = [1, 1, 1], ca_ratio = 1.633):
    """
    Create a basic crystal structure and return it as a list of `Atom` objects
    and box dimensions.

    Parameters
    ----------
    structure : string , bcc, fcc or hcp
        type of the crystal structure
    lattice_constant : float, default 1
        lattice constant of the crystal structure
    repetitions : list of ints of len 3
        of type [nx, ny, nz], repetions of the unit cell in x, y and z directions.
    ca_ratio : float, default 1.633
        ratio of c/a for hcp structures.
    """
    nx = repetitions[0]
    ny = repetitions[1]
    nz = repetitions[2]

    if structure == 'bcc':
        
        coord_no = 2
        natoms = coord_no*nx*ny*nz

        xfact = 1.
        yfact = 1.
        zfact = 1.

        unitcellx = np.zeros(coord_no)
        unitcelly = np.zeros(coord_no)
        unitcellz = np.zeros(coord_no)
        unitcellx[1] = 0.5*lattice_constant
        unitcelly[1] = 0.5*lattice_constant
        unitcellz[1] = 0.5*lattice_constant

    elif structure == 'fcc':

        coord_no = 4
        natoms = coord_no*nx*ny*nz

        xfact = 1.
        yfact = 1.
        zfact = 1.

        unitcellx = np.zeros(coord_no)
        unitcelly = np.zeros(coord_no)
        unitcellz = np.zeros(coord_no)
        unitcellx[1] = 0.5*lattice_constant
        unitcellz[1] = 0.5*lattice_constant
        unitcelly[2] = 0.5*lattice_constant
        unitcellz[2] = 0.5*lattice_constant
        unitcellx[3] = 0.5*lattice_constant
        unitcelly[3] = 0.5*lattice_constant

    elif structure == 'hcp':

        coord_no = 4
        natoms = coord_no*nx*ny*nz

        xfact = 1.
        yfact = np.sqrt(3)
        zfact = ca_ratio

        unitcellx = np.zeros(coord_no)
        unitcelly = np.zeros(coord_no)
        unitcellz = np.zeros(coord_no)
        unitcellx[1] = 0.5*lattice_constant
        unitcelly[1] = 0.5*lattice_constant*yfact
        unitcellx[2] = 0.5*lattice_constant
        unitcelly[2] = lattice_constant*(1.0/6.0)*yfact
        unitcellz[2] = 0.5*lattice_constant*zfact
        unitcelly[3] = 2.0*lattice_constant*(1.0/yfact)
        unitcellz[3] = 0.5*lattice_constant*zfact

    m = 0
    co = 1
    atoms = []
    xh = nx*lattice_constant*xfact
    yh = ny*lattice_constant*yfact
    zh = nz*lattice_constant*zfact
    boxdims = [[0, xh], [0, yh], [0, zh]]

    #create structure
    for i in range(1, nx+1):
        for j in range(1, ny+1):
            for k in range(1, nz+1):
                for l in range(1, coord_no+1):
                    m += 1
                    posx = (unitcellx[l-1]+(lattice_constant*xfact*(float(i)-1)))
                    posy = (unitcelly[l-1]+(lattice_constant*yfact*(float(j)-1)))
                    posz = (unitcellz[l-1]+(lattice_constant*zfact*(float(k)-1)))
                    atom = pc.Atom()
                    atom.set_x([posx, posy, posz])
                    atom.set_id(co)
                    atom.set_type(1)
                    atoms.append(atom)
                    co += 1

    return atoms, boxdims





"""
pyscal module for creating crystal structures.
"""

import pyscal.catom as pc
import numpy as np
import warnings

def make_crystal(structure, lattice_constant = 1.00, repetitions = None, ca_ratio = 1.633, noise = 0):
    """
    Create a basic crystal structure and return it as a list of `Atom` objects
    and box dimensions.

    Parameters
    ----------
    structure : {'bcc', 'fcc', 'hcp', 'diamond' or 'l12'}
        type of the crystal structure

    lattice_constant : float, optional
        lattice constant of the crystal structure, default 1

    repetitions : list of ints of len 3, optional
        of type `[nx, ny, nz]`, repetions of the unit cell in x, y and z directions.
        default `[1, 1, 1]`.

    ca_ratio : float, optional
        ratio of c/a for hcp structures, default 1.633

    noise : float, optional
        If provided add normally distributed noise with standard deviation `noise` to the atomic positions.

    Returns
    -------
    atoms : list of `Atom` objects
        list of all atoms as created by user input

    box : list of list of floats
        list of the type `[[xlow, xhigh], [ylow, yhigh], [zlow, zhigh]]` where each of them are the lower
        and upper limits of the simulation box in x, y and z directions respectively.

    Examples
    --------
    >>> atoms, box = make_crystal('bcc', lattice_constant=3.48, repetitions=[2,2,2])
    >>> sys = System()
    >>> sys.assign_atoms(atoms, box)

    """
    if repetitions == None:
        nx = 1
        ny = 1
        nz = 1
    else:
        nx = repetitions[0]
        ny = repetitions[1]
        nz = repetitions[2]

    #if noise > 0.1:
    #    warnings.warn("Value of noise is rather high. Atom positions might overlap")

    if structure == 'bcc':

        coord_no = 2
        atomtype = [1, 1]

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
        atomtype = [1, 1, 1, 1]

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
        atomtype = [1, 1, 1, 1]

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

    elif structure == 'diamond':

        coord_no = 8
        atomtype = [1, 1, 1, 1, 1, 1, 1, 1]

        natoms = coord_no*nx*ny*nz

        xfact = 1.
        yfact = 1.
        zfact = 1.

        unitcellx = np.zeros(coord_no)
        unitcelly = np.zeros(coord_no)
        unitcellz = np.zeros(coord_no)
        unitcellx[1]=0.25*lattice_constant
        unitcelly[1]=0.25*lattice_constant
        unitcellz[1]=0.25*lattice_constant
        unitcellx[2]=0.50*lattice_constant
        unitcelly[2]=0.50*lattice_constant
        unitcellz[2]=0.00*lattice_constant
        unitcellx[3]=0.75*lattice_constant
        unitcelly[3]=0.75*lattice_constant
        unitcellz[3]=0.25*lattice_constant
        unitcellx[4]=0.50*lattice_constant
        unitcelly[4]=0.00*lattice_constant
        unitcellz[4]=0.50*lattice_constant
        unitcellx[5]=0.00*lattice_constant
        unitcelly[5]=0.50*lattice_constant
        unitcellz[5]=0.50*lattice_constant
        unitcellx[6]=0.75*lattice_constant
        unitcelly[6]=0.25*lattice_constant
        unitcellz[6]=0.75*lattice_constant
        unitcellx[7]=0.25*lattice_constant
        unitcelly[7]=0.75*lattice_constant
        unitcellz[7]=0.75*lattice_constant

    elif structure == 'l12':

        coord_no = 4
        atomtype = [1, 2, 2, 2]

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

    else:
        raise ValueError("Unknown crystal structure")

    m = 0
    co = 1
    atoms = []
    xh = nx*lattice_constant*xfact
    yh = ny*lattice_constant*yfact
    zh = nz*lattice_constant*zfact
    box = [[xh, 0, 0], [0, yh, 0], [0, 0, zh]]

    #create structure
    for i in range(1, nx+1):
        for j in range(1, ny+1):
            for k in range(1, nz+1):
                for l in range(1, coord_no+1):
                    m += 1
                    posx = (unitcellx[l-1]+(lattice_constant*xfact*(float(i)-1)))
                    posy = (unitcelly[l-1]+(lattice_constant*yfact*(float(j)-1)))
                    posz = (unitcellz[l-1]+(lattice_constant*zfact*(float(k)-1)))
                    if noise > 0:
                        posx = np.random.normal(loc=posx, scale=noise)
                        posy = np.random.normal(loc=posy, scale=noise)
                        posz = np.random.normal(loc=posz, scale=noise)
                    atom = pc.Atom()
                    atom.pos = [posx, posy, posz]
                    atom.id = co
                    atom.type = atomtype[l-1]
                    atom.loc = co-1
                    atoms.append(atom)
                    co += 1

    return atoms, box

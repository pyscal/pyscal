import numpy as np
import gzip
import pyscal.catom as pca
from ase import Atom, Atoms
import gzip
import io


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

    raw = infile.split('.')
    if raw[-1] == 'gz' or  compressed:
        f = gzip.open(infile,'rt')
    else:
        f = open(infile,'r')

    data = []
    for line in f:
        data.append(line)

    no_atoms = data[5].split()
    nlev = 5

    try:
        no_atoms = np.array(no_atoms)
        no_atoms = no_atoms.astype(int)
    except ValueError:
        no_atoms = data[6].split()
        nlev = 6
        no_atoms = np.array(no_atoms)
        no_atoms = no_atoms.astype(int)
    except:
        raise ValueError("Unknown no of atoms")

    natoms = np.sum(no_atoms)
    atom_list = no_atoms

    scaling_factor = np.array(data[1].strip()).astype(float)
    xvector = np.array(data[2].strip().split()).astype(float)
    yvector = np.array(data[3].strip().split()).astype(float)
    zvector = np.array(data[4].strip().split()).astype(float)
    boxvecs = [xvector, yvector, zvector]

    xlow = 0
    xhigh = scaling_factor*xvector[0]
    ylow = 0
    yhigh = scaling_factor*yvector[1]
    zlow = 0
    zhigh = scaling_factor*zvector[2]

    if (data[nlev+1].strip().split()[0][0]=='s' or data[nlev+1].strip().split()[0][0]=='S'):
        selective_dynamics=True
        cord_system=data[nlev+2].strip()
        atom_start = nlev+3
    else:
        cord_system=data[nlev+1].strip()
        atom_start = nlev+2

    if cord_system in ['Cartesian', 'cartesian']:
        xscale = 1
        yscale = 1
        zscale = 1
    else:
        xscale = xhigh
        yscale = yhigh
        zscale = zhigh

    species = 1
    count = 0

    cum_list = np.cumsum(atom_list)
    i = atom_start
    atoms = []

    while i in range(atom_start,atom_start+natoms):
        if (count<cum_list[species-1]):
            raw = np.array(data[i].strip().split()[:3]).astype(float)
            typ = species
            x = float(raw[0])*xscale
            y = float(raw[1])*yscale
            z = float(raw[2])*zscale
            #if x,y,z are out of the box, they need to be put in
            if (x < xlow):
                x = x + (xhigh - xlow)
            elif (x > xhigh):
                x = x - (xhigh - xlow)
            if (y < ylow):
                y = y + (yhigh - ylow)
            elif (y > yhigh):
                y = y - (yhigh - ylow)
            if (z < zlow):
                z = z + (zhigh - zlow)
            elif (z > zhigh):
                z = z - (zhigh - zlow)

            count+=1
            idd = count
            atom = pca.Atom()
            atom.pos = [x, y, z]
            atom.id = idd
            atom.type = typ
            atom.loc = i-atom_start
            #atom = pc.Atom(pos=, id=idd, type=typ)
            atoms.append(atom)
            i+=1
        else:
            species+=1

    return atoms, boxvecs

def write_snap(sys, outfile, comments="pyscal"):
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

    atoms = sys.atoms
    atypes = [atom.type for atom in atoms]
    
    tt, cc  = np.unique(atypes, return_counts=True)
    atomgroups = [[] for x in range(len(tt))]
    
    for t in tt:
        for atom in atoms:
            if int(atom.type) == t:
                atomgroups[t-1].append(atom)

    fout.write("  ")
    for c in cc:
        fout.write("%d   "%int(c))
    fout.write("\n")

    fout.write("Cartesian\n")

    for i in range(len(atomgroups)):
        for atom in atomgroups[i]:
            pos = atom.pos
            fout.write(" %1.14f %1.14f %1.14f\n"%(pos[0], pos[1], pos[2]))

    fout.close()

def split_snaps(**kwargs):
	raise NotImplementedError("split method for mdtraj is not implemented")

def convert_snap(**kwargs):
	raise NotImplementedError("convert method for mdtraj is not implemented")

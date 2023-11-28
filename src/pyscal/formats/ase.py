import numpy as np
import gzip
import pyscal.catom as pca
from ase.io import read
from ase import Atom, Atoms
import gzip
import io
import os

#new function to wrap over ase objects
def read_snap(filename, check_triclinic=False):
    """
    Function to read from a ASE atoms objects

    Parameters
    ----------
    filename : str
	path of the ASE file or ASE object

    triclinic : bool, optional
        True if the configuration is triclinic

    """
    #We have to process atoms and atomic objects from ase
    #Known issues lammps -dump modified format
    #first get box
    if os.path.exists(filename):
        aseobject = read(filename)
    else:
        aseobject = filename
        
    a = np.array(aseobject.cell[0])
    b = np.array(aseobject.cell[1])
    c = np.array(aseobject.cell[2])

    box = np.array([a, b, c])

    #box and box dims are set. Now handle atoms
    chems = np.array(aseobject.get_chemical_symbols())
    atomsymbols = np.unique(aseobject.get_chemical_symbols())
    atomtypes = np.array(range(1, len(atomsymbols)+1))
    typedict = dict(zip(atomsymbols, atomtypes))

    #now start parsing atoms
    atoms = []
    positions = aseobject.positions
    for count, position in enumerate(positions):
        atom = pca.Atom()
        atom.pos = list(position)
        atom.id = (count+1)
        atom.type = typedict[chems[count]]
        atom.loc = count

        customdict = {'species': chems[count]}
        atom.custom = customdict
        atoms.append(atom)

    return atoms, box

def write_snap(**kwargs):
    raise NotImplementedError("write method for mdtraj is not implemented")

def split_snaps(**kwargs):
    raise NotImplementedError("split method for mdtraj is not implemented")
    
def convert_snap(sys, species=None):
    """
    Convert a given pyscal structure to ase object

    Parameters
    ----------
    sys : System object
        the system object to be converted

    species : list of str
        a list of species in the system

    Returns
    -------
    aseobject: ASE atoms object

    Notes
    -----
    ASE needs the species of atoms. If a property called `species`
    exist in :attr:`~pyscal.catom.Atom.custom`, this value is used for
    species. However if the value is not present, the keyword `species`
    is required. This should contain a mapping between :attr:`~pyscal.catom.Atom.type`
    and species name. For example, if `species` is `['Au', 'Ge']`, all atoms
    of type 1 are assigned as Au and those of type 2 are assigned as Ge.
    Note that ase is required to run this method.

    """
    #we only do a local import of ASE, this is not super nice
    #we can change this later depending on if ASE is to be treated
    #as a full dependency
    

    atoms = sys.atoms
    #get element strings
    if 'species' not in atoms[0].custom.keys():
        if species is None:
            raise ValueError("Species was not known! To convert to ase, species need to be provided using the species keyword")
        #otherwise we know the species
        types = [atom.type for atom in atoms]
        unique_types = np.unique(types)
        if not (len(unique_types) == len(species)):
            raise ValueError("Length of species and number of types found in system are different. Maybe you specified \"Au\" instead of [\"Au\"]")
        #now assign the species to custom
        
        for atom in atoms:
            custom = atom.custom
            custom['species'] = species[int(atom.type-1)]
        #we should also get the unique species key
        specieskey = "".join(species)
    else:
        #now if species are already there in custom
        #we can safely ignore any input
        types = [atom.type for atom in atoms]
        unique_types = np.unique(types)
        #now we know how many types are there
        species = []
        for ut in unique_types:
            for atom in atoms:
                if ut == atom.type:
                    species.append(atom.custom['species'])
                    break
        specieskey = "".join(species)
      
    cell = sys.box
    pbc = [1, 1, 1]

    #create ASE Atoms and assign everything
    aseobject = Atoms()
    aseobject.cell = cell
    aseobject.pbc = pbc
    
    #thats everything pretty much
    #now create ase Atom
    for atom in atoms:
        aseatom = Atom(atom.custom['species'], atom.pos)
        aseobject.append(aseatom)
    #done
    return aseobject

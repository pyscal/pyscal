import numpy as np
import gzip
from ase import Atom, Atoms
import gzip
import io
import os

#new function to wrap over ase objects
def read_snap(aseobject, check_triclinic=False):
    """
    Function to read from a ASE atoms objects

    Parameters
    ----------
    aseobject : ASE Atoms object
        name of the ASE atoms object

    triclinic : bool, optional
        True if the configuration is triclinic

    """
    #We have to process atoms and atomic objects from ase
    #Known issues lammps -dump modified format
    #first get box
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
    positions = aseobject.positions
    ids = []
    types = []
    species = []

    for count, position in enumerate(positions):
        ids.append(count+1)
        types.append(typedict[chems[count]])
        species.append(chems[count])

    atoms = {}
    atoms['positions'] = positions
    atoms['ids'] = ids
    atoms['types'] = types
    atoms['species'] = species
    atoms['ghost'] = [False for x in range(len(types))]

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
    
    #get element strings
    if 'species' not in sys.atoms.keys():
        sys.atoms["species"] = [None for x in range(sys.atoms.ntotal)]

    if sys.atoms.species[0] is None:
        if species is None:
            raise ValueError("Species was not known! To convert to ase, species need to be provided using the species keyword")
        #otherwise we know the species
        types = sys.atoms.types
        unique_types = np.unique(types)
        if not (len(unique_types) == len(species)):
            raise ValueError("Length of species and number of types found in system are different. Maybe you specified \"Au\" instead of [\"Au\"]")
        #now assign the species to custom
        atomspecies = []        
        for cc, typ in enumerate(types):
            atomspecies.append(species[int(typ-1)])
    else:
        atomspecies = sys.atoms.species
      
    cell = sys.box
    pbc = [1, 1, 1]

    #create ASE Atoms and assign everything
    aseobject = Atoms()
    aseobject.cell = cell
    aseobject.pbc = pbc
    
    #thats everything pretty much
    #now create ase Atom
    for count, pos in enumerate(sys.atoms.positions):
        aseatom = Atom(atomspecies[count], pos)
        aseobject.append(aseatom)
    
    #done
    return aseobject

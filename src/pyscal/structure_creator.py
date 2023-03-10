
import numpy as np
import warnings
import os

from pyscal.atoms import Atoms
from pyscal.attributes import read_yaml
structures = read_yaml(os.path.join(os.path.dirname(__file__), "data/structure_data.yaml"))

def make_crystal(structure, lattice_constant = 1.00, repetitions = None, ca_ratio = 1.633, noise = 0, element=None,
    return_structure_dict=False, structures=structures):
    """
    Create a basic crystal structure and return it as a list of `Atom` objects
    and box dimensions.

    Parameters
    ----------
    structure : {'sc', 'bcc', 'fcc', 'hcp', 'diamond', 'a15' or 'l12'}
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
    elif isinstance(repetitions, int):
        nx = repetitions
        ny = repetitions
        nz = repetitions
    else:
        nx = repetitions[0]
        ny = repetitions[1]
        nz = repetitions[2]

    if structure in structures.keys():
        sdict = structures[structure]
    else:
        raise ValueError("Unknown crystal structure")
    
    unique_types = np.unique(sdict["species"])
    print(element)
    if element is not None:
        if isinstance(element, str):
            element = [element]
        if not (len(element) == len(unique_types)):
            raise ValueError("Elements should equal number of species in the system")
        
        element_dict = dict([x for x in zip(unique_types, element)])
    else:
        element = [None for x in range(len(unique_types))]
    element_dict = dict([x for x in zip(unique_types, element)])
        
    m = 0
    co = 1
    natoms = sdict["natoms"]*nx*ny*nz
    positions = []
    types = []
    ids = []
    species = []

    xh = nx*lattice_constant*sdict["scaling_factors"][0]
    yh = ny*lattice_constant*sdict["scaling_factors"][1]
    zh = nz*lattice_constant*sdict["scaling_factors"][2]
    box = [[xh, 0, 0], [0, yh, 0], [0, 0, zh]]

    #create structure
    for i in range(1, nx+1):
        for j in range(1, ny+1):
            for k in range(1, nz+1):
                for l in range(1, sdict["natoms"]+1):
                    m += 1
                    posx = ((lattice_constant*sdict["positions"][l-1][0])+(lattice_constant*sdict["scaling_factors"][0]*(float(i)-1)))
                    posy = ((lattice_constant*sdict["positions"][l-1][1])+(lattice_constant*sdict["scaling_factors"][1]*(float(j)-1)))
                    posz = ((lattice_constant*sdict["positions"][l-1][2])+(lattice_constant*sdict["scaling_factors"][2]*(float(k)-1)))
                    if noise > 0:
                        posx = np.random.normal(loc=posx, scale=noise)
                        posy = np.random.normal(loc=posy, scale=noise)
                        posz = np.random.normal(loc=posz, scale=noise)
                    
                    positions.append([posx, posy, posz])
                    ids.append(co)
                    types.append(sdict["species"][l-1])
                    species.append(element_dict[sdict["species"][l-1]])
                    co += 1
    atoms = {}
    atoms['positions'] = positions
    atoms['ids'] = ids
    atoms['types'] = types
    atoms['species'] = species
    atoms['ghost'] = [False for x in range(len(types))]

    patoms = Atoms()
    patoms.from_dict(atoms)
    if return_structure_dict:
        return patoms, box, sdict

    return patoms, box

def _update_list_of_elements():
    """
    Only run when needed to update database
    """
    import mendeleev
    from mendeleev import element
    el_list = dir(mendeleev)
    el_dict = {}

    for el in el_list:
        if len(el) == 2:
            if el=="db":
                break
            chem = element(el)
            struct = chem.lattice_structure
            lc = chem.lattice_constant
            if (struct is not None) and (lc is not None):
                if struct.lower() in ['sc', 'bcc', 'fcc', 'hex', 'dia']:
                    el_dict[el] = {}
                    if struct.lower() == "dia":
                        struct = "diamond"
                    elif struct.lower() == "hex":
                        struct = "hcp"
                    el_dict[el]["structure"] = str(struct.lower())
                    el_dict[el]["lattice_constant"] = float(lc)
    return el_dict
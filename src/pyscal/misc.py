import pyscal.core as pc
import numpy as np
import warnings

def compare_atomic_env(infile, atomtype=2, precision=3, format="poscar", print_results=True):
    """
    Compare the atomic environment of given types of atoms
    in the inputfile. The comparison is made in terms of Voronoi
    volume and Voronoi fingerprint.

    Parameters
    ----------
    infile : string
        name of the inputfile

    atomtype: int, optional
        type of the atom
        default 2

    precision: float, optional
        precision for comparing Voronoi volumes
        default 3

    format: string, optional
        format of the input file
        default poscar

    print_results: bool, optional
        if True, print the results. If False, return the data
        instead. default True

    Returns
    -------
    vvx : list of floats
        unique Voronoi volumes. Returned only if print results is False

    vrx : list of strings
        unique Voronoi polyhedra. Returned only if print results is False

    vvc : list of ints
        number of unique quantities specified above. Returned only if print results is False
    """
    sys = pc.System()
    sys.read_inputfile(infile, format=format)
    sys.find_neighbors(method="voronoi")
    sys.calculate_vorovector()
    atoms = sys.atoms
    vols = []
    vors = []
    for atom in atoms:
        if atom.type == atomtype:
            vols.append(atom.volume)
            vors.append(" ".join(np.array(atom.vorovector).astype(str)))
    vols = np.array(vols)
    vols = np.round(vols, decimals=precision)
    vvx, vvc = np.unique(vols, return_counts=True)
    vrx, vrc = np.unique(vors, return_counts=True)
    if (len(vvx) != len(vrx)):
        warnings.warn("Different voronoi polyhedra with same volume!")

    if print_results:
        print("%d clusters found"%len(vvx))
        for i in range(len(vvx)):
            print("voro fingerprint = <%s>, vol = %.3f, counts = %d"%(vrx[i], vvx[i], vvc[i]))
    else:
        return vvx, vrx, vvc

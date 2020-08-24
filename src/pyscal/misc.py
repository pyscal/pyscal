import pyscal.core as pc
import numpy as np
import warnings

def compare_atomic_env(infile, atomtype=2, precision=2, format="poscar", print_results=True, return_system=False):
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

    return_system: bool, optional
        if True, return the system object.
        default False

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

    length_mismatch = False
    if (len(vvx) != len(vrx)):
        warnings.warn("Different voronoi polyhedra with same volume! Fingerprint wont be printed. Maybe change precision?")
        length_mismatch = True

    if precision > 3:
        warnings.warn("More than 3 digits of precision selected!")

    if print_results:
        print("%d clusters found"%len(vvx))
        if not length_mismatch:
            for i in range(len(vvx)):
                print("voro fingerprint = <%s>, vol = %.3f, counts = %d"%(vrx[i], vvx[i], vvc[i]))
        else:
            for i in range(len(vvx)):
                print("voro fingerprint = <x x x x>, vol = %.3f, counts = %d"%(vvx[i], vvc[i]))
        if return_system:
            return sys
    else:
        if return_system:
            return np.array([vvx, vrx, vvc]), sys
        else:
            return np.array([vvx, vrx, vvc])


def find_tetrahedral_voids(infile, format="poscar", print_results=True, return_system=False, 
                                direct_coordinates=True, precision=0.1):
    """
    Check for tetrahedral voids in the system

    Parameters
    ----------
    infile : string
        name of the input file

    format : string
        format of the input file, optional
        default poscar

    print_results: bool, optional
        if True, print the results. If False, return the data
        instead. default True

    return_system: bool, optional
        if True, return the system object.
        default False

    direct_coordinates: bool, optional
        if True, results are provided in direct coordinates
        default False

    precision: int, optional
        the number of digits to check for distances.
        default 1
    
    Returns
    -------
    types : list of atom types
    volumes : list of atom volumes
    pos : list of atom positions
    sys : system object, returns only if return_sys is True

    """

    sys = pc.System()
    sys.read_inputfile(infile, format=format)
    sys.find_neighbors(method="voronoi")
    sys.calculate_vorovector()
    atoms = sys.atoms

    volumes = []
    pos = []
    types = []

    box = sys.box
    boxx = box[0][1] - box[0][0]
    boxy = box[1][1] - box[1][0]
    boxz = box[2][1] - box[2][0]

    if len(atoms) < 21:
        warnings.warn("Very less number of atoms, results maybe wrong. See https://github.com/srmnitc/pyscal/issues/63 ")

    for atom in atoms:
        if atom.vorovector[1] == 4:
            #first to prevent small cells we need to filter
            neighs = np.unique(atom.neighbors)
            dists = np.array([sys.get_distance(atom, atoms[x]) for x in neighs])
            mindist = dists[np.argsort(dists)][0]
            mindistcount = 0
            for i in range(len(dists)):
                if mindist*(1.0-precision) <= dists[i] <= mindist*(1.0+precision):
                    mindistcount += 1

            if mindistcount == 4:
                volumes.append(atom.volume)
                if direct_coordinates:
                    p = atom.pos
                    pos.append([p[0]/boxx, p[1]/boxy, p[2]/boxz])
                else:
                    pos.append(atom.pos)
                types.append(atom.type)

    if print_results:
        if len(volumes) > 0:
            print("%d atoms found in tetrahedral position"%len(volumes))
            for i in range(len(volumes)):
                print("%d  volume %f type %d at %f, %f, %f"%(i+1, volumes[i], types[i], pos[i][0], pos[i][1], pos[i][2]))
        else:
            print("no atoms found in tetrahedral position")

    if return_system:
        return types, volumes, pos, sys
    else:
        return types, volumes, pos
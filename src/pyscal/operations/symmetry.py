import numpy as np
from pyscal.attributes import DocumentedKeywords
import spglib

def get_symmetry(system, angle_tolerance=-1.0,
        symprec=1e-5, hall_number=0):
    """
    Calculate space group number and other symmetry properties

    Parameters
    ----------
    system: pyscal System object
        the input system

    angle_tolerance: float, optional
        An experimental argument that controls angle tolerance between basis vectors. 
        Normally it is not recommended to use this argument. Default -1.0
        https://spglib.github.io/spglib/variable.html#angle-tolerance

    symprec: float, optional
        Distance tolerance in Cartesian coordinates to find crystal symmetry. 
        Default 1e-5
        https://spglib.github.io/spglib/variable.html#symprec

    hall_number: int, optional
        The argument to constrain the space-group-type search only for the Hall symbol corresponding to it. 
        Default 0
        https://spglib.github.io/spglib/dataset.html#hall-number

    Returns
    -------
    results: dict
        dict containing results. The keys are explained below:    
        international_space_group_number: 
            The space group type number defined in International Tables for Crystallography (ITA).
            https://it.iucr.org/Ac/contents/ 
        international_symbol:
            The (full) Hermann–Mauguin notation of space group type
            https://www.chemeurope.com/en/encyclopedia/Hermann-Mauguin_notation.html
        point_group:
            Symbol of the crystallographic point group in the Hermann–Mauguin notation 
            https://www.chemeurope.com/en/encyclopedia/Hermann-Mauguin_notation.html
    
    """
    res = spglib.get_symmetry_dataset((system.box,
        system.direct_coordinates, system.atoms.types))

    results['international_space_group_number'] = results["number"]
    results['international_symbol'] = results["international"]
    results['point_group'] = results["pointgroup"]
    return results

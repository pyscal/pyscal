import numpy as np
import os
import pyscal.csl as csl
import pyscal.crystal_structures as pcs
from pyscal.core import System, Atoms

class GrainBoundary:
    def __init__(self):
        self._axis = None
        self.limit = None
        self._sigma = None
        self._theta = None
        self.gb_index_limit = None
        self.gb_plane = None
        
        
    @property
    def axis(self):
        return self._axis
    
    @axis.setter
    def axis(self, axis):
        self._axis = np.array(axis)
    
    @property
    def sigma(self):
        return self._sigma
    
    @sigma.setter
    def sigma(self, val):
        self._sigma = val
    
    def get_possible_sigma(self):
        if self._axis is not None:
            sigmas, thetas = csl.get_sigma_list(self._axis, self.limit)
            return sigmas
        
    #theta is calculated; not set
    @property
    def theta(self):
        return np.degrees(self._theta)
    
    def calculate_theta(self):
        theta, m, n = csl.get_theta_m_n_list(self._axis, self._sigma)[0]
        R = csl.rot(self._axis, theta)
        M1, M2 = csl.create_minimal_cell_method_1(self._sigma, self._axis, R)
        self._theta = theta
        self._m = m
        self._n = n
        self._M1 = M1
        self._M2 = M2
        self._R = R
        return self._theta
    
    def calculate_possible_gb_planes(self):
        V1, V2, M, Gb = csl.create_possible_gb_plane_list(self._axis, 
                            self._m, self._n, self.gb_index_limit)
        return V1, V2, M, Gb

    def check_if_gb_is_valid(self, gb_vector):
        if self._theta is None:
            self.calculate_theta()
        validity, ortho_1, ortho_2 = csl.find_orthogonal_cell(self._axis, self._sigma, self._theta, 
                                                                   self._R, self._m, self._n, gb_vector,
                                                                   self._M1, self._M2, tol = 0.001)
        self._ortho_1 = ortho_1
        self._ortho_2 = ortho_2
        return validity
    
    def find_gb_character(self):
        if self.gb_plane is not None:
            self.gb_index_limit = max(np.abs(self.gb_plane)+1)
            outs = self.calculate_possible_gb_planes()
            main_str_array = np.array([" ".join(np.array(x).astype(str)) for x in outs[0]])
            search_array = " ".join(np.array(self.gb_plane).astype(str))
            arg = np.argwhere(main_str_array == search_array)
            return outs[3][arg[0][0]]

    def create_grain_boundary(self, axis, 
                              sigma, 
                              gb_plane,
                              limit=3, 
                              gb_index_limit=10,
                              ):
        self.axis = axis
        self.sigma = sigma
        self.limit = limit
        self.gb_plane = gb_plane
        self.gb_index_limit = gb_index_limit
        valid = self.check_if_gb_is_valid(gb_plane)
        if not valid:
            raise TypeError("GB cannot be created with the given input!")
        return valid
    
    def populate_grain_boundary(self, lattice, element=None, repetitions=(1,1,1), lattice_parameter=1, overlap=0.0):
        if lattice in pcs.structures.keys():
            structure = lattice
            element = element
            basis = pcs.structures[lattice]['positions']
            sdict = pcs.structures[lattice]
        elif lattice in pcs.elements.keys():
            structure = pcs.elements[lattice]['structure']
            element = lattice
            lattice_parameter = pcs.elements[lattice]['lattice_constant']
            basis = pcs.structures[structure]['positions']
            sdict = pcs.structures[structure]
        else:
            raise ValueError("Unknown lattice type")
        box, atoms1, atoms2 = csl.populate_gb(self._ortho_1, self._ortho_2, 
                        np.array(basis),
                        lattice_parameter,
                        dim=repetitions, 
                        overlap=overlap)
        atoms = Atoms()
        total_atoms = np.concatenate((atoms1, atoms2))
        if element is not None:
            atoms.from_dict({"positions": total_atoms, "species": [element for x in range(len(total_atoms))]})
        else:
            atoms.from_dict({"positions": total_atoms})
        sys = System()
        sys.box = box
        sys.atoms = atoms
        sys.atoms._lattice = structure
        sys.atoms._lattice_constant = lattice_parameter
        sys._structure_dict = sdict
        sys.remap_atoms_into_box()
        return sys
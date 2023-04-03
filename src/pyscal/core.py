"""


"""


import os
import numpy as np
import warnings
import itertools
from ase.io import write
import uuid
import gzip
import io
from scipy.special import sph_harm
import copy

from pyscal.atoms import Atoms, AttrSetter
import pyscal.csystem as pc
import pyscal.traj_process as ptp
from pyscal.formats.ase import convert_snap
import pyscal.structure_creator as pcs
import pyscal.operations.operations as po
#import pyscal.routines as routines
#import pyscal.visualization as pv


class System:
    """
    Python class for holding the properties of an atomic configuration 
    """    
    def __init__(self, filename=None, format="lammps-dump", 
                                            compressed = False, customkeys=None):
        self.initialized = True
        self.neighbors_found = False
        self.neighbor_method = None
        self.ghosts_created = False
        self.actual_box = None

        #box parameters
        self.rot = [[0,0,0], [0,0,0], [0,0,0]]
        self.rotinv = [[0,0,0], [0,0,0], [0,0,0]]
        self.boxdims = [0,0,0]
        self.triclinic = 0
        self._structure_dict = None
        self._atoms = Atoms()
        
        if filename is not None:
            self.read_inputfile(filename, format=format, 
                                            compressed = compressed, customkeys=customkeys)
    
    @classmethod
    def from_structure(cls, structure, lattice_constant = 1.00, repetitions = None, ca_ratio = 1.633, noise = 0, element=None, chemical_symbol=None):
        atoms, box, sdict = pcs.make_crystal(structure, lattice_constant=lattice_constant,
             repetitions=repetitions, ca_ratio=ca_ratio,
             noise=noise, element=element, return_structure_dict=True)
        obj = cls()
        obj.box = box
        obj.atoms = atoms
        obj.atoms._lattice = structure
        obj.atoms._lattice_constant = lattice_constant
        obj._structure_dict = sdict
        return obj

    def iter_atoms(self):
        return self.atoms.iter_atoms()

    @property
    def natoms(self):
        return self.atoms.natoms

    @property
    def concentration(self):
        return self.atoms.composition

    @property
    def composition(self):
        return self.atoms.composition

    @property
    def box(self):
        """
        Wrap for inbuilt box
        """
        if self.actual_box is not None:
            return self.actual_box
        else:
            return self._box

    @box.setter
    def box(self, userbox):
        """
        Box setter 
        """
        #we should automatically check for triclinic cells here
        summ = 0
        for i in range(3):
            box1 = np.array(userbox[i-1])
            box2 = np.array(userbox[i])
            summ += np.dot(box1, box2)/(np.linalg.norm(box1)*np.linalg.norm(box2))

        #check if the summ is zero
        if (np.abs(summ) > 1E-6):
            #this is a triclinic box
            rot = np.array(userbox).T
            rotinv = np.linalg.inv(rot)
            self.triclinic = 1
            self.rot = rot
            self.rotinv = rotinv

        self.boxdims[0] = np.sum(np.array(userbox[0])**2)**0.5
        self.boxdims[1] = np.sum(np.array(userbox[1])**2)**0.5
        self.boxdims[2] = np.sum(np.array(userbox[2])**2)**0.5
        self._box = userbox

    @property
    def box_dimensions(self):
        return [np.linalg.norm(self.box[x]) for x in range(3)]

    @property
    def direct_coordinates(self):
        dims = self.box_dimensions
        xfrac = np.array(self.atoms.positions)[:,0]/dims[0]
        yfrac = np.array(self.atoms.positions)[:,1]/dims[1]
        zfrac = np.array(self.atoms.positions)[:,2]/dims[2]
        coords = np.column_stack((xfrac, yfrac, zfrac))
        return coords    
    
    @property
    def volume(self):
        """
        Volume of box
        """
        vol = np.dot(np.cross(self._box[0], self._box[1]), self._box[2])
        return vol

    @property
    def atoms(self):
        return self._atoms
    
    @atoms.setter
    def atoms(self, atoms):
        """
        Set atoms
        """
        if(len(atoms['positions']) < 200):
            #we need to estimate a rough idea
            needed_atoms = 200 - len(atoms)
            #get a rough cell
            needed_cells = np.ceil(needed_atoms/len(atoms))
            nx = int(needed_cells**(1/3))
            nx = int(np.ceil(nx/2))

            if np.sum(self.box) == 0:
                raise ValueError("Simulation box should be initialized before atoms")
            atoms, box = self.repeat((nx, nx, nx), atoms=atoms, ghost=True, scale_box=True, assign=False, return_atoms=True)
            self.actual_box = self.box.copy()
            self.box = box

        self._atoms = atoms
        


    def add_atoms(self, atoms):
        """
        Cleanly add a given list of atoms

        Parameters
        ----------
        atoms : dict

        Returns
        -------
        None
        """ 
        ## MOVE TO ATOMS
        self._atoms.add_atoms(atoms)

    def repeat(self, repetitions, atoms=None, ghost=False, scale_box=True, assign=False, return_atoms=False):
        """
        """
        return po.repeat(self, repetitions, 
            atoms=atoms, ghost=ghost, 
            scale_box=scale_box, 
            return_atoms=return_atoms)

    def remap_atoms_into_box(self):
        """
        Go through atoms in the list and remap them into the bix
        """
        rot = np.array(self.box).T
        rotinv = np.linalg.inv(rot)

        for x in range(self.natoms):
            pos = pc.remap_atom_into_box(self.atoms["positions"][x], 
                self.triclinic,
                rot, 
                rotinv, 
                self.box_dimensions)
            self.atoms["positions"][x] = pos        

    def apply_mask(self, mask_type="primary", ids=None, indices=None, condition=None, selection=False):
        """

        Notes
        -----
        Masks can be used to exclude atoms from neighbor calculations. An atom for which
        mask is set to True is excluded from the calculation. There are two types of masks,
        `primary` or `secondary`. For example, neighbors are being calculated for a central
        atom `i`. The neighbor atom is denoted as `j`. If `primary` mask of `i` is True, no neighbor
        calculation is carried out for `i`. If it is False, `i` is considered. Now if `secondary`
        mask of `j` is True, it will not included in the list of neighbors of `i` even if it is within
        the cutoff distance. The `primary` mask of `j` has no effect in this situation.

        An example situation can be to calculate the local concentration around Ni atoms in a NiAl
        structure. In this case, the `primary` mask of all Al atoms can be set to True so that
        only `Ni` atoms are considered. Now, in a second case, the task is to count the number of Al
        atoms around each Ni atom. For this case, the `primary` mask of all Al atoms can be set to True,
        and the `secondary` mask of all Ni atoms can be set to True.

        The masks for ghost atoms are copied from the corresponding mask for real atoms.
        """
        #check if length of mask is equal to length of real atoms
        self.atoms.apply_mask(mask_type=mask_type, ids=ids, 
            indices=indices, condition=condition, selection=selection)

    def remove_mask(self, mask_type="primary", ids=None, indices=None, condition=None, selection=False):
        """
        Remove applied masks

        Parameters
        ----------
        mask_type: string, optional
            type of mask to be applied, either `primary`, `secondary` or `all`

        Returns
        -------
        None
        """
        self._atoms.remove_mask(mask_type=mask_type, ids=ids, 
            indices=indices, condition=condition, selection=selection)

    def apply_selection(self, ids=None, indices=None, condition=None):
        self._atoms.apply_selection(ids=ids, indices=indices, condition=condition)    
    
    def remove_selection(self, ids=None, indices=None, condition=None):
        self._atoms.remove_selection(ids=ids, indices=indices, condition=condition)
    
    def delete(self, ids=None, indices=None, condition=None, selection=False):
        self._atoms.delete(ids=ids, indices=indices, condition=condition, selection=selection)


    def embed_in_cubic_box(self, input_box=None, return_box=False):
        """
        Embedded the triclinic box in a cubic box
        """
        return po.embed_in_cubic_box(self, input_box=input_box,
            return_box=return_box) 

    def get_distance(self, pos1, pos2, vector=False):
        """
        Get the distance between two atoms.

        Parameters
        ----------
        pos1 : list
                first atom position
        pos2 : list
                second atom position
        vector: bool, optional
            If True, return the vector between two atoms

        Returns
        -------
        distance : double
                distance between the first and second atom.

        Notes
        -----
        Periodic boundary conditions are assumed by default.
        """
        
        if vector:
            return dist, diff
        else:
            return dist

    def get_concentration(self):
        """
        Return a dict containing the concentration of the system

        Parameters
        ----------
        None

        Returns
        -------
        condict : dict
            dict of concentration values
        """
        return self.concentration


    def read_inputfile(self, filename, format="lammps-dump", 
                                            compressed = False, customkeys=None):
        """

        Read input file that contains the information of system configuration.

        Parameters
        ----------
        filename : string
            name of the input file.

        format : {'lammps-dump', 'poscar', 'ase', 'mdtraj'}
            format of the input file, in case of `ase` the ASE Atoms object

        compressed : bool, optional
            If True, force to read a `gz` compressed format, default False.

        customkeys : list
            A list containing names of headers of extra data that needs to be read in from the
            input file.

        Returns
        -------
        None

        Notes
        -----
        `format` keyword specifies the format of the input file. Currently only
        a `lammps-dump` and `poscar` files are supported.  Additionaly, the widely
        use Atomic Simulation environment (https://wiki.fysik.dtu.dk/ase/ase/ase.html).
        mdtraj objects (http://mdtraj.org/1.9.3/) are also supported by using the keyword
        `'mdtraj'` for format. Please note that triclinic boxes are not yet supported for
        mdtraj format.
        Atoms object can also be used directly. This function uses the
        :func:`~pyscal.traj_process` module to process a file which is then assigned to system.

        `compressed` keyword is not required if a file ends with `.gz` extension, it is
        automatically treated as a compressed file.

        Triclinic simulation boxes can also be read in.

        If `custom_keys` are provided, this extra information is read in from input files if
        available. This information is can be accessed directly as `self.atoms['customkey']`


        """
        self.neighbors_found = False
        self.neighbor_method = None
        self.ghosts_created = False
        self.actual_box = None

        #box parameters
        self.rot = [[0,0,0], [0,0,0], [0,0,0]]
        self.rotinv = [[0,0,0], [0,0,0], [0,0,0]]
        self.boxdims = [0,0,0]
        self.triclinic = 0

        atoms, box = ptp.read_file(filename, format=format, 
                                    compressed=compressed, customkeys=customkeys,)
        self.box = box
        self.atoms = atoms

    def to_file(self, outfile, format='lammps-dump', customkeys=None, customvals=None,
                compressed=False, timestep=0, species=None):
        """
        Save the system instance to a trajectory file.

        Parameters
        ----------
        outfile : string
            name of the output file

        format : string, {'lammps-dump', 'lammps-data', 'poscar'}
            format of the output file, default `lammps-dump`
            Currently only `lammps-dump` format is supported.

        customkeys : list of strings, optional
            a list of extra atom wise values to be written in the output file.

        customvals : list or list of lists, optional
            If `customkey` is specified, `customvals` take an array of the same length
            as number of atoms, which contains the values to be written out.

        compressed : bool, optional
            If true, the output is written as a compressed file.

        timestep : int, optional
            timestep to be written to file. default 0

        species : None, optional
            species of the atoms. Required if any format other than 'lammps-dump' is used. Required
            for convertion to ase object.

        Returns
        -------
        None

        Notes
        -----
        `to_file` method can handle a number of file formats. The most customizable format is the
        `lammps-dump` which can take a custom options using customkeys and customvals. customkeys
        will be the header written to the dump file. It can be any Atom attribute, any property
        stored in custom variable of the Atom, or calculated q values which can be given by `q4`,
        `aq4` etc. External values can also be provided using `customvals` option. `customvals` array
        should be of the same length as the number of atoms in the system.

        For all other formats, ASE is used to write out the file, and hence the `species` keyword
        needs to be specified. If initially, an ASE object was used to create the System, `species`
        keyword will already be saved, and need not be specified. In other cases, `species` should
        be a list of atomic species in the System. For example `["Cu"]` or `["Cu", "Al"]`, depending
        on the number of species in the System. In the above case, atoms of type 1 will be mapped to
        Cu and of type 2 will be mapped to Al. For a complete list of formats that ASE can handle,
        see `here <https://wiki.fysik.dtu.dk/ase/ase/io/io.html>`_ . 
        """

        ptp.write_file(self, outfile, format = format,
            compressed = compressed, customkeys = customkeys, customvals = customvals,
            timestep = timestep, species = species)

    def to_ase(self, species=None):
        """
        Convert system to an ASE Atoms object

        Parameters
        ----------
        species : list of string
            The chemical species

        Returns
        -------
        None
        """
        return convert_snap(self, species=species)

    def reset_neighbors(self):
        """
        Reset the neighbors of all atoms in the system.

        Parameters
        ----------
        None

        Returns
        -------
        None

        Notes
        -----
        It is used automatically when neighbors are recalculated.

        """
        self.atoms["neighbors"] = []
        self.atoms["neighbordist"] = []
        self.atoms["temp_neighbors"] = []
        self.atoms["temp_neighbordist"] = []
        self.atoms["neighborweight"] = []
        self.atoms["diff"] = []
        self.atoms["r"] = []
        self.atoms["theta"] = []
        self.atoms["phi"] = []
        self.atoms["cutoff"] = []
        self.neighbors_found = False

        
        mapdict = {}
        mapdict["neighbors"] = {}
        mapdict["neighbors"]["index"] = "neighbors"
        mapdict["neighbors"]["distance"] = "neighbordist"
        mapdict["neighbors"]["weight"] = "neighborweight"
        mapdict["neighbors"]["displacement"] = "diff"
        mapdict["neighbors"]["cutoff"] = "cutoff"

        mapdict["neighbors"]["angle"] = {}
        mapdict["neighbors"]["angle"]["polar"] = "theta"
        mapdict["neighbors"]["angle"]["azimuthal"] = "phi"

        mapdict["neighbors"]["temporary"] = {}
        mapdict["neighbors"]["temporary"]["index"] = "temp_neighbors"
        mapdict["neighbors"]["temporary"]["distance"] = "temp_neighbordist"

        self.atoms._add_attribute(mapdict)

    def _check_neighbors(self):
        """
        Check if neighbors are calculated
        """
        if not self.neighbors_found:
            raise ValueError("This calculation needs neighbors to be calculated")

    def average_over_neighbors(self, key, include_self=True):
        """
        Perform a simple average over neighbor atoms

        Parameters
        ----------
        key: string
            atom property

        include_self: bool, optional
            If True, include the host atom in the calculation

        Returns
        -------

        """
        print(key)

        self._check_neighbors()

        if not key in self.atoms.keys():
            raise KeyError("required property not found!")

        test = self.atoms[key][0]

        if isinstance(test, list):
            raise TypeError("Averaging can only be done over 1D quantities")

        avgarr = []
        for i in range(len(self.atoms["positions"])):
            arr = []
            if include_self:
                arr.append(self.atoms[key][i])
            for j in self.atoms["neighbors"][i]:
                arr.append(self.atoms[key][j])
            avgarr.append(np.mean(arr))
        
        return avgarr 

    def find_neighbors(self, method='cutoff', cutoff=0, threshold=2, 
            voroexp=1, padding=1.2, nlimit=6, 
            cells=None, nmax=12, assign_neighbor=True):
        """

        Find neighbors of all atoms in the :class:`~pyscal.core.System`.

        Parameters
        ----------
        method : {'cutoff', 'voronoi', 'number'}
            `cutoff` method finds neighbors of an atom within a specified or adaptive cutoff distance from the atom.
            `voronoi` method finds atoms that share a Voronoi polyhedra face with the atom. Default, `cutoff`
            `number` method finds a specified number of closest neighbors to the given atom. Number only populates
            

        cutoff : { float, 'sann', 'adaptive'}
            the cutoff distance to be used for the `cutoff` based neighbor calculation method described above.
            If the value is specified as 0 or `adaptive`, adaptive method is used.
            If the value is specified as `sann`, sann algorithm is used.

        threshold : float, optional
            only used if ``cutoff=adaptive``. A threshold which is used as safe limit for calculation of cutoff.

        voroexp : int, optional
            only used if ``method=voronoi``. Power of the neighbor weight used to weight the contribution of each atom towards
            Steinhardt parameter values. Default 1.

        padding : double, optional
            only used if ``cutoff=adaptive`` or ``cutoff=number``. A safe padding value used after an adaptive cutoff is found. Default 1.2.

        nlimit : int, optional
            only used if ``cutoff=adaptive``. The number of particles to be considered for the calculation of adaptive cutoff.
            Default 6.
        
        cells : bool, optional
            If True, always use cell lists. Default None.

        nmax : int, optional
            only used if ``cutoff=number``. The number of closest neighbors to be found for each atom. Default 12
        

        Returns
        -------
        None

        Raises
        ------
        RuntimeWarning
            raised when `threshold` value is too low. A low threshold value will lead to 'sann' algorithm not converging
            when finding a neighbor. This function will try to automatically increase `threshold` and check again.

        RuntimeError
            raised when neighbor search was unsuccessful. This is due to a low `threshold` value.

        Notes
        -----
        This function calculates the neighbors of each particle. There are several ways to do this. A complete description of
        the methods can be `found here <https://pyscal.readthedocs.io/en/latest/nearestneighbormethods.html>`_.

        Method cutoff and specifying a cutoff radius uses the traditional approach being the one in which the neighbors of an atom
        are the ones that lie in the cutoff distance around it.

        In order to reduce time during the distance sorting during thefind_neighbors adaptive methods, pyscal sets an initial guess for a cutoff distance.
        This is calculated as,

        .. math:: r_{initial} = threshold * (simulation~box~volume/ number~of~particles)^{(1/3)}

        threshold is a safe multiplier used for the guess value and can be set using the `threshold` keyword.

        In Method cutoff, if ``cutoff='adaptive'``, an adaptive cutoff is found during runtime for each atom [1].
        Setting the cutoff radius to 0 also uses this algorithm. The cutoff for an atom i is found using,

        .. math:: r_c(i) = padding * ((1/nlimit) * \sum_{j=1}^{nlimit}(r_{ij}))

        padding is a safe multiplier to the cutoff distance that can be set through the keyword `padding`. `nlimit` keyword sets the
        limit for the top nlimit atoms to be taken into account to calculate the cutoff radius.

        In Method cutoff, if ``cutoff='sann'``, sann algorithm is used [2]. There are no parameters to tune sann algorithm.

        The second approach is using Voronoi polyhedra which also assigns a weight to each neighbor in the ratio of the face area between the two atoms.
        Higher powers of this weight can also be used [3]. The keyword `voroexp`
        can be used to set this weight.
        
        If method is `number`, instead of using a cutoff value for finding neighbors, a specified number of closest atoms are
        found. This number can be set through the argument `nmax`.

        If `cells` is None, cell lists are used if number of atoms are higher than 2500. If True, cell lists are always used.

        .. warning::

            Adaptive and number cutoff uses a padding over the intial guessed "neighbor distance". By default it is 2. In case
            of a warning that ``threshold`` is inadequate, this parameter should be further increased. High/low value
            of this parameter will correspond to the time taken for finding neighbors.

        References
        ----------
        .. [1] Stukowski, A, Model Simul Mater SC 20, 2012
        .. [2] van Meel, JA, Filion, L, Valeriani, C, Frenkel, D, J Chem Phys 234107, 2012
        .. [3] Haeberle, J, Sperl, M, Born, P, arxiv 2019

        """
        #first reset all neighbors
        self.reset_neighbors()
        self.filter = 0

        if threshold < 1:
            raise ValueError("value of threshold should be at least 1.00")

        if cells is None:
            cells = (self.natoms > 250)

        if method == 'cutoff':            
            if cutoff=='sann':    
                finished = False
                for i in range(1, 10):
                    finished = pc.get_all_neighbors_sann(self.atoms, 0.0, 
                        self.triclinic, self.rot, self.rotinv,
                        self.boxdims, threshold*i, cells)
                    if finished:
                        if i>1:
                            warnings.warn("found neighbors with higher threshold than default/user input")
                        break
                    warnings.warn("Could not find sann cutoff. trying with a higher threshold", RuntimeWarning)
                else:
                    raise RuntimeError("sann cutoff could not be converged. This is most likely, \
                        due to a low threshold value. Try increasing it.")
            
            elif cutoff=='adaptive' or cutoff==0:
                finished = pc.get_all_neighbors_adaptive(self.atoms, 0.0,
                    self.triclinic, self.rot, self.rotinv,
                    self.boxdims, threshold, nlimit, padding, cells)
                if not finished:
                    raise RuntimeError("Could not find adaptive cutoff")
            else:
                if cells:
                    pc.get_all_neighbors_cells(self.atoms, cutoff,
                        self.triclinic, self.rot, self.rotinv, self.boxdims)
                else:
                    pc.get_all_neighbors_normal(self.atoms, cutoff,
                        self.triclinic, self.rot, self.rotinv, self.boxdims)

        elif method == 'number':
            finished = pc.get_all_neighbors_bynumber(self.atoms, 0.0, 
                self.triclinic, self.rot, self.rotinv,
                self.boxdims, threshold, nmax, cells, assign_neighbor)
            if not finished:
                raise RuntimeError("Could not find enough neighbors - try increasing threshold")

        
        elif method == 'voronoi':
            clean_vertices = (cutoff>0)
            #CLEANING is TURNED OFF
            #if not clean_vertices:
                #copy the simulation cell
            backupbox = self._box.copy()
            if self.triclinic:
                if not self.ghosts_created:
                    atoms, box = self.repeat((1, 1, 1), ghost=True, scale_box=True, assign=False, return_atoms=True)
                    self._atoms = atoms
                    self = self.embed_in_cubic_box()
            pc.get_all_neighbors_voronoi(self.atoms, 0.0,
                self.triclinic, self.rot, self.rotinv,
                self.boxdims, voroexp)

            if self.triclinic:
                self._box = backupbox

            #now clean up
            #else:
            #    real_atomdict = {"positions":copy.copy(self.atoms.positions_for_all), 
            #     "ghost":copy.copy(self.atoms.positions_for_ghost)}
                #we need to call the method
                #this means alles good
            #    if self.actual_box is None:
            #        if self.triclinic:
            #            new_box = self.embed_in_cubic_box(inputbox=self._box, return_box=True)
            #            rot = np.array(new_box).T
            #            rotinv = np.linalg.inv(rot)
            #        else:
            #            new_box = self._box
            #            rot = [[0,0,0], [0,0,0], [0,0,0]]
            #            rotinv = [[0,0,0], [0,0,0], [0,0,0]]
            #    #ghosts are present
            #    else:
            #        if self.triclinic:
            #            new_box = self.embed_in_cubic_box(inputbox=self.actual_box, return_box=True)
            #            rot = np.array(new_box).T
            #            rotinv = np.linalg.inv(rot)
            #        else:
            #            new_box = self.actual_box
            #            rot = [[0,0,0], [0,0,0], [0,0,0]]
            #            rotinv = [[0,0,0], [0,0,0], [0,0,0]]

            #    boxdims = [0,0,0]
            #    boxdims[0] = np.sum(np.array(new_box[0])**2)**0.5
            #    boxdims[1] = np.sum(np.array(new_box[1])**2)**0.5
            #    boxdims[2] = np.sum(np.array(new_box[2])**2)**0.5
                
            #    pc.get_all_neighbors_voronoi(real_atomdict, 0.0,
            #        self.triclinic, rot, rotinv,
            #        boxdims, 1)                
                
            #    pc.clean_voronoi_vertices(real_atomdict, 
            #        self.atoms, 0.0,
            #        self.triclinic, rot, rotinv,
            #        boxdims, cutoff)

                #unique_vertices = []
                #for i in range(len(self.vertex_is_unique)):
                #    for j in range(len(self.vertex_is_unique[i])):
                #        if self.vertex_is_unique[i][j]:
                #            unique_vertices.append(self.vertex_positions[i][j])

                #self.atoms["vertex_positions_unique_skipcheck"] = unique_vertices

            #assign extra options
            mapdict = {}
            mapdict["voronoi"] = {}
            mapdict["voronoi"]["volume"] = "voronoi_volume"
            mapdict["voronoi"]["face"] = {}
            mapdict["voronoi"]["face"]["vertices"] = "face_vertices"
            mapdict["voronoi"]["face"]["perimeters"] = "face_perimeters"
            mapdict["voronoi"]["vertex"] = {}
            mapdict["voronoi"]["vertex"]["vectors"] = "vertex_vectors"
            mapdict["voronoi"]["vertex"]["numbers"] = "vertex_numbers"
            mapdict["voronoi"]["vertex"]["positions"] = "vertex_positions"
            #mapdict["voronoi"]["vertex"]["unique_positions"] = "vertex_positions_unique_nofilter"
            self.atoms._add_attribute(mapdict)

        self.neighbors_found = True

    def calculate_q(self, q, averaged=False, continuous_algorithm=False):
        """
        Find the Steinhardt parameter q_l for all atoms.

        Parameters
        ----------
        q : int or list of ints
            A list of all Steinhardt parameters to be found.

        averaged : bool, optional
            If True, return the averaged q values, default False
        
        continuous_algorithm: bool, optional
            See Notes for description.

        Returns
        -------
        q : list of floats
            calculated q values

        Notes
        -----
        Enables calculation of the Steinhardt parameters [1] q. The type of
        q values depend on the method used to calculate neighbors. See the description
        :func:`~pyscal.core.System.find_neighbors` for more details. 

        The option `continuous_algorithm` specifies which algorithm to use for calculations. If False, 
        an algorithm [3] is used. The C++ algorithm is faster is a large, consecutive number of q values (> 200)
        are to be calculated.

        This function creates three new attributes for this class: `qx`, `qx_real` and `qx_imag`,
        where `stands` for the q number.   

        References
        ----------
        .. [1] Steinhardt, PJ, Nelson, DR, Ronchetti, M. Phys Rev B 28, 1983
        .. [2] Lechner, W, Dellago, C, J Chem Phys, 2013
        """
        if isinstance(q, int):
            qq = [q]
        else:
            qq = q

        self._check_neighbors()

        if averaged:
            self._calculate_aq(qq)
            qvals = [self.atoms["avg_q%d"%x] for x in qq]
        else:    
            if continuous_algorithm:
                lm = max(qq)
                pc.calculate_q(self.atoms, lm)
            else:
                self._calculate_q(qq)
            qvals = [self.atoms["q%d"%x] for x in qq]
        return qvals

    def _calculate_q(self, qq):
        """
        Private method for calculation of qvals
        """
        for val in qq:
            pc.calculate_q_single(self.atoms, val)
  
        mapdict = {}
        mapdict["steinhardt"] = {}
        mapdict["steinhardt"]["generic"] = {}
        for val in qq:
            key1a = "q%d_norm"%val
            key1b = "q%d"%val
            key2 = "q%d_real"%val
            key3 = "q%d_imag"%val
            mapdict["steinhardt"]["generic"][key1a] = key1b
            mapdict["steinhardt"]["generic"][key2] = key2
            mapdict["steinhardt"]["generic"][key3] = key3
        self.atoms._add_attribute(mapdict)


    def _calculate_aq(self, qq):
        """
        Private method for calculation of avged qvals
        """

        todo_q = []
        for q in qq:
            keys = ["q%d"%q, "q%d_real"%q, "q%d_imag"%q]
            prod = []
            for key in keys:
                if key in self.atoms.keys():
                    prod.append(True)
                else:
                    prod.append(False)
            prod = np.prod(prod)
            if not prod:
                todo_q.append(q)

        _ = self._calculate_q(todo_q)

        #loop over atoms
        for val in qq:
            pc.calculate_aq_single(self.atoms, val)

        mapdict = {}
        mapdict["steinhardt"] = {}
        mapdict["steinhardt"]["average"] = {}
        for val in qq:
            key1a = "q%d_norm"%val
            key1b = "q%d"%val
            key2 = "q%d_real"%val
            key3 = "q%d_imag"%val
            mapdict["steinhardt"]["average"][key1a] = key1b
            mapdict["steinhardt"]["average"][key2] = key2
            mapdict["steinhardt"]["average"][key3] = key3
        self.atoms._add_attribute(mapdict)

    def calculate_disorder(self, averaged=False, q=6):
        """
        Calculate the disorder criteria for each atom
        
        Parameters
        ----------
        averaged : bool, optional
            If True, calculate the averaged disorder. Default False.
        q : int, optional
            The Steinhardt parameter value over which the bonds have to be calculated.
            Default 6.
        
        Returns
        -------
        None
        
        Notes
        -----
        Calculate the disorder criteria as introduced in [1]. The disorder criteria value for each atom is defined by,
        .. math::
            D_j = \\frac{1}{N_b^j} \sum_{i=1}^{N_b} [ S_{jj} + S_{kk} -2S_{jk}]
        where .. math:: S_{ij} = \sum_{m=-6}^6 q_{6m}(i) q_{6m}^*(i)
        
        Any q value other than six can also be used. This can be specified using the `q` argument.

        The keyword `averaged` is True, the disorder value is averaged over the atom and its neighbors. 
        For ordered systems, the value of disorder would be zero which would increase
        and reach one for disordered systems.

        This function creates two new attributes for this class: `disorder` and `avg_disorder`.
        
        References
        ----------
        .. [1] Kawasaki, T, Onuki, A, J. Chem. Phys. 135, 2011
        """
        #now routine for calculation of disorder

        keys = ["q%d_real"%q, "q%d_imag"%q]
        prod = []
        for key in keys:
            if key in self.atoms.keys():
                prod.append(True)
            else:
                prod.append(False)
        prod = np.prod(prod)
        if not prod:
            self.calculate_q(q)

        pc.calculate_disorder(self.atoms, q)

        mapdict = {}
        mapdict["steinhardt"] = {}
        mapdict["steinhardt"]["disorder"] = {}
        mapdict["steinhardt"]["disorder"]["norm"] = "disorder"

        if averaged:
            #average the disorder
            avg_arr = self.average_over_neighbors("disorder")
            self.atoms["avg_disorder"] = avg_arr
            mapdict["steinhardt"]["disorder"]["average"] = "avg_disorder"
        self.atoms._add_attribute(mapdict)


    def find_solids(self, bonds=0.5, threshold=0.5, avgthreshold=0.6, 
                          cluster=True, q=6, cutoff=0, right=True):
        """
        Distinguish solid and liquid atoms in the system.
        Parameters
        ----------
        bonds : int or float, optional
            Minimum number of solid bonds for an atom to be identified as
            a solid if the value is an integer. Minimum fraction of neighbors
            of an atom that should be solid for an atom to be solid if the
            value is float between 0-1. Default 0.5.
        
        threshold : double, optional
            Solid bond cutoff value. Default 0.5.
        
        avgthreshold : double, optional
            Value required for Averaged solid bond cutoff for an atom to be identified
            as solid. Default 0.6.
        
        cluster : bool, optional
            If True, cluster the solid atoms and return the number of atoms in the largest
            cluster.
        
        q : int, optional
            The Steinhardt parameter value over which the bonds have to be calculated.
            Default 6.
        
        cutoff : double, optional
            Separate value used for cluster classification. If not specified, cutoff used
            for finding neighbors is used.
        
        right: bool, optional
            If true, greater than comparison is to be used for finding solid particles. 
            default True.
        
        Returns
        -------
        solid : int
            Size of the largest solid cluster. Returned only if `cluster=True`.
        
        Notes
        -----
        The neighbors should be calculated before running this function.
        Check :func:`~pyscal.core.System.find_neighbors` method.
        
        `bonds` define the number of solid bonds of an atom to be identified as solid.
        Two particles are said to be 'bonded' if [1],
        .. math:: s_{ij} = \sum_{m=-6}^6 q_{6m}(i) q_{6m}^*(i) \geq threshold
        where `threshold` values is also an optional parameter.
        
        If the value of `bonds` is a fraction between 0 and 1, at least that much of an atom's neighbors
        should be solid for the atom to be solid.
        
        An additional parameter `avgthreshold` is an additional parameter to improve solid-liquid distinction.
        
        In addition to having a the specified number of `bonds`,
        
        .. math::  \langle s_{ij} \\rangle > avgthreshold
        
        also needs to be satisfied. In case another q value has to be used for calculation of S_ij, it can be
        set used the `q` attribute. In the above formulations, `>` comparison for `threshold` and `avgthreshold`
        can be changed to `<` by setting the keyword `right` to False.
        
        If `cluster` is True, a clustering is done for all solid particles. See :func:`~pyscal.csystem.find_clusters`
        for more details. 
        
        References
        ----------
        .. [1] Auer, S, Frenkel, D. Adv Polym Sci 173, 2005
        """
        #check if neighbors are found
        self._check_neighbors()

        if not isinstance(q, int):
            raise TypeError("q should be interger value")

        if not isinstance(threshold, (int, float)):
            raise TypeError("threshold should be a float value")
        else:
            if not ((threshold >= 0 ) and (threshold <= 1 )):
                raise ValueError("Value of threshold should be between 0 and 1")

        if not isinstance(avgthreshold, (int, float)):
            raise TypeError("avgthreshold should be a float value")
        else:
            if not ((avgthreshold >= 0 ) and (avgthreshold <= 1 )):
                raise ValueError("Value of avgthreshold should be between 0 and 1")

        #start identification routine
        #check the value of bonds and set criteria depending on that
        if isinstance(bonds, int):
            criteria = 0
        elif isinstance(bonds, float):
            if ((bonds>=0) and (bonds<=1.0)):
                criteria = 1
            else:
                raise TypeError("bonds if float should have value between 0-1")
        else:
             raise TypeError("bonds should be interger/float value")

        if right:
            compare_criteria = 0
        else:
            compare_criteria = 1

        self.calculate_q(q)

        #calculate solid neighs
        pc.calculate_bonds(self.atoms, q, 
            threshold, avgthreshold, bonds, 
            compare_criteria, criteria)

        mapdict = {}
        mapdict["steinhardt"] = {}
        mapdict["steinhardt"]["order"] = {}
        mapdict["steinhardt"]["order"]["bonds"] = "bonds"
        mapdict["steinhardt"]["order"]["sij"] = {}
        mapdict["steinhardt"]["order"]["sij"]["norm"] = "sij"
        mapdict["steinhardt"]["order"]["sij"]["average"] = "avg_sij"
        mapdict["steinhardt"]["order"]["sij"]["solid"] = "solid"
        self.atoms._add_attribute(mapdict)
        
        if cluster:
            lc = self.cluster_atoms(self.atoms.steinhardt.order.sij.solid, largest=True)
            return lc

    def find_largest_cluster(self):
        """
        Find largest cluster among given clusters
        
        Parameters
        ----------
        None

        Returns
        -------
        lc : int
            Size of the largest cluster.
        """
        if not "cluster" in self.atoms.keys():
            raise RuntimeError("cluster_atoms needs to be called first")

        clusterlist = [x for x in self.atoms["cluster"] if x != -1]
        xx, xxcounts = np.unique(clusterlist, return_counts=True)
        arg = np.argsort(xxcounts)[-1]
        largest_cluster_size = xxcounts[arg]
        largest_cluster_id = xx[arg]


        self.atoms["largest_cluster"] = [True if self.atoms["cluster"][x]==largest_cluster_id else False for x in range(len(self.atoms["cluster"]))]
        
        mapdict = {}
        mapdict["cluster"] = {}
        mapdict["cluster"]["largest"] = "largest_cluster"
        self.atoms._add_attribute(mapdict)

        return largest_cluster_size


    def cluster_atoms(self, condition, largest = True, cutoff=0):
        """
        Cluster atoms based on a property
        
        Parameters
        ----------
        condition : callable or atom property
            Either function which should take an :class:`~Atom` object, and give a True/False output
            or an attribute of atom class which has value or 1 or 0.
        
        largest : bool, optional
            If True returns the size of the largest cluster. Default False.
        
        cutoff : float, optional
            If specified, use this cutoff for calculation of clusters. By default uses the cutoff
            used for neighbor calculation.
        
        Returns
        -------
        lc : int
            Size of the largest cluster. Returned only if `largest` is True.
        
        Notes
        -----
        This function helps to cluster atoms based on a defined property. This property
        is defined by the user through the argument `condition` which is passed as a parameter.
        `condition` should be a boolean array the same length as number of atoms in the system.
        """
        
        self.apply_selection(condition=condition)
        pc.find_clusters(self.atoms, cutoff)

        mapdict = {}
        mapdict["cluster"] = {}
        mapdict["cluster"]["id"] = "cluster"
        self.atoms._add_attribute(mapdict)

        #done!
        lc = self.find_largest_cluster()
        #pcs.System.get_largest_cluster_atoms(self)
        self.remove_selection()
        if largest:
            return lc


    def calculate_rdf(self, rmin=0, rmax=5.00, bins=100):
        """
        Calculate the radial distribution function.
        
        Parameters
        ----------
        rmin : float, optional
            minimum value of the distance histogram. Default 0.0.
        
        rmax : float, optional
            maximum value of the distance histogram. Default 5.0

        bins : int
            number of bins in the histogram
                
        Returns
        -------
        rdf : array of ints
            Radial distribution function
        
        r : array of floats
            radius in distance units
        """
        self.find_neighbors(method="cutoff", cutoff=rmax)
        distances = list(itertools.chain(*self.atoms["neighbordist"]))

        hist, bin_edges = np.histogram(distances, bins=bins, 
            range=(rmin, rmax), density=True)
        edgewidth = np.abs(bin_edges[1]-bin_edges[0])
        hist = hist.astype(float)
        r = bin_edges[:-1]

        #get box density
        rho = self.natoms/self.volume
        shell_vols = (4./3.)*np.pi*((r+edgewidth)**3 - r**3)
        shell_rho = hist/shell_vols
        #now divide to get final value
        rdf = shell_rho/rho
        return rdf, r

    def calculate_angularcriteria(self):
        """
        Calculate the angular criteria for each atom
        
        Parameters
        ----------
        None
        
        Returns
        -------
        None
        
        Notes
        -----
        Calculates the angular criteria for each atom as defined in [1]_. Angular criteria is
        useful for identification of diamond cubic structures. Angular criteria is defined by,
        
        .. math::
            A = \sum_{i=1}^6 (\cos(\\theta_i) + \\frac{1}{3})^2
        
        where cos(theta) is the angle size suspended by each pair of neighbors of the central
        atom. A will have a value close to 0 for structures if the angles are close to 109 degrees.
        The calculated A parameter for each atom can be accessed by system.angular
        
        References
        ----------
        .. [1] Uttormark, MJ, Thompson, MO, Clancy, P, Phys. Rev. B 47, 1993
        """
        self._check_neighbors()
        angulars = []

        for count, pos1 in enumerate(self.atoms["positions"]):
            
            dists = []
            distneighs = []
            distvectors = []

            for count2, neigh in enumerate(self.atoms["neighbors"][count]):
                pos2 = self.atoms["positions"][neigh]
                dist = self.atoms["neighbordist"][count][count2]
                vectors = self.atoms["diff"][count][count2]
                dists.append(dist)
                distneighs.append(neigh)
                distvectors.append(vectors)

            args = np.argsort(dists)
            #find top four
            topfourargs = np.array(args)[:4]

            combos = list(itertools.combinations(topfourargs, 2))
            costhetasum = 0

            for combo in combos:
                vec1 = distvectors[combo[0]]
                vec2 = distvectors[combo[1]]
                modvec1 = np.sqrt(np.sum([x**2 for x in vec1]))
                modvec2 = np.sqrt(np.sum([x**2 for x in vec2]))
                costheta = np.dot(vec1, vec2)/(modvec1*modvec2)
                costhetasum += (costheta +(1./3.))**2
            angulars.append(costhetasum)

        self.atoms["angular"] = angulars
        
        mapdict = {}
        mapdict["angular_parameters"] = {}
        mapdict["angular_parameters"]["diamond_angle"] = "angular"
        self.atoms._add_attribute(mapdict)

    def calculate_chiparams(self, angles=False):
        """
        Calculate the chi param vector for each atom
        
        Parameters
        ----------
        angles : bool, optional
            If True, return the list of cosines of all neighbor pairs
        
        Returns
        -------
        angles : array of floats
            list of all cosine values, returned only if `angles` is True.
        
        Notes
        -----
        This method tries to distinguish between crystal structures by finding the cosines of angles
        formed by an atom with its neighbors. These cosines are then historgrammed with bins
        `[-1.0, -0.945, -0.915, -0.755, -0.705, -0.195, 0.195, 0.245, 0.795, 1.0]` to find a vector for
        each atom that is indicative of its local coordination. Compared to chi parameters from chi_0 to
        chi_7 in the associated publication, the vector here is from chi_0 to chi_8. This is due to an additional
        chi parameter which measures the number of neighbors between cosines -0.705 to -0.195.
        Parameter `nlimit` specifies the number of nearest neighbors to be included in the analysis to find the cutoff.
        If parameter `angles` is true, an array of all cosine values is returned. The publication further provides
        combinations of chi parameters for structural identification which is not implemented here. The calculated
        chi params can be accessed using :attr:`~pyscal.catom.chiparams`.
        
        References
        ----------
        .. [1] Ackland, Jones, Phys. Rev. B 73, 2006
        """

        self._check_neighbors()

        bins = [-1.0, -0.945, -0.915, -0.755, -0.705, -0.195, 0.195, 0.245, 0.795, 1.0]
        chiparams = []
        cosines = []

        for count, pos in enumerate(self.atoms["positions"]):

            dists = self.atoms["neighbordist"][count]
            neighs = self.atoms["neighbors"][count]

            args = range(len(dists))
            combos = list(itertools.combinations(args, 2))
            costhetas = []
            
            for combo in combos:
                vec1 = self.atoms["diff"][count][combo[0]]
                vec2 = self.atoms["diff"][count][combo[1]]
                modvec1 = np.linalg.norm(vec1)
                modvec2 = np.linalg.norm(vec2)
                costheta = np.dot(vec1, vec2)/(modvec1*modvec2)
                #found costheta
                costhetas.append(costheta)


            #now add according to classification in paper
            chivector = np.histogram(costhetas, bins=bins)
            chiparams.append(chivector[0])
            if angles:
                cosines.append(costhetas)
        
        self.atoms["chiparams"] = chiparams
        
        mapdict = {}
        mapdict["angular_parameters"] = {}
        mapdict["angular_parameters"]["chi_params"] = "chiparams"
        
        if angles:
            self.atoms["cosines"] = cosines
            mapdict["angular_parameters"]["cosines"] = "cosines"

        self.atoms._add_attribute(mapdict)
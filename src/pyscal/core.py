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

import pyscal.csystem as pc
import pyscal.traj_process as ptp
from pyscal.formats.ase import convert_snap

#import pyscal.routines as routines
#import pyscal.visualization as pv


class System:
    """
    Python class for holding the properties of an atomic configuration 
    """
    def __init__(self):
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
        self._atoms = {}

    #overload get methods
    def __getattr__(self, name):
        if name in self.atoms.keys():
            res = [self.atoms[name][x] for x in range(len(self.atoms[name])) if self.atoms["ghost"][x]==False]
            return res            
        else:
            raise AttributeError("Attribute %s not found"%name)

    @property
    def natoms(self):
        if 'positions' in self.atoms.keys():
            nop = np.sum([1 for x in range(len(self.atoms["positions"])) if self.atoms["ghost"][x]==False])
            return nop
        else:
            return 0

    @property
    def concentration(self):
        return self.get_concentration()

    @property
    def composition(self):
        return self.get_concentration()

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
    def atoms(self):
        """
        Atom access
        """
        return self._atoms

    @atoms.setter
    def atoms(self, atoms):
        """
        Set atoms
        """
        #we need to check atoms and add necessary keys
        if not 'positions' in atoms.keys():
            raise ValueError('positions is a necessary key in atoms')
        nop = len(atoms["positions"])
        for key, val in atoms.items():
            if not (len(val)==nop):
                raise ValueError("All times in the atoms dict should have same length as positions")
        #now add necessary keys-ids, types, ghost
        if not 'ids' in atoms.keys():
            atoms['ids'] = [x+1 for x in range(nop)]
        if not 'types' in atoms.keys():
            atoms['types'] = [1 for x in range(nop)]
        if not 'ghost' in atoms.keys():
            atoms['ghost'] = [False for x in range(nop)]

        if(len(atoms['positions']) < 200):
            #we need to estimate a rough idea
            needed_atoms = 200 - len(atoms)
            #get a rough cell
            needed_cells = np.ceil(needed_atoms/len(atoms))
            nx = int(needed_cells**(1/3))
            nx = int(np.ceil(nx/2))

            if np.sum(self.box) == 0:
                raise ValueError("Simulation box should be initialized before atoms")
            atoms = self.repeat((nx, nx, nx), atoms=atoms, ghost=True, scale_box=True)

        self._atoms = atoms

    #iterator for atoms
    def iter_atoms(self):
        """
        Iter over atoms
        """
        for i in range(len(self.atoms["positions"])):
            if not self.atoms["ghost"][i]:
                rdict = {}
                for key in self.atoms.keys():
                    rdict[key] = self.atoms[key][i]
            yield rdict 

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
        if not 'positions' in atoms.keys():
            raise ValueError('positions is a necessary key in atoms')
        nop = len(atoms["positions"])
        for key, val in atoms.items():
            if not (len(val)==nop):
                raise ValueError("All times in the atoms dict should have same length as positions")
        
        #now add necessary keys-ids, types, ghost
        maxid = max(self.atoms["ids"])
        if not 'ids' in atoms.keys():
            atoms['ids'] = [maxid+x+1 for x in range(nop)]
        else:
            for i in atoms['ids']:
                if i in self.atoms['ids']:
                    raise ValueError("Atom id already exists, unique ID is required")

        if not 'types' in atoms.keys():
            atoms['types'] = [1 for x in range(nop)]
        if not 'ghost' in atoms.keys():
            atoms['ghost'] = [False for x in range(nop)]

        for key in self.atoms.keys():
            self.atoms[key] = [*self.atoms[key], *atoms[key]]

    def repeat(self, reps, atoms=None, ghost=False, scale_box=True):
        """
        Replicate simulation cell
        
        Parameters
        ----------
        reps : list of ints of size 3
            repetitions in each direction
        
        atoms : list of atoms, optional
            if not provided, use atoms that are assigned

        ghost : bool, optional
            If True, assign the new atoms as ghost instead of actual atoms
        """
        
        box = self.box        
        self.actual_box = box.copy()

        if atoms is None:
            atoms = self.atoms

        newatoms = []
        idstart = len(atoms) + 1

        x1 = -reps[0]
        x2 = reps[0]+1
        y1 = -reps[1]
        y2 = reps[1]+1
        z1 = -reps[2]
        z2 = reps[2]+1
        xs = 2*reps[0] + 1
        ys = 2*reps[1] + 1
        zs = 2*reps[2] + 1
        tp = atoms['types'][0]

        positions = []
        ids = []
        types = []
        ghosts = []

        for i in range(x1, x2):
            for j in range(y1, y2):
                for k in range(z1, z2):
                    if (i==j==k==0):
                        continue
                    for pos in atoms['positions']:
                        pos = (pos + i*np.array(box[0]) + j*np.array(box[1]) + k*np.array(box[2]))
                        positions.append(pos)
                        ids.append(idstart)
                        idstart += 1
                        types.append(tp)
                        ghosts.append(ghost)

        if scale_box:
            box[0] = xs*np.array(box[0])
            box[1] = ys*np.array(box[1])
            box[2] = zs*np.array(box[2])
            self.box = box
        if ghost:
            self.ghosts_created = True

        atoms['positions'] = [*atoms['positions'], *positions]
        atoms['ids'] = [*atoms['ids'], *ids]
        atoms['types'] = [*atoms['types'], *types]
        atoms['ghost'] = [*atoms['ghost'], *ghosts]

        return atoms


    def embed_in_cubic_box(self,):
        """
        Embedded the triclinic box in a cubic box
        """
        #first task is to create a box representation
        
        box = self._box
        backupbox = box.copy()
        
        a = np.array(box[0])
        b = np.array(box[1])
        c = np.array(box[2])

        cosa = np.dot(b, c)/(np.linalg.norm(b)*np.linalg.norm(c))
        cosb = np.dot(c, a)/(np.linalg.norm(c)*np.linalg.norm(a))
        cosc = np.dot(a, b)/(np.linalg.norm(a)*np.linalg.norm(b))

        lx = np.linalg.norm(a)
        xy = np.linalg.norm(b)*cosc
        xz = np.linalg.norm(c)*cosb
        ly = np.sqrt(np.linalg.norm(b)*np.linalg.norm(b) - xy*xy)
        yz = (np.linalg.norm(b)*np.linalg.norm(c)*cosa - xy*xz)/ly
        lz = np.sqrt(np.linalg.norm(c)*np.linalg.norm(c) - xz*xz - yz*yz)

        xlo = ylo = zlo = 0
        xhi = lx
        yhi = ly
        zhi = lz

        xlo_bound = xlo + min(0.0,xy,xz,xy+xz)
        xhi_bound = xhi + max(0.0,xy,xz,xy+xz)
        ylo_bound = ylo + min(0.0,yz)
        yhi_bound = yhi + max(0.0,yz)
        zlo_bound = zlo
        zhi_bound = zhi

        newbox = np.array([[xhi_bound-xlo_bound, 0, 0], [0, yhi_bound-ylo_bound, 0], [0, 0, zhi_bound-zlo_bound]])
        
        self.newbox = newbox
        self.box = newbox
        self.actual_box = None

    def get_distance(self, pos1, pos2):
        """
        Get the distance between two atoms.

        Parameters
        ----------
        pos1 : list
                first atom position
        pos2 : list
                second atom position

        Returns
        -------
        distance : double
                distance between the first and second atom.

        Notes
        -----
        Periodic boundary conditions are assumed by default.
        """
        return pc.get_abs_distance(pos1, pos2, self.triclinic, 
            self.rot, self.rotinv, self.boxdims, 0.0, 0.0, 0.0)

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
        typelist = [self.atoms["types"][x] for x in range(len(self.atoms["positions"])) if self.atoms["ghost"][x]==False]
        types, typecounts = np.unique(typelist, return_counts=True)
        concdict = {}
        for c, t in enumerate(types):
            concdict[str(t)] = typecounts[c]
        return concdict


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
        self._atoms = {}

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

    def to_ase(self, species):
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
        self.atoms["neighbor_dist"] = []
        self.atoms["temp_neighbors"] = []
        self.atoms["temp_neighbordist"] = []
        self.atoms["neighborweight"] = []
        self.atoms["diff"] = []
        self.atoms["r"] = []
        self.atoms["theta"] = []
        self.atoms["phi"] = []
        self.atoms["cutoff"] = []

    def find_neighbors(self, method='cutoff', cutoff=None, threshold=2, 
            filter=None, voroexp=1, padding=1.2, nlimit=6, 
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

        filter : {'None', 'type', 'type_r'}, optional
            apply a filter to nearest neighbor calculation. If the `filter` keyword is set to
            `type`, only atoms of the same type would be included in the neighbor calculations. 
            If `type_r`, only atoms of a different type will be included in the calculation. Default None.

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

        In order to reduce time during the distance sorting during the adaptive methods, pyscal sets an initial guess for a cutoff distance.
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
            cells = (len(self.positions) > 250)

        if filter == 'type':
            self.filter = 1
        elif filter == 'type_r':
            self.filter = 2

        if method == 'cutoff':            
            if cutoff=='sann':    
                finished = False
                for i in range(1, 10):
                    finished = pc.get_all_neighbors_sann(self.atoms, 0.0, 
                        self.triclinic, self.filter, self.rot, self.rotinv,
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
                    self.triclinic, self.filter, self.rot, self.rotinv,
                    self.boxdims, threshold, nlimit, padding, cells)
                if not finished:
                    raise RuntimeError("Could not find adaptive cutoff")
            else:
                if cells:
                    pc.get_all_neighbors_cells(self.atoms, cutoff,
                        self.triclinic, self.filter, self.rot, self.rotinv, self.boxdims)
                else:
                    pc.get_all_neighbors_normal(self.atoms, cutoff,
                        self.triclinic, self.filter, self.rot, self.rotinv, self.boxdims)

        elif method == 'number':
            finished = pc.get_all_neighbors_bynumber(self.atoms, 0.0, 
                self.triclinic, self.filter, self.rot, self.rotinv,
                self.boxdims, threshold, nmax, cells, assign_neighbor)
            if not finished:
                raise RuntimeError("Could not find enough neighbors - try increasing threshold")

        '''
        elif method == 'voronoi':
            
            self.voroexp = int(voroexp)

            #copy the simulation cell
            backupbox = self._box.copy()
            if self.triclinic:
                if not self.ghosts_created:
                    atoms  = self.get_atoms()
                    atoms = self.repeat((1, 1, 1), atoms=atoms, ghost=True, scale_box=True)
                    self.set_atoms(atoms)
                    self.embed_in_cubic_box()
            #self.embed_in_cubic_box()
            self.get_all_neighbors_voronoi()

            #replace box
            self.box = backupbox
        '''
        self.neighbors_found = True

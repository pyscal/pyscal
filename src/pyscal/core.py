"""


"""


import pyscal.traj_process as ptp
import os
import numpy as np
import warnings
import pyscal.csystem as pc
from pyscal.catom import Atom
import itertools

#------------------------------------------------------------------------------------------------------------
"""
System class definitions
"""
#------------------------------------------------------------------------------------------------------------

def test():
    """
    A simple function to test if the module works

    Parameters
    ----------
    None

    Returns
    -------
    works : bool
        True if the module works and could create a System and Atom object
        False otherwise.

    """
    try:
        s = System()
        a = Atom()
        return True
    except:
        return False

class System(pc.System):
    """
    A python/pybind11 hybrid class for holding the properties of a system.

    Attributes
    ----------
    box : list of list of floats
        A list containing the dimensions of the simulation box in the format
        `[[x_low, x_high], [y_low, y_high], [z_low, z_high]]`

    atoms : list of :class:`~pyscal.catom.Atom` objects

    Notes
    -----
    A `System` consists of two
    major components - the simulation box and the atoms. All the associated variables
    are then calculated using this class.

    .. note::

        atoms can be accessed or set as :attr:`~pyscal.core.System.atoms`. However, due to
        technical reasons individual atoms should be accessed using the
        :func:`~pyscal.core.System.get_atom` method. An atom can be assigned
        to the atom using the :func:`~pyscal.core.System.set_atom` method.

    Examples
    --------
    >>> sys = System()
    >>> sys.read_inputfile('atoms.dat')
    """
    def __init__(self):

        self.initialized = True
        self.neighbors_found = False
        self.neighbor_method = None
        pc.System.__init__(self)

    def read_inputfile(self, filename, format="lammps-dump", frame=-1, compressed = False, customkeys=None, is_triclinic = False):
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

        frame : int
            If the trajectory contains more than one time step, the slice can be specified
            using the `frame` option.

            .. note::

                works only with `lammps-dump` format.

        customkeys : list
            A list containing names of headers of extra data that needs to be read in from the
            input file.

        is_triclinc : bool, optional
            Only used in the case of `format='ase'`. If the read ase object is triclinic, this
            options should be set to True.

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

        `frame` keyword allows to read in a particular slice from a long trajectory. If all slices
        need to analysed, this is a very inefficient way. For handling multiple time slices,
        the :func:`~pyscal.traj_process` module offers a better set of tools.

        Triclinic simulation boxes can also be read in for `lammps-dump`. No special keyword is
        necessary.

        If `custom_keys` are provided, this extra information is read in from input files if
        available. This information is not passed to the C++ instance of atom, and is stored
        as a dictionary. It can be accessed directly as `atom.custom['customval']`


        """
        if customkeys == None:
            customkeys = []

        if format == 'lammps-dump':
            #check customkeys and assign a variable
            customread = (len(customkeys) > 0)

            if frame != -1:
                #split the traj and returns set of filenames
                filenames = ptp.split_traj_lammps_dump(filename, compressed=compressed)

                #reassign filename
                try:
                    filename = filenames[frame]
                except:
                    raise IOError("frame %d is not found in the trajectory"%frame)

                #now if file exists
                if os.path.exists(filename):
                    atoms, boxdims, box, triclinic = ptp.read_lammps_dump(filename, compressed=compressed, check_triclinic=True, box_vectors=True, customkeys=customkeys)
                    self.atoms = atoms
                    self.box = boxdims

                    if triclinic:
                        #we have to input rotation matrix and the inverse rotation matrix
                        rot = box.T
                        rotinv = np.linalg.inv(rot)
                        self.assign_triclinic_params(rot, rotinv)
                else:
                    raise IOError("input file %s not found"%filename)

                #now remove filenames
                for file in filenames:
                    os.remove(file)

            elif os.path.exists(filename):
                atoms, boxdims, box, triclinic = ptp.read_lammps_dump(filename, compressed=compressed, check_triclinic=True, box_vectors=True, customkeys=customkeys)
                self.atoms = atoms
                self.box = boxdims

                if triclinic:
                    rot = box.T
                    rotinv = np.linalg.inv(rot)
                    self.assign_triclinic_params(rot, rotinv)
            else:
                raise IOError("input file %s not found"%filename)


        elif format == 'poscar':
            if os.path.exists(filename):
                atoms, boxdims, box = ptp.read_poscar(filename, compressed=compressed, box_vectors = True)
                self.atoms = atoms
                self.box = boxdims
                if is_triclinic:
                    rot = box.T
                    rotinv = np.linalg.inv(rot)
                    self.assign_triclinic_params(rot, rotinv)
            else:
                raise IOError("input file %s not found"%filename)

        elif format == 'ase':
            atoms, boxdims, box = ptp.read_ase(filename, box_vectors = True)
            self.atoms = atoms
            self.box = boxdims
            if is_triclinic:
                rot = box.T
                rotinv = np.linalg.inv(rot)
                self.assign_triclinic_params(rot, rotinv)

        elif format == 'mdtraj':
            atoms, boxdims, box = ptp.read_mdtraj(filename, box_vectors = True)
            self.atoms = atoms
            self.box = boxdims
            if is_triclinic:
                rot = box.T
                rotinv = np.linalg.inv(rot)
                self.assign_triclinic_params(rot, rotinv)

        else:
            raise TypeError("format recieved an unknown option %s"%format)

    def get_atom(self, index):
        """

        Get the :class:`~pyscal.catom.Atom` object at the queried position in the list of all atoms
        in the :class:`~pyscal.core.System`.

        Parameters
        ----------
        index : int
            index of required atom in the list of all atoms.

        Returns
        -------
        atom : Atom object
            atom object at the queried position.
        """
        atom = self.cget_atom(index)
        return atom

    def set_atom(self, atom):
        """
        Return the atom to its original location after modification.

        Parameters
        ----------
        atom : Atom
                atom to be replaced

        Returns
        -------
        None

        Notes
        -----
        For example, an :class:`~pyscal.catom.Atom` at location `i` in the list of all atoms in
        :class:`~pyscal.core.System` can be queried by,
        ``atom = System.get_atom(i)``, then any kind of modification, for example, the
        position of the `Atom` can done by, ``atom.pos = [2.3, 4.5, 4.5]``. After
        modification, the `Atom` can be set back to its position in `System` by
        :func:`~pyscal.core.System.set_atom`.

        Although the complete list of atoms can be accessed or set using ``atoms = sys.atoms``,
        `get_atom` and `set_atom` functions should be used for accessing individual atoms.
        If an atom already exists at that index in the list, it will be overwritten and will
        lead to loss of information.

        """
        self.cset_atom(atom)


    def calculate_rdf(self, histobins=100, histomin=0.0, histomax=None):
        """
        Calculate the radial distribution function.

        Parameters
        ----------
        histobins : int
            number of bins in the histogram
        histomin : float, optional
            minimum value of the distance histogram. Default 0.0.

        histomax : float, optional
            maximum value of the distance histogram. Default, the maximum value
            in all pair distances is used.

        Returns
        -------
        rdf : array of ints
            Radial distribution function
        r : array of floats
            radius in distance units

        """
        distances = self.get_pairdistances()

        if histomax == None:
            histomax = max(distances)

        hist, bin_edges = np.histogram(distances, bins=histobins, range=(histomin, histomax))
        edgewidth = np.abs(bin_edges[1]-bin_edges[0])
        hist = hist.astype(float)
        r = bin_edges[:-1]

        #get box density
        boxvecs = self.get_boxvecs()
        vol = np.dot(np.cross(boxvecs[0], boxvecs[1]), boxvecs[2])
        natoms = self.nop
        rho = natoms/vol

        shell_vols = (4./3.)*np.pi*((r+edgewidth)**3 - r**3)
        shell_rho = hist/shell_vols
        #now divide to get final value
        rdf = shell_rho/rho

        return rdf, r


    def get_qvals(self, q, averaged = False):
        """
        Get the required q_l (Steinhardt parameter) values of all atoms.

        Parameters
        ----------
        q_l : int or list of ints
            required q_l value with l from 2-12

        averaged : bool, optional
            If True, return the averaged q values, default False

        Returns
        -------
        qvals : list of floats
            list of q_l of all atoms.

        Notes
        -----
        The function returns a list of
        q_l values in the same order as the list of the atoms in the system.

        """
        if isinstance(q, int):
            if q in range(2, 13):
                if averaged:
                    rq = self.cget_aqvals(q)
                else:
                    rq = self.cget_qvals(q)
                return rq
            else:
                raise ValueError("the value of q should be between 2 and 12")

        else:
            for qq in q:
                if not qq in range(2, 13):
                    raise ValueError("the value of q should be between 2 and 12")
            if averaged:
                rq = [ self.cget_aqvals(qq) for qq in q ]
            else:
                rq = [ self.cget_qvals(qq) for qq in q ]
            return rq


    def get_distance(self, atom1, atom2, vector=False):
        """
        Get the distance between two atoms.

        Parameters
        ----------
        atom1 : `Atom` object
                first atom
        atom2 : `Atom` object
                second atom
        vector : bool, optional
            If True, the displacement vector connecting the atoms
            is also returned. default false.

        Returns
        -------
        distance : double
                distance between the first and second atom.

        Notes
        -----
        Periodic boundary conditions are assumed by default.
        """
        if vector:
            displacement_vector = self.get_absdistance_vector(atom1, atom2)
            return self.get_absdistance(atom1, atom2), displacement_vector
        else:
            return self.get_absdistance(atom1, atom2)


    def find_neighbors(self, method='cutoff', cutoff=None, threshold=2, filter=None,
                                            voroexp=1, padding=1.2, nlimit=6, cells=False,
                                                nmax=12):
        """

        Find neighbors of all atoms in the :class:`~pyscal.core.System`.

        Parameters
        ----------
        method : {'cutoff', 'voronoi', 'number'}
            `cutoff` method finds neighbors of an atom within a specified or adaptive cutoff distance from the atom.
            `voronoi` method finds atoms that share a Voronoi polyhedra face with the atom. Default, `cutoff`
            `number` method finds a specified number of closest neighbors to the given atom.

        cutoff : { float, 'sann', 'adaptive'}
            the cutoff distance to be used for the `cutoff` based neighbor calculation method described above.
            If the value is specified as 0 or `adaptive`, adaptive method is used.
            If the value is specified as `sann`, sann algorithm is used.

        threshold : float, optional
            only used if ``cutoff=adaptive``. A threshold which is used as safe limit for calculation of cutoff.

        filter : {'None', 'type'}, optional
            apply a filter to nearest neighbor calculation. If the `filter` keyword is set to
            `type`, only atoms of the same type would be included in the neighbor calculations. Default None.

        voroexp : int, optional
            only used if ``method=voronoi``. Power of the neighbor weight used to weight the contribution of each atom towards
            Steinhardt parameter values. Default 1.

        padding : double, optional
            only used if ``cutoff=adaptive`` or ``cutoff=number``. A safe padding value used after an adaptive cutoff is found. Default 1.2.

        nlimit : int, optional
            only used if ``cutoff=adaptive``. The number of particles to be considered for the calculation of adaptive cutoff.
            Default 6.

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
        
        If method os `number`, instead of using a cutoff value for finding neighbors, a specified number of closest atoms are
        found. This number can be set through the argument `nmax`.

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
        self.reset_allneighbors()
        self.filter = 0

        if filter == 'type':
            # type corresponds to 1
            self.filter = 1

        if method == 'cutoff':
            if cutoff=='sann':
                if threshold < 1:
                    raise ValueError("value of threshold should be at least 1.00")
                self.usecells = (len(self.atoms) > 4000)
                finished = self.get_all_neighbors_sann(threshold)
                #if it finished without finding neighbors
                if not finished:
                    finallydone = False
                    for i in range(1,10):
                        #threshold value is probably too low
                        #try increasing threshold
                        warnings.warn("Could not find sann cutoff. trying with a higher threshold", RuntimeWarning)
                        self.reset_allneighbors()
                        newfinished = self.get_all_neighbors_sann(threshold*i)
                        if newfinished:
                            finallydone = True
                            warnings.warn("found neighbors with higher threshold than default/user input")
                            break

                    if not finallydone:
                        raise RuntimeError("sann cutoff could not be converged. This is most likely, \
                        due to a low threshold value. Try increasing it.")
            elif cutoff=='adaptive' or cutoff==0:
                if threshold < 1:
                    raise ValueError("value of threshold should be at least 1.00")

                self.usecells =  (len(self.atoms) > 4000)
                finished = self.get_all_neighbors_adaptive(threshold, nlimit, padding)
                if not finished:
                    raise RuntimeError("Could not find adaptive cutoff")
            else:
                #warnings.warn("THIS RAN")
                self.set_neighbordistance(cutoff)
                if len(self.atoms) > 2300:
                #if cells:
                    self.get_all_neighbors_cells()
                else:
                    self.get_all_neighbors_normal()

        elif method == 'number':
            if threshold < 1:
                raise ValueError("value of threshold should be at least 1.00")

            self.usecells =  (len(self.atoms) > 4000)
            finished = self.get_all_neighbors_bynumber(threshold, nmax)
            if not finished:
                raise RuntimeError("Could not find enough neighbors - try increasing threshold")

        elif method == 'voronoi':
            self.voroexp = int(voroexp)
            self.get_all_neighbors_voronoi()

        self.neighbors_found = True



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
        self.reset_allneighbors()

    def calculate_vorovector(self, edge_cutoff=0.05, area_cutoff=0.01, edge_length=False):
        """
        get the voronoi structure identification vector.

        Parameters
        ----------

        edge_cutoff : float, optional
            cutoff for edge length. Default 0.05.


        area_cutoff : float, optional
            cutoff for face area. Default 0.01.

        edge_length : bool, optional
            if True, a list of unrefined edge lengths are returned. Default false.

        Returns
        -------
        vorovector : array like, int
            array of the form (n3, n4, n5, n6)

        Notes
        -----
        Returns a vector of the form `(n3, n4, n5, n6)`, where `n3` is the number
        of faces with 3 vertices, `n4` is the number of faces with 4
        vertices and so on. This can be used to identify structures [1] [2].

        The keywords `edge_cutoff` and `area_cutoff` can be used to tune the values to minimise
        the effect of thermal distortions. Edges are only considered in the analysis if the
        `edge_length/sum(edge_lengths)` is at least `edge_cutoff`. Similarly, faces are only
        considered in the analysis if the  `face_area/sum(face_areas)` is at least `face_cutoff`.

        References
        ----------
        .. [1] Finney, JL, Proc. Royal Soc. Lond. A 319, 1970
        .. [2] Tanemura, M, Hiwatari, Y, Matsuda, H,Ogawa, T, Ogita, N, Ueda, A. Prog. Theor. Phys. 58, 1977

        """
        atoms = self.atoms
        for atom in atoms:
            #start looping over and eliminating short edges
            st = 1
            refined_edges = []
            edge_lengths = []

            for vno in atom.face_vertices:

                vphase = atom.vertex_numbers[st:st+vno]
                edgecount = 0
                dummy_edge_lengths = []

                #now calculate the length f each edge
                for i in range(-1, len(vphase)-1):
                    #get pairs of indices
                    #verts are i, i+1
                    ipos = atom.vertex_vectors[vphase[i]*3:vphase[i]*3+3]
                    jpos = atom.vertex_vectors[vphase[i+1]*3:vphase[i+1]*3+3]

                    #now calculate edge length
                    edgeln = np.sqrt((ipos[0]-jpos[0])**2 + (ipos[1]-jpos[1])**2 + (ipos[2]-jpos[2])**2)
                    dummy_edge_lengths.append(edgeln)

                edge_lengths.append(dummy_edge_lengths)
                st += (vno+1)

            #now all the edge lengths are saved
            for c, ed in enumerate(edge_lengths):
                #normalise the edge lengths
                norm = (ed/np.sum(ed))
                #apply face area cutoff
                if (atom.neighbor_weights[c] > area_cutoff):
                    #check for edge length cutoff
                    edgecount = len([cc for cc,x in enumerate(norm) if x > edge_cutoff])
                    refined_edges.append(edgecount)

            #now loop over refined edges and collect n3, n4, n5, n6
            vorovector = [0, 0, 0, 0]

            for ed in refined_edges:
                if ed == 3:
                    vorovector[0] += 1
                elif ed == 4:
                    vorovector[1] += 1
                elif ed == 5:
                    vorovector[2] += 1
                elif ed == 6:
                    vorovector[3] += 1

            atom.edge_lengths = edge_lengths
            atom.vorovector = vorovector
        self.atoms = atoms


    def calculate_q(self, q, averaged = False):
        """
        Find the Steinhardt parameter q_l for all atoms.

        Parameters
        ----------
        q_l : int or list of ints
            A list of all Steinhardt parameters to be found from 2-12.

        averaged : bool, optional
            If True, return the averaged q values, default False

        Returns
        -------
        None

        Notes
        -----
        Enables calculation of the Steinhardt parameters [1] q from 2-12. The type of
        q values depend on the method used to calculate neighbors. See the description
        :func:`~pyscal.core.System.find_neighbors` for more details. If the keyword `average` is set to True,
        the averaged versions of the bond order parameter [2] is returned.

        References
        ----------
        .. [1] Steinhardt, PJ, Nelson, DR, Ronchetti, M. Phys Rev B 28, 1983
        .. [2] Lechner, W, Dellago, C, J Chem Phys, 2013

        """
        if isinstance(q, int):
            qq = [q]
        else:
            qq = q

        for ql in qq:
            if not ql in range(2,13):
                raise ValueError("value of q should be between 2 and 13")

        self.ccalculate_q(qq)

        if averaged:
            self.ccalculate_aq(qq)


    def find_solids(self, bonds=0.5, threshold=0.5, avgthreshold=0.6, cluster=True, q=6, new_algo=False, cutoff=0):
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
        set used the `q` attribute.


        References
        ----------
        .. [1] Auer, S, Frenkel, D. Adv Polym Sci 173, 2005

        """
        #check if neighbors are found
        if not self.neighbors_found:
            raise RuntimeError("neighbors should be calculated before finding solid atoms. Run System.find_neighbors.")

        if not isinstance(q, int):
            raise TypeError("q should be interger value")
        else:
            if not ((q >= 2 ) and (q <= 12 )):
                raise ValueError("Value of q should be between 2 and 12")

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
            self.criteria = 0
        elif isinstance(bonds, float):
            if ((bonds>=0) and (bonds<=1.0)):
                self.criteria = 1
            else:
                raise TypeError("bonds if float should have value between 0-1")
        else:
             raise TypeError("bonds should be interger/float value")

        #Set the vlaue of q
        self.solidq = q
        #first calculate q
        self.ccalculate_q([q])
        #self.calculate_q(6)
        #calculate solid neighs
        self.set_nucsize_parameters(bonds, threshold, avgthreshold)
        self.calculate_frenkelnumbers()
        #now find solids
        self.find_solid_atoms()

        if cluster:
            def ccondition(atom):
                return atom.solid

            lc = self.cluster_atoms(ccondition, largest=True, new_algo=new_algo, cutoff=cutoff)
            return lc


    def cluster_atoms(self, condition, largest = True, new_algo = False, cutoff=0):
        """
        Cluster atoms based on a property

        Parameters
        ----------
        condition : callable or atom property
            Either function which should take an :class:`~Atom` object, and give a True/False output
            or an attribute of atom class which has value or 1 or 0.

        largest : bool, optional
            If True returns the size of the largest cluster. Default False.

        Returns
        -------
        lc : int
            Size of the largest cluster. Returned only if `largest` is True.


        Notes
        -----
        This function helps to cluster atoms based on a defined property. This property
        is defined by the user through the argument `condition` which is passed as a parameter.
        `condition` can be of two types. The first type is a function which takes an 
        :class:`~Atom` object and should give a True/False value. `condition` can also be an
        :class:`~Atom` attribute or a value from `custom` values stored in an atom.


        When clustering, the code loops over each atom and its neighbors. If the
        `condition` is true for both host atom and the neighbor, they are assigned to
        the same cluster. For example, a condition to cluster solid atoms would be,

        .. code:: python

            def condition(atom):
                #if both atom is solid
                return (atom1.solid)

        The same can be done by passing `"solid"` as the condition argument instead of the above
        function. Passing a function allows to evaluate complex conditions, but is slower than
        passing an attribute.

        """
        testatom = self.atoms[0]

        #test the condition
        isatomattr = False

        try:
            out = condition(testatom)
            if out not in [True, False, 0, 1]:
                raise RuntimeError("The output of condition should be either True or False. Received %s"%str(out))

        except:
            try:
                out = self.get_custom(testatom, [condition])[0]
                if out not in [True, False, 0, 1]:
                    raise RuntimeError("The output of condition should be either True or False. Received %s"%str(out))
                isatomattr = True        
            except:
                raise RuntimeError("condition did not work")
        
        #now loop
        atoms = self.atoms

        if isatomattr:
            for atom in atoms:
                atom.condition = self.get_custom(atom, [condition])[0]
        else:
            for atom in atoms:
                cval = condition(atom)
                atom.condition = cval
        
        self.atoms = atoms
        self.cfind_clusters_recursive(cutoff)

        #done!
        lc = self.find_largest_cluster()
        #pcs.System.get_largest_cluster_atoms(self)

        if largest:
            return lc


    def calculate_solidneighbors(self):
        """
        Find Solid neighbors of all atoms in the system.

        Parameters
        ----------
        None

        Returns
        -------
        None

        Notes
        -----
        A solid bond is considered between two atoms if the `connection <https://pyscal.readthedocs.io/en/latest/solidliquid.html>`_
        between them is greater than 0.6.


        """
        self.calculate_frenkelnumbers()

    def find_clusters(self, recursive = True, largest = True):
        """
        Find the clusters of all atoms in the system.

        Parameters
        ----------
        recursive : Bool, optional
            If True, use a recursive clustering algorithm, otherwise use an id based clustering.
            The difference in values between two methods can be upto 3 particles. Default True.

        largest : Bool, optional
            If True, return the number of particles in the largest cluster. Default True.

        Returns
        -------
        cluster : int
            The size of the largest cluster in the system. Only returned if `largest` is set to True.

        Notes
        -----
        Go through all the atoms in the system and cluster them together based on the `issolid` parameter of the atom.
        To cluster based on any user defined criteria, you can use `set_solid` method of `Atom` to explicitely
        set the `issolid` value.

        .. warning::

            This function is deprecated and will be removed in a future release. Please use
            :func:`~System.cluster_atoms` instead.

        """
        warnings.simplefilter('always', DeprecationWarning)
        warnings.warn("This function is deprecated - use cluster_atoms instead", DeprecationWarning)

        if recursive:
            self.cfind_clusters_recursive()
        else:
            self.cfind_clusters()

        if largest:
            cluster = self.find_largest_cluster()
            return cluster

    def find_largestcluster(self):
        """
        Find the largest solid cluster of atoms in the system from all the clusters.

        Parameters
        ----------
        None

        Returns
        -------
        cluster : int
            the size of the largest cluster

        Notes
        -----
        :func:`pyscal.core.System.find_clusters` has to be used before using this function.

        """
        return self.find_largest_cluster()


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
        The calculated A parameter for each atom is stored in :attr:`~pyscal.catom.Atom.angular`.

        References
        ----------
        .. [1] Uttormark, MJ, Thompson, MO, Clancy, P, Phys. Rev. B 47, 1993
        """

        atoms = self.atoms

        for atom in atoms:
            dists = []
            distneighs = []
            distvectors = []

            neighs = atom.neighbors

            for neigh in neighs:
                dist, vectors = self.get_distance(atom, atoms[neigh], vector=True)
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
            atom.angular = costhetasum

        self.atoms = atoms

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

        bins = [-1.0, -0.945, -0.915, -0.755, -0.705, -0.195, 0.195, 0.245, 0.795, 1.0]
        atoms = self.atoms

        for atom in atoms:
            dists = []
            distneighs = []
            distvectors = []

            #dists = atom.neighbor_distance
            #distneighs = atom.neighbors
            #distvectors = atom.

            for neigh in neighs:
                dist, vectors = self.get_distance(atom, atoms[neigh], vector=True)
                dists.append(dist)
                distneighs.append(neigh)
                distvectors.append(vectors)

            args = np.argsort(dists)
            topargs = np.array(args)

            combos = list(itertools.combinations(topargs, 2))
            costhetas = []
            for combo in combos:
                vec1 = distvectors[combo[0]]
                vec2 = distvectors[combo[1]]
                modvec1 = np.sqrt(np.sum([x**2 for x in vec1]))
                modvec2 = np.sqrt(np.sum([x**2 for x in vec2]))
                costheta = np.dot(vec1, vec2)/(modvec1*modvec2)
                #found costheta
                costhetas.append(costheta)


            #now add according to classification in paper
            chivector = np.histogram(costhetas, bins=bins)
            atom.chiparams = chivector[0]
            if angles:
                atom.custom['cosines'] = costhetas


        self.atoms = atoms

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

        The keyword `averaged` is True, the disorder value is averaged over the atom and its neighbors. The disorder value
        can be accessed using :attr:`~pyscal.catom.disorder` and the averaged version can be accessed using
        :attr:`~pyscal.catom.avg_disorder`. For ordered systems, the value of disorder would be zero which would increase
        and reach one for disordered systems.

        References
        ----------
        .. [1] Kawasaki, T, Onuki, A, J. Chem. Phys. 135, 2011

        """
        #now routine for calculation of disorder
        if q in range(2,13):
            self.solidq = q
        else:
            raise ValueError("q values should be between 2-12")

        self.ccalculate_disorder()

        if averaged:
            self.ccalculate_avg_disorder()


    def calculate_sro(self, reference_type=1, average=True, shells=2):
        """
        Calculate short range order

        Parameters
        ----------
        reference_type: int, optional
            type of the atom to be used a reference. default 1

        average: bool, optional
            if True, average over all atoms of the reference type in the system.
            default True.

        Returns
        -------
        vec: list of float
            The short range order averaged over the whole system for atom of
            the reference type. Only returned if `average` is True. First value is SRO
            of the first neighbor shell and the second value corresponds to the second
            nearest neighbor shell.

        Notes
        -----
        Calculates the short range order for an AB alloy using the approach by
        Cowley [1]. Short range order is calculated as,

        .. math::

            \\alpha_i = 1 - \\frac{n_i}{m_A c_i}

        where n_i is the number of atoms of the non reference type among the c_i atoms
        in the ith shell. m_A is the concentration of the non reference atom. Please
        note that the value is calculated for shells 1 and 2 by default. In order for
        this to be possible, neighbors have to be found first using the :func:`~pyscal.core.System.find_neighbors`
        method. The selected neighbor method should include the second shell as well. For this
        purpose `method=cutoff` can be chosen with a cutoff long enough to include the second
        shell. In order to estimate this cutoff, one can use the :func:`~pyscal.core.System.calculate_rdf`
        method.

        References
        ----------
        .. [1] Cowley J. M., PR 77(5), 1950.

        """
        if not self.neighbors_found:
            raise RuntimeError("Neighbors not found, please find neighbors using the cutoff method")
        if not reference_type in [1,2]:
            raise ValueError("reference atom type should be either 1 or 2")

        atoms = self.atoms

        try:
            type1 = len([1 for atom in atoms if atom.type == 1])
            type2 = len([1 for atom in atoms if atom.type == 2])
            concv = [type1, type2]
        except:
            raise RuntimeError("There should be two atom types")
        if not ((type1>0) and (type2>0)):
            raise RuntimeError("There should be two atom types")

        mref = concv[reference_type-1]/float(np.sum(concv))
        motr = 1 - mref

        for atom in atoms:
            if atom.type == reference_type:
                neighs = atom.neighbors
                #get all neighbor distances
                distances = [self.get_distance(atom, atoms[n]) for n in neighs]
                avgdistance = np.mean(distances)
                distsplit = avgdistance/2.00
                #find two shells of atoms
                if shells == 2:
                    set1 = [neighs[n] for n, dist in enumerate(distances) if dist <= avgdistance]
                    set2 = [neighs[n] for n, dist in enumerate(distances) if dist > avgdistance]
                    #now evaluate types of atoms
                    set1ref = np.sum([1 for n in set1 if atoms[n].type == reference_type])
                    set1otr = len(set1) - set1ref
                    set2ref = np.sum([1 for n in set2 if atoms[n].type == reference_type])
                    set2otr = len(set2) - set2ref
                    #now calculate values
                    shell1 = 1 - (set1otr/(len(set1)*motr))
                    shell2 = 1 - (set2otr/(len(set2)*motr))
                    atom.sro = [shell1, shell2]
                elif shells == 1:
                    set1ref = np.sum([1 for n in neighs if atoms[n].type == reference_type])
                    set1otr = len(neighs) - set1ref
                    shell1 = 1 - (set1otr/(len(neighs)*motr))
                    atom.sro = [shell1]

        #add atoms
        self.atoms = atoms
        #now if avg is reqd, find it
        if average:
            vec = np.zeros(2)
            count = 0
            for atom in atoms:
                if atom.type == reference_type:
                    vec += np.array(atom.sro)
                    count += 1
            return vec/float(count)

    def get_custom(self, atom, customkeys):
        """
        Get a custom attribute from Atom

        Parameters
        ----------
        atom : Atom object

        customkeys : list of strings
            the list of keys to be found

        Returns
        -------
        vals : list
            array of custom values

        """
        #first option - maybe it appears
        vals = []
        for ckey in customkeys:
            #if the key is there - ignore it
            if not ckey in atom.custom.keys():
                #try to get from attribute
                try:
                    val = getattr(atom, ckey)
                    vals.append(val)

                except AttributeError:
                    #since attr failed, check if they are q or aq values
                    if ckey[0] == 'q':
                        qkey = ckey[1:]
                        #try to acess this value
                        val = atom.get_q(int(qkey))
                        #add this pair to dict - as string vals
                        vals.append(val)

                    elif ckey[:2] == 'aq':
                        qkey = ckey[2:]
                        val = atom.get_q(int(qkey), averaged=True)
                        vals.append(val)

                    else:
                        raise AttributeError("custom key was not found")
            else:
                val = atom.custom[ckey]
                vals.append(val)
        return vals

    def to_file(self, outfile, format='lammps-dump', customkeys=None, compressed=False, timestep=0):
        """
        Save the system instance to a trajectory file.

         Parameters
        ----------
        outfile : string
            name of the output file

        format : string, optional
            format of the output file, default `lammps-dump`
            Currently only `lammps-dump` format is supported.

        customkeys : list of strings, optional
            a list of extra atom wise values to be written in the output file.

        compressed : bool, optional
            If true, the output is written as a compressed file.

        timestep : int, optional
            timestep to be written to file. default 0

        Returns
        -------
        None

        Notes
        -----

        """
        if customkeys == None:
            customkeys = []



        boxdims = self.box
        atoms = self.atoms

        if len(customkeys) > 0:
            cvals = [self.get_custom(atom, customkeys) for atom in atoms]

        #open files for writing
        if compressed:
            gz = gzip.open(outfile,'w')
            dump = io.BufferedReader(gz)
        else:
            gz = open(outfile,'w')
            dump = gz

        #now write
        dump.write("ITEM: TIMESTEP\n")
        dump.write("%d\n" % timestep)
        dump.write("ITEM: NUMBER OF ATOMS\n")
        dump.write("%d\n" % len(atoms))
        dump.write("ITEM: BOX BOUNDS\n")
        dump.write("%f %f\n" % (boxdims[0][0], boxdims[0][1]))
        dump.write("%f %f\n" % (boxdims[1][0], boxdims[1][1]))
        dump.write("%f %f\n" % (boxdims[2][0], boxdims[2][1]))

        #now write header
        if len(customkeys) > 0:
            ckey = " ".join(customkeys)
            title_str = "ITEM: ATOMS id type x y z %s\n"% ckey
        else:
            title_str = "ITEM: ATOMS id type x y z\n"

        dump.write(title_str)

        for cc, atom in enumerate(atoms):
            pos = atom.pos
            if len(customkeys) > 0:
                cval_atom = " ".join(np.array(list(cvals[cc])).astype(str))
                atomline = ("%d %d %f %f %f %s\n")%(atom.id, atom.type, pos[0], pos[1], pos[2], cval_atom)
            else:
                atomline = ("%d %d %f %f %f\n")%(atom.id, atom.type, pos[0], pos[1], pos[2])

            dump.write(atomline)

        dump.close()

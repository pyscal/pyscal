"""
Main module of pyscal. This module contains definitions of the two major
classes in pyscal - the :class:`~System` and :class:`~Atom`. Both classes
inherit the C++ classes of the same name and provide additional functionality
by combining various methods, validation etc.

"""


import pyscal.traj_process as ptp
import pyscal.pickle as pp
import os
import numpy as np
import warnings
import pyscal.csystem as pc
from pyscal.catom import Atom

#------------------------------------------------------------------------------------------------------------
"""
System class definitions
"""
#------------------------------------------------------------------------------------------------------------


class System(pc.System):
    """
    A c++ class for holding the properties of a system.

    Attributes
    ----------
    box : list of list of floats
        A list containing the dimensions of the simulation box in the format
        `[[x_low, x_high], [y_low, y_high], [z_low, z_high]]`

    atoms : list of :class:`~pyscal.core.Atom` objects

        .. note::

            atoms can be accessed or set as `System.atoms`. However, due to
            technical reasons individual atoms should be accessed using the
            :func:`~pyscal.core.System.get_atom` method. An atom can be assigned
            to the atom using the :func:`~pyscal.core.System.set_atom` method.

    Notes
    -----
    A `System` consists of two
    major components - the simulation box and the atoms. All the associated variables
    are then calculated using this class.

    Examples
    --------
    >>> sys = System()
    >>> sys.read_inputfile('atoms.dat')
    """
    def __init__(self):

        self.initialized = True
        self.neighbors_found = False
        pc.System.__init__(self)

    def read_inputfile(self, filename, format="lammps-dump", frame=-1, compressed = False, customkeys=[]):
        """

        Read input file that contains the information of system configuration.

        Parameters
        ----------
        filename : string
            name of the input file.

        format : {'lammps-dump', 'poscar'}
            format of the input file

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

        Returns
        -------
        None

        Notes
        -----
        `format` keyword specifies the format of the input file. Currently only
        a `lammps-dump` and `poscar` files are supported.  This function uses the
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
                atoms, boxdims = ptp.read_poscar(filename, compressed=compressed)
                self.atoms = atoms
                self.box = boxdims
            else:
                raise IOError("input file %s not found"%filename)
        else:
            raise TypeError("format recieved an unknown option %s"%format)

    def get_atom(self, index):
        """

        Get the :class:`~pyscal.core.Atom` object at the queried position in the list of all atoms
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
        For example, an :class:`~pyscal.core.Atom` at location `i` in the list of all atoms in
        :class:`~pyscal.core.System` can be queried by,
        ``atom = System.get_atom(i)``, then any kind of modification, for example, the
        position of the `Atom` can done by, ``atom.set_pos([2.3, 4.5, 4.5])``. After
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


    def get_distance(self, atom1, atom2):
        """
        Get the distance between two atoms.

        Parameters
        ----------
        atom1 : `Atom` object
                first atom
        atom2 : `Atom` object
                second atom

        Returns
        -------
        distance : double
                distance between the first and second atom.

        Notes
        -----
        Periodic boundary conditions are assumed by default.
        """
        return self.get_absdistance(atom1, atom2)


    def find_neighbors(self, method='cutoff', cutoff=None, threshold=2, filter=None,
                                            voroexp=1, face_cutoff=0.002, padding=1.2, nlimit=6):
        """

        Find neighbors of all atoms in the :class:`~pyscal.core.System`.

        Parameters
        ----------
        method : {'cutoff', 'voronoi'}
            `cutoff` method finds neighbors of an atom within a specified or adaptive cutoff distance from the atom.
            `voronoi` method finds atoms that share a Voronoi polyhedra face with the atom. Default, `cutoff`

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

        face_cutoff : double, optional
            only used if ``method=voronoi``. The minimum fraction of total voronoi face area a single phase should have in order to
            include it in the analysis of voronoi polyhedra to find `(n_3, n_4, n_5, n_6)` vector. Default 0.002

        padding : double, optional
            only used if ``cutoff=adaptive``. A safe padding value used after an adaptive cutoff is found. Default 1.2.

        nlimit : int, optional
            only used if ``cutoff=adaptive``. The number of particles to be considered for the calculation of adaptive cutoff.
            Default 6.

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

        In Method cutoff, if ``cutoff='adaptive'``, an adaptive cutoff is found during runtime for each atom [1]_.
        Setting the cutoff radius to 0 also uses this algorithm. The cutoff for an atom i is found using,

        .. math:: r_c(i) = padding * ((1/nlimit) * \sum_{j=1}^{nlimit}(r_{ij}))

        padding is a safe multiplier to the cutoff distance that can be set through the keyword `padding`. `nlimit` keyword sets the
        limit for the top nlimit atoms to be taken into account to calculate the cutoff radius.

        In Method cutoff, if ``cutoff='sann'``, sann algorithm is used [2]_. There are no parameters to tune sann algorithm.

        The second approach is using Voronoi polyhedra which also assigns a weight to each neighbor in the ratio of the face area between the two atoms.
        Higher powers of this weight can also be used [3]_. The keyword `voroexp`
        can be used to set this weight.

        .. warning::

            Adaptive cutoff uses a padding over the intial guessed "neighbor distance". By default it is 2. In case
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
                finished = self.get_all_neighbors_adaptive(threshold, nlimit, padding)
                if not finished:
                    raise RuntimeError("Could not find adaptive cutoff")
            else:
                #warnings.warn("THIS RAN")
                self.set_neighbordistance(cutoff)
                self.get_all_neighbors_normal()

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
        vertices and so on. This can be used to identify structures [1]_ [2]_.

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
        Enables calculation of the Steinhardt parameters [1]_ q from 2-12. The type of
        q values depend on the method used to calculate neighbors. See the description
        :func:`~pyscal.core.System.find_neighbors` for more details. If the keyword `average` is set to True,
        the averaged versions of the bond order parameter [2]_ is returned.

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


    def find_solids(self, bonds=7, threshold=0.5, avgthreshold=0.6, cluster=True):
        """
        Distinguish solid and liquid atoms in the system.

        Parameters
        ----------
        bonds : int, optional
            Minimum number of solid bonds for an atom to be identified as
            a solid. Default 7.

        threshold : double, optional
            Solid bond cutoff value. Default 0.5.

        avgthreshold : double, optional
            Value required for Averaged solid bond cutoff for an atom to be identified
            as solid. Default 0.6.

        cluster : bool, optional
            If True, cluster the solid atoms and return the number of atoms in the largest
            cluster.

        Returns
        -------
        solid : int
            Size of the largest solid cluster. Returned only if `cluster=True`.

        Notes
        -----
        The neighbors should be calculated before running this function.
        Check :func:`~pyscal.core.System.find_neighbors` method.

        `bonds` define the number of solid bonds of an atom to be identified as solid.
        Two particles are said to be 'bonded' if [1]_,

        .. math:: s_{ij} = \sum_{m=-6}^6 q_{6m}(i) q_{6m}^*(i) \geq threshold

        where `threshold` values is also an optional parameter.

        An additional parameter `avgthreshold` is an additional parameter to improve solid-liquid distinction.
        In addition to having a the specified number of `bonds`,

        .. math::  \langle s_{ij} \\rangle > avgthreshold

        also needs to be satisfied.


        References
        ----------
        .. [1] Auer, S, Frenkel, D. Adv Polym Sci 173, 2005

        """
        #check if neighbors are found
        if not self.neighbors_found:
            raise RuntimeError("neighbors should be calculated before finding solid atoms. Run System.find_neighbors.")

        if not isinstance(bonds, int):
            raise TypeError("bonds should be interger value")

        if not isinstance(threshold, float):
            raise TypeError("threshold should be a float value")
        else:
            if not ((threshold >= 0 ) and (threshold <= 1 )):
                raise ValueError("Value of threshold should be between 0 and 1")

        if not isinstance(avgthreshold, float):
            raise TypeError("avgthreshold should be a float value")
        else:
            if not ((avgthreshold >= 0 ) and (avgthreshold <= 1 )):
                raise ValueError("Value of avgthreshold should be between 0 and 1")

        #start identification routine
        #first calculate q
        self.ccalculate_q([6])
        #self.calculate_q(6)
        #calculate solid neighs
        self.set_nucsize_parameters(bonds, threshold, avgthreshold)
        self.calculate_frenkelnumbers()
        #now find solids
        self.find_solid_atoms()

        if cluster:
            def ccondition(atom):
                return atom.solid

            lc = self.cluster_atoms(ccondition, largest=True)
            return lc


    def cluster_atoms(self, condition, largest = True):
        """
        Cluster atoms based on a property

        Parameters
        ----------
        condition : callable
            function which should take an :class:`~Atom` object, and give a True/False output

        largest : bool, optional
            If True returns the size of the largest cluster. Default False.

        Returns
        -------
        lc : int
            Size of the largest cluster. Returned only if `largest` is True.


        Notes
        -----
        This function helps to cluster atoms based on a defined property. This property
        is defined by the user through the function `condition` which is passed as a parameter.
        For each atom in the system, the `condition` should give a True/False values.

        When clustering, the code loops over each atom and its neighbors. If the
        `condition` is true for both host atom and the neighbor, they are assigned to
        the same cluster. For example, a condition to cluster solid atoms would be,

        .. code:: python

            def condition(atom):
                #if both atom is solid
                return (atom1.solid)

        Check examples for more details.

        """
        testatom = self.atoms[0]

        #test the condition
        try:
            out = condition(testatom)
            if out not in [True, False]:
                raise RuntimeError("The output of condition should be either True or False. Received %s"%str(out))
        except:
            raise RuntimeError("condition did not work")

        #now loop
        atoms = self.atoms
        for atom in atoms:
            cval = condition(atom)
            atom.condition = cval
        self.atoms = atoms

        self.cfind_clusters_recursive()

        #done!
        lc = self.find_largest_cluster()
        #pcs.System.get_largest_cluster_atoms(self)

        if largest:
            return lc




    def calculate_nucsize(self, frenkelnums, threshold, avgthreshold):
        """
        Calculate the size of the largest cluster in the given system.

        Parameters
        ----------
        None

        Returns
        -------
        cluster size : int
            size of the largest solid cluster in liquid (number of atoms)

        Notes
        -----
        Calculation of the the size of the largest solid cluster needs various prerequisites that can be set
        by the functions `set_nucsize_parameters`.

        .. warning::

            This function is deprecated and will be removed in a future release. Please use
            :func:`System.find_solids` instead.

        """
        #print("this raaan")
        #warnings.simplefilter('always', DeprecationWarning)
        #warnings.warn("This function is deprecated - use find_solids instead", DeprecationWarning)

        self.ccalculate_q([6])
        self.set_nucsize_parameters(frenkelnums, threshold, avgthreshold)
        return self.ccalculate_nucsize()

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
        A solid bond is considered between two atoms if the connection
        betweem them is greater than 0.6.

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


    def prepare_pickle(self):
        """
        Prepare the system for pickling and create a picklable system

        Parameters
        ----------
        None

        Returns
        -------
        psys : picklable system object

        Notes
        -----
        This function prepares the system object for pickling. From
        a user perspective, the :func:`~pyscal.core.System.to_pickle` method should be used
        directly.

        See also
        --------
        to_file()

        """

        #get the basic system indicators
        indicators = self.get_indicators()
        #get box dims and triclinic params if triclinic
        box = self.box
        if indicators[6] == 1:
            rot = self.get_triclinic_params()
        else:
            rot = 0

        #now finally get atoms
        atoms = self.atoms
        #convert them to picklabale atoms
        patoms = [pp.pickle_atom(atom) for atom in atoms]

        #create System instance and assign things
        psys = pp.pickleSystem()
        psys.indicators = indicators
        psys.atoms = patoms
        psys.box = box
        psys.rot = rot

        return psys

    def to_pickle(self, file):
        """
        Save a system to file

        Parameters
        ----------
        file : string
            name of output file

        Returns
        -------
        None

        Notes
        -----
        This function can be used to save a :class:`~pyscal.core.System` object directly to
        file. This retains all the calculated quantities of the system,
        including the atoms and their properties. This can be useful to
        restart the calculation. The function uses `numpy.save` method
        to save the information. Hence pickling between different versions
        of python could lead to issues.

        .. warning::

            Pickling between different versions of numpy or python could be incompatible.

        """
        psys = self.prepare_pickle()
        np.save(file, psys, allow_pickle=True)


    def from_pickle(self, file):
        """
        Read the contents of :class:`~pyscal.core.System` object from a
        pickle file.

        Parameters
        ----------
        file : string
            name of input file

        Returns
        -------
        None

        Notes
        -----
        This function can be used to set up a system
        from a file. A system object needs to be created first.

        Examples
        --------

        >>> sys = System()
        >>> sys.from_pickle(filename)

        """
        if os.path.exists(file):
            psys = np.load(file, allow_pickle=True).flatten()[0]
        else:
            raise IOError("file does not exist")
        #set up indicators
        self.set_indicators(psys.indicators)
        #unpickle atoms
        self.atoms = [pp.unpickle_atom(atom) for atom in psys.atoms]
        self.box = psys.box


        #if triclinic, get those
        if psys.indicators[6] == 1:
            rot = psys.rot
            rotinv = np.linalg.inv(rot)
            self.assign_triclinic_params(rot, rotinv)

    def to_file(self, outfile, format='lammps-dump', custom=[], compressed=False):
        """
        Save the system instance to a trajectory file.

        Parameters
        ----------
        outfile : string
            name of the output file

        format : string, optional
            format of the output file, default `lammps-dump`
            Currently only `lammps-dump` format is supported.

        custom : list of strings, optional
            a list of extra atom wise values to be written in the output file.

        compressed : bool, optional
            If true, the output is written as a compressed file.

        Returns
        -------
        None

        Notes
        -----

        """

        def get_custom(atom, customkeys):
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


        boxdims = self.box
        atoms = self.atoms

        if len(custom) > 0:
            cvals = [get_custom(atom, custom) for atom in atoms]

        #open files for writing
        if compressed:
            gz = gzip.open(outfile,'w')
            dump = io.BufferedReader(gz)
        else:
            gz = open(outfile,'w')
            dump = gz

        #now write
        dump.write("ITEM: TIMESTEP\n")
        dump.write("0\n")
        dump.write("ITEM: NUMBER OF ATOMS\n")
        dump.write("%d\n" % len(atoms))
        dump.write("ITEM: BOX BOUNDS\n")
        dump.write("%f %f\n" % (boxdims[0][0], boxdims[0][1]))
        dump.write("%f %f\n" % (boxdims[1][0], boxdims[1][1]))
        dump.write("%f %f\n" % (boxdims[2][0], boxdims[2][1]))

        #now write header
        if len(custom) > 0:
            ckey = " ".join(custom)
            title_str = "ITEM: ATOMS id type x y z %s\n"% ckey
        else:
            title_str = "ITEM: ATOMS id type x y z\n"

        dump.write(title_str)

        for cc, atom in enumerate(atoms):
            pos = atom.pos
            if len(custom) > 0:
                cval_atom = " ".join(np.array(list(cvals[cc])).astype(str))
                atomline = ("%d %d %f %f %f %s\n")%(atom.id, atom.type, pos[0], pos[1], pos[2], cval_atom)
            else:
                atomline = ("%d %d %f %f %f\n")%(atom.id, atom.type, pos[0], pos[1], pos[2])

            dump.write(atomline)

        dump.close()

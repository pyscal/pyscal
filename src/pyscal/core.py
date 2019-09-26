"""
Main module of pyscal. This module contains definitions of the two major
classes in pyscal - the :class:`~System` and :class:`~Atom`. Both classes
inherit the C++ classes of the same name and provide additional functionality
by combining various methods, validation etc.

"""

import pyscal.catom as pca
import pyscal.csystem as pcs
import pyscal.traj_process as ptp
import pyscal.pickle as pp
import os
import numpy as np
import warnings

#------------------------------------------------------------------------------------------------------------
"""
System class definitions
"""
#------------------------------------------------------------------------------------------------------------


class System(pcs.System):
    """
    A c++ class for holding the properties of a system.

    Notes
    -----
    A `System` consists of two
    major components - the simulation box and the atoms. All the associated variables
    are then calculated over these.

    Examples
    --------
    >>> sys = System()
    >>> sys.read_inputfile()
    """
    def __init__(self):
        self.initialized = True
        #self.solid_params_set = False
        self.neighbors_found = False
        self.atoms = []
        self.nop = 0
        #this method can be done more
        #we can remove checks on the cpp side
        pcs.System.__init__(self)

    #try a get attribute
    def __getattr__(self, name):
        #now check if the first letter is q

        if name == 'box':
            #then extract the next two parts
            box6dim = pcs.System.get_box(self)
            pbox = [[box6dim[0], box6dim[1]], [box6dim[2], box6dim[3]], [box6dim[4], box6dim[5]]]
            return pbox

        elif name == 'box_vectors':
            pbox = pcs.System.get_boxvecs(self)
            return pbox

        else:
            raise AttributeError(name)

    #overload the setattr function to overload
    def __setattr__(self, variable, value):
        if variable == 'box':
            #check for length
            if not len(value) == 3:
                raise TypeError('box should be of the form [[box_x_low, box_x_high], [box_y_low, box_y_high], [box_z_low, box_z_high]]')

            if not all(all(isinstance(x, (int, float)) for x in xx) for xx in value):
                raise TypeError("all values must be float")

            pcs.System.set_box(self, value)

        elif variable == 'box_vectors':
            raise AttributeError("box_vectors are calculated from box. If using triclinic, read in a file instead.")

        #finally assign the variables
        super(System, self).__setattr__(variable, value)


    def read_inputfile(self, filename, format="lammps-dump", frame=-1, compressed = False, customkeys=[]):
        """

        Read input file containing the information of a time slice.

        Parameters
        ----------
        filename : string
            name of the input file to be read in

        format : {'lammps-dump', 'poscar'}
            format of the input file

        compressed : bool, optional
            If True, force to read a `gz` compressed format, default False.

        frame : int
            If the trajectory contains more than one time slice, the slice can be specified
            using the `frame` option.
            Alert: works only with `lammps-dump` format.

        customkeys : list
            A list containing names of headers of extra data that needs to be read in from the
            input file.

        Returns
        -------
        None

        Notes
        -----
        `format` keyword specifies the format of the input file. Currently only
        a `lammps-dump` and `poscar` files are supported. However, this restriction can easily
        be overcome using the :func:`~System.assign_atoms` method from system where a list of atoms
        and box vectors are directly provided to the system. This function itself uses the
        :func:`~pyscal.traj_process` module to process a file which is then assigned to system
        using :func:`~System.assign_atoms`.

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

        See Also
        --------
        assign_atoms
        """
        self.customkeys = customkeys

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
                    raise FileNotFoundError("frame %d is not found in the trajectory"%frame)

                #now if file exists
                if os.path.exists(filename):
                    tatoms, boxdims, box, triclinic = ptp.read_lammps_dump(filename, compressed=compressed, check_triclinic=True, box_vectors=True, customkeys=customkeys)
                    #self.tatoms = tatoms
                    #convert these atoms to normal atoms and save them
                    self.atoms = [self.copy_atom(tatom) for tatom in tatoms]
                    #now get catoms
                    catoms = [self.copy_atom_to_catom(atom) for atom in self.atoms]
                    pcs.System.assign_particles(self, catoms, boxdims)
                    self.nop = len(catoms)

                    if triclinic:
                        #we have to input rotation matrix and the inverse rotation matrix
                        rot = box.T
                        rotinv = np.linalg.inv(rot)
                        pcs.System.assign_triclinic_params(self, rot, rotinv)
                else:
                    raise FileNotFoundError("input file %s not found"%filename)

                #now remove filenames
                for file in filenames:
                    os.remove(file)

            elif os.path.exists(filename):
                tatoms, boxdims, box, triclinic = ptp.read_lammps_dump(filename, compressed=compressed, check_triclinic=True, box_vectors=True, customkeys=customkeys)
                #convert these atoms to normal atoms and save them
                self.atoms = [self.copy_atom(tatom) for tatom in tatoms]
                #now get catoms
                catoms = [self.copy_atom_to_catom(atom) for atom in self.atoms]

                pcs.System.assign_particles(self, catoms, boxdims)
                self.nop = len(catoms)

                if triclinic:
                    #we have to input rotation matrix and the inverse rotation matrix
                    rot = box.T
                    rotinv = np.linalg.inv(rot)
                    pcs.System.assign_triclinic_params(self, rot, rotinv)
            else:
                raise FileNotFoundError("input file %s not found"%filename)


        elif format == 'poscar':
            if os.path.exists(filename):
                tatoms, boxdims = ptp.read_poscar(filename, compressed=compressed)
                #convert these atoms to normal atoms and save them
                self.atoms = [self.copy_atom(tatom) for tatom in tatoms]
                #now get catoms
                catoms = [self.copy_atom_to_catom(atom) for atom in self.atoms]

                pcs.System.assign_particles(self, catoms, boxdims)
                self.nop = len(catoms)

            else:
                raise FileNotFoundError("input file %s not found"%filename)
        else:
            raise TypeError("format recieved an unknown option %s"%format)

    def copy_atom(self, atom1):
        """
        Copies the essential infomation from one atom to another

        Parameters
        ----------
        atom1 : atom object
            source atom

        Returns
        -------
        atom2 : copied atom

        Notes
        -----
        This method only copies the essential information. Positions, type and id
        and customvals
        """
        atom2 = Atom()
        atom2.pos = atom1.pos
        atom2.loc = atom1.loc
        atom2.id = atom1.id
        atom2.type = atom1.type
        atom2.custom = atom1.custom
        return atom2

    def assign_atoms(self, atoms, box):
        """

        Assign atoms and box vectors to :class:`~System`.

        Parameters
        ----------
        atoms : list of `Atom` objects
            list consisting of all atoms
        box   : list of list of floats
            list which consists of the box dimensions in the format-
            `[[box_x_low, box_x_high], [box_y_low, box_y_high], [box_z_low, box_z_high]]`

        Returns
        -------
        None

        Notes
        -----
        Receive a vector of atom objects which is stored instead
        of reading in the input file. If this method is used, there is no need of using the
        :func:`~System.read_inputfile` method. Also using this function allows for reading of multiple
        file formats which are not supported by the inbuilt :func:`~System.read_inputfile` method.

        See Also
        --------
        read_inputfile
        """
        #self.no_of_atoms = len(atoms)
        pcs.System.assign_particles(self, atoms, box)

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
        distances = pcs.System.get_pairdistances(self)

        if histomax == None:
            histomax = max(distances)

        hist, bin_edges = np.histogram(distances, bins=histobins, range=(histomin, histomax))
        edgewidth = np.abs(bin_edges[1]-bin_edges[0])
        hist = hist.astype(float)
        r = bin_edges[:-1]

        #get box density
        boxvecs = pcs.System.get_boxvecs(self)
        vol = np.dot(np.cross(boxvecs[0], boxvecs[1]), boxvecs[2])
        natoms = pcs.System.get_nop(self)
        rho = natoms/vol

        shell_vols = (4./3.)*np.pi*((r+edgewidth)**3 - r**3)
        shell_rho = hist/shell_vols
        #now divide to get final value
        rdf = shell_rho/rho

        return rdf, r






    def get_atoms(self):
        """

        Get a list of all :class:`~Atom` objects that belong to the system.

        Parameters
        ----------
        None

        Returns
        -------
        allatoms : list of `Atom` objects
            all atoms in the system
        """
        atomcs = pcs.System.get_allatoms(self)
        atoms = [self.copy_catom_to_atom(xx) for xx in atomcs]
        return atoms

    def get_qvals(self, q, averaged = False):
        """
        Get the required q values of all atoms.

        Parameters
        ----------
        q : int or list of ints
            required q value from 2-12

        averaged : bool, optional
            If True, return the averaged q values,
            If False, return the non averaged ones
            default False

        Returns
        -------
        qvals : list of floats
            list of qvalue of all atoms.

        Notes
        -----
        The function returns a list of
        q values in the same order as that of the atoms in the system.

        """
        if isinstance(q, int):
            if q in range(2, 13):
                if averaged:
                    rq = pcs.System.get_aqvals(self, q)
                else:
                    rq = pcs.System.get_qvals(self, q)
                return rq
            else:
                raise ValueError("the value of q should be between 2 and 12")

        else:
            for qq in q:
                if not qq in range(2, 13):
                    raise ValueError("the value of q should be between 2 and 12")
            if averaged:
                rq = [ pcs.System.get_aqvals(self, qq) for qq in q ]
            else:
                rq = [ pcs.System.get_qvals(self, qq) for qq in q ]
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

        """

        atom1c = self.copy_atom_to_catom(atom1)
        atom2c = self.copy_atom_to_catom(atom2)
        return pcs.System.get_absdistance(self, atom1c, atom2c)


    def find_neighbors(self, method="cutoff", cutoff=None, threshold=2, filter=None,
                                            voroexp=1, face_cutoff=0.002, padding=1.2, nlimit=6):
        """

        Find neighbors of all atoms in the :class:`~System`.

        Parameters
        ----------
        method : {'cutoff', 'voronoi'}
            `cutoff` method finds atoms within a specified or adaptive cutoff distance of the host atom
            `voronoi` method finds atoms that share a Voronoi polyhedra face with the host atom. Default, `cutoff`

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
            only used if ``method=voronoi``. Power of the neighbor weight used to weight the contribution of each atom towards the q
            values. Default 1.

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
            raised when `threshold` value is too less. A low threshold value will lead to 'sann' algorithm not converging
            in finding a neighbor. The function will try to automatically increase `threshold` and check again.

        RuntimeError
            raised when neighbor search was unsuccessful. This is due to a low `threshold` value.

        Notes
        -----
        This function calculates the neighbors of each particle. There are several ways to do this. A complete description of
        the methods can be `found here <https://pyscal.readthedocs.io/en/latest/nearestneighbormethods.html>`_.

        Method cutoff and specifying a cutoff radius uses the traditional approach being the one in which the neighbors of an atom
        are the ones that lie in the cutoff distance around it.

        In order to reduce time during sorting of distances during the adaptive methods, pyscal guesses an initial large cutoff radius.
        This is calculated as,

        .. math:: r_{initial} = threshold * (simulation~box~volume/ number~of~particles)^{(1/3)}

        threshold is a safe multiplier used for the guess value and can be set using the `threshold` keyword.

        In Method cutoff, if ``cutoff='adaptive'``, an adaptive cutoff is decided during runtime for each atom to find its neighbors [1]_.
        Setting the cutoff radius to 0 also triggers this algorithm. The cutoff for an atom i is found using,

        .. math:: r_c(i) = padding * ((1/nlimit) * \sum_{j=1}^{nlimit}(r_{ij}))

        padding is a safe multiplier to the cutoff distance that can be set through the keyword `padding`. `nlimit` keyword sets the
        limit for the top nlimit atoms to be taken into account to calculate the cutoff radius.

        In Method cutoff, if ``cutoff='sann'``, sann algorithm is used [2]_. There are no parameters to customise sann behaviour.

        The second approach is using Voronoi polyhedra. All the atoms that share a Voronoi polyhedra face with the host atoms are considered
        its neighbors. A corresponding neighborweight is also assigned to each neighbor in the ratio of the face area between the two atoms.
        This weight can later be used to weight steinhardt parameters. Higher powers of this weight can also be used [3]_. The keyword `voroexp`
        can be used to set this weight. If `voroexp` is set to 0, the neighbors would be calculated using Voronoi method, but Steinhardts
        parameters could be calculated normally.

        Keyword `filter` can be used to filter the neighbors based on a condition. Choosing ``filter='type'`` only considers an atom as
        a neighbor if both the neighbor atom and host atom are of the same type.

        .. warning::

            Adaptive cutoff uses a padding over the intial guessed "neighbor distance". By default it is 2. In case
            of a warning that ``threshold`` is inadequate, it should be further increased. High/low value
            of this parameter will correspond to the time taken for finding neighbors.

        References
        ----------
        .. [1] Stukowski, A, Model Simul Mater SC 20, 2012
        .. [2] van Meel, JA, Filion, L, Valeriani, C, Frenkel, D, J Chem Phys 234107, 2012
        .. [3] Haeberle, J, Sperl, M, Born, P, arxiv 2019

        """
        #first reset all neighbors
        pcs.System.reset_allneighbors(self)
        pcs.System.set_filter(self, 0)

        if filter == 'type':
            # type corresponds to 1
            pcs.System.set_filter(self, 1)

        if method == 'cutoff':
            if cutoff=='sann':
                if threshold < 1:
                    raise ValueError("value of threshold should be at least 1.00")
                finished = pcs.System.get_all_neighbors_sann(self, threshold)
                #if it finished without finding neighbors
                if not finished:
                    finallydone = False
                    for i in range(1,10):
                        #threshold value is probably too low
                        #try increasing threshold
                        warnings.warn("Could not find sann cutoff. trying with a higher threshold", RuntimeWarning)
                        pcs.System.reset_allneighbors(self)
                        newfinished = pcs.System.get_all_neighbors_sann(self, threshold*i)
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
                finished = pcs.System.get_all_neighbors_adaptive(self, threshold, nlimit, padding)
                if not finished:
                    raise RuntimeError("Could not find adaptive cutoff")
            else:
                #warnings.warn("THIS RAN")
                pcs.System.set_neighbordistance(self, cutoff)
                pcs.System.get_all_neighbors_normal(self)

        elif method == 'voronoi':
            pcs.System.set_face_cutoff(self, face_cutoff)
            pcs.System.set_alpha(self, int(voroexp))
            pcs.System.get_all_neighbors_voronoi(self)

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
        pcs.System.reset_allneighbors(self)

    def calculate_vorovector(self, edge_cutoff=0.05, area_cutoff=0.01, edge_length=False):
        """
        get the voronoi structure identification vector.

        Parameters
        ----------
        edge_cutoff : float, optional
            cutoff for edge length. Default 0.05.


        area_cutoff : float, optional        self.initialized = 23
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
        atoms = pcs.System.get_allatoms(self)

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


    def calculate_q(self, q, averaged = False):
        """
        Find the bond order parameter q for all atoms.

        Parameters
        ----------
        qs : int or list of ints
            A list of all q params to be found from 2-12.

        averaged : bool, optional
            If True, return the averaged q values,
            If False, return the non averaged ones
            default False

        Returns
        -------
        None

        Notes
        -----
        Enables calculation of the Steinhardt parameters q from 2-12. The type of
        q values depend on the method used to calculate neighbors. See the description
        :func:`~System.find_neighbors` for more details. If the keyword `average` is set to True,
        the averaged versions of the bond order parameter [1]_ is returned.

        References
        ----------
        .. [1] Lechner, W, Dellago, C, J Chem Phys, 2013

        """
        if isinstance(q, int):
            qq = [q]
        else:
            qq = q

        for ql in qq:
            if not ql in range(2,13):
                raise ValueError("value of q should be between 2 and 13")

        pcs.System.calculate_q(self, qq)

        if averaged:
            pcs.System.calculate_aq(self, qq)

    def get_largestcluster(self):
        """
        Get id of the the largest cluster.

        Parameters
        ----------
        None

        Returns
        -------
        clusterid : int
            id of the largest cluster

        Notes
        -----
        id is only available if the largest cluster has already
        been found using the :func:`~System.find_solids` with keyword `cluster=True`. Otherwise it returns the default values.

        """
        return pcs.System.get_largestcluster(self)

    def find_solids(self, bonds=7, threshold=0.5, avgthreshold=0.6, cluster=True):
        """
        Distinguish solid and liquid atoms in the system.

        Parameters
        ----------
        bonds : int, optional
            Minimum number of solid bonds for an atom to be identified as
            a solid. Default 7.

        threshold : double, optional
            The cutoff value of connection between two atoms for them to be def
            ined as having a bond. Default 0.5.

        avgthreshold : double, optional
            Averaged value of connection between an atom and its neighbors for
            an atom to be solid. This threshold is known to improve the solid-liquid
            distinction in interfaces between solid and liquid. Default 0.6.

        cluster : bool, optional
            If True, cluster the solid atoms and return the number of atoms in the largest
            cluster.

        Returns
        -------
        lc : int
            Size of the largest cluster of solid atoms. Returned only if `cluster=True`.

        Notes
        -----
        The number of solid atoms in liquid using method found in [1]_ and example
        usage (not with this module) can be found in [2]_. The neighbors should be calculated
        before running this function. Check :func:`~System.find_neighbors` method.

        `bonds` define the number of solid bond required for an atom to be considered solid.
        Two particles are said to be 'bonded' if,

        .. math:: s_{ij} = \sum_{m=-6}^6 q_{6m}(i) q_{6m}^*(i) \geq threshold

        where `threshold` values is also an optional parameter.

        `avgthreshold` is an additional parameter to improve solid-liquid distinction.
        In addition to having a the specified number of `bonds`,

        .. math::  \langle s_{ij} \\rangle > avgthreshold

        also needs to be satisfied.


        References
        ----------
        .. [1] Auer, S, Frenkel, D. Adv Polym Sci 173, 2005
        .. [2] Diaz Leines, G, Drautz, R, Rogal, J, J Chem Phys 146, 2017

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
        pcs.System.calculate_q(self, [6])
        #self.calculate_q(6)
        #calculate solid neighs
        pcs.System.set_nucsize_parameters(self, bonds, threshold, avgthreshold)
        pcs.System.calculate_frenkelnumbers(self)
        #now find solids
        pcs.System.find_solid_atoms(self)

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

        When clustering happens, the code loops over each atom and its neighbors. If the
        host atom and the neighbor both have True value for `condition`, they are put
        in the same cluster. For example, if the atoms need to be clustered over if they
        are solid or not, corresponding condition would be,

        .. code:: python

            def condition(atom):
                #if both atom is solid
                return (atom1.solid)

        Check examples for more details.

        """
        testatom = self.get_atom(0)

        #test the condition
        try:
            out = condition(testatom)
            if out not in [True, False]:
                raise RuntimeError("The output of condition should be either True or False. Received %s"%str(out))
        except:
            raise RuntimeError("condition did not work")

        #now loop
        nop = pcs.System.get_nop(self)
        for i in range(nop):
            atom = self.get_atom(i)
            cval = condition(atom)
            atom.condition = cval
            self.set_atom(atom)

        #atom conditions are set
        #now can clustering function
        pcs.System.find_clusters_recursive(self)

        #done!
        lc = pcs.System.find_largest_cluster(self)
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
        warnings.simplefilter('always', DeprecationWarning)
        warnings.warn("This function is deprecated - use find_solids instead", DeprecationWarning)

        pcs.System.calculate_q(self, [6])
        pcs.System.set_nucsize_parameters(self, frenkelnums, threshold, avgthreshold)
        return pcs.System.calculate_nucsize(self)

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
        pcs.System.calculate_frenkelnumbers(self)

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
            pcs.System.find_clusters_recursive(self)
        else:
            pcs.System.find_clusters(self)

        if largest:
            cluster = pcs.System.find_largest_cluster(self)
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
        :func:`System.find_clusters` has to be used before using this function.

        """
        return pcs.System.find_largest_cluster(self)

    def copy_catom_to_atom(self, atomc, destination=None):
        """
        Used to copy a C++ `Atom` object to python `Atom` object.

        Parameters
        ----------
        atomc : C++ `Atom` object
            the input atom

        Returns
        -------
        atom : python `Atom` object
            output atom

        Notes
        -----
        This function is used to make sure that the user gets a python
        object rather than a C++ one.



        """
        if destination is None:
            atom = Atom()
        else:
            atom = destination

        atom.pos = atomc.get_x()
        atom.solid = bool(atomc.get_solid())
        atom.structure = atomc.get_structure()

        cinfo = atomc.get_cluster()
        atom.surface = bool(cinfo[1])
        atom.largest_cluster = bool(cinfo[2])
        atom.cluster = cinfo[3]

        atom.neighbors = atomc.get_neighbors()
        atom.neighbor_weights = atomc.get_neighborweights()
        atom.allq = atomc.get_allq()
        atom.allaq = atomc.get_allaq()
        atom.id = atomc.get_id()
        atom.loc = atomc.get_loc()
        atom.type = atomc.get_type()

        #atom.set_vorovector(atomc.get_vorovector())
        atom.volume = atomc.get_volume()
        atom.avg_volume = atomc.get_avgvolume()
        atom.face_vertices = atomc.get_facevertices()
        atom.face_perimeters = atomc.get_faceperimeters()
        atom.vertex_numbers = atomc.get_vertexnumbers()
        atom.vertex_vectors = atomc.get_vertexvectors()
        atom.condition = atomc.get_condition()
        atom.avg_connection = atomc.get_avgconnection()
        atom.bonds = atomc.get_bonds()
        #atom.set_allqcomps(atomc.get_allqcomps())
        return atom

    def copy_atom_to_catom(self, atom):
        """
        Used to copy a python `Atom` object to C++ `Atom` object.

        Parameters
        ----------
        atom : python `Atom` object
            output atom

        Returns
        -------
        atomc : C++ `Atom` object
            the input atom

        Notes
        -----
        A python object is converted to a C++ object which is then
        passed on the system for calculation.

        """
        atomc = pca.Atom()
        atomc.set_x(atom.pos)
        atomc.set_solid(atom.solid)
        atomc.set_structure(atom.structure)

        #prep cluster values
        cinfo = [int(atom.solid), int(atom.surface), int(atom.largest_cluster), atom.cluster]
        atomc.set_cluster(cinfo)
        atomc.set_neighbors(atom.neighbors)
        atomc.set_neighborweights(atom.neighbor_weights)
        atomc.set_allq(atom.allq)
        atomc.set_allaq(atom.allaq)
        atomc.set_id(atom.id)
        atomc.set_loc(atom.loc)
        atomc.set_type(atom.type)
        #atomc.set_vorovector(atom.get_vorovector())
        atomc.set_volume(atom.volume)
        atomc.set_avgvolume(atom.avg_volume)
        atomc.set_facevertices(atom.face_vertices)
        atomc.set_faceperimeters(atom.face_perimeters)
        atomc.set_vertexnumbers(atom.vertex_numbers)
        atomc.set_vertexvectors(atom.vertex_vectors)
        atomc.set_condition(atom.condition)
        atomc.set_avgconnection(atom.avg_connection)
        atomc.set_bonds(atom.bonds)
        #atomc.set_allqcomps(atom.get_allqcomps())
        return atomc


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
        a user perspective, the :func:`~System.to_file` method should be used
        directly.

        See also
        --------
        to_file()

        """

        #get the basic system indicators
        indicators = pcs.System.get_indicators(self)

        #get box dims and triclinic params if triclinic
        boxdims = self.box
        if indicators[6] == 1:
            rot = pcs.System.get_triclinic_params(self)
        else:
            rot = 0

        #now finally get atoms
        atoms = pcs.System.get_allatoms(self)
        #convert them to picklabale atoms
        patoms = [self.copy_catom_to_atom(atom) for atom in atoms]

        #create System instance and assign things
        psys = pp.System()
        psys.indicators = indicators
        psys.atoms = patoms
        psys.boxdims = boxdims
        psys.rot = rot

        return psys

    def to_file(self, file):
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
        This function can be used to save a :class:`~System` object directly to
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


    def from_file(self, file):
        """
        Populate the empty system from file

        Parameters
        ----------
        file : string
            name of output file

        Returns
        -------
        None

        Notes
        -----
        This function can be used to set up a system
        from a file. A system needs to be created first.

        Examples
        --------

        >>> sys = System()
        >>> sys.from_file(filename)

        """
        if os.path.exists(file):
            psys = np.load(file, allow_pickle=True).flatten()[0]
        else:
            raise FileNotFoundError("file does not exist")
        #set up indicators
        pcs.System.set_indicators(self, psys.indicators)
        #unpickle atoms
        catoms = [self.copy_atom_to_catom(atom) for atom in psys.atoms]
        boxdims = psys.boxdims


        #if triclinic, get those
        if psys.indicators[6] == 1:
            rot = psys.rot
            rotinv = np.linalg.inv(rot)
            pcs.System.assign_triclinic_params(self, rot, rotinv)

        #assign atoms and box
        pcs.System.reassign_particles(self, catoms, boxdims)

import pybop.ccore as pc
import pybop.traj_process as ptp
import os
import numpy as np

"""
Definitions of class Atom.

"""
class Atom(pc.Atom):
    """
    A c++ class for holding the properties of a single atom. Various properties of the atom
    can be accessed through member functions which are described below in detail. Atoms can
    be created individually or directly by reading in a file. Check the examples for more 
    details on how atoms are created. For creating atoms directly from an input file check
    the documentation of `System` class.

    Although an `Atom` object can be created independently, Atoms should be thought of 
    inherently as members of the `System` class. All the properties that define an `atom` are
    relative to the `System` class. `System` has a list of all atoms using which the neighbors
    of an atom, if its solid and so on can be calculated. All the calculated properties of an
    atom which depend on any other atom hence should be calculated through `System`. Please
    check the examples section of the documentation for more details. 

    Parameters
    ----------
    pos : list of floats of length 3, default [0,0,0]
        position of the `Atom`

    id : int, default 0
        id of the `Atom`

    Examples
    --------
    >>> #method 1 - individually
    >>> atom = Atom()
    >>> #now set positions of the atoms
    >>> atom.set_x([23.0, 45.2, 34.2])
    >>> #now set id
    >>> atom.set_id(23)

    See also
    --------
    get_x
    set_x
    set_id
    get_id
    System
    
    """
    def __init__(self, pos=[0,0,0], id=0, type=1):
        """
        Deafults args
        """
        pc.Atom.__init__(self)     
        pc.Atom.set_x(self, pos)
        pc.Atom.set_id(self, id)
        pc.Atom.set_type(self, type)

    #now wrapping for other normal functions
    def get_x(self):
        """
        
        Get the position of the atom. Meaningful values are only returned if the atoms are
        set before using this function.

        Parameters
        ----------
        None
        
        Returns
        -------
        x : array of float
            contains the position of the atom in the form [posx, posy, posz], where
            posx is the x coordinate of the atom, posy is the y coordinate and posz 
            is the z coordinate. 

        Examples
        --------
        >>> atom = Atom()
        >>> x = atom.get_x()

        See also
        --------
        set_x
        Atom
        System      
        """
        x = pc.Atom.get_x(self)
        return x

    def set_x(self, pos):
        """
        
        Set the position of the atom. 

        Parameters
        ----------
        x : list of floats of length 3
            list contains three values which are the position coordinates of the atom with
            respect to the simulation box.

        Returns
        -------
        None

        Examples
        --------
        >>> atom = Atom()
        >>> x = atom.set_x([23.0, 45.2, 34.2])

        See also
        --------
        get_x

        """
        if len(pos) == 3:
            pc.Atom.set_x(self, pos)

    def get_solid(self):
        """
        Find if an atom is solid or not.

        Parameters
        ----------
        None

        Returns
        -------
        issolid : 1 if solid, 0 otherwise
        """
        return pc.Atom.get_solid(self)

    def get_structure(self):
        """
        Get the structure of an atom.

        Parameters
        ----------
        None

        Returns
        -------
        structure : structural value
        """
        return pc.Atom.get_structure(self)

    def set_solid(self, issolid):
        """
        Find if an atom is solid or not.

        Parameters
        ----------
        issolid : int, 0 or 1 
            1 if the atom is set to solid, 0 otherwise

        Returns
        -------
        None
        """
        if int(issolid) in [0, 1]: 
            pc.Atom.set_solid(self, int(issolid))
        else:
            raise ValueError("Value of issolid should be either 0 or 1")

    def set_structure(self, structure):
        """
        Set the structure of an atom.

        Parameters
        ----------
        structure : int
            structure of the atom

        Returns
        -------
        None
        """
        return pc.Atom.set_structure(self, structure)


    def get_volume(self, averaged = False):
        """
        
        Get the voronoi colume of the atom. Meaningful values are only returned if the neighbors
        are calculated using voronoi method.

        Parameters
        ----------
        averaged : bool, default False
            If True, averaged version of the volume is returned.
        
        Returns
        -------
        volume : float
            voronoi volume of the atom.

        Examples
        --------
        >>> volume = atom.get_volume()

        See also
        --------
        find_neighbors
        Atom
        System      
        """
        if averaged:
            vol = pc.Atom.get_avgvolume(self)
        else:    
            vol = pc.Atom.get_volume(self)
        
        return vol

    def get_cluster(self):
        """
        Get the cluster properties of the atom. The cluster properties of the atom
        include four different properties as listed below. The properties are only
        returned if they are calculated before using 'calculate_nucsize' function 
        before.

        Parameters
        ----------
        None
        
        Returns
        -------
        cluster : list of int of length 4
            cluster is a vector of four values. they are described below-
                issolid - which is 1 if the atom is solid, 0 otherwise
                issurface - 1 if the atom has liquid neighbors, 0 otherwise
                lcluster - 1 if the atom belongs to the largest cluster, 0 otherwise
                belongsto - which gives the id of the cluster that the atom belongs to.
        
        Examples
        --------
        >>> cinfo = atom.get_cluster()

        See also
        --------
        set_nucsize_parameters
        calculate_nucsize

        """
        x = pc.Atom.get_cluster(self)
        return x

    def get_neighbors(self):
        """
        
        Returns the neighbors indices of the atom. The list returned consistes of the indices
        of neighbor atom which indicate their position in the list of all atoms. The neighbors
        of an atom can be calculated from the `System` object that it belongs to.

        Parameters
        ----------
        None
        
        Returns
        -------
        x : list of int
            list of neighbor indices of the atom.

        Examples
        --------
        neighbors = atom.get_neighbors()

        See also
        --------
        set_neighbors
        set_neighborweights
        get_neighborweights

        """
        return pc.Atom.get_neighbors(self)

    def set_neighbors(self, neighs):
        """
        
        Set the neighbors of an atom manually.

        Parameters
        ----------
        neighs : list of ints
            index of the neighbor atoms 
        
        Returns
        -------
        None

        Examples
        --------
        atom.set_neighbors([0,23,11,22,334,112,11])

        See also
        --------
        get_neighbors
        set_neighborweights
        get_neighborweights
        """
        pc.Atom.set_neighbors(self, neighs)


    def get_coordination(self):
        """
        
        Returns the coordination number of the atom. `get_allneighbors` function of the `System` class
        has to be used before accessing coordination numbers. 

        Parameters
        ----------
        None
        
        Returns
        -------
        cn : int
            coordination number of the atom.

        Examples
        --------
        neighbors = atom.get_neighbors()

        See also
        --------
        set_neighbors
        set_neighborweights
        get_neighborweights

        """
        return pc.Atom.get_coordination(self)

    def get_neighborweights(self):
        """
        
        Get the neighbor weights of the atom. The neighbor weights are used weight the 
        contribution of each neighboring atom towards the q value of the host atom. By 
        default, each neighbor has a weight of 1 each. However, if the neighbors are calculated
        using the `System.get_allneighbors(method='voronoi')`, each neighbor atom gets a 
        weight proportional to the face area shared between the neighboring atom and the 
        host atom. This can sometimes be helpful in controlling the contribution of atoms
        with low face areas due to the thermal vibrations at high temperature.

        Parameters
        ----------
        None
        
        Returns
        -------
        x : list of float
            neighbor weights

        Examples
        --------
        >>> weights = atom.get_neighborweights()

        See also
        --------
        get_neighbors
        set_neighbors
        set_neighborweights
        """

        return pc.Atom.get_neighborweights(self)

    def set_neighborweights(self, weights):
        """
        
        Set the neighbor weights of an atom.

        Parameters
        ----------
        weights : list of floats 
            weights of the neighbor atoms 
        
        Returns
        -------
        None

        Examples
        --------
        >>> atom.set_neighborweights([0.1, 0.2, 0.2, 0.4, 0.1])

        See also
        --------
        get_neighbors
        set_neighbors
        set_neighborweights
        """
        pc.Atom.set_neighborweights(self, weights)

    def get_q(self, q, averaged = False):
        """
        get q value of the atom. The q value can be either normal or can
        be averaged according to Lechner, W, Dellago, C. JCP 129, 2008.
        The averaged version can be obtained by using keyword
        `averaged = True`.

        Parameters
        ----------
        q : int or list of int
            number of the required q - from 2-12

        averaged : bool, default False
            If True, return the averaged q values,
            If False, return the non averaged ones

        Returns
        -------
        q : float or list of floats
            The queried q value

        Examples
        --------
        >>> q2 = atom.get_q(2, averaged = True)
        >>> q24 = atom.get_q([2, 4])

        See also
        --------
        set_q
        get_aq
        set_aq      
        """
        if isinstance(q, int):
            if averaged:
                rq = pc.Atom.get_aq(self, q)
            else:
                rq = pc.Atom.get_q(self, q)
            return rq

        else:
            if averaged:
                rq = [ pc.Atom.get_aq(self, qq) for qq in q ]
            else:
                rq = [ pc.Atom.get_q(self, qq) for qq in q ]
            return rq

    def set_q(self, q, d, averaged = False):
        """
        set the q value of the atom. If `averaged = True`, sets the averaged versions of
        the q parameters.

        Parameters
        ----------
        q : int or list of ints
            number of the required q - from 2-12
        d : float or list of floats
            the q value to set
        averaged : bool, default False
            If True, return the averaged q values,
            If False, return the non averaged ones            

        Returns
        -------
        None

        Examples
        --------
        >>> atom.set_q(2, 0.24, averaged = True)
        >>> atom.set_q([2,4], [0.24, 0.05])

        See also
        --------
        set_aq
        get_aq
        get_q
        """
        if isinstance(q, int):
            if averaged:
                pc.Atom.set_aq(self, q, d)
            else:
                pc.Atom.set_q(self, q, d)
        else:
            if averaged:
                for count, qq in enumerate(q):
                    pc.Atom.set_aq(self, qq, d[count])
            else:
                for count, qq in enumerate(q):
                    pc.Atom.set_q(self, qq, d[count])


    def get_id(self):
        """
        
        get  the id of the atom.

        Parameters
        ----------
        None

        Returns
        -------
        id : int
            id of the atom

        Examples
        --------
        >>> id = atom.get_id()

        See also
        --------
        set_id
        """
        return pc.Atom.get_id(self)

    def set_id(self, idd):
        """
        
        set  the id of the atom.

        Parameters
        ----------
        idd : int
            id of the atom

        Returns
        -------
        None

        Examples
        --------
        >>> atom.set_id(2)

        See also
        --------
        get_id

        """
        pc.Atom.set_id(self, idd)

    def get_loc(self):
        """
        
        get  the location of the atom.

        Parameters
        ----------
        None

        Returns
        -------
        loc : int
            loc of the atom

        Examples
        --------
        >>> loc = atom.get_loc()

        See also
        --------
        set_loc
        """
        return pc.Atom.get_loc(self)

    def set_loc(self, idd):
        """
        
        set  the loc of the atom.

        Parameters
        ----------
        idd : int
            loc of the atom

        Returns
        -------
        None

        Examples
        --------
        >>> atom.set_loc(2)

        See also
        --------
        get_loc

        """
        pc.Atom.set_loc(self, idd)

    def get_type(self):
        """
        
        get  the type(species) of the atom.

        Parameters
        ----------
        None

        Returns
        -------
        type : int
            species/type of the atom

        Examples
        --------
        >>> t = atom.get_type()

        See also
        --------
        set_type
        """
        return pc.Atom.get_type(self)

    def set_type(self, tt):
        """
        
        set  the type of the atom.

        Parameters
        ----------
        tt : int
            type of the atom

        Returns
        -------
        None

        Examples
        --------
        >>> atom.set_type(2)

        See also
        --------
        get_type

        """
        pc.Atom.set_type(self, tt)


    def get_vorovector(self):
        """
        get the voronoi structure identification vector. Returns a
        vector of the form (n3, n4, n5, n6), where n3 is the number
        of faces with 3 vertices, n4 is the number of faces with 4
        vertices and so on. This can be used to identify structures.

        Parameters
        ----------
        None

        Returns
        -------
        vorovector : array like, int
            array of the form (n3, n4, n5, n6)
        """
        return pc.Atom.get_vorovector(self)

"""
System class definitions
"""
class System(pc.System):
    """
    A c++ class for holding the properties of a system. A `System` consists of two
    major components - the simulation box and the atoms. All the associated variables
    are then calculated over these.

    A `System` can be set and populated by reading an input file in lammps dump format.
    This enables for automatic reading of all atomic positions and the simulation box.

    Examples
    --------
    >>> sys = System()
    >>> sys.read_inputfile()
    """
    def __init__(self):
        self.initialized = True
        pc.System.__init__(self)

    def read_inputfile(self, filename, format="lammps-dump", frame=-1, compressed = False):
        """
        
        Read input file containing the information of a time slice from a molecular dynamics
        simulation. 
        
        The format of the input file is specified using the `format` keyword. Currently only
        a `lammps-dump` file is supported. However, this restriction can easily be overcome 
        using the `assign_particles` method from system where a list of atoms and box vectors 
        are directly provided to the system.

        `use_c` is no deprecated and no longer used to ease transition to python 3. 
        
        Parameters
        ----------
        filename : string
            name of the input file to be read in

        format : string, `lammps-dump` or `poscar`
            format of the input file

        compressed : bool, default False
            If True, force to read a `gz` compressed format. However, if a file ends with `.gz`
            extension, it is automatically treated as a compressed file and this keyword is not
            necessary

        frame : int
            If the trajectory contains more than one time slice, the slice can be specified
            using the `frame` option. Alert: works only with `lammps-dump` format. 

        use_c : bool, default False, deprecated
            If True, use the `read_particle_file` method from c++ module. This might be faster
            but only accepts file format of `lammps-dump` type with a particular header layout.
            Also `compressed` keyword doesnt work anymore. This keyword is deprecated and only
            kept for compatibility reasons. Use of this keyword is not recommended.


        Returns
        -------
        None

        See Also
        --------
        assign_particles
        """
        if format == 'lammps-dump':
            if frame != -1:
                #split the traj and returns set of filenames
                filenames = ptp.split_traj_lammps_dump(filename, compressed=compressed)
                #reassign filename
                filename = filenames[frame]
                if os.path.exists(filename):
                    atoms, boxdims, box, triclinic = ptp.read_lammps_dump(filename, compressed=compressed, check_triclinic=True, box_vectors=True)
                    pc.System.assign_particles(self, atoms, boxdims)
                    if triclinic:
                        #we have to input rotation matrix and the inverse rotation matrix
                        rot = box.T
                        rotinv = np.linalg.inv(rot)
                        pc.System.assign_triclinic_params(self, rot, rotinv)
                else:
                    raise FileNotFoundError("input file not found")
                #now remove filenames
                for file in filenames:
                    os.remove(file)

            elif os.path.exists(filename):
                atoms, boxdims, box, triclinic = ptp.read_lammps_dump(filename, compressed=compressed, check_triclinic=True, box_vectors=True)
                pc.System.assign_particles(self, atoms, boxdims)
                if triclinic:
                    #we have to input rotation matrix and the inverse rotation matrix
                    rot = box.T
                    rotinv = np.linalg.inv(rot)
                    pc.System.assign_triclinic_params(self, rot, rotinv)
            else:
                raise FileNotFoundError("input file not found")

        elif format == 'poscar':
            if os.path.exists(filename):
                atoms, boxdims = ptp.read_poscar(filename, compressed=compressed)
                pc.System.assign_particles(self, atoms, boxdims)
            else:
                raise FileNotFoundError("input file not found")


    def assign_atoms(self, atoms, box):
        """
        
        Assign atoms directly. Receive a vector of atom objects which is stored instead
        of reading in the input file. If this method is used, there is no need of using
        `read_inputfile` method. Also using this function allows for reading of multiple
        file formats which are not supported by the inbuilt `read_inputfile` method.

        Parameters
        ----------
        atoms : list of `Atom` objects
            list consisting of all atoms
        box   : list of list of floats
            list which consists of the box dimensions in the format-
            [[box_x_low, box_x_high], [box_y_low, box_y_high], [box_z_low, box_z_high]]

        Returns
        -------
        None

        See Also
        --------
        read_inputfile
        """
        pc.System.assign_particles(self, atoms, box)

    def calculate_rdf(self, histobins=100, histomin=0.0, histomax=None):
        """
        Calculate the radial distribution function. It is calculated by finding distances
        between all pairs of atoms and then creating a histogram from it.

        Parameters
        ----------
        histobins : int
            number of bins in the histogram
        histomin : float, optional
            minimum value of the distance histogram, if not specified, 0.0 is taken as the
            minimum.
        histomax : float, optional
            maximum value of the distance histogram. If not specified, the maximum value
            in all pair distances is used.

        Returns
        -------
        rdf : array of ints
            Radial distribution function
        r : array of floats
            radius in distance units

        """
        distances = pc.System.get_pairdistances(self)

        if histomax == None:
            histomax = max(distances)

        hist, bin_edges = np.histogram(distances, bins=histobins, range=(histomin, histomax))
        edgewidth = np.abs(bin_edges[1]-bin_edges[0])
        hist = hist.astype(float)
        r = bin_edges[:-1]

        #get box density
        boxvecs = pc.System.get_boxvecs(self)
        vol = np.dot(np.cross(boxvecs[0], boxvecs[1]), boxvecs[2])
        natoms = pc.System.get_nop(self)
        rho = natoms/vol

        shell_vols = (4./3.)*np.pi*((r+edgewidth)**3 - r**3)
        shell_rho = hist/shell_vols
        #now divide to get final value
        rdf = shell_rho/rho

        return rdf, r



    def get_largestcluster(self):
        """
        Get id of the the largest cluster. id is only available if the largest cluster has already
        been found. Otherwise it returns the default values.

        Parameters
        ----------
        None

        Returns
        -------
        clusterid : int 
            id of the largest cluster

        See Also
        --------
        get_allneighbors
        calculate_nucsize
        set_nucsize_parameters
        """
        return pc.System.get_largestcluster(self)

    def set_nucsize_parameters(self, cutoff, minfrenkel, threshold, avgthreshold):
        """
        Set the value of parameters for calculating the largest solid cluster in the
        liquid, a detailed description of the order parameter can be found in  
        Diaz Leines et al, JCP 146(2017). http://doi.org/10.1063/1.4980082.

        The number of atoms in the largest solid cluster in liquid is often used as an
        order parameter in the study of nucleation during solidification. In order to
        actually calculate the largest solid cluster, `calculate_nucsize` has to be 
        called after setting the parameters.

        Parameters
        ----------
        cutoff : float
            cutoff distance for calculating neighbors

        minfrenkel : int
            Minimum number of solid bonds for an atom to be identified as
            a solid.
        
        threshold : double
            The cutoff value of connection between two atoms for them to be def
            ined as having a bond.
        
        avgthreshold : double
            Averaged value of connection between an atom and its neighbors for 
            an atom to be solid. This threshold is known to improve the solid-liquid
            distinction in interfaces between solid and liquid. 

        Returns
        -------
        None

        See Also
        --------
        calculate_nucsize

        Examples
        --------
        >>> st.set_nucsize_parameters(7,0.5,0.5)
        """
        pc.System.set_nucsize_parameters(self, cutoff, minfrenkel, threshold, avgthreshold)

    def calculate_nucsize(self):
        """
        Calculate the size of the largest cluster in the given system. Calculation
        the size of the largest cluster needs various prerequisites that can be set
        by the functions `set_nucsize_parameters`. A detailed description of the order 
        parameter can be found in Diaz Leines et al, JCP 146(2017). 
        http://doi.org/10.1063/1.4980082.

        The number of atoms in the largest solid cluster in liquid is often used as an
        order parameter in the study of nucleation during solidification.

        Parameters
        ----------
        None

        Returns
        -------
        cluster size : int
            size of the largest solid cluster in liquid (number of atoms)
        """
        return pc.System.calculate_nucsize(self)

    def get_atom(self, index):
        """
        
        Get the `Atom` object at the queried position in the list of all atoms
        in the `System`.

        Parameters
        ----------
        index : int
            index of required atom in the list of all atoms.

        Returns
        -------
        atom : Atom object
            atom object at the queried position.
        """
        atomc = pc.System.get_atom(self, index)
        atom = self.copy_catom_to_atom(atomc)
        return atom

    def set_atom(self, atom):
        """
        
        Return the atom to its original location after modification. For example, an
        `Atom` at location `i` in the list of all atoms in `System` can be queried by,
        `atom = System.get_atom(i)`, then any kind of modification, for example, the 
        position of the `Atom` can done by, `atom.set_x([2.3, 4.5, 4.5])`. After 
        modification, the `Atom` can be set back to its position in `System` by
        `System.set_atom(atom)`.

        Parameters
        ----------
        atom : Atom
                atom to be replaced

        Returns
        -------
        None
        """
        atomc = self.copy_atom_to_catom(atom)
        pc.System.set_atom(self, atomc)

    def get_allatoms(self):
        """
        
        Get a list of all `Atom` objects that belong to the system.

        Parameters
        ----------
        None

        Returns
        -------
        allatoms : list of `Atom` objects
            all atoms in the system
        """
        atomcs = pc.System.get_allatoms(self)
        atoms = [self.copy_catom_to_atom(xx) for xx in atomcs]
        return atoms

    def get_box(self, box_vectors=False):
        """
        
        Get the dimensions of the simulation box.

        Parameters
        ----------
        box_vectors : bool, default False
            If True, return the whole box dimesions.

        Returns
        -------
        boxdims : list of box dimensions
            If `box_vectors` is false:
            the return value consists of the vector of values in the form-
            [[box_x_low, box_x_high], [box_y_low, box_y_high], [box_z_low, box_z_high]]
            If `box_vectors` is true:
            return the box vectors of the form
            [[x1, x2, x3], [y1, y2, y3], [z1, z2, z3]]
        """
        if not box_vectors:
            box6dim = pc.System.get_box(self)
            pbox = [[box6dim[0], box6dim[1]], [box6dim[2], box6dim[3]], [box6dim[4], box6dim[5]]]
        else:
            pbox = pc.System.get_boxvecs(self)
        return pbox

    def set_box(self, box):
        """
        
        Set the dimensions of the simulation box.

        Parameters
        ----------
        boxdims : list of box dimensions of length 6
            the return value consists of the vector of values in the form-
            [[box_x_low, box_x_high], [box_y_low, box_y_high], [box_z_low, box_z_high]]

        Returns
        -------
        None
        """
        pc.System.set_box(self, box)

    def get_qvals(self, q, averaged = False):
        """
        Get the required q values of all atoms. The function returns a list of 
        q values in the same order as that of the atoms in the system.

        Parameters
        ----------
        q : int or list of ints
            required q value from 2-12
        averaged : bool, default False
            If True, return the averaged q values,
            If False, return the non averaged ones 

        Returns
        -------
        qvals : list of floats
            list of qvalue of all atoms.
        """
        if isinstance(q, int):
            if averaged:
                rq = pc.System.get_aqvals(self, q)
            else:
                rq = pc.System.get_qvals(self, q)
            return rq

        else:
            if averaged:
                rq = [ pc.System.get_aqvals(self, qq) for qq in q ]
            else:
                rq = [ pc.System.get_qvals(self, qq) for qq in q ]
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
        return pc.System.get_absdistance(self, atom1c, atom2c)

    def get_neighbors(self, method="cutoff", cutoff=None, nlimit=6, threshold=1.21, filter=None):
        """
        
        Find neighbors of all atoms in the `System`. There are few methods to do this, the 
        traditional approach being the one in which the neighbors of an atom are the ones that lie
        in a cutoff distance around it. The second approach is using Voronoi polyhedra. All the atoms
        that share a Voronoi polyhedra face with the host atoms are considered its neighbors.

        Finally there is also an adaptive cutoff method if `cutoff` is specified as `'adaptive'` or 0. 
        In this method, the six closest atoms to the
        host atom are found and a cutoff distance is determined as threshold*(1/nlimit)*sum(r_ij).

        Parameters
        ----------
        method : `cutoff` or `voronoi` or `adaptive-cutoff`, default: `cutoff`
            `cutoff` method finds atoms within a specified cutoff distance of the host atom
            `voronoi` method finds atoms that share a Voronoi polyhedra face with the host atom.

        cutoff : float or `adaptive`
            the cutoff distance to be used for the `cutoff` based neighbor calculation method
            described above.
            If the value is specified as 0 or `'adaptive'`, the `adaptive-cutoff` method is used.
        
        nlimit : int
            only used if `cutoff='adaptive'`. The number of atoms used to calculate the adaptive
            cutoff.

        threshold : float
            only used if `cutoff='adaptive'`. The threshold for the adaptive cutoff. If the average
            distance between host atom and `nlimit` neighbors of it is `a`, and `threshold` is 2, then
            host atom will have neighbors between `a` and `2a`.

        filter : string - `None` or `type`, default None
            apply a filter to nearest neighbor calculation. If the `filter` keyword is set to
            `type`, only atoms of the same type would be included in the neighbor calculations. 
        
        Returns
        -------
        None

        See also
        --------
        reset_allneighbors
        """
        #first reset all neighbors
        pc.System.reset_allneighbors(self)
        pc.System.set_filter(self, 0)

        if filter == 'type':
            # type corresponds to 1
            pc.System.set_filter(self, 1)

        if method == 'cutoff':
            if cutoff=='adaptive' or cutoff==0:
                pc.System.get_all_neighbors_adaptive(self, nlimit, threshold)
            else:    
                pc.System.set_neighbordistance(self, cutoff)
                pc.System.get_all_neighbors_normal(self)

        elif method == 'voronoi':
            pc.System.get_all_neighbors_voronoi(self)
            

    def reset_neighbors(self):
        """
        
        Reset the neighbors of all atoms in the system. This should be used before recalculating neighbors
        with two different approaches.

        Parameters
        ----------
        None
        
        Returns
        -------
        None

        See also
        --------
        get_allneighbors
        """
        pc.System.reset_allneighbors(self)

    def calculate_q(self, q, averaged = False):
        """
        Find the bond order parameter q for all atoms. Any of the q parameters from 2-12 can be provided.
        If the averaged versions are to be calculated, the keyword `averaged = True` should be set.
        See Lechner, Dellago, JCP 129, 2008. for a description of the
        averaged bond order parameters.

        Parameters
        ----------
        qs : int or list of ints
            A list of all q params to be found from 2-12. 
        averaged : bool, default False
            If True, return the averaged q values,
            If False, return the non averaged ones 

        Returns
        -------
        None
        """
        if isinstance(q, int):
            qq = [q]
        else:
            qq = q
            
        pc.System.calculate_q(self, qq)

        if averaged:
            pc.System.calculate_aq(self, qq)

    def calculate_frenkelnumbers(self):
        """
        Find frenkel numbers of all atoms in the system. Frenkel number is the number of solid
        neighbors that an atom has. A solid bond is considered between two atoms if the connection
        betweem them is greater than 0.6.

        Parameters
        ----------
        None
        
        Returns
        -------
        None
        
        See also
        --------
        get_number_from_bond
        find_clusters
        find_largest_cluster
        set_nucsize_parameters
        """
        pc.System.calculate_frenkelnumbers(self)

    def find_clusters(self, recursive = True, largest = True):
        """
        Find the clusters of all atoms in the system. Go through all the atoms and cluster them
        together based on the `issolid` parameter of the atom. To cluster based on any user defined criteria,
        you can use `set_solid` method of `Atom` to explicitely set the `issolid` value.

        Parameters
        ----------
        recursive : Bool, default True
            If True, use a recursive clustering algorithm, otherwise use an id based clustering.
            The difference in values between two methods can be upto 3 particles.
        largest : Bool, default True
            If True, return the number of particles in the largest cluster.
        
        Returns
        -------
        cluster : int
            The size of the largest cluster in the system. Only returned if `largest` is set to True.
        

        See also
        --------
        calculate_frenkelnumbers
        get_number_from_bond
        find_largest_cluster
        set_nucsize_parameters
        """
        if recursive:
            pc.System.find_clusters_recursive(self)
        else:    
            pc.System.find_clusters(self)

        if largest:
            cluster = pc.System.find_largest_cluster(self)
            return cluster            

    def find_largest_cluster(self):
        """
        Find the largest solid cluster of atoms in the system from all the clusters. `find_clusters`
        has to be used before using this function.

        Parameters
        ----------
        None
        
        Returns
        -------
        cluster : int
            the size of the largest cluster

        See also
        --------
        calculate_frenkelnumbers
        find_clusters
        get_number_from_bond
        set_nucsize_parameters
        """
        return pc.System.find_largest_cluster(self)

    def copy_catom_to_atom(self, atomc):
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
        """
        atom = Atom()
        atom.set_x(atomc.get_x())
        atom.set_solid(atomc.get_solid())
        atom.set_structure(atomc.get_structure())
        atom.set_cluster(atomc.get_cluster())
        atom.set_neighbors(atomc.get_neighbors())
        atom.set_neighborweights(atomc.get_neighborweights())
        atom.set_allq(atomc.get_allq())
        atom.set_allaq(atomc.get_allaq())
        atom.set_id(atomc.get_id())
        atom.set_loc(atomc.get_loc())
        atom.set_type(atomc.get_type())
        atom.set_vorovector(atomc.get_vorovector())
        atom.set_volume(atomc.get_volume())
        atom.set_avgvolume(atomc.get_avgvolume())
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

        """
        atomc = pc.Atom()
        atomc.set_x(atom.get_x())
        atomc.set_solid(atom.get_solid())
        atomc.set_structure(atom.get_structure())
        atomc.set_cluster(atom.get_cluster())
        atomc.set_neighbors(atom.get_neighbors())
        atomc.set_neighborweights(atom.get_neighborweights())
        atomc.set_allq(atom.get_allq())
        atomc.set_allaq(atom.get_allaq())
        atomc.set_id(atom.get_id())
        atomc.set_loc(atom.get_loc())
        atomc.set_type(atom.get_type())
        atomc.set_vorovector(atom.get_vorovector())
        atomc.set_volume(atom.get_volume())
        atomc.set_avgvolume(atom.get_avgvolume())
        return atomc
































    








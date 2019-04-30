import pybop.ccore as pc
import os
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
    def __init__(self, pos=[0,0,0], id=0):
        """
        Deafults args
        """     
        pc.Atom.set_x(self, pos)
        pc.Atom.set_id(self, id)

    #now wrapping for other normal functions
    def get_x(self):
        """
        TD
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
        TD
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
        TD
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
        TD
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
        TD
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
        TD
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
        TD
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
                pc.Atom.set_aq(self, q)
            else:
                pc.Atom.set_q(self, q)
        else:
            if averaged:
                for count, qq in enumerate(q):
                    pc.Atom.set_aq(self, qq, d[count])
            else:
                for count, qq in enumerate(q):
                    pc.Atom.set_q(self, qq, d[count])


    def get_id(self):
        """
        TD
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
        TD
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

    def get_qlm(self, q, averaged = False):
        """
        Get the real and imaginary qlm values of the atom. If `averaged = True`,
        return the averaged qlm values.

        Parameters
        ----------
        q : int or list of ints
            number of the required q - from 2-12

        Returns
        -------
        qlms : 2D array of 2q+1 values or n 2D arrays
            the first part of the array is the 2q+1 real values
            second part is the 2q+1 imaginary values.
        averaged : bool, default False
            If True, return the averaged q values,
            If False, return the non averaged ones    
        """
        if isinstance(q, int):
            if averaged:
                rq = pc.Atom.get_aqlm(self, q)
            else:
                rq = pc.Atom.get_qlm(self, q)
            return rq

        else:
            if averaged:
                rq = [ pc.Atom.get_aqlm(self, qq) for qq in q ]
            else:
                rq = [ pc.Atom.get_qlm(self, qq) for qq in q ]
            return rq

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
        return pc.Atom.get_vorovector()

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

    def read_inputfile(self, filename, format="lammps-dump"):
        """
        TD
        Read input file containing the information of a time slice from a molecular dynamics
        simulation. 
        

        As of now, the file should be a lammps dump format and can only have a
        specific header format. That is-
        id type mass x y z vx vy vz
        However, this restriction can easily be overcome using the `assign_particles` method
        from system where a list of atoms and box vectors are directly provided to the system.
        
        Parameters
        ----------
        filename : string
            name of the input file to be read in

        format : `lammps-dump`
            format of the input file

        Returns
        -------
        None

        See Also
        --------
        assign_particles
        """
        if format == 'lammps-dump':
            if os.path.exists(filename):
                filename = unicode(filename, "utf-8")
                pc.System.read_inputfile(self, filename)

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

    def set_nucsize_parameters(self, minfrenkel, threshold, avgthreshold):
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
        pc.System.set_nucsize_parameters(self, minfrenkel, threshold, avgthreshold)

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
        pc.System.calculate_nucsize(self)

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
        return pc.System.get_atom(self, index)

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
        pc.System.set_atom(self, atom)

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
        return pc.System.get_allatoms(self)

    def get_box(self):
        """
        Get the dimensions of the simulation box.

        Parameters
        ----------
        None

        Returns
        -------
        boxdims : list of box dimensions of length 2
            the return value consists of the vector of values in the form-
            [[box_x_low, box_x_high], [box_y_low, box_y_high], [box_z_low, box_z_high]]
        """
        box6dim = pc.System.get_box(self)
        pbox = [[box6dim[0], box6dim[1]], [box6dim[2], box6dim[3]], [box6dim[4], box6dim[5]]]
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
        return pc.System.get_absdistance(self, atom1, atom2)

    def get_neighbors(self, method="cutoff", cutoff=None):
        """
        Find neighbors of all atoms in the `System`. There are two methods to do this, the 
        traditional approach being the one in which the neighbors of an atom are the ones that lie
        in a cutoff distance around it. The second approach is using Voronoi polyhedra. All the atoms
        that share a Voronoi polyhedra face with the host atoms are considered its neighbors.

        Parameters
        ----------
        method : `cutoff` or `voronoi`, default: `cutoff`
            `cutoff` method finds atoms within a specified cutoff distance of the host atom
            `voronoi` method finds atoms that share a Voronoi polyhedra face with the host atom.

        cutoff : float
            the cutoff distance to be used for the `cutoff` based neighbor calculation method
            described above.
        
        Returns
        -------
        None

        See also
        --------
        reset_allneighbors
        """
        #first reset all neighbors
        pc.System.reset_allneighbors(self)

        if method == 'cutoff':
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
            q = [q]

        pc.System.calculate_q(self, q)

        if averaged:
            pc.System.calculate_aq(self, q)


    def get_connection(self, atom1, atom2):
        """
        Get the connection between two atoms. Connection is defined by Qlm(i).Qlm(j). Normally,
        a connection of more than 0.6 is considered a solid bond.

        Parameters
        ----------
        atom1 : `Atom` object
                first atom
        atom2 : `Atom` object
                second atom
        
        Returns
        -------
        connection : double
            connection between the first and second atom.

        See also
        --------
        calculate_frenkelnumbers
        find_clusters
        find_largest_cluster
        set_nucsize_parameters
        """
        return pc.System.get_number_from_bond(self, atom1, atom2)

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

    def find_clusters(self):
        """
        Find the clusters of all atoms in the system. Go through all the atoms and cluster them
        together.

        Parameters
        ----------
        None
        
        Returns
        -------
        None

        See also
        --------
        calculate_frenkelnumbers
        get_number_from_bond
        find_largest_cluster
        set_nucsize_parameters
        """
        pc.System.find_clusters(self)

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


















    








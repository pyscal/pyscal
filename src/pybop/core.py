import pybop.ccore as pc

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

		return pc.Atom.get_neighborweights()

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

	def get_q(self, q):
		"""
        get q value of the atom. 

        Parameters
        ----------
        q : int or list of int
            number of the required q - from 2-12

        Returns
        -------
        q : float or list of floats
            The queried q value

        Examples
        --------
        >>> q2 = atom.get_q(2)
        >>> q24 = atom.get_q([2, 4])

        See also
        --------
        set_q
        get_aq
        set_aq		
		"""
		if isinstance(q, int):
			rq = pc.Atom.get_q(self, q)
			return rq

		else:
			rq = [ pc.Atom.get_q(self, qq) for qq in q ]
			return rq

	def set_q(self, q, d):
		"""
        set the q value of the atom.

        Parameters
        ----------
        q : int or list of ints
            number of the required q - from 2-12
        d : float or list of floats
            the q value to set

        Returns
        -------
        None

        Examples
        --------
        >>> atom.set_q(2, 0.24)
        >>> atom.set_q([2,4], [0.24, 0.05])

        See also
        --------
        set_aq
        get_aq
        get_q
		"""
		if isinstance(q, int):
			pc.Atom.set_q(self, q)

		else:
			for count, qq in enumerate(q):
				pc.Atom.set_q(self, qq, d[count])

	def get_aq(self, q):
		"""
        get averaged q value of the atom. 

        Parameters
        ----------
        q : int or list of int
            number of the required q - from 2-12

        Returns
        -------
        q : float or list of floats
            The queried q value

        Examples
        --------
        >>> q2 = atom.get_aq(2)
        >>> q24 = atom.get_aq([2, 4])

        See also
        --------
        set_aq
        get_q
        set_q		
		"""
		if isinstance(q, int):
			rq = pc.Atom.get_aq(self, q)
			return rq

		else:
			rq = [ pc.Atom.get_q(self, qq) for qq in q ]
			return rq

	def set_aq(self, q, d):
		"""
        set the averaged q value of the atom.

        Parameters
        ----------
        q : int or list of ints
            number of the required q - from 2-12
        d : float or list of floats
            the q value to set

        Returns
        -------
        None

        Examples
        --------
        >>> atom.set_aq(2, 0.24)
        >>> atom.set_aq([2,4], [0.24, 0.05])

        See also
        --------
        set_q
        get_q
        get_aq
		"""
		if isinstance(q, int):
			pc.Atom.set_aq(self, q)

		else:
			for count, qq in enumerate(q):
				pc.Atom.set_aq(self, qq, d[count])

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

	def get_qlm(self, q):
		"""
        Get the real and imaginary qlm values of the atom.

        Parameters
        ----------
        q : int or list of ints
            number of the required q - from 2-12

        Returns
        -------
        qlms : 2D array of 2q+1 values or n 2D arrays
            the first part of the array is the 2q+1 real values
            second part is the 2q+1 imaginary values.

		"""
		if isinstance(q, int):
			rq = pc.Atom.get_qlm(self, q)
			return rq

		else:
			rq = [ pc.Atom.get_qlm(self, qq) for qq in q ]
			return rq

	def get_aqlm(self, q):
		"""
                Get the real and imaginary aqlm values of the atom.

                Parameters
                ----------
                q : int
                    number of the required q - from 2-12

                Returns
                -------
                qlms : 2D array of 2q+1 values
                    the first part of the array is the 2q+1 real values
                    second part is the 2q+1 imaginary values.
		"""
		if isinstance(q, int):
			rq = pc.Atom.get_aqlm(self, q)
			return rq

		else:
			rq = [ pc.Atom.get_aqlm(self, qq) for qq in q ]
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

	








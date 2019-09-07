import pyscal.ccore as pc
import pyscal.traj_process as ptp
import pyscal.pickle as pp
import os
import numpy as np
import warnings

"""
Definitions of class Atom.

"""
class Atom(pc.Atom):
    """
    Class to store atom details.

    Parameters
    ----------
    pos : list of floats of length 3
        position of the `Atom`, default [0,0,0]

    id : int
        id of the `Atom`, default 0

    type : int
        type of the `Atom`, default 1

    Notes
    -----
    A c++ class for holding the properties of a single atom. Various properties of the atom
    can be accessed through member functions which are described below in detail. Atoms can
    be created individually or directly by reading in a file. Check the examples for more 
    details on how atoms are created. For creating atoms directly from an input file check
    the documentation of `System` class.

    Although an `Atom` object can be created independently, `Atom` should be thought of 
    inherently as members of the `System` class. All the properties that define an atom are
    relative to the `System` class. `System` has a list of all atoms using which the neighbors
    of an `Atom`, if its solid and so on can be calculated. All the calculated properties of an
    atom which depend on any other atom, hence should be calculated through `System`. Please
    check the examples section of the documentation for more details. 


    Examples
    --------
    >>> #method 1 - individually
    >>> atom = Atom()
    >>> #now set positions of the atoms
    >>> atom.set_x([23.0, 45.2, 34.2])
    >>> #now set id
    >>> atom.set_id(23)
    
    """
    def __init__(self, pos=[0,0,0], id=0, type=1):
        """
        Deafults args
        """
        a=1
        pc.Atom.__init__(self)     
        pc.Atom.set_x(self, pos)
        pc.Atom.set_id(self, id)
        pc.Atom.set_type(self, type)

    #now wrapping for other normal functions
    def get_pos(self):
        """
        
        Get the position of the `Atom`. 

        Parameters
        ----------
        None
        
        Returns
        -------
        pos : array of float
            contains the position of the atom in the form `[posx, posy, posz]`, where
            `posx` is the x coordinate of the atom, `posy` is the y coordinate and `posz` 
            is the z coordinate. 

        Notes
        -----
        Meaningful values are only returned if the atoms are
        set before using this function.

        Examples
        --------
        >>> atom = Atom()
        >>> x = atom.get_pos()
      
        """
        x = pc.Atom.get_x(self)
        return x

    def set_pos(self, pos):
        """
        
        Set the position of the ``Atom``. 

        Parameters
        ----------
        pos : list of floats of length 3
            list contains three values which are the position coordinates of the `Atom` with
            respect to the simulation box.

        Returns
        -------
        None

        Examples
        --------
        >>> atom = Atom()
        >>> x = atom.set_pos([23.0, 45.2, 34.2])

        """
        if len(pos) == 3:
            try:
                pos = np.array(pos).astype(float)
                pc.Atom.set_x(self, pos)
            except:
                raise TypeError("Position values should be float.")
        else:
            raise ValueError("Length of position array should be 3.")

    def get_solid(self):
        """
        Find if an atom is solid or not.

        Parameters
        ----------
        None

        Returns
        -------
        issolid : int
            1 if solid, 0 otherwise
        
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
        structure : int  
            structural value
        
        Notes
        -----
        As of now it is not calculated using any inbuilt function. This can be used to store
        structure values that are calculated using other methods.

        """
        return pc.Atom.get_structure(self)

    def set_solid(self, issolid):
        """
        Set the solidity of atom

        Parameters
        ----------
        issolid : {0, 1} 
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
        pc.Atom.set_structure(self, structure)


    def get_volume(self, averaged = False):
        """
        
        Get the voronoi volume of the atom. 

        Parameters
        ----------
        averaged : bool  
            If True, averaged version of the volume is returned, default False
        
        Returns
        -------
        volume : float
            voronoi volume of the atom.

        Notes
        -----
        Calculation of volume happens during the neighbor calculation method. Meaningful 
        values are only returned if the neighbors are calculated using voronoi method.
        If keyword `averaged` is set to True, the volume of an atom is calculated as an
        average over itself and its neighbors. 

        Examples
        --------
        >>> volume = atom.get_volume()

  
        """
        if averaged:
            vol = pc.Atom.get_avgvolume(self)
        else:    
            vol = pc.Atom.get_volume(self)
        
        return vol

    def get_cluster(self):
        """
        Get the cluster properties of the atom. 

        Parameters
        ----------
        None
        
        Returns
        -------
        cluster : list of int of length 4  
            `cluster` is a vector of four values. they are described below-  
                `issolid` : which is 1 if the atom is solid, 0 otherwise  
                `issurface` : 1 if the atom has liquid neighbors, 0 otherwise  
                `lcluster` : 1 if the atom belongs to the largest cluster, 0 otherwise  
                `belongsto` : which gives the id of the cluster that the atom belongs to.  
        
        
        Notes
        -----
        The cluster properties of the atom include four different properties. The properties are only
        returned if they are calculated before using ``calculate_nucsize`` function before.

        Examples
        --------
        >>> cinfo = atom.get_cluster()

        """
        x = pc.Atom.get_cluster(self)
        return x

    def get_neighbors(self):
        """
        
        Returns the neighbors indices of the atom. 

        Parameters
        ----------
        None
        
        Returns
        -------
        x : list of int
            list of neighbor indices of the atom.

        Notes
        -----
        The list returned consistes of the indices
        of neighbor atom which indicate their position in the list of all atoms. The neighbors
        of an atom can be calculated from the ``System`` object that it belongs to.

        Examples
        --------
        neighbors = atom.get_neighbors()

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
        
        Returns the coordination number of the atom.  

        Parameters
        ----------
        None
        
        Returns
        -------
        cn : int
            coordination number of the atom.

        Notes
        -----
        ``System.get_neighbors`` function has to be used before accessing coordination numbers.

        Examples
        --------
        neighbors = atom.get_neighbors()

        """
        return pc.Atom.get_coordination(self)

    def get_neighborweights(self):
        """
        
        Get the neighbor weights of the atom. 
        
        Parameters
        ----------
        None
        
        Returns
        -------
        x : list of float
            neighbor weights

        Notes
        -----
        The neighbor weights are used to weight the contribution of each neighboring atom towards the 
        q value of the host atom [1]_. By default, each neighbor has a weight of 1 each. 
        However, if the neighbors are calculated using the  `System.get_neighbors(method='voronoi')`, each neighbor 
        atom gets a weight proportional to the face area shared between the neighboring atom and the 
        host atom (or higher powers [2]_). This can sometimes be helpful in controlling the contribution 
        of atoms with low face areas due to the thermal vibrations at high temperature.
        
        References
        ----------
        .. [1] Mickel, W, Kapfer, SC, Schroder-Turk, GE, Mecke, K, J Chem Phys 138, 2013
        .. [2] Haeberle, J, Sperl, M, Born, P 2019

        Examples
        --------
        >>> weights = atom.get_neighborweights()

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

        """
        pc.Atom.set_neighborweights(self, weights)

    def get_q(self, q, averaged = False):
        """
        get q value of the atom. 

        Parameters
        ----------
        q : int or list of int
            number of the required q - from 2-12

        averaged : bool
            If True, return the averaged q values,
            If False, return the non averaged ones
            default False

        Returns
        -------
        q : float or list of floats
            The queried q value

        Notes
        -----
        The q value can be either normal or can be averaged [1]_
        The averaged version can be obtained by using keyword
        `averaged = True`.

        References
        ----------
        .. [1] Lechner, W, Dellago, C, J Chem Phys, 2013

        Examples
        --------
        >>> q2 = atom.get_q(2, averaged = True)
        >>> q24 = atom.get_q([2, 4])

        """
        if isinstance(q, int):
            if not q in range(2, 13):
                raise ValueError("q values should be in range 2-12")
            if averaged:
                rq = pc.Atom.get_aq(self, q)
            else:
                rq = pc.Atom.get_q(self, q)
            return rq

        else:
            for qq in q:
                if not qq in range(2, 13):
                    raise ValueError("q values should be in range 2-12")
            if averaged:
                rq = [ pc.Atom.get_aq(self, qq) for qq in q ]
            else:
                rq = [ pc.Atom.get_q(self, qq) for qq in q ]
            return rq

    def set_q(self, q, d, averaged = False):
        """
        set the q value of the atom. 

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

        Notes
        -----
        The q value can be either normal or can be averaged [1]_
        The averaged version can be obtained by using keyword
        `averaged = True`.

        References
        ----------
        .. [1] Lechner, W, Dellago, C, J Chem Phys, 2013

        Examples
        --------
        >>> atom.set_q(2, 0.24, averaged = True)
        >>> atom.set_q([2,4], [0.24, 0.05])

        """
        if isinstance(q, int):
            if not q in range(2, 13):
                raise ValueError("q values should be in range 2-12")
            if averaged:
                pc.Atom.set_aq(self, q, d)
            else:
                pc.Atom.set_q(self, q, d)
        else:
            for qq in q:
                if not qq in range(2, 13):
                    raise ValueError("q values should be in range 2-12")
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

        """
        pc.Atom.set_id(self, idd)

    def get_loc(self):
        """

        get  the location of the atom in the array of all atoms in ``System``

        Parameters
        ----------
        None

        Returns
        -------
        loc : int
            loc of the atom

        Notes
        -----
        `loc` and `id` are different. The id of the atom is read in from the 
        input file. `loc` however is the position of the atom in the list of
        all atoms of the system.

        Examples
        --------
        >>> loc = atom.get_loc()

        """
        return pc.Atom.get_loc(self)

    def set_loc(self, idd):
        """
        
        set  the ``loc`` of the atom. 

        Parameters
        ----------
        idd : int
            loc of the atom

        Returns
        -------
        None

        Notes
        -----
        When an atom is put back in the ``System``, it will
        be reset if another atom exists in that ``loc``

        Examples
        --------
        >>> atom.set_loc(2)

        """
        #warnings.warn("If the loc of atom is changed and set to system, it will overwrite the existing data, if any.")
        pc.Atom.set_loc(self, idd)

    def get_type(self):
        """
        
        get  the ``type`` (species) of the atom.

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

        """
        if isinstance(tt, int):
            pc.Atom.set_type(self, tt)
        else:
            raise ValueError("type value should be integer")


    def get_vorovector(self):
        """
        get the voronoi structure identification vector. 

        Parameters
        ----------
        None

        Returns
        -------
        vorovector : array like, int  
            array of the form (n3, n4, n5, n6)

        Notes
        -----
        Returns a vector of the form ``(n3, n4, n5, n6)``, where ``n3`` is the number
        of faces with 3 vertices, ``n4`` is the number of faces with 4
        vertices and so on. This can be used to identify structures [1]_ [2]_.

        References
        ----------
        .. [1] Finney, JL, Proc. Royal Soc. Lond. A 319, 1970
        .. [2] Tanemura, M, Hiwatari, Y, Matsuda, H,Ogawa, T, Ogita, N, Ueda, A. Prog. Theor. Phys. 58, 1977

        """
        return pc.Atom.get_vorovector(self)

    def get_facevertices(self):
        """
        get the number of vertices of the voronoi face shared between an atom and its neighbors. 
        
        Parameters
        ----------
        None

        Returns
        -------
        facevertices : array like, int
            array of the vertices

        Notes
        -----
        Returns a vector with number of entries equal to the number of neighbors. 
        The corresponding atom indices can be obtained through ``Atom.get_neighbors``
        A shorter version of this vector
        in a condensed form is available through ``Atom.get_vorovector``.
  
        """
        return pc.Atom.get_facevertices(self)

#------------------------------------------------------------------------------------------------------------
"""
System class definitions
"""
#------------------------------------------------------------------------------------------------------------


class System(pc.System):
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
        pc.System.__init__(self)

    def read_inputfile(self, filename, format="lammps-dump", frame=-1, compressed = False):
        """
        
        Read input file containing the information of a time slice.
        
        Parameters
        ----------
        filename : string
            name of the input file to be read in

        format : {'lammps-dump', 'poscar'}
            format of the input file

        compressed : bool
            If True, force to read a `gz` compressed format, default False.

        frame : int
            If the trajectory contains more than one time slice, the slice can be specified
            using the `frame` option. 
            Alert: works only with `lammps-dump` format. 

        Returns
        -------
        None

        Notes
        -----
        `format` keyword specifies the format of the input file. Currently only
        a `lammps-dump` and `poscar` files are supported. However, this restriction can easily 
        be overcome using the `assign_particles` method from system where a list of atoms 
        and box vectors are directly provided to the system. This function itself uses the
        `pyscal.traj_process` module to process a file which is then assigned to system
        using `pyscal.core.assign_atoms`.

        `compressed` keyword is not required if a file ends with `.gz` extension, it is 
        automatically treated as a compressed file and this keyword is not necessary.

        `frame` keyword allows to read in a particular slice from a long trajectory. If all slices
        need to analysed, this is a very inefficient way. For handling multiple time slices,
        the `pyscal.traj_process` module offers a better set of tools.

        Triclinic simulation boxes can also be read in for lammps-dump. No special keyword is
        necessary.
        

        See Also
        --------
        assign_particles
        """
        if format == 'lammps-dump':
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
                    atoms, boxdims, box, triclinic = ptp.read_lammps_dump(filename, compressed=compressed, check_triclinic=True, box_vectors=True)
                    pc.System.assign_particles(self, atoms, boxdims)
                    if triclinic:
                        #we have to input rotation matrix and the inverse rotation matrix
                        rot = box.T
                        rotinv = np.linalg.inv(rot)
                        pc.System.assign_triclinic_params(self, rot, rotinv)
                else:
                    raise FileNotFoundError("input file %s not found"%filename)
                
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
                raise FileNotFoundError("input file %s not found"%filename)

        elif format == 'poscar':
            if os.path.exists(filename):
                atoms, boxdims = ptp.read_poscar(filename, compressed=compressed)
                pc.System.assign_particles(self, atoms, boxdims)
            else:
                raise FileNotFoundError("input file %s not found"%filename)
        else:
            raise TypeError("format recieved an unknown option %s"%format)


    def assign_atoms(self, atoms, box):
        """
        
        Assign atoms and box vectors to `System`. 

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
        of reading in the input file. If this method is used, there is no need of using
        `read_inputfile` method. Also using this function allows for reading of multiple
        file formats which are not supported by the inbuilt `read_inputfile` method.

        See Also
        --------
        read_inputfile
        """
        pc.System.assign_particles(self, atoms, box)

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
        been found. Otherwise it returns the default values.

        """
        return pc.System.get_largestcluster(self)

    
    def set_nucsize_parameters(self, cutoff, minfrenkel, threshold, avgthreshold):
        """
        Set the value of parameters for distinguishing solid and liquid atoms.

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

        Notes
        -----
        This function sets the parameters that can then be used to distinguish solid
        and liquid atoms. The number of atoms in the largest solid cluster in liquid 
        is often used as an order parameter in the study of nucleation during solidification. 
        A detailed description of the order parameter can be found in [1]_ and example 
        usage (not with this module) can be found in [2]_. This function merely sets the
        different parameters.  In order to actually complete the calculation, 
        `calculate_nucsize` has to be called after setting the parameters.

        References
        ----------
        .. [1] Auer, S, Frenkel, D. Adv Polym Sci 173, 2005
        .. [2] Diaz Leines, G, Drautz, R, Rogal, J, J Chem Phys 146, 2017

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
        For example, an `Atom` at location `i` in the list of all atoms in `System` can be queried by,
        ``atom = System.get_atom(i)``, then any kind of modification, for example, the 
        position of the `Atom` can done by, ``atom.set_pos([2.3, 4.5, 4.5])``. After 
        modification, the `Atom` can be set back to its position in `System` by
        ``System.set_atom(atom)``.

        If an atom already exists at that index in the list, it will be overwritten.
        
        """
        atomc = self.copy_atom_to_catom(atom)
        pc.System.set_atom(self, atomc)

    def get_atoms(self):
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
        box_vectors : bool, optional
            If True, return the whole box dimesions, default False

        Returns
        -------
        boxdims : list of box dimensions
            If `box_vectors` is false:
            the return value consists of the vector of values in the form-
            `[[box_x_low, box_x_high], [box_y_low, box_y_high], [box_z_low, box_z_high]]`
            If `box_vectors` is true:
            return the box vectors of the form
            `[[x1, x2, x3], [y1, y2, y3], [z1, z2, z3]]`
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
            `[[box_x_low, box_x_high], [box_y_low, box_y_high], [box_z_low, box_z_high]]`

        Returns
        -------
        None
        """
        pc.System.set_box(self, box)

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
                    rq = pc.System.get_aqvals(self, q)
                else:
                    rq = pc.System.get_qvals(self, q)
                return rq
            else:
                raise ValueError("the value of q should be between 2 and 12")

        else:
            for qq in q:
                if not qq in range(2, 13):
                    raise ValueError("the value of q should be between 2 and 12")
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


    def get_neighbors(self, method="cutoff", cutoff=None, threshold=2, filter=None, 
                                            voroexp=1, face_cutoff=0.002, padding=1.2, nlimit=6):
        """
        
        Find neighbors of all atoms in the `System`. 

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

            Adaptive cutoff
            uses a padding over the intial guessed "neighbor distance". By default it is 2. In case
            of a warning that ``threshold`` is inadequate, it should be further increased. High/low value
            of this parameter will correspond to the time taken for finding neighbors.

        filter : {'None', 'type'}, optional
            apply a filter to nearest neighbor calculation. If the `filter` keyword is set to
            `type`, only atoms of the same type would be included in the neighbor calculations. Default None.

        voroexp : int, optional 
            only used if ``method=voronoi``. Power of the neighbor weight used to weight the contribution of each atom towards the q 
            values. Default 1. 

            Higher powers can sometimes lead to better resolution. Works only with ``voronoi``
            neighbour method. See  arXiv:1906.08111v1 for more details.

        face_cutoff : double, optional
            only used if ``method=voronoi``. The minimum fraction of total voronoi face area a single phase should have in order to
            include it in the analysis of voronoi polyhedra to find (n_3, n_4, n_5, n_6) vector. Default 0.002

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
        This weight can later be used to weight steinhardt parameters. Higher powers of this weight can also be used. The keyword `voroexp`
        can be used to set this weight. If `voroexp` is set to 0, the neighbors would be calculated using Voronoi method, but Steinhardts
        parameters could be calculated normally.

        Keyword `filter` can be used to filter the neighbors based on a condition. Choosing ``filter='type'`` only considers an atom as
        a neighbor if both the neighbor atom and host atom are of the same type.

        References
        ----------
        .. [1] Stukowski, A, Model Simul Mater SC 20, 2012
        .. [2] van Meel, JA, Filion, L, Valeriani, C, Frenkel, D, J Chem Phys 234107, 2012

        """
        #first reset all neighbors
        pc.System.reset_allneighbors(self)
        pc.System.set_filter(self, 0)

        if filter == 'type':
            # type corresponds to 1
            pc.System.set_filter(self, 1)

        if method == 'cutoff':
            if cutoff=='sann':
                if threshold < 1:
                    raise ValueError("value of threshold should be at least 1.00")
                finished = pc.System.get_all_neighbors_sann(self, threshold)
                #if it finished without finding neighbors
                if not finished:
                    finallydone = False
                    for i in range(1,10):
                        #threshold value is probably too low
                        #try increasing threshold
                        warnings.warn("Could not find sann cutoff. trying with a higher threshold", RuntimeWarning)
                        pc.System.reset_allneighbors(self)
                        newfinished = pc.System.get_all_neighbors_sann(self, threshold*i)
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
                finished = pc.System.get_all_neighbors_adaptive(self, threshold, nlimit, padding)
                if not finished:
                    raise RuntimeError("Could not find adaptive cutoff")
            else:
                #warnings.warn("THIS RAN")    
                pc.System.set_neighbordistance(self, cutoff)
                pc.System.get_all_neighbors_normal(self)

        elif method == 'voronoi':
            pc.System.set_face_cutoff(self, face_cutoff)
            pc.System.set_alpha(self, int(voroexp))
            pc.System.get_all_neighbors_voronoi(self)
            

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
        pc.System.reset_allneighbors(self)

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
        `System.get_neighbors` for more details. If the keyword `average` is set to True,
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
        
        pc.System.calculate_q(self, qq)

        if averaged:
            pc.System.calculate_aq(self, qq)

    def calculate_frenkelnumbers(self):
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
        pc.System.calculate_frenkelnumbers(self)

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
        `find_clusters` has to be used before using this function.

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

        Notes
        -----
        This function is used to make sure that the user gets a python
        object rather than a C++ one.

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
        atom.set_facevertices(atomc.get_facevertices())
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
        atomc.set_facevertices(atom.get_facevertices())
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

        """

        #get the basic system indicators
        indicators = pc.System.get_indicators(self)

        #get box dims and triclinic params if triclinic
        boxdims = self.get_box()
        if indicators[6] == 1:
            rot = pc.System.get_triclinic_params(self)
        else:
            rot = 0

        #now finally get atoms
        atoms = pc.System.get_allatoms(self)
        #convert them to picklabale atoms
        patoms = [pp.pickle_atom(atom) for atom in atoms]

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

        """
        psys = self.prepare_pickle()
        np.save(file, psys)

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

        """

        psys = np.load(file, allow_pickle=True).flatten()[0]
        #set up indicators
        pc.System.set_indicators(self, psys.indicators)
        #unpickle atoms
        hatoms = [pp.unpickle_atom(atom) for atom in psys.atoms]
        catoms = [self.copy_atom_to_catom(atom) for atom in hatoms]
        boxdims = psys.boxdims


        #if triclinic, get those
        if psys.indicators[6] == 1:
            rot = psys.rot
            rotinv = np.linalg.inv(rot)
            pc.System.assign_triclinic_params(self, rot, rotinv)

        #assign atoms and box
        pc.System.reassign_particles(self, catoms, boxdims)

































    








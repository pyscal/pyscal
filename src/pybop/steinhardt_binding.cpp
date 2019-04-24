/*
This is the main file of steinhardt module where the bindings of the classes System and Atom are
provided. Additionally, the docstrings for the functions are also provided here. The docstrings should
be elaborate and ideally follow pep-8 conventions.
 */

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include <vector>
#include <string>
#include "steinhardt.h"


namespace py = pybind11;
using namespace std;


PYBIND11_MODULE(core, m) {
//bindings for Atom class
//------------------------------------------------------------------
    py::class_<Atom>(m,"Atom",             
        R"doc(
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


        )doc"
        )
        .def(py::init< >())
        
        .def("get_x",&Atom::gx,
            R"doc(
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

                )doc")

        .def("set_x",&Atom::sx,
            R"doc(
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

                )doc")


        .def("get_cluster",&Atom::gcluster,
            R"doc(
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

                )doc")

        .def("get_neighbors",&Atom::gneighbors,
            R"doc(
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

                )doc")

        .def("set_neighbors",&Atom::sneighbors,
            R"doc(
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

                )doc")

        .def("get_neighborweights",&Atom::gneighborweights,
            R"doc(
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

                )doc")

        .def("set_neighborweights",&Atom::sneighborweights,
            R"doc(
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

                )doc")

        .def("set_custom",&Atom::scustom,
            R"doc(
                NOT CLEAR
                Set custom values of an Atom

                Parameters
                ----------
                vector double : list of custom values
                
                Returns
                -------
                None

                )doc")

        .def("get_custom",&Atom::gcustom,
            R"doc(
                NOT CLEAR
                returns the custom values of an atom.

                Parameters
                ----------
                None
                
                Returns
                -------
                custom : vector double
                    custom values of the atom.

                )doc")

        .def("get_q",&Atom::gq,
            R"doc(
                get q value of the atom. 

                Parameters
                ----------
                q : int
                    number of the required q - from 2-12

                Returns
                -------
                q : float
                    The queried q value

                Examples
                --------
                >>> q2 = atom.get_q(2)

                See also
                --------
                set_q
                get_aq
                set_aq

                )doc")

        .def("get_id",&Atom::gid,
            R"doc(
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

                )doc")

        .def("set_id",&Atom::sid,
            R"doc(
                set  the id of the atom.

                Parameters
                ----------
                id : int
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

                )doc")


        .def("set_q",&Atom::sq,
            R"doc(
                set the q value of the atom.

                Parameters
                ----------
                q : int
                    number of the required q - from 2-12
                d : double
                    the q value to set

                Returns
                -------
                None

                Examples
                --------
                >>> atom.set_q(2, 0.24)

                See also
                --------
                set_aq
                get_aq
                get_q

                )doc")

        .def("get_aq",&Atom::gaq,
            R"doc(
                get averaged q value of the atom. 

                Parameters
                ----------
                q : int
                    number of the required q - from 2-12

                Returns
                -------
                q : float
                    The queried q value

                Examples
                --------
                >>> q2 = atom.get_q(2)

                See also
                --------
                set_q
                get_q
                set_aq

                )doc")


        .def("set_aq",&Atom::saq,
            R"doc(
                set the averaged q value of the atom.

                Parameters
                ----------
                q : int
                    number of the required q - from 2-12
                d : double
                    the q value to set

                Returns
                -------
                None

                Examples
                --------
                >>> atom.set_aq(2, 0.24)

                See also
                --------
                set_q
                get_q
                get_aq

                )doc")

        .def("get_qlm",&Atom::gqlm,
            R"doc(
                
                Get the real and imaginary qlm values of the atom.

                Parameters
                ----------
                q : int
                    number of the required q - from 2-12

                Returns
                -------
                qlms : 2D array of 2q+1 values
                    the first part of the array is the 2q+1 real values
                    second part is the 2q+1 imaginary values.
                
                )doc")


        .def("get_aqlm",&Atom::gaqlm,
            R"doc(
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
                
                )doc")

        .def("get_vorovector",&Atom::gvorovector,
            R"doc(
                
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
                
                )doc")

    ; 

    //bindings and documentation for individual functions
    py::class_<System>(m,"System",R"doc(
                
                A c++ class for holding the properties of a system. A `System` consists of two
                major components - the simulation box and the atoms. All the associated variables
                are then calculated over these.

                A `System` can be set and populated by reading an input file in lammps dump format.
                This enables for automatic reading of all atomic positions and the simulation box.

                Examples
                --------
                >>> sys = System()
                >>> sys.set_inputfile("atoms.dat")
                >>> sys.read_inputfile()
                
        )doc"
        )

        .def(py::init< >())
        
        .def("set_inputfile",&System::set_inputfile,
            R"doc(
                Set input file

                Parameters
                ----------
                None

                Returns
                -------
                None

                See Also
                --------
                read_inputfile - read an input file
                
               )doc")        

        .def("get_largestcluster",&System::glargestclusterid,
            R"doc(
                get id of the the largest cluster

                Parameters
                ----------
                None

                Returns
                -------
                clusterid : int 
                    id of the largest cluster

                See Also
                --------
                read_inputfile - read an input file

                )doc")        

        .def("set_cutoff",&System::set_neighbordistance,
            R"doc(
                Set cutoff for neighbor calculation.

                Parameters
                ----------
                cutoff : double
                    cutoff distance for calculating neighbors
                
                Returns
                -------
                None
                
                See Also
                --------
                read_inputfile - read an input file

                )doc")        

        .def("set_nucsize_parameters",&System::set_nucsize_parameters,
            R"doc(
                Set the value of parameters for calculating the largest cluster in the
                liquid, a detailed description of the order parameter can be found in  
                Diaz Leines et al, JCP 146(2017). http://doi.org/10.1063/1.4980082 .

                Parameters
                ----------
                minfrenkel : int
                    Minimum number of solid connections for an atom to be identified as
                    a solid.
                threshold : double
                    The cutoff value of connection between two atoms for them to be def
                    ined as having a bond.
                avgthreshold : double
                    Averaged value of connection between an atom and its neighbors for 
                    an atom to be solid.

                Returns
                -------
                None

                See Also
                --------
                set_inputfile - sets the input file for reading inc
                set_neighbordistance - sets the cutoff distance for neighbors of an atom.

                Examples
                --------
                >>> st.set_nucsize_parameters(7,0.5,0.5)

                )doc")
        

        .def("assign_particles",&System::assign_particles,
            R"doc(
                Assign atoms directly. Receive a vector of atom objects which is stored instead
                of reading in the input file. If this method is used, there is no need of using
                read_inputfile method.

                Parameters
                ----------
                atoms : vector Atoms
                    vector of Atom class instances
                box   : vector double
                    vector of box dimensions

                Returns
                -------
                None

                See Also
                --------
                read_inputfile - read an input file

                )doc")

        .def("calculate_nucsize",&System::calculate_nucsize,
            R"doc(
                Calculate the size of the largest cluster in the given system. Calculation
                the size of the largest cluster needs various prerequisites that can be set
                by the functions set_neighbordistance and set_nucsize_parameters. 
                For a detailed description of how the calculation works see-
                Diaz Leines et al, JCP 146(2017). http://doi.org/10.1063/1.4980082

                Parameters
                ----------
                None

                Returns
                -------
                cluster size : int
                    size of the largest cluster in number of atoms

                    )doc"
            )


        .def("get_atom",&System::gatom,
            R"doc(
                Access function that returns the Atom object at the queried position.

                Parameters
                ----------
                index : int
                        index of required atom

                Returns
                -------
                atom : Atom object
                    atom object at the queried position

                    )doc"
            )

       .def("set_atom",&System::satom,
            R"doc(
                return the atom to its original location after modification.

                Parameters
                ----------
                atom : Atom
                        atom to be replaced

                Returns
                -------
                None

                    )doc"
            )


        .def("get_allatoms",&System::gallatoms,
            R"doc(
                Access function that returns the a vector of Atom objects.

                Parameters
                ----------
                None

                Returns
                -------
                allatoms : vctor of Atom objects
                    all atoms in the system

                    )doc"
            )


        .def("get_box",&System::gboxdims,
            R"doc(
                Access function that returns the low and high values of sim box.

                Parameters
                ----------
                None

                Returns
                -------
                boxdims : vector of box dimensions

                    )doc"
            )


        .def("get_qvals",&System::gqvals,
            R"doc(
                return the required q values of all atoms.

                Parameters
                ----------
                q : int
                    required q value from 2-12

                Returns
                -------
                qvals : list of double
                    list of qth qvalue of all atoms.

                    )doc"
            )


        .def("get_aqvals",&System::gaqvals,
            R"doc(
                return the required aq values of all atoms.

                Parameters
                ----------
                q : int
                    required q value from 2-12

                Returns
                -------
                qvals : list of double
                    list of qth aqvalue of all atoms.

                    )doc"
            )


        .def("assign_cluster_info",&System::get_largest_cluster_atoms,
            R"doc(
                Assigns parameters such as if the atom belongs to the largest cluster,
                if it is on the surface.

                Parameters
                ----------
                None
                
                Returns
                -------
                None

                    )doc"
            )

        .def("read_inputfile",&System::read_particle_file,
            R"doc(
                Read a single snapshot of the lammps dump file and assign the positions
                and ids to an array of Atom objects stored in the parent class.

                Parameters
                ----------
                None
                
                Returns
                -------
                None

                    )doc"
            )

        .def("read_particleinstance",&System::read_particle_instance,
            R"doc(
                Read a single snapshot of the lammps dump file and assign the positions
                and ids to an array of Atom objects stored in the parent class.

                Parameters
                ----------
                None
                
                Returns
                -------
                None

                    )doc"
            )

        .def("get_absdistance", (double (System::*) (Atom, Atom))  &System::get_abs_distance,
            R"doc(
                Get the distance between two atoms.

                Parameters
                ----------
                atom1 : Atom object
                        first atom
                atom2 : Atom object
                        second atom
                
                Returns
                -------
                distance : double
                        distance between the first and second atom.
                
                    )doc"
            )

        .def("get_allneighbors", (void (System::*) (string &))  &System::get_all_neighbors, py::arg("method"), 
            R"doc(
                Find neighbors of all atoms in the system.

                Parameters
                ----------
                method
                
                Returns
                -------
                None

                    )doc"
            )

        .def("reset_allneighbors", (void (System::*) ()) &System::reset_all_neighbors,
            R"doc(
                Reset the neighbors of all atoms in the system.

                Parameters
                ----------
                None
                
                Returns
                -------
                None

                    )doc"
            )

        .def("calculate_complexQLM6",&System::calculate_complexQLM_6,
            R"doc(
                Find complex qlm 6 values for all atoms.

                Parameters
                ----------
                None
                
                Returns
                -------
                None

                    )doc"
            )


        .def("calculate_q",&System::calculate_q,
            R"doc(
                Find the q  value for all atoms. The required q indices(2-12) should be set
                using the set_reqd_qs function.

                Parameters
                ----------
                None
                
                Returns
                -------
                None

                    )doc"
            )

        .def("calculate_aq",&System::calculate_aq,
            R"doc(
                Find the averaged q  value for all atoms. The required q indices(2-12) should be set
                using the set_reqd_qs function.

                Parameters
                ----------
                None
                
                Returns
                -------
                None

                    )doc"
            )

        .def("set_reqdqs",&System::set_reqd_qs,
            R"doc(
                Set the list of qvalues to be calculated which will be done through 
                calculate_q function.

                Parameters
                ----------
                qlist : array of int
                    a list of q values to be calculated from 2-12
                
                Returns
                -------
                None

                    )doc"
            )


        .def("get_number_from_bond", (double (System::*) (Atom, Atom))  &System::get_number_from_bond,
            R"doc(
                Get the connection between two atoms. Connection is defined by Qlm(i).Qlm(j)

                Parameters
                ----------
                atom1 : Atom object
                        first atom
                atom2 : Atom object
                        second atom
                
                Returns
                -------
                connection : double
                        connection between the first and second atom.
                
                    )doc"
            )

        .def("calculate_frenkelnumbers",&System::calculate_frenkel_numbers,
            R"doc(
                Find frenkel numbers of all atoms in the system.

                Parameters
                ----------
                None
                
                Returns
                -------
                None

                    )doc"
            )

        .def("find_clusters",&System::find_clusters,
            R"doc(
                Find he clusters of all atoms in the system.

                Parameters
                ----------
                None
                
                Returns
                -------
                None

                    )doc"
            )

        .def("find_largest_cluster",&System::largest_cluster,
            R"doc(
                Find the largest in the system.

                Parameters
                ----------
                None
                
                Returns
                -------
                cluster : int
                    the size of the largest cluster

                    )doc"
            )


        ;


#ifdef VERSION_INFO
    m.attr("__version__") = VERSION_INFO;
#else
    m.attr("__version__") = "dev";
#endif
}

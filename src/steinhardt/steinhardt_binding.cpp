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

//module name would be steinhardt
/*
What we could think of is to compile this module as something like _steinhardt and provide the 
functions through a python module - the advantage of this is that the error checking and so on
can be done in the wrapping python code, which is great.
 */

PYBIND11_PLUGIN(steinhardt) {
//bindings for Atom class
//------------------------------------------------------------------
    py::module m("steinhardt");
    py::class_<Atom>(m,"Atom",             
        R"doc(
        A c++ class for holding the properties of a single atom. Various properties of the atom
        can be accessed through member functions which are described below.

        Examples
        --------
        atom = Atom()

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

                )doc")

        .def("get_cluster",&Atom::gcluster,
            R"doc(
                Get the cluster properties of the atom. The cluster properties of the atom
                include

                Parameters
                ----------
                None
                
                Returns
                -------
                cluster : vector int
                    cluster is a vector of four values. they are described below-
                        issolid - which is 1 if the atom is solid, 0 otherwise
                        issurface - 1 if the atom has liquid neighbors, 0 otherwise
                        lcluster - 1 if the atom belongs to the largest cluster, 0 otherwise
                        belongsto - which gives the id of the cluster that the atom belongs to.
                )doc")

        .def("get_neighbors",&Atom::gneighbors,
            R"doc(
                Returns the neighbors of the atom.

                Parameters
                ----------
                None
                
                Returns
                -------
                x : vector int
                    neighbor indices
                )doc")

        .def("set_neighbors",&Atom::sneighbors,
            R"doc(
                Set the neighbors of an atom

                Parameters
                ----------
                array : index of the neighbor atoms 
                
                Returns
                -------
                None

                )doc")

        .def("get_neighborweights",&Atom::gneighborweights,
            R"doc(
                Returns the neighbor weights of the atom.

                Parameters
                ----------
                None
                
                Returns
                -------
                x : vector float
                    neighbor weights
                )doc")

        .def("set_neighborweights",&Atom::sneighborweights,
            R"doc(
                Set the neighbor weights of an atom

                Parameters
                ----------
                array like float: weights of the neighbor atoms 
                
                Returns
                -------
                None

                )doc")

        .def("set_custom",&Atom::scustom,
            R"doc(
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
                get  q value of the atom.

                Parameters
                ----------
                q : int
                    number of the required q - from 2-12

                Returns
                -------
                q : float
                    The queried q value
                )doc")

        .def("get_id",&Atom::gid,
            R"doc(
                get  q value of the atom.

                Parameters
                ----------
                q : int
                    number of the required q - from 2-12

                Returns
                -------
                q : float
                    The queried q value
                )doc")


        .def("set_q",&Atom::sq,
            R"doc(
                set the q values of the atom.

                Parameters
                ----------
                q : int
                    number of the required q - from 2-12
                d : double
                    the q value to set

                Returns
                -------
                None
                )doc")

        .def("get_aq",&Atom::gaq,
            R"doc(
                get  avg q value of the atom.

                Parameters
                ----------
                q : int
                    number of the required q - from 2-12

                Returns
                -------
                q : float
                    The queried q value
                )doc")


        .def("set_aq",&Atom::saq,
            R"doc(
                set the avg q values of the atom.

                Parameters
                ----------
                q : int
                    number of the required q - from 2-12
                d : double
                    the q value to set

                Returns
                -------
                None
                )doc")

        .def("get_qlm",&Atom::gqlm,
            R"doc(
                get the real and imaginary qlm values of the atom.

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
                get the real and imaginary aqlm values of the atom.

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
        Class to hold a steinhardt system. It includes all the atoms and
        other system properties.

        A list of variables that can be set directly is provided.

        Attributes
        ----------
        File operations
        ---------------
        inputfile : string
            Name of the input file to read the atom information
        
        Simulation box
        --------------
        nop : int
            Number of atoms in the system.
        minfrenkel : int
            minimum number of frenkel connections to be identified as a 
            solid.
        boxx : float
            x dimension of the box
        boxy : float
            y dimension of the box
        boxz : float
            z dimension of the box
        neighbordistance : float
            cutoff distance to be used for neighbor calculations.
            accessible from python module as cutoff.

        Calculation of largest cluster
        ------------------------------
        threshold : float
            The cutoff value of connection between two atoms for them to be 
            defined as having a bond.
        avgthreshold : float
            Averaged value of connection between an atom and its neighbors for 
            an atom to be solid.
        maxclusterid : int
            id of the biggest cluster.

        )doc")

        .def(py::init< >())
        //.def_readwrite("inputfile", &System::inputfile)
        //.def_readwrite("nop", &System::nop)
        //.def_readwrite("minfrenkel", &System::minfrenkel)
        //.def_readwrite("boxx", &System::boxx)
        //.def_readwrite("boxy", &System::boxy)
        //.def_readwrite("boxz", &System::boxz)
        //.def_readwrite("cutoff", &System::neighbordistance)
        //.def_readwrite("threshold", &System::threshold)
        //.def_readwrite("avgthreshold", &System::avgthreshold)
        //.def_readwrite("maxclusterid", &System::maxclusterid)
        //minfrenkel function
        .def("set_inputfile",&System::set_inputfile,
            R"doc(
                Set input file

                Parameters
                ----------
                
                Returns
                -------
                
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

                Parametersb
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

        .def("get_allneighbors", (void (System::*) ()) &System::get_all_neighbors,
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

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include <vector>
#include <string>
#include "steinhardt.h"


namespace py = pybind11;
using namespace std;

PYBIND11_PLUGIN(steinhardt) {

//bindings for Atom class
//------------------------------------------------------------------
    py::module m("steinhardt");
    py::class_<Atom>(m,"Atom",             
        R"doc(
            Class to hold the details of an atom. The various 
            variables that can be accessed from python module
            are mentioned in attributes.

            Attributes
            ----------
            Basic atom properties
            ---------------------
            x : float
                x coordinate of atom
            y : float
                y coordinate of atom
            z : float
                z coordinate of atom
            id : int
                id of the atom

            Neighbor related properties
            ---------------------------
            n_neighbors : int
                number of neighbors of the atom
                Note that to access the list of neighbors, 
                gneighbors() method needs to be used.

            Cluster related properties
            --------------------------
            frenkelnumber : int
                frenkelnumber of the atom.
            issolid : int
                0 or 1. 1 if atom is solid.
            structure : int
                structure of the atom.
            belongsto : int
                id of the cluster to which atom belongs to.

       )doc"
        )
        .def(py::init< >())
        //.def_readwrite("id", &Atom::id)
        //.def_readwrite("x", &Atom::posx)
        //.def_readwrite("y", &Atom::posy)
        //.def_readwrite("z", &Atom::posz)
        //.def_readwrite("n_neighbors", &Atom::n_neighbors)
        //.def_readwrite("frenkelnumber", &Atom::frenkelnumber)
        //.def_readwrite("issolid", &Atom::issolid)
        //.def_readwrite("structure", &Atom::structure)
        //.def_readwrite("belongsto", &Atom::belongsto)

        .def("get_x",&Atom::gx,
            R"doc(
                Returns the position of the atom.

                Parameters
                ----------
                None
                
                Returns
                -------
                x : vector double
                    vector of x coordinates [posx, posy, posz]
                )doc")

        .def("get_cluster",&Atom::gcluster,
            R"doc(
                Returns the cluster properties of the atom.

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

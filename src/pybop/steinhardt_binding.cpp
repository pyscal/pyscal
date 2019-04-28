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


PYBIND11_MODULE(ccore, m) {
//bindings for Atom class
//------------------------------------------------------------------
    py::class_<Atom>(m,"Atom")
        .def(py::init< >())
        .def("get_x",&Atom::gx)
        .def("set_x",&Atom::sx)
        .def("get_cluster",&Atom::gcluster)
        .def("get_neighbors",&Atom::gneighbors)
        .def("get_coordination",&Atom::gnneighbors)
        .def("set_neighbors",&Atom::sneighbors)
        .def("get_neighborweights",&Atom::gneighborweights)
        .def("set_neighborweights",&Atom::sneighborweights)
        .def("set_custom",&Atom::scustom)
        .def("get_custom",&Atom::gcustom)
        .def("get_q",&Atom::gq)
        .def("get_id",&Atom::gid)
        .def("set_id",&Atom::sid)
        .def("set_q",&Atom::sq)
        .def("get_aq",&Atom::gaq)
        .def("set_aq",&Atom::saq)
        .def("get_qlm",&Atom::gqlm)
        .def("get_aqlm",&Atom::gaqlm)
        .def("get_vorovector",&Atom::gvorovector); 

    //bindings and documentation for individual functions
    py::class_<System>(m,"System")
        .def(py::init< >())
        .def("read_inputfile",&System::read_particle_file,
            R"doc(
                Read input file containing the information of a time slice from a molecular dynamics
                simulation. As of now, the file should be a lammps dump format and can only have a
                specific header format. That is-
                id type mass x y z vx vy vz
                However, this restriction can easily be overcome using the `assign_particles` method
                from system where a list of atoms and box vectors are directly provided to the system.
                This functions only sets the input file
                Parameters
                ----------
                None

                Returns
                -------
                None

                See Also
                --------
                assign_particles
                
               )doc")        

        .def("get_largestcluster",&System::glargestclusterid,
            R"doc(
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

                )doc")             

        .def("set_nucsize_parameters",&System::set_nucsize_parameters,
            R"doc(
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

                )doc")
        

        .def("assign_particles", &System::assign_particles,
            R"doc(
                Assign atoms directly. Receive a vector of atom objects which is stored instead
                of reading in the input file. If this method is used, there is no need of using
                `read_inputfile` method. Also using this function allows for reading of multiple
                file formats which are not supported by the inbuilt `read_inputfile` method.

                Parameters
                ----------
                atoms : list of `Atom` objects
                    list consisting of all atoms
                box   : list of float of length 6
                    list which consists of the box dimensions in the format-
                    [box_x_low, box_x_high, box_y_low, box_y_high, box_z_low, box_z_high]

                Returns
                -------
                None

                See Also
                --------
                read_inputfile


                )doc")

        .def("calculate_nucsize",&System::calculate_nucsize,
            R"doc(
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

                    )doc"
            )


        .def("get_atom",  &System::gatom,
            R"doc(
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

                    )doc"
            )

       .def("set_atom", &System::satom,
            R"doc(
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

                    )doc"
            )


        .def("get_allatoms",&System::gallatoms,
            R"doc(
                Get a list of all `Atom` objects that belong to the system.

                Parameters
                ----------
                None

                Returns
                -------
                allatoms : list of `Atom` objects
                    all atoms in the system

                    )doc"
            )


        .def("get_box",&System::gboxdims,
            R"doc(
                Get the dimensions of the simulation box.

                Parameters
                ----------
                None

                Returns
                -------
                boxdims : list of box dimensions of length 2
                    the return value consists of the vector of values in the form-
                    [[box_x_low, box_x_high], [box_y_low, box_y_high], [box_z_low, box_z_high]]

                    )doc"
            )

        .def("set_box",&System::sbox,
            R"doc(
                Set the dimensions of the simulation box.

                Parameters
                ----------
                boxdims : list of box dimensions of length 6
                    the return value consists of the vector of values in the form-
                    [[box_x_low, box_x_high], [box_y_low, box_y_high], [box_z_low, box_z_high]]

                Returns
                -------
                None

                    )doc"
            )

        .def("get_qvals",&System::gqvals,
            R"doc(
                Get the required q values of all atoms. The function returns a list of 
                q values in the same order as that of the atoms in the system.

                Parameters
                ----------
                q : int
                    required q value from 2-12

                Returns
                -------
                qvals : list of floats
                    list of qvalue of all atoms.

                    )doc"
            )


        .def("get_aqvals",&System::gaqvals,
            R"doc(
                Get the required averaged q values of all atoms. The function returns a list of 
                q values in the same order as that of the atoms in the system.

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


        .def("get_absdistance", (double (System::*) (Atom, Atom))  &System::get_abs_distance,
            R"doc(
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
                
                    )doc"
            )

        .def("get_allneighbors", &System::get_all_neighbors, py::arg("method"), py::arg("cutoff") = 3.0, 
            R"doc(
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

                    )doc"
            )

        .def("reset_allneighbors", &System::reset_all_neighbors,
            R"doc(
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

                    )doc"
            )

        .def("calculate_q",&System::calculate_q,
            R"doc(
                Find the bond order parameter - q  for all atoms. Any of the q parameters from 2-12 can be provided. 

                Parameters
                ----------
                qs : list of ints
                    A list of all q params to be found from 2-12. Even if a single value is required, it
                    has be enclosed in a list. For example, [2] or [2,4,6] and so on.
                
                Returns
                -------
                None

                See also
                --------
                calculate_aq

                    )doc"
            )

        .def("calculate_aq",&System::calculate_aq,
            R"doc(
                Find the averaged bond order parameter q for all atoms. Any of the q parameters from 2-12 can be provided.
                If `calculate_aq` function is used before `calculate_q`, q parameters are first calculated before
                the averaged versions are calculated. See Lechner, Dellago, JCP 129, 2008. for a description of the
                averaged bond order parameters.

                Parameters
                ----------
                qs : list of ints
                    A list of all q params to be found from 2-12. Even if a single value is required, it
                    has be enclosed in a list. For example, [2] or [2,4,6] and so on.
                
                Returns
                -------
                None

                See also
                --------
                calculate_q

                    )doc"
            )

        .def("get_number_from_bond", (double (System::*) (Atom, Atom))  &System::get_number_from_bond,
            R"doc(
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
                
                    )doc"
            )

        .def("calculate_frenkelnumbers",&System::calculate_frenkel_numbers,
            R"doc(
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
                    
                    )doc"
            )

        .def("find_clusters",&System::find_clusters,
            R"doc(
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

                    )doc"
            )

        .def("find_largest_cluster",&System::largest_cluster,
            R"doc(
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

                    )doc"
            )


        ;


#ifdef VERSION_INFO
    m.attr("__version__") = VERSION_INFO;
#else
    m.attr("__version__") = "dev";
#endif
}

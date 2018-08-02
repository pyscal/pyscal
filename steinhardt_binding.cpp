#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include <vector>
#include "steinhardt.h"

namespace py = pybind11;
using namespace std;

PYBIND11_MODULE(steinhardt, m) {
    

    //bindings and documentation for individual functions
    py::class_<System>(m,"System")
        .def(py::init< >())

        //minfrenkel function
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
        
        .def("set_inputfile",&System::set_inputfile,
            R"doc(
                Set the inputfile for reading in for calculations. Currently, only a lammps
                dump file can be used.

                Parameters
                ----------
                inputfile : string
                    filename of the file to be read

                Returns
                -------
                None
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
                )doc")

        .def("set_neighbordistance",&System::set_neighbordistance,
            R"doc(
                Set the cutoff distance for determining the neighbours of an atom.

                Parameters
                ----------
                cutoff : double
                    neighbor distance

                Returns
                -------
                None
                    )doc"
            )
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

        .def("gatom",&System::gatom,
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

        .def("glargestclusterid",&System::glargestclusterid,
            R"doc(
                Access function that returns the id of largest cluster. This can be used in 
                combination with gid() method of Atom to find if an atom belongs to the 
                largest cluster. eg - if( atom.gbelongsto()==system.glargestclusterid() )

                Parameters
                ----------
                None
                
                Returns
                -------
                cluster id  : int
                    id of the largest cluster.

                    )doc"
            )

        ;

//bindings for Atom class
//------------------------------------------------------------------

    py::class_<Atom>(m,"Atom")
        .def(py::init< >())     
        .def("gx",&Atom::gx,
            R"doc(
                Returns the coordinates of the atom.

                Parameters
                ----------
                None
                
                Returns
                -------
                x : vector double
                    position coordinates
                )doc")

        .def("gneighbors",&Atom::gneighbors,
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

        .def("gn_neighbors",&Atom::gn_neighbors,
            R"doc(
                Returns the number of neighbors of the atom.

                Parameters
                ----------
                None
                
                Returns
                -------
                n_neighbors : int 
                    number of neighbors
                )doc")

        .def("gfrenkelnumber",&Atom::gfrenkelnumber,
            R"doc(
                Returns the frenkelnumber of the atom.

                Parameters
                ----------
                None
                
                Returns
                -------
                frenkelnumber : int
                    frenkelnumber
                )doc")

        .def("gissolid",&Atom::gissolid,
            R"doc(
                returns the solidity of an atom.

                Parameters
                ----------
                None
                
                Returns
                -------
                issolid : int
                    A value of 1 is returned if the atom is solid, 0 otherwise.
                )doc")

        .def("gid",&Atom::gid,
            R"doc(
                Returns the id of the atom.

                Parameters
                ----------
                None
                
                Returns
                -------
                id : int
                    id of the atom
                )doc")

        .def("gbelongsto",&Atom::gbelongsto,
            R"doc(
                Returns the cluster id of the atom.

                Parameters
                ----------
                None
                
                Returns
                -------
                id : int
                    cluster id of the atom. A value of -1 means that the atom is liquid.
                )doc")

        .def("sx",&Atom::sx,
            R"doc(
                Set the coordinates of the atom.

                Parameters
                ----------
                pos : vector double
                    a vector containing x,y and z coordinates of the atom
                
                Returns
                -------
                None

                )doc")

        .def("sid",&Atom::sid,
            R"doc(
                sets the id of the atom.

                Parameters
                ----------
                n : int
                    id of the atom

                Returns
                -------
                None.
                )doc")

    ;
/*
    m.def("get_absDistance",&get_absDistance,
        R"pbdoc(Calculate distance between two atoms with PBC
        
        Parameters
        ------------
        atom1 : Atom object
                first atom
        atom2 : Atom object
                second atom
        box   : simulation box

        Returns
        ------------
        abs   : Distance between the atoms 

        )pbdoc"
        );

    m.def("convert_SphericalCoordinates",&convert_SphericalCoordinates,
        R"pbdoc(Find spherical coordinates
        
        Parameters
        ------------
        diff : double vector
               diff in coordinates between two atoms

        Returns
        ------------
        res   : The spherical coordinates 

        )pbdoc"
        );

    m.def("get_AllNeighborsandDistances",&get_AllNeighborsandDistances,
        R"pbdoc(get all neighbors, neighbor distances, theta, phi and r values
        
        Parameters
        ------------
        atoms : atom vector
                vector of all atoms in the system

        Returns
        ------------
        atoms : atom vector
                vector of all atoms in the system

        )pbdoc"
        );

    m.def("calculate_complexQLM",&calculate_complexQLM,
        R"pbdoc(calculate q component value for all atoms
        
        Parameters
        ------------
        atoms : atom vector
                vector of all atoms in the system
        box   : Simbox class
        param : params class

        Returns
        ------------
        atoms : atom vector
                vector of all atoms in the system

        )pbdoc"
        );

    m.def("calculate_aQ",&calculate_aQ,
        R"pbdoc(calculate average q values
        
        Parameters
        ------------
        atoms : atom vector
                vector of all atoms in the system
        box   : Simbox class
        param : params class

        Returns
        ------------
        atoms : atom vector
                vector of all atoms in the system

        )pbdoc"
        );

    m.def("QLM",&QLM,
        R"pbdoc(calculate QLM values
        
        Parameters
        ------------
        q       :int
                        qvalue
        m       :int
                        m value
        theta   :double - angle
        phi     :double - angle
        
        Returns
        ------------
        rqlm : double vector
                real and imaginary values of qlm

        )pbdoc"
        );

    m.def("YLM",&YLM,
        R"pbdoc(calculate YLM values
        
        Parameters
        ------------
        q       :int
                        qvalue
        m       :int
                        m value
        theta   :double - angle
        phi     :double - angle
        
        Returns
        ------------
        rqlm : double vector
                real and imaginary values of ylm

        )pbdoc"
        );
    
    m.def("PLM",&PLM,
        R"pbdoc(calculate PLM values
        
        Parameters
        ------------
        q       :int
                        qvalue
        m       :int
                        m value
        cos(theta)   :double - cos of angle
        
        Returns
        ------------
        mplm : double 
                       
        )pbdoc"
        );
    
    m.def("find_solids",&find_solids,
        R"pbdoc(calculate average q values
        
        Parameters
        ------------
        atoms : atom vector
                vector of all atoms in the system
        box   : Simbox class
        param : params class

        Returns
        ------------
        atoms : atom vector
                vector of all atoms in the system

        )pbdoc"
        );    
    m.def("find_clusters",&find_clusters,
        R"pbdoc(calculate average q values
        
        Parameters
        ------------
        atoms : atom vector
                vector of all atoms in the system
        box   : Simbox class
        param : params class

        Returns
        ------------
        atoms : atom vector
                vector of all atoms in the system

        )pbdoc"
        );  
    
    m.def("largest_cluster",&largest_cluster,
        R"pbdoc(calculate average q values
        
        Parameters
        ------------
        atoms : atom vector
                vector of all atoms in the system
        box   : Simbox class
        param : params class

        Returns
        ------------
        atoms : atom vector
                vector of all atoms in the system

        )pbdoc"
        );  

    m.def("nucsize",&nucsize,
        R"pbdoc(calculate average q values
        
        Parameters
        ------------
        atoms : atom vector
                vector of all atoms in the system
        box   : Simbox class
        param : params class

        Returns
        ------------
        atoms : atom vector
                vector of all atoms in the system

        )pbdoc"
        );

*/
#ifdef VERSION_INFO
    m.attr("__version__") = VERSION_INFO;
#else
    m.attr("__version__") = "dev";
#endif
}
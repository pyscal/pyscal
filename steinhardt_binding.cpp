#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include <vector>
#include "steinhardt.h"

namespace py = pybind11;
using namespace std;

PYBIND11_MODULE(steinhardt, m) {
    
    //bindings for Atom class
    //nothing to bind right now;
    //its internal to system
    /*
    py::class_<Atom>(m,"Atom")
        .def(py::init<vector<double>,vector<double>,int>())
        .def("gx",&Atom::gx)
        .def("gv",&Atom::gv)
        .def("gid",&Atom::gid)
        .def("gr",&Atom::gr)
        .def("gtheta",&Atom::gtheta)
        .def("gphi",&Atom::gphi)
        .def("gaq",&Atom::gaq)
        .def("gneighdist",&Atom::gneighdist)
        .def("gneighbors",&Atom::gneighbors)
        .def("gneighborcount",&Atom::gneighborcount)
        .def("gneighborcount",&Atom::gneighborcount)
        .def("gfrenkelcons",&Atom::gfrenkelcons)
        .def("gissolid",&Atom::gissolid)
        .def("gcluster",&Atom::gcluster)
        ;
    */
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
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
        .def("set_minfrenkel",&System::set_minfrenkel)
        .def("set_inputfile",&System::set_inputfile)
        .def("set_neighbordistance",&System::set_neighbordistance)
        .def("calculate_largestClusterparameter_Full",&System::calculate_largestClusterparameter_Full)
        ;
    /*  
    m.def("get_diff",&get_diff,
        R"pbdoc(Calculate the diff in cordinate between two atoms with PBC
        
        Parameters
        ------------
        atom1 : Atom object
                first atom
        atom2 : Atom object
                second atom
        box   : simulation box

        Returns
        ------------
        diff   : diff between coordinates of the atoms 

        )pbdoc"
        );


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
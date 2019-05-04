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
        .def("set_cluster",&Atom::scluster)
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
        .def("get_loc",&Atom::gloc)
        .def("set_loc",&Atom::sloc)
        .def("get_allq",&Atom::gallq)
        .def("get_allaq",&Atom::gallaq)
        .def("set_allq",&Atom::sallq)
        .def("set_allaq",&Atom::sallaq)
        .def("get_type",&Atom::gtype)
        .def("set_type",&Atom::stype)
        .def("set_q",&Atom::sq)
        .def("get_aq",&Atom::gaq)
        .def("set_aq",&Atom::saq)
        .def("get_qlm",&Atom::gqlm)
        .def("get_aqlm",&Atom::gaqlm)
        .def("get_vorovector",&Atom::gvorovector)
        .def("set_vorovector",&Atom::svorovector); 

    //bindings and documentation for individual functions
    py::class_<System>(m,"System")
        .def(py::init< >())
        .def("read_inputfile",&System::read_particle_file)        
        .def("get_largestcluster",&System::glargestclusterid)             
        .def("set_nucsize_parameters",&System::set_nucsize_parameters)
        .def("assign_particles", &System::assign_particles)
        .def("calculate_nucsize",&System::calculate_nucsize)
        .def("get_atom",  &System::gatom)
        .def("set_atom", &System::satom)
        .def("get_allatoms",&System::gallatoms)
        .def("get_box",&System::gboxdims)
        .def("set_box",&System::sbox)
        .def("get_qvals",&System::gqvals)
        .def("get_aqvals",&System::gaqvals)
        .def("get_absdistance", (double (System::*) (Atom, Atom))  &System::get_abs_distance)
        .def("get_all_neighbors_normal",&System::get_all_neighbors_normal)
        .def("get_all_neighbors_voronoi",&System::get_all_neighbors_voronoi)
        .def("set_neighbordistance", &System::set_neighbordistance)
        .def("reset_allneighbors", &System::reset_all_neighbors)
        .def("calculate_q",&System::calculate_q)
        .def("calculate_aq",&System::calculate_aq)
        .def("get_number_from_bond", (double (System::*) (Atom, Atom))  &System::get_number_from_bond)
        .def("calculate_frenkelnumbers",&System::calculate_frenkel_numbers)
        .def("find_clusters",&System::find_clusters)
        .def("find_largest_cluster",&System::largest_cluster)
        ;


#ifdef VERSION_INFO
    m.attr("__version__") = VERSION_INFO;
#else
    m.attr("__version__") = "dev";
#endif
}

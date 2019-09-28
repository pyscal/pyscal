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
#include "system.h"


namespace py = pybind11;
using namespace std;

class PySystem : public System {
public:
  using System::System;
};



PYBIND11_MODULE(csystem, m) {
    py::options options;
    options.disable_function_signatures();


    //bindings and documentation for individual functions
    py::class_<System>(m,"System")
        .def(py::init< >())
        //.def("assign_particles", &System::assign_particles, py::keep_alive<1, 2>())
        .def_property("largest_clusterid", &System::glargestclusterid, &System::slargestclusterid)
        .def("set_nucsize_parameters",&System::set_nucsize_parameters)
        .def_property("box", &System::gbox, &System::sbox )
        .def_property("atoms", &System::get_atoms, &System::set_atoms)
        //.def("get_atoms", &System::get_atoms)
        //.def("set_atoms", &System::set_atoms)
        .def("calculate_nucsize",&System::calculate_nucsize)
        .def("cget_atom",  &System::gatom)
        .def("cset_atom", &System::satom)
        .def_property("voroexp", &System::galpha, &System::salpha)
        .def_property("nop",&System::gnop,&System::snop)
        .def("cget_qvals",&System::gqvals)
        .def("cget_aqvals",&System::gaqvals)
        .def("get_absdistance", (double (System::*) (Atom, Atom))  &System::get_abs_distance)
        .def("get_all_neighbors_normal",&System::get_all_neighbors_normal)
        .def("get_all_neighbors_sann",&System::get_all_neighbors_sann)
        .def("get_all_neighbors_adaptive",&System::get_all_neighbors_adaptive)
        .def("get_all_neighbors_voronoi",&System::get_all_neighbors_voronoi)
        .def("set_neighbordistance", &System::set_neighbordistance)
        .def("reset_allneighbors", &System::reset_all_neighbors)
        .def("ccalculate_q",&System::calculate_q)
        .def("ccalculate_aq",&System::calculate_aq)
        .def("get_number_from_bond", (double (System::*) (Atom, Atom))  &System::get_number_from_bond)
        .def("calculate_frenkelnumbers",&System::calculate_frenkel_numbers)
        .def("cfind_clusters",&System::find_clusters)
        .def("cfind_clusters_recursive",&System::find_clusters_recursive)
        .def("find_largest_cluster",&System::largest_cluster)
        .def("get_largest_cluster_atoms",&System::get_largest_cluster_atoms)
        .def_property("filter",&System::gfilter,&System::sfilter)
        .def("assign_triclinic_params",&System::assign_triclinic_params)
        .def("get_triclinic_params",&System::get_triclinic_params)
        .def("get_boxvecs",&System::gboxvecs)
        .def("get_pairdistances",&System::get_pairdistances)
        .def("find_average_volume",&System::find_average_volume)
        .def("get_indicators",&System::get_indicators)
        .def("set_indicators",&System::set_indicators)
        .def("find_solid_atoms",&System::find_solid_atoms)
        ;


#ifdef VERSION_INFO
    m.attr("__version__") = VERSION_INFO;
#else
    m.attr("__version__") = "dev";
#endif
}

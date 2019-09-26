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


PYBIND11_MODULE(csystem, m) {
    py::options options;
    options.disable_function_signatures();

    //bindings and documentation for individual functions
    py::class_<System>(m,"System")
        .def(py::init< >())
        .def("read_inputfile",&System::read_particle_file)
        .def("get_largestcluster",&System::glargestclusterid)
        .def("set_nucsize_parameters",&System::set_nucsize_parameters)
        .def("assign_particles", &System::assign_particles)
        .def("reassign_particles", &System::reassign_particles)
        .def("calculate_nucsize",&System::calculate_nucsize)
        .def("get_atom",  &System::gatom)
        .def("set_atom", &System::satom)
        .def("set_alpha", &System::salpha)
        .def("get_allatoms",&System::gallatoms)
        .def("get_box",&System::gboxdims)
        .def("get_nop",&System::gnop)
        .def("set_box",&System::sbox)
        .def("get_qvals",&System::gqvals)
        .def("get_aqvals",&System::gaqvals)
        .def("get_absdistance", (double (System::*) (Atom, Atom))  &System::get_abs_distance)
        .def("get_all_neighbors_normal",&System::get_all_neighbors_normal)
        .def("get_all_neighbors_sann",&System::get_all_neighbors_sann)
        .def("get_all_neighbors_adaptive",&System::get_all_neighbors_adaptive)
        .def("get_all_neighbors_voronoi",&System::get_all_neighbors_voronoi)
        .def("set_neighbordistance", &System::set_neighbordistance)
        .def("reset_allneighbors", &System::reset_all_neighbors)
        .def("calculate_q",&System::calculate_q)
        .def("calculate_aq",&System::calculate_aq)
        .def("get_number_from_bond", (double (System::*) (Atom, Atom))  &System::get_number_from_bond)
        .def("calculate_frenkelnumbers",&System::calculate_frenkel_numbers)
        .def("find_clusters",&System::find_clusters)
        .def("find_clusters_recursive",&System::find_clusters_recursive)
        .def("find_largest_cluster",&System::largest_cluster)
        .def("get_largest_cluster_atoms",&System::get_largest_cluster_atoms)
        .def("set_filter",&System::sfilter)
        .def("assign_triclinic_params",&System::assign_triclinic_params)
        .def("get_triclinic_params",&System::get_triclinic_params)
        .def("get_boxvecs",&System::gboxvecs)
        .def("get_pairdistances",&System::get_pairdistances)
        .def("find_average_volume",&System::find_average_volume)
        .def("set_face_cutoff",&System::set_face_cutoff)
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

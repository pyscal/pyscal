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
        //no binding - delete
        .def("read_inputfile",&System::read_particle_file)
        //change to property
        .def("get_largestcluster",&System::glargestclusterid)
        //needed
        .def("set_nucsize_parameters",&System::set_nucsize_parameters)
        //make atoms and box property
        .def("assign_particles", &System::assign_particles)
        //needed
        .def("reassign_particles", &System::reassign_particles)
        //old fn
        .def("calculate_nucsize",&System::calculate_nucsize)
        //needed
        .def("get_atom",  &System::gatom)
        //needed
        .def("set_atom", &System::satom)
        //property
        .def("set_alpha", &System::salpha)
        //make to peorpety
        .def("get_allatoms",&System::gallatoms)
        //property
        .def("get_box",&System::gboxdims)
        //property
        .def("get_nop",&System::gnop)
        //property
        .def("set_box",&System::sbox)
        //needed
        .def("get_qvals",&System::gqvals)
        //needed
        .def("get_aqvals",&System::gaqvals)
        //needed
        .def("get_absdistance", (double (System::*) (Atom, Atom))  &System::get_abs_distance)
        //all of them as it is
        .def("get_all_neighbors_normal",&System::get_all_neighbors_normal)
        .def("get_all_neighbors_sann",&System::get_all_neighbors_sann)
        .def("get_all_neighbors_adaptive",&System::get_all_neighbors_adaptive)
        .def("get_all_neighbors_voronoi",&System::get_all_neighbors_voronoi)
        //needed
        .def("set_neighbordistance", &System::set_neighbordistance)
        //needed
        .def("reset_allneighbors", &System::reset_all_neighbors)
        //needed - but rename the binding
        .def("calculate_q",&System::calculate_q)
        .def("calculate_aq",&System::calculate_aq)
        //needed
        .def("get_number_from_bond", (double (System::*) (Atom, Atom))  &System::get_number_from_bond)
        //needed
        .def("calculate_frenkelnumbers",&System::calculate_frenkel_numbers)
        //needed
        .def("find_clusters",&System::find_clusters)
        .def("find_clusters_recursive",&System::find_clusters_recursive)
        //property
        .def("find_largest_cluster",&System::largest_cluster)
        //internally?
        .def("get_largest_cluster_atoms",&System::get_largest_cluster_atoms)
        //property
        .def("set_filter",&System::sfilter)
        //needed?
        .def("assign_triclinic_params",&System::assign_triclinic_params)
        //needed?
        .def("get_triclinic_params",&System::get_triclinic_params)
        //property
        .def("get_boxvecs",&System::gboxvecs)
        //needed
        .def("get_pairdistances",&System::get_pairdistances)
        //needed
        .def("find_average_volume",&System::find_average_volume)
        //remove
        .def("set_face_cutoff",&System::set_face_cutoff)
        //property
        .def("get_indicators",&System::get_indicators)
        //property
        .def("set_indicators",&System::set_indicators)
        //needed
        .def("find_solid_atoms",&System::find_solid_atoms)
        ;


#ifdef VERSION_INFO
    m.attr("__version__") = VERSION_INFO;
#else
    m.attr("__version__") = "dev";
#endif
}

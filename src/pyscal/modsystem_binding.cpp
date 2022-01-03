#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>
#include <vector>
#include <string>
#include "modsystem.h"
#include <map>
#include <string>
#include <any>

namespace py = pybind11;
using namespace std;


PYBIND11_MODULE(csystem, m) {
    py::options options;
    options.disable_function_signatures();
    m.def("get_abs_distance", &get_abs_distance);
    m.def("get_distance_vector", &get_distance_vector);
    m.def("reset_all_neighbors", &reset_all_neighbors);
    m.def("get_all_neighbors_normal", &get_all_neighbors_normal);
    m.def("get_all_neighbors_cells", &get_all_neighbors_cells);
    m.def("get_all_neighbors_bynumber", &get_all_neighbors_bynumber);
    m.def("get_all_neighbors_sann", &get_all_neighbors_sann);
    m.def("get_all_neighbors_adaptive", &get_all_neighbors_adaptive);
    m.def("calculate_q", &calculate_q);
    m.def("calculate_q_atom", &calculate_q_atom);
    m.def("calculate_q_single", &calculate_q_single);
    m.def("calculate_aq_single", &calculate_aq_single);
    m.def("calculate_disorder", &calculate_disorder);
    m.def("calculate_bonds", &calculate_bonds);
    m.def("find_clusters", &find_clusters);

#ifdef VERSION_INFO
    m.attr("__version__") = VERSION_INFO;
#else
    m.attr("__version__") = "dev";
#endif
}
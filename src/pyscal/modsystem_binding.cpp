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

PYBIND11_MAKE_OPAQUE(std::vector<Atom, std::allocator<Atom>>);
using AtomList = std::vector<Atom, std::allocator<Atom>>;

PYBIND11_MODULE(cmodsystem, m) {
    py::options options;
    options.disable_function_signatures();
    m.def("get_abs_distance", &get_abs_distance);
    m.def("reset_all_neighbors", &reset_all_neighbors);
    m.def("get_all_neighbors_normal", &get_all_neighbors_normal);
    m.def("get_all_neighbors_cells", &get_all_neighbors_cells);

#ifdef VERSION_INFO
    m.attr("__version__") = VERSION_INFO;
#else
    m.attr("__version__") = "dev";
#endif
}
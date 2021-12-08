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
    m.def("perturb_atom", &perturb_atom);

    //py::class_<AtomList>(m, "AtomList")
    //    .def(py::init<>())
    //    .def("pop_back", &AtomList::pop_back)
    //    /* There are multiple versions of push_back(), etc. Select the right ones. */
    //    .def("push_back", (void (AtomList::*)(Atom &)) &AtomList::push_back)
    //    .def("back", (Atom &(AtomList::*)()) &AtomList::back)
    //    .def("__len__", [](const AtomList &v) { return v.size(); })
    //    .def("__iter__", [](AtomList &v) {
    //       return py::make_iterator(v.begin(), v.end());
    //    }, py::keep_alive<0, 1>());

#ifdef VERSION_INFO
    m.attr("__version__") = VERSION_INFO;
#else
    m.attr("__version__") = "dev";
#endif
}
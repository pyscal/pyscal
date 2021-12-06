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

PYBIND11_MODULE(cmodsystem, m) {
    py::options options;
    options.disable_function_signatures();
    m.def("get_abs_distance", &get_abs_distance);

#ifdef VERSION_INFO
    m.attr("__version__") = VERSION_INFO;
#else
    m.attr("__version__") = "dev";
#endif
}
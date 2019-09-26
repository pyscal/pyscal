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
#include "atom.h"


namespace py = pybind11;
using namespace std;

class PyAtom : public Atom {
public:
  using Atom::Atom;
};

vector<double> vv{0,0,0};

PYBIND11_MODULE(catom, m) {
    py::options options;
    options.disable_function_signatures();

//bindings for Atom class
//------------------------------------------------------------------
    py::class_<Atom, PyAtom>(m,"Atom")
        .def(py::init < vector<double>, int , int >(), py::arg("pos")=vv, py::arg("id")=0, py::arg("type")=0)
        .def_property("pos", &Atom::gx, &Atom::sx)
        .def_property("cluster",&Atom::gcluster, &Atom::scluster)
        .def_property("neighbors",&Atom::gneighbors, &Atom::sneighbors)
        .def("get_coordination",&Atom::gnneighbors)
        .def_property("neighborweights",&Atom::gneighborweights, &Atom::sneighborweights)
        .def_property("bonds",&Atom::gfrenkelnumber, &Atom::sfrenkelnumber)
        .def_property("id",&Atom::gid,&Atom::sid)
        .def_property("condition",&Atom::gcondition, &Atom::scondition)
        .def_property("solid",&Atom::gsolid,&Atom::ssolid)
        .def_property("structure",&Atom::gstructure,&Atom::sstructure)
        .def_property("loc",&Atom::gloc,&Atom::sloc)
        .def_property("allq",&Atom::gallq,&Atom::sallq)
        .def_property("allaq",&Atom::gallaq,&Atom::sallaq)
        .def_property("type",&Atom::gtype,&Atom::stype)
        .def("get_q",&Atom::gq, R"mydelimiter(
            Test documentation

            Parameters
            ----------
            None
        )mydelimiter")
        .def("set_q",&Atom::sq)
        .def("get_aq",&Atom::gaq)
        .def("set_aq",&Atom::saq)
        .def("get_qlm",&Atom::gqlm)
        .def("get_aqlm",&Atom::gaqlm)
        .def_property("volume",&Atom::gvolume,&Atom::svolume)
        .def_property("avg_volume",&Atom::gavgvolume,&Atom::savgvolume)
        .def_property("face_vertices", &Atom::gfacevertices, &Atom::sfacevertices)
        .def_property("face_perimeters", &Atom::gfaceperimeters, &Atom::sfaceperimeters)
        .def_property("vertex_numbers", &Atom::gvertex_numbers, &Atom::svertex_numbers)
        .def_property("vertex_vectors", &Atom::gvertex_vectors, &Atom::svertex_vectors)
        .def_property("avg_connection", &Atom::gasij, &Atom::sasij)
        ;


#ifdef VERSION_INFO
    m.attr("__version__") = VERSION_INFO;
#else
    m.attr("__version__") = "dev";
#endif
}

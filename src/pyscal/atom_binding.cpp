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

        .def_property("pos", &Atom::gx, &Atom::sx, R"mydelimiter(
            List of floats of the type [x, y, z], default [0, 0, 0]
            Position of the atom.
        )mydelimiter")

        .def_property("cluster",&Atom::gcluster, &Atom::scluster, R"mydelimiter(
            int
            Id number of the cluster that the atom belongs to.
        )mydelimiter")

        .def_property("neighbors",&Atom::gneighbors, &Atom::sneighbors, R"mydelimiter(
            List of ints
            List of neighbors of the atom. The list contains indices of neighbor
            atoms which indicate their position in the list of all atoms.
        )mydelimiter")

        .def("get_coordination",&Atom::gnneighbors, R"mydelimiter(
            Get the coordination number of the atom

            Parameters
            ----------
            None

            Returns
            -------
            coordination : int
                coordination number of the atom.

        )mydelimiter")

        .def_property("neighborweights",&Atom::gneighborweights, &Atom::sneighborweights, R"mydelimiter(
            List of floats
            Used to weight the contribution of each neighbor atom towards the value of
            Steinhardt's parameters. By default, each atom has a weight of 1 each. However,
            if :func:`~pyscal.core.System.find_neighbors` is used with `method='voronoi'`,
            each neighbor gets a weight proportional to the area shared between the neighboring
            atom and host atom.
        )mydelimiter")

        .def_property("bonds",&Atom::gfrenkelnumber, &Atom::sfrenkelnumber, R"mydelimiter(
            Int
            The number of solid connections of an atom.
        )mydelimiter")
        .def_property("id",&Atom::gid,&Atom::sid, R"mydelimiter(

        )mydelimiter")

        .def_property("condition",&Atom::gcondition, &Atom::scondition, R"mydelimiter(

        )mydelimiter")
        .def_property("solid",&Atom::gsolid,&Atom::ssolid, R"mydelimiter(

        )mydelimiter")
        .def_property("surface",&Atom::gsurface,&Atom::ssurface, R"mydelimiter(

        )mydelimiter")
        .def_property("largest_cluster",&Atom::glcluster,&Atom::slcluster, R"mydelimiter(

        )mydelimiter")
        .def_property("structure",&Atom::gstructure,&Atom::sstructure, R"mydelimiter(

        )mydelimiter")
        .def_property("loc",&Atom::gloc,&Atom::sloc, R"mydelimiter(

        )mydelimiter")
        .def_property("allq",&Atom::gallq,&Atom::sallq, R"mydelimiter(

        )mydelimiter")
        .def_property("allaq",&Atom::gallaq,&Atom::sallaq, R"mydelimiter(

        )mydelimiter")
        .def_property("type",&Atom::gtype,&Atom::stype, R"mydelimiter(

        )mydelimiter")
        .def_property("custom",&Atom::gcustom,&Atom::scustom, R"mydelimiter(

        )mydelimiter")
        .def("get_q",&Atom::gq, R"mydelimiter(
            Test documentation

            Parameters
            ----------
            None
        )mydelimiter")
        .def("set_q",&Atom::sq, R"mydelimiter(

        )mydelimiter")
        .def("get_aq",&Atom::gaq, R"mydelimiter(

        )mydelimiter")
        .def("set_aq",&Atom::saq, R"mydelimiter(

        )mydelimiter")
        .def("get_qlm",&Atom::gqlm, R"mydelimiter(

        )mydelimiter")
        .def("get_aqlm",&Atom::gaqlm, R"mydelimiter(

        )mydelimiter")
        .def_property("volume",&Atom::gvolume,&Atom::svolume, R"mydelimiter(

        )mydelimiter")
        .def_property("avg_volume",&Atom::gavgvolume,&Atom::savgvolume, R"mydelimiter(

        )mydelimiter")
        .def_property("face_vertices", &Atom::gfacevertices, &Atom::sfacevertices, R"mydelimiter(

        )mydelimiter")
        .def_property("face_perimeters", &Atom::gfaceperimeters, &Atom::sfaceperimeters, R"mydelimiter(

        )mydelimiter")
        .def_property("vertex_numbers", &Atom::gvertex_numbers, &Atom::svertex_numbers, R"mydelimiter(

        )mydelimiter")
        .def_property("vertex_vectors", &Atom::gvertex_vectors, &Atom::svertex_vectors, R"mydelimiter(

        )mydelimiter")
        .def_property("edge_lengths", &Atom::gedgelengths, &Atom::sedgelengths, R"mydelimiter(

        )mydelimiter")
        .def_property("vorovector", &Atom::gvorovector, &Atom::svorovector, R"mydelimiter(

        )mydelimiter")

        .def_property("avg_connection", &Atom::gasij, &Atom::sasij, R"mydelimiter(

        )mydelimiter")
        .def("get_q", (double (Atom::*) (bool, int))  &Atom::gq_big, py::arg("averaged"), py::arg("q"), R"mydelimiter(

        )mydelimiter")
        .def("get_q", (vector<double> (Atom::*) (bool, vector<int>))  &Atom::gq_big, py::arg("averaged"), py::arg("q"), R"mydelimiter(

        )mydelimiter")
        .def("set_q", (void (Atom::*) (bool, int, double))  &Atom::sq_big, py::arg("averaged"), py::arg("q"), py::arg("vals"), R"mydelimiter(

        )mydelimiter")
        .def("set_q", (void (Atom::*) (bool, vector<int>, vector<double>))  &Atom::sq_big, py::arg("averaged"), py::arg("q"), py::arg("vals"), R"mydelimiter(

        )mydelimiter")

        ;


#ifdef VERSION_INFO
    m.attr("__version__") = VERSION_INFO;
#else
    m.attr("__version__") = "dev";
#endif
}

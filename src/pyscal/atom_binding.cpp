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


vector<double> vv{0,0,0};

PYBIND11_MODULE(catom, m) {
    py::options options;
    options.disable_function_signatures();

//bindings for Atom class
//------------------------------------------------------------------
py::class_<Atom>(m,"Atom", R"mydelimiter(
        Class to store atom details.

        Parameters
        ----------
        pos : list of floats of length 3
            position of the `Atom`, default [0,0,0]

        id : int
            id of the `Atom`, default 0

        type : int
            type of the `Atom`, default 1

        Notes
        -----
        A pybind11 class for holding the properties of a single atom. Various properties of the atom
        can be accessed through  the attributes and member functions which are described below in detail. Atoms can
        be created individually or directly by reading a file. Check the examples for more
        details on how atoms are created. For creating atoms directly from an input file check
        the documentation of :class:`~pyscal.core.System` class.

        Although an `Atom` object can be created independently, `Atom` should be thought of
        inherently as members of the :class:`~pyscal.core.System` class. All the properties that define an atom are
        relative to the parent class. :class:`~pyscal.core.System` has a list of all atoms. All the properties of an
        atom, hence should be calculated through :class:`~pyscal.core.System`.

        Examples
        --------
        >>> #method 1 - individually
        >>> atom = Atom()
        >>> #now set positions of the atoms
        >>> atom.pos = [23.0, 45.2, 34.2]
        >>> #now set id
        >>> atom.id = 23
        >>> #now set type
        >>> atom.type = 1
        >>> #Setting through constructor
        >>> atom = Atom([23.0, 45.2, 34.2], 23, 1)

        References
        ----------
        `Creation of atoms <https://pyscal.readthedocs.io/en/latest/examples.html>`_.

    )mydelimiter")

    //-------------------------------------------------------
    // Constructor, Destructor
    //-------------------------------------------------------
    .def(py::init < vector<double>, int , int >(), py::arg("pos")=vv, py::arg("id")=0, py::arg("type")=0)

    //-------------------------------------------------------
    // Basic Atom properties
    //-------------------------------------------------------
    .def_readwrite("pos", &Atom::pos, R"mydelimiter(
        *List of floats of the type [x, y, z], default [0, 0, 0]*.
        Position of the atom.
    )mydelimiter")

    .def_readwrite("id", &Atom::id, R"mydelimiter(
        *int*.
        Id of the atom.
    )mydelimiter")

    .def_readwrite("condition", &Atom::condition, R"mydelimiter(
        *int*.
        condition that specifies if an atom is included in the clustering algorithm or not.
        Only atoms with the value of condition=1 will be used for clustering in
        :func:`~pyscal.core.System.cluster_atoms`.
    )mydelimiter")

    .def_readwrite("mask", &Atom::mask, R"mydelimiter(
        *bool*.
        Mask variable for atom. If mask is true, the atom is ignored from calculations.
    )mydelimiter")

    .def_readwrite("loc", &Atom::loc, R"mydelimiter(
        *int*.
        indicates the position of the atom in the list of all atoms.
    )mydelimiter")

    .def_readwrite("type", &Atom::type, R"mydelimiter(
        *int*.
        int specifying type of the atom.
    )mydelimiter")

    .def_readwrite("ghost", &Atom::ghost, R"mydelimiter(
        *int*.
        int specifying ghost status of the atom.
    )mydelimiter")

    .def_readwrite("custom", &Atom::custom, R"mydelimiter(
        *dict*.
        dictionary specfying custom values for an atom. The module only stores the id, type and
        position of the atom. If any extra values need to be stored, they can be stored in custom
        using `atom.custom = {"velocity":12}`. :func:`~pyscal.core.System.read_inputfile` can also
        read in extra atom information. By default, custom values are treated as string.
    )mydelimiter")

    //-------------------------------------------------------
    // Neighbor related properties
    //-------------------------------------------------------
    .def_readwrite("neighbors",&Atom::neighbors, R"mydelimiter(
        *List of ints*.
        List of neighbors of the atom. The list contains indices of neighbor
        atoms which indicate their position in the list of all atoms.
    )mydelimiter")

    .def_readwrite("neighbor_distance", &Atom::neighbordist, R"mydelimiter(
        *List of floats*.
        List of neighbor distances of the atom.
    )mydelimiter")

    .def_readwrite("coordination", &Atom::n_neighbors, R"mydelimiter(
        *int*.
        coordination number of the atom. Coordination will only be updated
        after neighbors are calculated using :func:`~pyscal.core.System.find_neighbors`.
    )mydelimiter")

    .def_readwrite("neighbor_weights",&Atom::neighborweight, R"mydelimiter(
        *List of floats*.
        Used to weight the contribution of each neighbor atom towards the value of
        Steinhardt's parameters. By default, each atom has a weight of 1 each. However,
        if :func:`~pyscal.core.System.find_neighbors` is used with `method='voronoi'`,
        each neighbor gets a weight proportional to the area shared between the neighboring
        atom and host atom.
    )mydelimiter")

   .def_readwrite("cutoff", &Atom::cutoff, R"mydelimiter(
        *double*.
        cutoff used for finding neighbors for each atom.
    )mydelimiter")
   
    ;


#ifdef VERSION_INFO
    m.attr("__version__") = VERSION_INFO;
#else
    m.attr("__version__") = "dev";
#endif
}

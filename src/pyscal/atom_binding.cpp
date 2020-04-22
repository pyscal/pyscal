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

    .def(py::init < vector<double>, int , int >(), py::arg("pos")=vv, py::arg("id")=0, py::arg("type")=0)

    .def_property("pos", &Atom::gx, &Atom::sx, R"mydelimiter(
        *List of floats of the type [x, y, z], default [0, 0, 0]*.
        Position of the atom.
    )mydelimiter")

    .def_property("cluster",&Atom::gcluster, &Atom::scluster, R"mydelimiter(
        *int*.
        identification number of the cluster that the atom belongs to.
    )mydelimiter")

    .def_property("neighbors",&Atom::gneighbors, &Atom::sneighbors, R"mydelimiter(
        *List of ints*.
        List of neighbors of the atom. The list contains indices of neighbor
        atoms which indicate their position in the list of all atoms.
    )mydelimiter")

    .def_property("neighbor_distance",&Atom::gneighdist, &Atom::sneighdist, R"mydelimiter(
        *List of floats*.
        List of neighbor distances of the atom.
    )mydelimiter")

    .def_property("coordination",&Atom::gnneighbors, &Atom::snneighbors,R"mydelimiter(
        *int*.
        coordination number of the atom. Coordination will only be updated
        after neighbors are calculated using :func:`~pyscal.core.System.find_neighbors`.
    )mydelimiter")

    .def_property("neighbor_weights",&Atom::gneighborweights, &Atom::sneighborweights, R"mydelimiter(
        *List of floats*.
        Used to weight the contribution of each neighbor atom towards the value of
        Steinhardt's parameters. By default, each atom has a weight of 1 each. However,
        if :func:`~pyscal.core.System.find_neighbors` is used with `method='voronoi'`,
        each neighbor gets a weight proportional to the area shared between the neighboring
        atom and host atom.
    )mydelimiter")

    .def_property("bonds",&Atom::gfrenkelnumber, &Atom::sfrenkelnumber, R"mydelimiter(
        *Int*.
        The number of solid bonds of an atom.
    )mydelimiter")

    .def_property("id",&Atom::gid,&Atom::sid, R"mydelimiter(
        *int*.
        Id of the atom.
    )mydelimiter")

    .def_property("condition",&Atom::gcondition, &Atom::scondition, R"mydelimiter(
        *int*.
        condition that specifies if an atom is included in the clustering algorithm or not.
        Only atoms with the value of condition=1 will be used for clustering in
        :func:`~pyscal.core.System.cluster_atoms`.
    )mydelimiter")

    .def_property("solid",&Atom::gsolid,&Atom::ssolid, R"mydelimiter(
        *bool*.
        True if the atom is solid, False otherwise. Solid atoms are only identified
        after using the :func:`~pyscal.core.System.find_solids` function.
    )mydelimiter")

    .def_property("surface",&Atom::gsurface,&Atom::ssurface, R"mydelimiter(
        *bool*.
        True if the atom has at least one liquid neighbor, False otherwise. Surface atoms are only identified
        after using the :func:`~pyscal.core.System.find_solids` function.

    )mydelimiter")

    .def_property("largest_cluster",&Atom::glcluster,&Atom::slcluster, R"mydelimiter(
        *bool*.
        True if the atom belongs to the largest cluster, False otherwise. Largest cluster is only identified
        after using the :func:`~pyscal.core.System.cluster_atoms` function.
    )mydelimiter")

    .def_property("structure",&Atom::gstructure,&Atom::sstructure, R"mydelimiter(
        *int*.
        Indicates the structure of atom. Not used currently.
    )mydelimiter")

    .def_property("mask",&Atom::gmask, &Atom::smask, R"mydelimiter(
        *bool*.
        Mask variable for atom. If mask is true, the atom is ignored from calculations.
    )mydelimiter")

    .def_property("loc",&Atom::gloc,&Atom::sloc, R"mydelimiter(
        *int*.
        indicates the position of the atom in the list of all atoms.
    )mydelimiter")

    .def_property("allq",&Atom::gallq,&Atom::sallq, R"mydelimiter(
        *list of floats*.
        list of all q values of the atom.
    )mydelimiter")

    .def_property("allaq",&Atom::gallaq,&Atom::sallaq, R"mydelimiter(
        *list of floats*.
        list of all averaged q values of the atom.
    )mydelimiter")

    .def_property("type",&Atom::gtype,&Atom::stype, R"mydelimiter(
        *int*.
        int specifying type of the atom.
    )mydelimiter")

    .def_property("cutoff",&Atom::gcutoff,&Atom::scutoff, R"mydelimiter(
        *double*.
        cutoff used for finding neighbors for each atom.
    )mydelimiter")

    .def_property("custom",&Atom::gcustom,&Atom::scustom, R"mydelimiter(
        *dict*.
        dictionary specfying custom values for an atom. The module only stores the id, type and
        position of the atom. If any extra values need to be stored, they can be stored in custom
        using `atom.custom = {"velocity":12}`. :func:`~pyscal.core.System.read_inputfile` can also
        read in extra atom information. By default, custom values are treated as string.
    )mydelimiter")

    .def_property("volume",&Atom::gvolume,&Atom::svolume, R"mydelimiter(
        *float*. Voronoi volume of the atom. The Voronoi volume is only calculated if neighbors
        are found using the :func:`~pyscal.core.System.find_neighbors` using the `method='voronoi'`
        option.

    )mydelimiter")
    .def_property("avg_volume",&Atom::gavgvolume,&Atom::savgvolume, R"mydelimiter(
        *float*. Averaged version of the Voronoi volume which is calculated as an average over
        itself and its neighbors. Only calculated when the :func:`~pyscal.core.System.find_neighbors`
        using the `method='voronoi'` option is used.
    )mydelimiter")

    .def_property("face_vertices", &Atom::gfacevertices, &Atom::sfacevertices, R"mydelimiter(
          *list of floats*. A list of the number of vertices shared between an atom and its
          neighbors.  Only calculated when the :func:`~pyscal.core.System.find_neighbors`
          using the `method='voronoi'` option is used.
    )mydelimiter")

    .def_property("face_perimeters", &Atom::gfaceperimeters, &Atom::sfaceperimeters, R"mydelimiter(
          *list of floats*. List consisting of the perimeters of each Voronoi face of an atom.
          Only calculated when the :func:`~pyscal.core.System.find_neighbors`
          using the `method='voronoi'` option is used.
    )mydelimiter")

    .def_property("vertex_numbers", &Atom::gvertex_numbers, &Atom::svertex_numbers, R"mydelimiter(
          *list of floats*. For each Voronoi face of the atom, this values includes a List
          of vertices that constitute the face.  Only calculated when the :func:`~pyscal.core.System.find_neighbors`
          using the `method='voronoi'` option is used.
    )mydelimiter")

    .def_property("vertex_vectors", &Atom::gvertex_vectors, &Atom::svertex_vectors, R"mydelimiter(
          *list of floats*. A list of positions of each vertex of the Voronoi polyhedra of
          the atom.  Only calculated when the :func:`~pyscal.core.System.find_neighbors`
          using the `method='voronoi'` option is used.
    )mydelimiter")

    .def_property("edge_lengths", &Atom::gedgelengths, &Atom::sedgelengths, R"mydelimiter(
          *list of floats*. For each face, this vector contains the lengths of edges
          that make up the Voronoi polyhedra of the atom.  Only calculated when the :func:`~pyscal.core.System.find_neighbors`
          using the `method='voronoi'` option is used.
    )mydelimiter")

    .def_property("vorovector", &Atom::gvorovector, &Atom::svorovector, R"mydelimiter(
          *list of ints*. A vector of the form `(n3, n4, n5, n6)` where n3 is the number of faces with 3 vertices,
          n4 is the number of faces with 4 vertices and so on. This can be used to identify structures [1][2].
          Vorovector is calculated if the :func:`~pyscal.core.System.calculate_vorovector` method is used.

          References
          ----------
          .. [1] Finney, JL, Proc. Royal Soc. Lond. A 319, 1970
          .. [2] Tanemura, M, Hiwatari, Y, Matsuda, H,Ogawa, T, Ogita, N, Ueda, A. Prog. Theor. Phys. 58, 1977

    )mydelimiter")

    .def_property("sij", &Atom::gsij, &Atom::ssij, R"mydelimiter(
          *float*. Value of s_ij which is used for identification of solid atoms. s_ij is defined by

          .. math:: s_{ij} = \sum_{m=-l}^l q_{lm}(i) q_{lm}^*(i)

    )mydelimiter")

    .def_property("avg_sij", &Atom::gasij, &Atom::sasij, R"mydelimiter(
          *float*. Value of averaged s_ij which is used for identification of solid atoms. s_ij is defined by

          .. math:: s_{ij} = \sum_{m=-l}^l q_{lm}(i) q_{lm}^*(i)

    )mydelimiter")

    .def("get_q", (double (Atom::*) (int q, bool))  &Atom::gq_big,  py::arg(), py::arg("averaged")=false, R"mydelimiter(
          Calculate the steinhardt parameter q_l value.

          Parameters
          ----------
          q : int or list of ints
              number of the required q_l - from 2-12

          averaged : bool, optional
              If True, return the averaged q values,
              If False, return the non averaged ones
              default False

          Returns
          -------
          q_l : float or list of floats
              the value(s) of the queried Steinhardt parameter(s).

          Notes
          -----
          Please check this `link <https://pyscal.readthedocs.io/en/latest/steinhardtparameters.html>`_
          for more details about Steinhardts parameters and the averaged versions.

          Meaningful values are only returned if :func:`~pyscal.core.System.calculate_q` is used.
    )mydelimiter")

    .def("get_q", (vector<double> (Atom::*) (vector<int>, bool))  &Atom::gq_big, py::arg(), py::arg("averaged")=false )

    .def("set_q", (void (Atom::*) (int, double, bool))  &Atom::sq_big, py::arg(), py::arg(), py::arg("averaged")=false, R"mydelimiter(

          Set the value of steinhardt parameter q_l.

          Parameters
          ----------
          q : int or list of ints
              number of the required q_l - from 2-12

          val : float or list of floats
              value(s) of Steinhardt parameter(s).

          averaged : bool, optional
              If True, return the averaged q values,
              If False, return the non averaged ones
              default False

          Returns
          -------
          None

    )mydelimiter")

    .def("set_q", (void (Atom::*) (vector<int>, vector<double>, bool))  &Atom::sq_big, py::arg(), py::arg(), py::arg("averaged")=false)


    .def_property("angular",&Atom::gangular, &Atom::sangular, R"mydelimiter(
        *Float*.
        The value of angular parameter A of an atom. The angular parameter measures the tetrahedral coordination of an atom.
        Meaningful values are only returned if the property is calculated using :func:`~pyscal.core.System.calculate_angularcriteria`.
    )mydelimiter")

    .def_property("chiparams",&Atom::gchiparams, &Atom::schiparams, R"mydelimiter(
        *Float*.
        The value of chiparameter of an atom. The return value is a vector of length 8.
        Meaningful values are only returned if chi params are calculated using :func:`~pyscal.core.System.calculate_chiparams`.
    )mydelimiter")

    .def_property("avg_angular",&Atom::gavgangular, &Atom::savgangular, R"mydelimiter(
        *Float*.
        The average angular parameter value. Not used currently.
    )mydelimiter")

    .def_property("disorder",&Atom::gdisorder, &Atom::sdisorder, R"mydelimiter(
        *Float*.
        The value of disorder parameter.
    )mydelimiter")

    .def_property("avg_disorder",&Atom::gavgdisorder, &Atom::savgdisorder, R"mydelimiter(
        *Float*.
        The value of averaged disorder parameter.
    )mydelimiter")

    .def_property("sro",&Atom::gsro, &Atom::ssro, R"mydelimiter(
        *Float*.
        The value of short range order parameter.
    )mydelimiter")

    .def("get_qlm", &Atom::get_qcomps, py::arg(), py::arg("averaged")=false, R"mydelimiter(
          Get the q_lm values.

          Parameters
          ----------
          q : int
              number of the required q_l - from 2-12

          averaged : bool, optional
              If True, return the averaged qlm values,
              If False, return the non averaged ones
              default False

          Returns
          -------
          q_lm : complex vector
              vector of complex numbers.

          Meaningful values are only returned if :func:`~pyscal.core.System.calculate_q` is used.
    )mydelimiter")

    ;


#ifdef VERSION_INFO
    m.attr("__version__") = VERSION_INFO;
#else
    m.attr("__version__") = "dev";
#endif
}

.. pyscal documentation master file, created by
   sphinx-quickstart on Wed Apr 24 14:11:20 2019.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

pyscal- A python module for structural analysis of atomic environments
======================================================================

.. image:: https://dev.azure.com/sarathrmenon/pyscal/_apis/build/status/srmnitc.pyscal?branchName=master
    :target: https://dev.azure.com/sarathrmenon/pyscal/_build/latest?definitionId=1&branchName=master
    :width: 20%

.. image:: https://codecov.io/gh/srmnitc/pyscal/branch/master/graph/badge.svg
  :target: https://codecov.io/gh/srmnitc/pyscal
  :width: 15 %

.. image:: https://mybinder.org/badge_logo.svg
   :target: https://mybinder.org/v2/gh/srmnitc/pybop/master?filepath=examples%2F
   :width: 15 %

.. image:: https://anaconda.org/pyscal/pyscal/badges/installer/conda.svg
   :target: https://anaconda.org/conda-forge/pyscal
   :width: 13 %

.. image:: https://img.shields.io/conda/dn/conda-forge/pyscal.svg
   :target: https://conda.anaconda.org/pyscal
   :width: 13 %

.. image:: https://joss.theoj.org/papers/168eca482155601dc517523899527a4e/status.svg
   :target: https://joss.theoj.org/papers/168eca482155601dc517523899527a4e
   :width: 20 %

.. image:: https://img.shields.io/conda/pn/conda-forge/pyscal.svg
  :target: https://anaconda.org/conda-forge/pyscal
  :width: 20 %

**pyscal** is a python module for the calculation of local atomic structural environments including Steinhardt's bond orientational order parameters [1]_ during post-processing
of atomistic simulation data. The core functionality of pyscal is written in C++ with python wrappers using
`pybind11 <https://pybind11.readthedocs.io/en/stable/intro.html>`_  which allows for fast calculations and
easy extensions in python.

Steinhardt's order parameters are widely used for the identification of crystal structures [3]_. They are also used to distinguish
if an atom is in a solid or liquid environment [4]_. pyscal is inspired by the
`BondOrderAnalysis <https://homepage.univie.ac.at/wolfgang.lechner/bondorderparameter.html>`_ code,
but has since incorporated many additional features and modifications. The pyscal module includes the following functionalities:

Highlights
----------

* calculation of Steinhardt's order parameters and their averaged version [2]_.
* links with the `Voro++ <http://math.lbl.gov/voro++/>`_ code, for the calculation of Steinhardt parameters weighted using the face areas of Voronoi polyhedra [3]_.
* classification of atoms as solid or liquid [4]_.
* clustering of particles based on a user defined property.
* methods for calculating radial distribution functions, Voronoi volumes of particles, number of vertices and face area of Voronoi polyhedra, and coordination numbers.
* calculation of angular parameters to identify diamond structure [5]_.

.. [1]  `Steinhardt, P. J., Nelson, D. R., & Ronchetti, M. (1983). Physical Review B, 28 <https://journals.aps.org/prb/abstract/10.1103/PhysRevB.28.784>`_.
.. [2]  `Lechner, W., & Dellago, C. (2008). The Journal of Chemical Physics, 129 <https://aip.scitation.org/doi/full/10.1063/1.2977970>`_.
.. [3]  `Mickel, W., Kapfer, S. C., Schröder-Turk, G. E., & Mecke, K. (2013). The Journal of Chemical Physics, 138 <https://aip.scitation.org/doi/full/10.1063/1.4774084>`_.
.. [4]  `Auer, S., & Frenkel, D. (2005). Advances in Polymer Science, 173 <https://link.springer.com/chapter/10.1007/b99429>`_.
.. [5]  `Uttormark, M. J., Thompson, M. O., Clancy, P. (1993). Physical Review B, 47 <https://journals.aps.org/prb/abstract/10.1103/PhysRevB.47.15717>`_.

Getting started
---------------

.. toctree::

    gettingstarted

Download
--------

.. toctree::

    download

Examples
--------

.. toctree::

    examples


Methods implemented in pyscal
-----------------------------

.. toctree::

    methods

News and updates
----------------

.. toctree::

    news

pyscal reference
----------------

.. toctree::

    modules


Support, contributing and extending
-----------------------------------

.. toctree::

    extending
    helpandsupport
    common_issues

Credits
-------

.. toctree::

    citing
    acknowledgements

License
-------

.. toctree::

    license

.. pyscal documentation master file, created by
   sphinx-quickstart on Wed Apr 24 14:11:20 2019.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

pyscal- python Structural Environment Calculator
================================================

.. image:: https://travis-ci.com/srmnitc/pyscal.svg?branch=master
    :target: https://travis-ci.com/srmnitc/pyscal
    :width: 13 %

.. image:: https://codecov.io/gh/srmnitc/pyscal/branch/master/graph/badge.svg
  :target: https://codecov.io/gh/srmnitc/pyscal
  :width: 15 %

.. image:: https://mybinder.org/badge_logo.svg
   :target: https://mybinder.org/v2/gh/srmnitc/pybop/master?filepath=examples%2F
   :width: 15 %

.. image:: https://img.shields.io/badge/License-GPLv3-blue.svg
   :target: https://www.gnu.org/licenses/gpl-3.0
   :width: 13 %

**pyscal** is a python module for the calculation of Steinhardt's bond orientational order parameters [1]_ during post-processing
of atomistic simulation data. The core functionality of pyscal is written in C++ with python wrappers using
`pybind11 <https://pybind11.readthedocs.io/en/stable/intro.html>`_  which allows for fast calculations, and
easy extension in python.

Steinhardt's order parameters are widely used for identification of crystal structures [3]_. They are also used to identify
if an atom is solid or liquid [4]_. pyscal is inspired by
`BondOrderAnalysis <https://homepage.univie.ac.at/wolfgang.lechner/bondorderparameter.html>`_ code,
but has since incorporated many additions and modifications. pyscal module includes the following functionality-

Highlights
----------

* calculation of Steinhardt's order parameters and their averaged version [2]_.
* links with `Voro++ <http://math.lbl.gov/voro++/>`_ code, for calculation of Steinhardt parameters weighted using face area of Voronoi polyhedra [3]_.
* classification of atoms as solid or liquid [4]_.
* clustering of particles based on a user defined property.
* methods for calculating radial distribution function, voronoi volume of particles, number of vertices and face area of voronoi polyhedra and coordination number.

.. [1]  `Steinhardt, PJ, Nelson, DR, Ronchetti, M. Phys Rev B 28, 1983 <https://journals.aps.org/prb/abstract/10.1103/PhysRevB.28.784>`_.
.. [2]  `Lechner, W, Dellago, C, J Chem Phys, 2013 <https://aip.scitation.org/doi/full/10.1063/1.2977970>`_.
.. [3]  `Mickel, W, Kapfer, SC, Schroder-Turk, GE, Mecke, K, J Chem Phys 138, 2013 <https://aip.scitation.org/doi/full/10.1063/1.4774084>`_.
.. [4]  `Auer, S, Frenkel, D. Adv Polym Sci 173, 2005 <https://link.springer.com/chapter/10.1007/b99429>`_.

Getting started
---------------

.. toctree::

    gettingstarted


Examples
--------

.. toctree::

    examples

Download
--------

.. toctree::

    download

Methods
-------

.. toctree::

    methods


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

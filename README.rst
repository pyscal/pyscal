.. image:: https://dev.azure.com/sarathrmenon/pyscal/_apis/build/status/srmnitc.pyscal?branchName=master
    :target: https://dev.azure.com/sarathrmenon/pyscal/_build/latest?definitionId=1&branchName=master
    :width: 13%

.. image:: https://codecov.io/gh/srmnitc/pyscal/branch/master/graph/badge.svg
  :target: https://codecov.io/gh/srmnitc/pyscal
  :width: 13 %

.. image:: https://readthedocs.org/projects/pyscal/badge/?version=latest
    :target: https://pyscal.readthedocs.io/en/latest/?badge=latest
    :width: 50%

.. image:: https://mybinder.org/badge_logo.svg
   :target: https://mybinder.org/v2/gh/srmnitc/pybop/master?filepath=examples%2F
   :width: 13 %

.. image:: https://img.shields.io/badge/License-GPLv3-blue.svg
   :target: https://www.gnu.org/licenses/gpl-3.0
   :width: 12 %

.. image:: https://anaconda.org/pyscal/pyscal/badges/installer/conda.svg
   :target: https://conda.anaconda.org/pyscal
   :width: 13 %

.. image:: https://joss.theoj.org/papers/168eca482155601dc517523899527a4e/status.svg
   :target: https://joss.theoj.org/papers/168eca482155601dc517523899527a4e
   :width: 20 %

.. image:: https://anaconda.org/pyscal/pyscal/badges/installer/conda.svg
  :target: https://anaconda.org/conda-forge/pyscal
  :width: 13 %

pyscal - python Structural Environment Calculator
=================================================

complete documentation with examples available `here <https://pyscal.com/>`_.

**pyscal** is a python module for the calculation of local atomic structural environments including Steinhardt's bond orientational order parameters [1]_ during post-processing
of atomistic simulation data. The core functionality of pyscal is written in C++ with python wrappers using `pybind11 <https://pybind11.readthedocs.io/en/stable/intro.html>`_  which allows for fast calculations with possibilities for easy expansion in python.

Steinhardt's order parameters are widely used for identification of crystal structures [2]_. They are also used to identify if an atom is solid or liquid [3]_. pyscal is inspired by `BondOrderAnalysis <https://homepage.univie.ac.at/wolfgang.lechner/bondorderparameter.html>`_ code, but has since incorporated many additions and modifications. pyscal module includes the following functionality-

Highlights
----------

* fast and efficient calculations using C++ and expansion using python.
* calculation of Steinhardt's order parameters and their averaged version [4]_.
* links with `Voro++ <http://math.lbl.gov/voro++/>`_ code, for calculation of Steinhardt parameters weighted using face area of Voronoi polyhedra [3]_.
* classification of atoms as solid or liquid [4]_.
* clustering of particles based on a user defined property.
* methods for calculating radial distribution function, voronoi volume of particles, number of vertices and face area of voronoi polyhedra and coordination number.

.. [1]  `Steinhardt, P. J., Nelson, D. R., & Ronchetti, M. (1983). Physical Review B, 28 <https://journals.aps.org/prb/abstract/10.1103/PhysRevB.28.784>`_.
.. [2]  `Lechner, W., & Dellago, C. (2008). The Journal of Chemical Physics, 129 <https://aip.scitation.org/doi/full/10.1063/1.2977970>`_.
.. [3]  `Mickel, W., Kapfer, S. C., Schröder-Turk, G. E., & Mecke, K. (2013). The Journal of Chemical Physics, 138 <https://aip.scitation.org/doi/full/10.1063/1.4774084>`_.
.. [4]  `Auer, S., & Frenkel, D. (2005). Advances in Polymer Science, 173 <https://link.springer.com/chapter/10.1007/b99429>`_.

Citing the work
---------------

If you use pyscal in your work, the citation of the `following article <https://joss.theoj.org/papers/10.21105/joss.01824>`_ will be greatly appreciated:

Sarath Menon, Grisell Díaz Leines and Jutta Rogal (2019). pyscal: A python module for structural analysis of atomic environments. Journal of Open Source Software, 4(43), 1824, https://doi.org/10.21105/joss.01824

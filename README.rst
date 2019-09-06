.. image:: https://travis-ci.com/srmnitc/pyscal.svg?branch=master
    :target: https://travis-ci.com/srmnitc/pyscal
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
   :width: 13 %


pyscal - python Structural Environment Calculator 
=================================================

complete documentation with examples available `here <https://pyscal.readthedocs.io/>`_.

**pyscal** is a python module for calculation of local structural environment  including Steinhardt's bond order parameters [#]_. The core functionality of **pyscal** is written in C++ with python wrappers using `pybind11 <https://pybind11.readthedocs.io/en/stable/intro.html>`_ . This allows for fast calculations with possibilities for seamless expansion in python. 

Stenhardt's order parameters or local bond order parameters are widely used in simulations to determine crystal structure of atoms. They are also used to identify if an atom is solid or liquid and the identify the solid atoms in a liquid, which is often used an order parameter to study progress of reaction during solidification. Additionally **pyscal** provides an easy environment for reading in snapshots from molecular dynamics trajectories, allowing for easy post-processing of simulation data. **pyscal** is inspired by `BondOrderAnalysis <https://homepage.univie.ac.at/wolfgang.lechner/bondorderparameter.html>`_ code, but has since incorporated many additions and modifications. The main highlights of **pyscal** are given below.

Highlights
----------

* Fast and efficient library using C++.
* Calculation of Steinhardts order parameters and their averaged version [#]_.
* Links with `Voro++ <http://math.lbl.gov/voro++/>`_ code, to enable weighted calculation of local bond order parameters [#]_.
* An environment for easy processing of simulation data.
* Distinction of solid and liquid atoms using q6 parameter [#]_.
* Cluster algorithms of solid atoms to find largest cluster of solid atoms in liquid [#]_.


.. [#]  Steinhardt, PJ, Nelson, DR, Ronchetti, M. PRB 28, 1983.
.. [#]  Lechner, W, Dellago, C. JCP 129, 2008.
.. [#]  Mickel, W, kapfer, SC, Shroder-Turk, GE, Mecke, K. JCP 138, 2013.
.. [#]  Auer, S, Frenkel, D. APS 173, 2005., Diaz Leines, G, Drautz, R, Rogal, J. JCP 146, 2017.
.. [#]  Diaz Leines, G, Drautz, R, Rogal, J. JCP 146, 2017.

Citing the work
---------------
We are currently preparing a publication. Until it is ready, if you want to use the code in your work, it would be great if you let `me <mailto:sarath.menon@rub.de>`_ know about it. 

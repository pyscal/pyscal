.. pybop documentation master file, created by
   sphinx-quickstart on Wed Apr 24 14:11:20 2019.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

pybop
=================================

**pybop** is a python module for calculation of Steinhardt's bond order parameters [#]_. The core functionality of **pybop** is written in C++ with python wrappers using `pybind11 <https://pybind11.readthedocs.io/en/stable/intro.html>`_ . This allows for fast calculations with possibilities for seamless expansion in python. 

Stenhardt's order parameters or local bond order parameters are widely used in simulations to determine crystal structure of atoms. They are also used to identify if an atom is solid or liquid and the identify the solid atoms in a liquid, which is often used an order parameter to study progress of reaction during solidification. Additionally **pybop** provides an easy environment for reading in snapshots from molecular dynamics trajectories, allowing for easy post-processing of simulation data. **pybop** is inspired by `BondOrderAnalysis <https://homepage.univie.ac.at/wolfgang.lechner/bondorderparameter.html>`_ code, but has since incorporated many additions and modifications. The main highlights of **pybop** are given below.

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

Documentation
-------------

.. toctree::
    
    gettingstarted
    pybop
    examples
    jupyter

Contributing and extending
--------------------------
**pybop** welcomes contribution and extension to the module. Rather than local modifications, we request that the modifications be submitted through a pull request. **pybop** follows the `pep8 <https://www.python.org/dev/peps/pep-0008/>`_ style. The quality of the code can be checked by `pylint <https://www.pylint.org/>`_ library by running ``pylint yourscript.py``. Rather than the style of code, what is more important is the documentation. We request that all improvements are documented in detail.

Help and support
----------------
In case of bugs and feature improvements, you are welcome to create a new issue on the `github repo <https://github.com/srmnitc/pybop>`_. In case of other questions, please contact `us <mailto:sarath.menon@rub.de>`_.

Citing the work
---------------
We are currently preparing a publication. Until it is ready, if you want to use the code in your work, it would be great if you let `me <mailto:sarath.menon@rub.de>`_ know about it. 
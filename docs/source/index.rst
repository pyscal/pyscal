.. pyscal documentation master file, created by
   sphinx-quickstart on Wed Apr 24 14:11:20 2019.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

pyscal- python Structural Environment Calculator 
================================================

.. image:: https://travis-ci.com/srmnitc/pyscal.svg?branch=master
    :target: https://travis-ci.com/srmnitc/pyscal
    :width: 20 %

.. image:: https://codecov.io/gh/srmnitc/pyscal/branch/master/graph/badge.svg
  :target: https://codecov.io/gh/srmnitc/pyscal
  :width: 20 %

.. image:: https://mybinder.org/badge_logo.svg
   :target: https://mybinder.org/v2/gh/srmnitc/pybop/master?filepath=examples%2F
   :width: 20 %

.. image:: https://img.shields.io/badge/License-GPLv3-blue.svg
   :target: https://www.gnu.org/licenses/gpl-3.0
   :width: 20 %

**pyscal** is a python module for calculation of Steinhardt's bond orientational order parameters [Steinhardt1983]_ during post-processing of atomistic simulation data. The core functionality of pyscal is written in C++ with python wrappers using `pybind11 <https://pybind11.readthedocs.io/en/stable/intro.html>`_  which allows for fast calculations with possibilities for easy expansion in python. 

Steinhardt's order parameters are widely used for identification of crystal structures [Mickel2013]_. They are also used to identify if an atom is solid or liquid [Auer2005]_. pyscal is inspired by `BondOrderAnalysis <https://homepage.univie.ac.at/wolfgang.lechner/bondorderparameter.html>`_ code, but has since incorporated many additions and modifications. pyscal module includes the following functionality-  

Highlights
----------

* Fast and efficient calculations using C++ and expansion using python.  
* Calculation of Steinhardt's order parameters and their averaged version [Lechner2008]_.
* Links with `Voro++ <http://math.lbl.gov/voro++/>`_ code, to enable weighted calculation of local bond order parameters [#]_.
* An environment for easy processing of simulation data.
* Distinction of solid and liquid atoms using q6 parameter [#]_.
* Cluster algorithms of solid atoms to find largest cluster of solid atoms in liquid [#]_.

.. [Steinhardt1983]  Steinhardt, PJ, Nelson, DR, Ronchetti, M. PRB 28, 1983.
.. [Lechner2008]  Lechner, W, Dellago, C. JCP 129, 2008.
.. [Mickel2013]  Mickel, W, kapfer, SC, Shroder-Turk, GE, Mecke, K. JCP 138, 2013.
.. [Auer2005]  Auer, S, Frenkel, D. APS 173, 2005.
.. [#]  Diaz Leines, G, Drautz, R, Rogal, J. JCP 146, 2017.

First steps
-----------

.. toctree::
    
    gettingstarted
    examples
    jupyter

Methods
-------

.. toctree::
    
    methods
    

API documentation
-----------------

.. toctree::
    
    pyscal
	
Support, contributing and extending
-----------------------------------

.. toctree::
    
    extending
    helpandsupport
    common_issues
    license
    citing
    acknowledgements


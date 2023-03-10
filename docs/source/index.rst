.. pyscal documentation master file, created by
   sphinx-quickstart on Wed Apr 24 14:11:20 2019.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

pyscal
======

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


**pyscal** is a python module for the calculation of local atomic
structural environments including `Steinhardt's bond orientational order
parameters <https://journals.aps.org/prb/abstract/10.1103/PhysRevB.28.784>`__
during post-processing of atomistic simulation data. The core
functionality of pyscal is written in C++ with python wrappers using
`pybind11 <https://pybind11.readthedocs.io/en/stable/intro.html>`__
which allows for fast calculations with possibilities for easy expansion
in python.

Steinhardt's order parameters are widely used for `identification of
crystal
structures <https://aip.scitation.org/doi/full/10.1063/1.4774084>`__.
They are also used to identify if an atom is `solid or
liquid <https://link.springer.com/chapter/10.1007/b99429>`__. pyscal is
inspired by
`BondOrderAnalysis <https://homepage.univie.ac.at/wolfgang.lechner/bondorderparameter.html>`__
code, but has since incorporated many additions and modifications.


.. toctree::
   :hidden:
   :maxdepth: 2

   Installation <gettingstarted>
   prologue/extending
   prologue/helpandsupport
   prologue/citing
   prologue/acknowledgements
   prologue/license


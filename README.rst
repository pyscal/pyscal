..  image:: https://dev.azure.com/sarathrmenon/pyscal/_apis/build/status/srmnitc.pyscal?branchName=master
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

.. image:: https://anaconda.org/pyscal/pyscal/badges/installer/conda.svg
   :target: https://anaconda.org/conda-forge/pyscal
   :width: 13 %

.. image:: https://joss.theoj.org/papers/168eca482155601dc517523899527a4e/status.svg
   :target: https://joss.theoj.org/papers/168eca482155601dc517523899527a4e
   :width: 20 %

.. image:: https://img.shields.io/conda/dn/conda-forge/pyscal.svg
   :target: https://anaconda.org/conda-forge/pyscal
   :width: 13 %

.. image:: https://img.shields.io/conda/pn/conda-forge/pyscal.svg
   :target: https://anaconda.org/conda-forge/pyscal
   :width: 20 %

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

Installation
------------

**For python 3**

pyscal can be installed directly using `Conda <https://docs.conda.io/en/latest/>`_ from the `conda-forge channel <https://conda-forge.org/>`_ by the following statement-

.. code:: console

    conda install -c conda-forge pyscal

**For python 2**

Python 2 version of pyscal can be installed using-

.. code:: console

    conda install -c pyscal pyscal

**From repository**

pyscal can be built from the repository by-

.. code:: console

    git clone https://github.com/srmnitc/pyscal.git
    cd pyscal
    python setup.py install --user

complete documentation with examples available `here <https://pyscal.com/>`_.

News
----

- November 30, 2019 : pyscal reaches over 1500 downloads on conda-forge channel. Thanks for the support!

- November 21, 2019 : pyscal is selected as the E-CAM module of the month. See the news `here <https://www.e-cam2020.eu/pyscal-a-python-module-for-structural-analysis-of-atomic-environments/>`_.

- November 3, 2019 : `Version 2.1.1 <https://github.com/srmnitc/pyscal/releases/tag/2.2.1>`_ of pyscal is released. This release
  introduces cell lists to speed up calculations for large number of atoms.

- November 1, 2019 : pyscal paper is accepted in the Journal of Open Source Software. See the paper `here <https://joss.theoj.org/papers/10.21105/joss.01824>`_.

- October 29, 2019 : `Version 2.0.1 <https://github.com/srmnitc/pyscal/releases/tag/2.0.1>`_ of pyscal is released.

- October 17, 2019 : Publication for pyscal submitted to the Journal of Open Source Software. See the review `here <https://github.com/openjournals/joss-reviews/issues/1824>`_.

- October 2, 2019 : `Version 2.0.0 <https://github.com/srmnitc/pyscal/releases/tag/2.0.0>`_ of pyscal is released. This version contains
  significant update to the code, increasing the speed and memory footprint.

- August 4, 2019 : `Version 1.0.3 <https://github.com/srmnitc/pyscal/releases/tag/1.0.3>`_ of pyscal is released.

- July 12, 2019 : `Version 1.0.1 <https://github.com/srmnitc/pyscal/releases/tag/v1.0.1>`_ of pyscal is released.

- July 12, 2019 : `Version 1.0.0 <https://github.com/srmnitc/pyscal/releases/tag/v1.0.0>`_ of pyscal is released.

Citing the work
---------------

If you use pyscal in your work, the citation of the `following article <https://joss.theoj.org/papers/10.21105/joss.01824>`_ will be greatly appreciated:

Sarath Menon, Grisell Díaz Leines and Jutta Rogal (2019). pyscal: A python module for structural analysis of atomic environments. Journal of Open Source Software, 4(43), 1824, https://doi.org/10.21105/joss.01824

Citation in bib format can be downloaded `here <https://rubde-my.sharepoint.com/:u:/g/personal/sarath_menon_rub_de/Ecfuz7X8__ZJiz73k-dvvpEBjjMU6VJvg0v-hDtsFd3Kkw?download=1>`_.

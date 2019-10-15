
Methods to calculate neighbors of a particle
--------------------------------------------

pyscal includes different methods to explore the local environment of a
particle that rely on the calculation of nearest neighbors. Various
approaches to compute the neighbors of particles are discussed here.

Fixed cutoff method
~~~~~~~~~~~~~~~~~~~

The most common method to calculate the nearest neighbors of an atom is using a
cutoff radius. The neighborhood of an atom for calculation of
Steinhardt's parameters [1]_ is often carried out using this method. Commonly, a cutoff is
selected as the first minimum of the radial distribution functions. Once a cutoff is
selected, the neighbors of an atom are those that fall within this
selected radius. The following code snippet will use the cutoff method
to calculate neighbors. Please check the `examples section <https://pyscal.readthedocs.io/en/latest/examples.html#basic-examples>`_ for basic use
of the module. In this example, ``conf.dump`` is assumed to be the input configuration
of the system. A cutoff radius of 3 is assumed for calculation of
neighbors.

.. code:: python

    import pyscal.core as pc
    sys = pc.System()
    sys.read_inputfile('conf.dump')
    sys.find_neighbors(method='cutoff', cutoff=3)

Adaptive cutoff methods
~~~~~~~~~~~~~~~~~~~~~~~

A fixed cutoff radius can introduce limitations to explore the local environment of the particle in
some cases:

-  At finite temperatures, when thermal fluctuations take place, the selection
   of a fixed cutoff may result in an inaccurate description of the local environment.

-  If there is more than one structure present in the system, for
   example, bcc and fcc, the selection of cutoff such that it includes
   the first shell of both structures can be difficult.

In order to achieve a more accurate description of the local environment, various adaptive approaches
have been proposed. Two of the methods implemented in the module are
discussed below.

Solid angle based nearest neighbor algorithm (SANN)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

SANN algorithm [2]_ determines the cutoff radius by counting the solid
angles around an atom and equating it to :math:`4\pi`. The algorithm
solves the following equation iteratively.


.. math:: R_i^{(m)} = \frac{\sum_{j=1}^m r_{i,j}}{m-2} < r_{i, m+1}

where :math:`i` is the host atom, :math:`j` are its neighbors with :math:`r_{ij}`
is the distance between atoms :math:`i` and :math:`j`.
:math:`R_i` is the cutoff radius for each particle :math:`i` which is
found by increasing the neighbor of neighbors :math:`m` iteratively. For a
description of the algorithm and more details, please check the reference [2]_. SANN
algorithm can be used to find the neighbors by,

.. code:: python

    import pyscal.core as pc
    sys = pc.System()
    sys.read_inputfile('conf.dump')
    sys.find_neighbors(method='cutoff', cutoff='sann')

Since SANN algorithm involves sorting, a sufficiently large cutoff is
used in the beginning to reduce the number entries to be sorted. This
parameter is calculated by,


  .. math::  r_{initial} = \mathrm{threshold} \times \bigg(\frac{\mathrm{Simulation~box~volume}}{\mathrm{Number~of~particles}}\bigg)^{\frac{1}{3}}

a tunable ``threshold`` parameter can be set through function arguments.

Adaptive cutoff method
^^^^^^^^^^^^^^^^^^^^^^

An adaptive cutoff specific for each atom can also be found using an
algorithm similar to adaptive common neighbor analysis [3]_. This adaptive
cutoff is calculated by first making a list of all neighbor distances
for each atom similar to SANN method. Once this list is available,
then the cutoff is calculated from,


  .. math::  r_{cut}(i) = \mathrm{padding}\times \bigg(\frac{1}{\mathrm{nlimit}} \sum_{j=1}^{\mathrm{nlimit}} r_{ij} \bigg)

This method can be chosen by,

.. code:: python

    import pyscal.core as pc
    sys = pc.System()
    sys.read_inputfile('conf.dump')
    sys.find_neighbors(method='cutoff', cutoff='adaptive')

The ``padding`` and ``nlimit`` parameters in the above equation can be
tuned using the respective keywords.

Either of the adaptive method can be used to find neighbors, which can
then be used to calculate Steinhardt's parameters or their averaged version.

Voronoi tessellation
~~~~~~~~~~~~~~~~~~~~

`Voronoi tessellation <https://en.wikipedia.org/wiki/Voronoi_diagram>`_
provides a completely parameter free geometric
approach for calculation of neighbors. `Voro++ <http://math.lbl.gov/voro++/>`_ code is used for
Voronoi tessellation. Neighbors can be calculated using this method by,

.. code:: python

    import pyscal.core as pc
    sys = pc.System()
    sys.read_inputfile('conf.dump')
    sys.find_neighbors(method='voronoi')

Finding neighbors using Voronoi tessellation also calculates a weight
for each neighbor. The weight of a neighbor :math:`j` towards a host
atom :math:`i` is given by,


  .. math::  W_{ij} = \frac{A_{ij}}{\sum_{j=1}^N A_{ij}}

where :math:`A_{ij}` is the area of Voronoi facet between atom :math:`i` and :math:`j`,
:math:`N` are all the neighbors identified through Voronoi
tessellation. This weight can be used later for calculation of
weighted Steinhardt's parameters. Optionally, it is possible to choose
the exponent for this weight. Option ``voroexp`` is used to set this
option. For example if ``voroexp=2``, the weight would be calculated as,

  .. math::  W_{ij} = \frac{A_{ij}^2}{\sum_{j=1}^N A_{ij}^2}


.. [1] `Steinhardt, PJ, Nelson, DR, Ronchetti, M. Phys Rev B 28, 1983 <https://journals.aps.org/prb/abstract/10.1103/PhysRevB.28.784>`_.
.. [2] `van Meel, JA, Filion, L, Valeriani, C, Frenkel, D, J Chem Phys 234107, 2012 <https://aip.scitation.org/doi/full/10.1063/1.4729313>`_.
.. [3] `Stukowski, A, Model Simul Mater SC 20, 2012 <https://iopscience.iop.org/article/10.1088/0965-0393/20/4/045021/meta>`_.

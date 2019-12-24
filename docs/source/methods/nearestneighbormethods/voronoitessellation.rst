
Voronoi tessellation
--------------------

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

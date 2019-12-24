
Voronoi weighted Steinhardt's parameters
----------------------------------------

In order to improve the resolution of crystal structures Mickel et al [1]_
proposed weighting the contribution of each neighbor to the Steinhardt
parameters by the ratio of the area of the Voronoi facet shared between
the neighbor and host atom. The weighted parameters are given by,

.. math::  q_{lm} (i) =  \frac{1}{N(i)} \sum_{j=1}^{N(i)} \frac{A_{ij}}{A} Y_{lm}(\pmb{r}_{ij})

where :math:`A_{ij}` is the area of the Voronoi facet between atoms
:math:`i` and :math:`j` and :math:`A` is the sum of the face areas of
atom :math:`i`. In pyscal, the area weights are already assigned
during the neighbor calculation phase when the Voronoi method is used to
calculate neighbors in the :func:`~pyscal.core.System.find_neighbors`.
The Voronoi weighted Steinhardt's parameters can be
calculated as follows,

.. code:: python

    sys.find_neighbors(method='voronoi')
    sys.calculate_q([4, 6])
    q = sys.get_qvals([4, 6])

The weighted Steinhardt's parameters can also be averaged as described
above. Once again, the keyword ``averaged=True`` can be used for this
purpose.

.. code:: python

    sys.find_neighbors(method='voronoi')
    sys.calculate_q([4, 6], averaged=True)
    q = sys.get_qvals([4, 6], averaged=True)

It was also proposed that higher powers of the weight [2]_
:math:`\frac{A_{ij}^{\alpha}}{A(\alpha)}` where :math:`\alpha = 2, 3` can also
be used, where :math:`A(\alpha) = \sum_{j=1}^{N(i)} A_{ij}^{\alpha}` The value of this can be set using the keyword ``voroexp``
during the neighbor calculation phase.

.. code:: python

    sys.find_neighbors(method='voronoi', voroexp=2)

If the value of ``voroexp`` is set to 0, the neighbors would be found
using Voronoi method, but the calculated Steinhardt's parameters will
not be weighted.

.. [1] `Mickel, W, Kapfer, SC, Schroder-Turk, GE, Mecke, K, J Chem Phys 138, 2013 <https://aip.scitation.org/doi/full/10.1063/1.4774084>`_.
.. [2] `Haeberle, J, Sperl, M, Born, P Arxiv 2019 <https://arxiv.org/abs/1906.08111>`_.


Steinhardt's bond orientational order parameters
------------------------------------------------

Steinhardt's parameters
~~~~~~~~~~~~~~~~~~~~~~~

Steinhardt's bond orientational order parameters [1]_ are a set of parameters
based on `spherical harmonics <https://en.wikipedia.org/wiki/Spherical_harmonics>`_
to explore the local atomic environment. These parameters have been used
extensively for various uses such as distinction of crystal structures,
identification of solid and liquid atoms and identification of defects [2]_.

These parameters, which are rotationally and translationally invariant
are defined by,

.. math::  q_l (i) =  \Big(  \frac{4\pi}{2l+1}  \sum_{m=-l}^l | q_{lm}(i) |^2 \Big )^{\frac{1}{2}}

where,

.. math::  q_{lm} (i) =  \frac{1}{N(i)} \sum_{j=1}^{N(i)} Y_{lm}(\pmb{r}_{ij})

in which :math:`Y_{lm}` are the spherical harmonics and :math:`N(i)` is
the number of neighbours of particle :math:`i`, :math:`\pmb{r}_{ij}` is
the vector connecting particles :math:`i` and :math:`j`, and :math:`l`
and :math:`m` are both intergers with :math:`m \in [-l,+l]`. Various
parameters have found specific uses, such as :math:`q_2` and :math:`q_6`
for identification of crystallinity, :math:`q_6` for identification of
solidity, and :math:`q_4` and :math:`q_6` for distinction of crystal
structures [2]_. Commonly this method uses a cutoff radius to identify the neighbors of an atom. The cutoff can be chosen
based on different methods available. Once the cutoff is chosen and
neighbors are calculated, the calculation of Steinhardt parameters is
straightforward.

.. code:: python

    sys.calculate_q([4, 6])
    q = sys.get_qvals([4, 6])

Averaged Steinhardt's parameters
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

At high temperatures, thermal vibrations affect the atomic positions.
This in turn leads to overlapping distributions of :math:`q_l`
parameters, which makes the identification of
crystal structures difficult. To address this problem, the averaged
version :math:`\bar{q}_l` of Steinhardt's parameters was introduced by Lechner
and Dellago [3]_. :math:`\bar{q}_l` is given by,

.. math::  \bar{q}_l (i) =  \Big(  \frac{4\pi}{2l+1}  \sum_{m=-l}^l \Big| \frac{1}{\tilde{N}(i)} \sum_{k=0}^{\tilde{N}(i)} q_{lm}(k) \Big|^2 \Big )^{\frac{1}{2}}

where the sum from :math:`k=0` to :math:`\tilde{N}(i)` is over all the
neighbors and the particle itself. The averaged parameters takes into
account the first neighbor shell and also information from the
neighboring atoms and thus reduces the overlap between the
distributions. Commonly :math:`\bar{q}_4` and :math:`\bar{q}_6` are used
in identification of crystal structures.
Averaged versions can be calculated by setting the
keyword ``averaged=True`` as follows.

.. code:: python

    sys.calculate_q([4, 6], averaged=True)
    q = sys.get_qvals([4, 6], averaged=True)


Voronoi weighted Steinhardt's parameters [2]_
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In order to improve the resolution of crystal structures Mickel et al [2]_
proposed weighting the contribution of each neighbor to the Steinhardt
parameters by the ratio of the area of the Voronoi facet shared between
the neighbor and host atom. The weighted parameters are given by,

.. math::  q_{lm} (i) =  \frac{1}{N(i)} \sum_{j=1}^{N(i)} \frac{A_{ij}}{A} Y_{lm}(\pmb{r}_{ij})

where :math:`A_{ij}` is the area of the Voronoi facet between atoms
:math:`i` and :math:`j` and :math:`A` is the sum of the face areas of
atom :math:`i`. In ``pyscal``, the area weights are already assigned
during the neighbor calculation phase when the Voronoi method is used to
calculate neighbors (see `documentation <https://pyscal.readthedocs.io/en/latest/pyscal.html#pyscal.core.System.find_neighbors>`_).
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

It was also proposed that higher powers of the weight [4]_
:math:`\frac{A_{ij}^{\alpha}}{A(\alpha)}` where :math:`\alpha = 2, 3` can also
be used. The value of this can be set using the keyword ``voroexp``
during the neighbor calculation phase.

.. code:: python

    sys.find_neighbors(method='voronoi', voroexp=2)

If the value of ``voroexp`` is set to 0, the neighbors would be found
using Voronoi method, but the calculated Steinhardt's parameters will
not be weighted.


.. [1] `Steinhardt, PJ, Nelson, DR, Ronchetti, M. Phys Rev B 28, 1983 <https://journals.aps.org/prb/abstract/10.1103/PhysRevB.28.784>`_.
.. [2] `Mickel, W, Kapfer, SC, Schroder-Turk, GE, Mecke, K, J Chem Phys 138, 2013 <https://aip.scitation.org/doi/full/10.1063/1.4774084>`_.
.. [3] `Lechner, W, Dellago, C, J Chem Phys, 2013 <https://aip.scitation.org/doi/full/10.1063/1.2977970>`_.
.. [4] `Haeberle, J, Sperl, M, Born, P Arxiv 2019 <https://arxiv.org/abs/1906.08111>`_.

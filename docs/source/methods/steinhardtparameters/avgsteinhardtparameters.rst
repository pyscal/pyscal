
Averaged Steinhardt's parameters
--------------------------------

At high temperatures, thermal vibrations affect the atomic positions.
This in turn leads to overlapping distributions of :math:`q_l`
parameters, which makes the identification of
crystal structures difficult. To address this problem, the averaged
version :math:`\bar{q}_l` of Steinhardt's parameters was introduced by Lechner
and Dellago [1]_. :math:`\bar{q}_l` is given by,

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

.. [1] `Lechner, W, Dellago, C, J Chem Phys, 2013 <https://aip.scitation.org/doi/full/10.1063/1.2977970>`_.

..  note:: Associated methods

    :func:`~pyscal.core.System.find_neighbors`
    :func:`~pyscal.core.System.calculate_q`
    :func:`~pyscal.core.System.get_qvals`
    :func:`~pyscal.catom.Atom.get_q`

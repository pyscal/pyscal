
Disorder parameter
-------------------

Kawasaki and Onuki [1]_ proposed a disorder variable based
on Steinhardt's order paramaters which can be used to distinguish between ordered
and disordered structures.

The disorder variable for an atom is defined as,

.. math::

    D_j = \frac{1}{n_b^j} \sum_{k \in neighbors } [S_{jj} + S_{kk} - 2S_{jk}]


where S is given by,

.. math::

    S_{jk} = \sum_{-l \leq m \leq l} q_{lm}^j (q_{lm}^k)^*

l = 6 was used in the original publication as it is a good indicator of crystallinity.
However, l = 4 can also be used for treating bcc structures.
An averaged disorder parameter for each atom can also be calculated in pyscal,

.. math::

    \bar{D}_j = \frac{1}{n_b^j} \sum_{k \in neighbors } D_j

In pyscal, disorder parameter can be calculated by the following code-block,

.. code:: python

    import pyscal.core as pc
    sys = pc.System()
    sys.read_inputfile('conf.dump')
    sys.find_neighbors(method='cutoff', cutoff=0)
    sys.calculate_q(6)
    sys.calculate_disorder(averaged=True, q=6)

The value of q can be replaced with whichever is required from 2-12. The calculated values can be accessed by,
:attr:`~pyscal.catom.Atom.disorder` and :attr:`~pyscal.catom.Atom.avg_disorder`.

.. [1] Kawasaki .T, Onuki .A, JCP 135, 174109 (2011)

..  note:: Associated methods

    :func:`~pyscal.core.System.calculate_q`
    :func:`~pyscal.core.System.calculate_disorder`
    :attr:`~pyscal.catom.Atom.disorder`
    :attr:`~pyscal.catom.Atom.avg_disorder`
    

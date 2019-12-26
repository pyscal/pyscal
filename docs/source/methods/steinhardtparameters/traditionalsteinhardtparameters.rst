
Steinhardt's parameters
-----------------------

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
based on `different methods available <https://pyscal.readthedocs.io/en/latest/nearestneighbormethods.html>`_. Once the cutoff is chosen and
neighbors are calculated, the calculation of Steinhardt's parameters is
straightforward.

.. code:: python

    sys.calculate_q([4, 6])
    q = sys.get_qvals([4, 6])


.. [1] `Steinhardt, PJ, Nelson, DR, Ronchetti, M. Phys Rev B 28, 1983 <https://journals.aps.org/prb/abstract/10.1103/PhysRevB.28.784>`_.
.. [2] `Mickel, W, Kapfer, SC, Schroder-Turk, GE, Mecke, K, J Chem Phys 138, 2013 <https://aip.scitation.org/doi/full/10.1063/1.4774084>`_.

..  note:: Associated methods

    :func:`~pyscal.core.System.find_neighbors`
    :func:`~pyscal.core.System.calculate_q`
    :func:`~pyscal.core.System.get_qvals`
    :func:`~pyscal.catom.Atom.get_q`
    `Example <http://pyscal.com/en/latest/examples/steinhardtparameters/calculateq.html>`_

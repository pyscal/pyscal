
Steinhardt's bond orientational order parameters
------------------------------------------------

Steinhardt's parameters
~~~~~~~~~~~~~~~~~~~~~~~

Steinhardt's bond orientational order parameters are a set of parameters
based on the local atomic environment. These parameters have been used
extensively for various uses such as distinction of crystal structures,
identification of solid and liquid atoms and identification of defects.

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
structures. The traditional method makes use of a cutoff radius which
takes into account the neighbors of an atom. The cutoff can be chosen
based on different methods available. Once the cutoff is chosen and
neighbors are calculated, the calculation of Steinhardt parameters are
straightforward.

.. code:: python

    sys.calculate_q([4, 6])
    q = sys.get_qvals([4, 6])

Averaged Steinhardt's parameters
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

At high temperatures, thermal vibrations affect the atomic positions.
This in turn leads to overlapping distributions of :math:`q_l`
parameters, which in turn leads to difficulty in identification of
crystal structures. In order to address this problem, the averaged
version :math:`\bar{q}_l` of Steinhardt's parameters was introduced by
Lechner and Dellago. :math:`\bar{q}_l` is given by,

.. math::  \bar{q}_l (i) =  \Big(  \frac{4\pi}{2l+1}  \sum_{m=-l}^l \Big| \frac{1}{\tilde{N}(i)} \sum_{k=0}^{\tilde{N}(i)} q_{lm}(k) \Big|^2 \Big )^{\frac{1}{2}} 

where the sum from :math:`k=0` to :math:`\tilde{N}(i)` is over all the
neighbors and the particle itself. The averaged parameters takes into
account the first neighbor shell and also information from the
neighboring atoms and thus reduces the overlap between the
distributions. Averaged versions can be calculated by setting the
keyword ``averaged=True`` as follows.

.. code:: python

    sys.calculate_q([4, 6], averaged=True)
    q = sys.get_qvals([4, 6], averaged=True)

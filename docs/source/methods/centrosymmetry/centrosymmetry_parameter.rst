Centrosymmetry parameter
========================

Centrosymmetry parameter (CSP) was introduced by `Kelchner et
al. <https://journals.aps.org/prb/abstract/10.1103/PhysRevB.58.11085>`__
to identify defects in crystals. The parameter measures the loss of
local symmetry. For an atom with :math:`N` nearest neighbors, the
parameter is given by,

.. math::


   \mathrm{CSP} = \sum_{i=1}^{N/2} \big | \textbf{r}_i + \textbf{r}_{i+N/2} \big |^2

:math:`\textbf{r}_i` and :math:`\textbf{r}_{i+N/2}` are vectors from the
central atom to two opposite pairs of neighbors. There are two main
methods to identify the opposite pairs of neighbors as described in
`this publication <https://arxiv.org/abs/2003.08879>`__. The first of
the approaches is called `Greedy Edge
Selection <https://iopscience.iop.org/article/10.1088/0965-0393/20/4/045021/meta>`__\ (GES)
and is implemented in `LAMMPS <https://lammps.sandia.gov/>`__ and
`Ovito <https://www.ovito.org/>`__. GES algorithm calculates a weight
:math:`w_{ij} = |\textbf{r}_i + \textbf{r}_j|` for all combinations of
neighbors around an atom and calculates CSP over the smallest
:math:`N/2` weights.

A centrosymmetry parameter calculation using GES algorithm can be
carried out as follows-

.. code:: python

    import pyscal.core as pc
    sys = pc.System()
    sys.read_inputfile('conf.dump')
    sys.find_neighbors(method='voronoi')
    sys.calculate_centrosymmetry(nmax = 12)

``nmax`` parameter specifies the number of nearest neighbors to be
considered for the calculation of CSP. The second algorithm is called
the `Greedy Vertex
Matching <https://dl.acm.org/doi/book/10.5555/1206879>`__ and is
implemented in `AtomEye <http://li.mit.edu/Archive/Graphics/A/>`__ and
`Atomsk <https://atomsk.univ-lille.fr/>`__. This algorithm orders the
neighbors atoms in order of increasing distance from the central atom.
From this list, the closest neighbor is paired with its lowest weight
partner and both atoms removed from the list. This process is continued
until no more atoms are remaining in the list. CSP calculation using
this algorithm can be carried out by,

.. code:: python

    import pyscal.core as pc
    sys = pc.System()
    sys.read_inputfile('conf.dump')
    sys.find_neighbors(method='voronoi')
    sys.calculate_centrosymmetry(nmax = 12, algorithm = "gvm")

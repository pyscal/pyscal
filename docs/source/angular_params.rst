Angular parameter
-----------------

This illustrates the use of angular parameter to identify diamond
structure. Angular parameter was introduced by `Uttormark et
al. <https://journals.aps.org/prb/abstract/10.1103/PhysRevB.47.15717>`__,
and measures the tetrahedrality of the local atomic structure. An atom
belonging to diamond structure has four nearest neighbors which gives
rise to six three body angles around the atom. The angular parameter
:math:`A` is then defined as,

:math:`A = \sum_{i=1}^6 (\cos(\theta_i)+\frac{1}{3})^2`

An atom belonging to diamond structure would show the value of angular
params close to 0. The following example illustrates the use of this
parameter.

.. code:: python

    import pyscal.core as pc
    import pyscal.crystal_structures as pcs
    import numpy as np
    import matplotlib.pyplot as plt

Create structures
~~~~~~~~~~~~~~~~~

The first step is to create some structures using the :mod:`~pyscal.crystal_structures` module and assign it to a System. This can be done as
follows-

.. code:: python

    atoms, box = pcs.make_crystal('diamond', lattice_constant=4, repetitions=[3,3,3])
    sys = pc.System()
    sys.atoms = atoms
    sys.box = box

Now we can find the neighbors of all atoms. In this case we will use an
adaptive method using the :func:`~pyscal.core.find_neighbors` which can find an individual cutoff for each atom.

.. code:: python

    sys.find_neighbors(method='cutoff', cutoff='adaptive')

Finally, the angular criteria can be calculated by,

.. code:: python

    sys.calculate_angularcriteria()

The above function assigns the angular value for each atom in the attribute :attr:`~pyscal.catom.angular` which can be
accessed using,

.. code:: python

    atoms = sys.atoms
    angular = [atom.angular for atom in atoms]

The angular values are zero for atoms that belong to diamond structure.

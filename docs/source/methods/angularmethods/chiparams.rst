
:math:`\chi` parameters for structural identification
---------------------------------------------------------

:math:`\chi` parameters introduced by Ackland and Jones [1]_ measures all local angles created
by an atom with its neighbors and creates a histogram of these angles to produce vector which can be
used to identify structures. After finding the neighbors of an atom, :math:`\cos \theta_{ijk}` for
atoms j and k which are neighbors of i is calculated for all combinations of i, j and k. The set of all
calculated cosine values are then added to a histogram with the following bins - [-1.0, -0.945, -0.915,
-0.755, -0.705, -0.195, 0.195, 0.245, 0.795, 1.0]. Compared to :math:`\chi` parameters from
:math:`\chi_0` to :math:`\chi_7` in the associated publication, the vector calculated in pyscal contains
values from :math:`\chi_0` to :math:`\chi_8` which is due to an additional :math:`\chi` parameter which
measures the number of neighbors between cosines -0.705 to -0.195. The :math:`\chi` vector is characteristic
of the local atomic environment and can be used to identify crystal structures, details of which can be found
in the publication [1]_.

:math:`\chi` parameters can be calculated in pyscal using,

.. code:: python

    import pyscal.core as pc
    sys = pc.System()
    sys.read_inputfile('conf.dump')
    sys.find_neighbors(method='cutoff', cutoff='adaptive')
    sys.calculate_chiparams()

The calculated values for each atom can be accessed using :attr:`~pyscal.catom.Atom.chiparams`.

.. [1] Ackland, Jones, Phys. Rev. B 73, 2006

..  note:: Associated methods

    :func:`~pyscal.core.System.calculate_chiparams`
    :attr:`~pyscal.catom.Atom.chiparams`

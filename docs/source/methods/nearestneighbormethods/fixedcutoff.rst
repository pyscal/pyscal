
Fixed cutoff method
-------------------

The most common method to calculate the nearest neighbors of an atom is using a
cutoff radius. The neighborhood of an atom for calculation of
Steinhardt's parameters [1]_ is often carried out using this method. Commonly, a cutoff is
selected as the first minimum of the radial distribution functions. Once a cutoff is
selected, the neighbors of an atom are those that fall within this
selected radius. The following code snippet will use the cutoff method
to calculate neighbors. Please check the `examples section <https://pyscal.readthedocs.io/en/latest/examples.html#basic-examples>`_ for basic use
of the module. In this example, ``conf.dump`` is assumed to be the input configuration
of the system. A cutoff radius of 3 is assumed for calculation of
neighbors.

.. code:: python

    import pyscal.core as pc
    sys = pc.System()
    sys.read_inputfile('conf.dump')
    sys.find_neighbors(method='cutoff', cutoff=3)


.. [1] `Steinhardt, PJ, Nelson, DR, Ronchetti, M. Phys Rev B 28, 1983 <https://journals.aps.org/prb/abstract/10.1103/PhysRevB.28.784>`_.

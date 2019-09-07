Saving systems and atoms to file
--------------------------------

pyscal offers tools to save ``System`` and ``Atom`` classes to file.
This makes use of ``numpy.save`` and ``numpy.load`` methods to save the
information. Information including calculated neighbors, q values,
solidity, voronoi volume etc are all saved. Saving this information is
beneficial in systems with a large number of atoms as it saves time
which would otherwise be spent recalculating information. This example
illustrates the use of ``pyscal.pickle`` module.

.. code:: python

    import pyscal.core as pc
    import pyscal.traj_process as ptp

First a system is set up using an input file

.. code:: python

    sys = pc.System()
    sys.read_inputfile('examples/conf.dump')

In the next step, the neighbors are calculated, and :math:`\bar{q}_4`
and :math:`\bar{q}_6` are also calculated.

.. code:: python

    sys.get_neighbors(method='cutoff', cutoff=3.6)
    sys.calculate_q([4, 6], averaged=True)

In order to prevent recalculation of neighbors, for example if later we
need to calculate :math:`\bar{q}_8` and :math:`\bar{q}_{10}`, we can
save the set of systems to a file.

Saving individual system
~~~~~~~~~~~~~~~~~~~~~~~~

Single System instances can be saved without ``pyscal.pickle`` module
directly, similar to how pandas DataFrames are saved to file.

.. code:: python

    sys.to_file('test_system.npy')

Thats it! The information is saved in the while. Once again,
``pyscal.pickle.read_systems`` can be used to read the System instance.
Alternatively, a new System can be created and the information can be
read in from a file.

.. code:: python

    new_sys = pc.System()
    new_sys.from_file('test_system.npy')

This system retains all the information and hence it can be used for
further calculations

.. code:: python

    new_sys.calculate_q([8, 10], averaged=True)

Here :math:`\bar{q}_8` and :math:`\bar{q}_{10}` were calculated without
having to find neighbors again.

Saving atoms
~~~~~~~~~~~~

Alternatively, a list of atoms can also be saved to file.

.. code:: python

    atoms = new_sys.get_atoms()

This list of atoms can be saved to file using,

.. code:: python

    pp.write_atoms('atoms.npy', atoms)

Similar to Systems, atoms can also be read in from the file,

.. code:: python

    natoms = pp.read_atoms('atoms.npy')

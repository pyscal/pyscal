Saving systems and atoms to file
--------------------------------

pyscal offers tools to save ``System`` and ``Atom`` classes to file.
This makes use of ``numpy.save`` and ``numpy.load`` methods to save the
information. Information including calculated neighbors, q values,
solidity, voronoi volume etc are all saved. Saving this information is
beneficial in systems with a large number of atoms as it saves time
which would otherwise be spent recalculating information. This example
illustrates the use of ``pyscal.pickle`` module.

First an input file - an md trajectory of 10 time slices, each
containing 500 atoms is used. The trajectory is split into individual
time slices using ``pyscal.traj_process`` module.

.. code:: python

    import pyscal.core as pc
    import pyscal.traj_process as ptp

Saving a list of systems
~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: python

    files = ptp.split_trajectory('traj.light')

In the next step, each of the individual files are read into a system,
the neighbors are calculated, and :math:`\bar{q}_4` and
:math:`\bar{q}_6` are also calculated. The systems are then stored in an
array.

.. code:: python

    systems = []
    for file in files:
        sys = pc.System()
        sys.read_inputfile(file)
        sys.get_neighbors(method='cutoff', cutoff=3.6)
        sys.calculate_q([4, 6], averaged=True)
        systems.append(sys)

In order to prevent recalculation of neighbors, for example if later we
need to calculate :math:`\bar{q}_8` and :math:`\bar{q}_{10}`, we can
save the set of systems to a file.

.. code:: python

    import pyscal.pickle as pp

.. code:: python

    pp.write_systems('systems.npy', systems)

All the information from the system are stored in this file. At a later
point, the systems can be read in using,

.. code:: python

    rsystems = pp.read_systems('systems.npy')

Now further calculations can be carried out without starting over again!

.. code:: python

    for sys in rsystems:
        sys.calculate_q([8, 10], averaged=True)

Saving individual system
~~~~~~~~~~~~~~~~~~~~~~~~

Instead of saving a whole list of systems, single System instances can
also be saved. This can be done without ``pyscal.pickle`` module,
similar to how pandas DataFrames are saved to file. First, select one
system for testing,

.. code:: python

    test_system = rsystems[0]

This test system can be saved to a file by,

.. code:: python

    test_system.to_file('test_system.npy')

Thats it! The information is saved in the while. Once again,
``pyscal.pickle.read_systems`` can be used to read the System instance.
Alternatively, a new System can be created and the information can be
read in from a file.

.. code:: python

    new_sys = pc.System()
    new_sys.from_file('test_system.npy')

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

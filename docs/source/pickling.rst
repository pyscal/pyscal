Saving System and Atom objects to file
--------------------------------------

pyscal offers tools to save :class:`~pyscal.core.System` and :class:`~pyscal.catom.Atom` classes to file.
There are two methods to do this. The :func:`~pyscal.core.System.to_file` method can save a
trajectory file which contains all the atom positions and system box
dimensions.

The second approach makes use of `numpy.save <https://docs.scipy.org/doc/numpy/reference/generated/numpy.save.html>`_
and `numpy.load <https://docs.scipy.org/doc/numpy/reference/generated/numpy.load.html>`_
methods to save the information. Information including calculated
neighbors, q values, solidity, voronoi volume etc are all saved. Saving
this information is beneficial in systems with a large number of atoms
as it saves time which would otherwise be spent recalculating
information. The :mod:`~pyscal.pickle`
module provides the base functions for pickling support.

.. code:: python

    import pyscal.core as pc
    import pyscal.traj_process as ptp

First a system is set up using an input file

.. code:: python

    sys = pc.System()
    sys.read_inputfile('conf.bcc')

In the next step, the neighbors are calculated, and :math:`\bar{q}_4`
and :math:`\bar{q}_6` are also calculated.

.. code:: python

    sys.find_neighbors(method='cutoff', cutoff=3.6)
    sys.calculate_q([4, 6], averaged=True)

In order to prevent recalculation of neighbors, for example if later we
need to calculate :math:`\bar{q}_8` and :math:`\bar{q}_{10}`, we can
save the set of systems to a file.

Saving individual system
~~~~~~~~~~~~~~~~~~~~~~~~

Single System instances can be saved without :mod:`~pyscal.pickle` module
directly, similar to how pandas `DataFrames <https://pandas.pydata.org/pandas-docs/stable/reference/api/pandas.DataFrame.html>`_ are saved to file.

.. code:: python

    sys.to_pickle('test_system.npy')

Thats it! The information is saved in the while. Once again,
:func:`~pyscal.pickle.read_systems` can be used to read the System instance.
Alternatively, a new System can be created and the information can be
read in from a file.

.. code:: python

    new_sys = pc.System()
    new_sys.from_pickle('test_system.npy')

This system retains all the information and hence it can be used for
further calculations

.. code:: python

    new_sys.calculate_q([8, 10], averaged=True)

Here :math:`\bar{q}_8` and :math:`\bar{q}_{10}` were calculated without
having to find neighbors again.

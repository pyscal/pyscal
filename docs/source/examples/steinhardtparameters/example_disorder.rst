Disorder variable
-----------------

In this example, `disorder variable <http://pyscal.com/en/latest/methods/steinhardtparameters/disorder.html>`__ which was
introduced to measure the disorder of a system is explored. We start by
importing the necessary modules. We will use
:mod:`~pyscal.crystal_structures` to create the necessary crystal
structures.

.. code:: python

    import pyscal.core as pc
    import pyscal.crystal_structures as pcs
    import matplotlib.pyplot as plt
    import numpy as np

First an fcc structure with a lattice constant of 4.00 is created.

.. code:: python

    fcc_atoms, fcc_box = pcs.make_crystal('fcc', lattice_constant=4, repetitions=[4,4,4])

The created atoms and box are assigned to a
:class:`~pyscal.core.System` object.

.. code:: python

    fcc = pc.System()
    fcc.atoms = fcc_atoms
    fcc.box = fcc_box

The next step is `find the neighbors <http://pyscal.com/en/latest/methods/nearestneighbormethods/nearestneighbormethods.html>`_, and the calculate the `Steinhardt
parameter <http://pyscal.com/en/latest/methods/steinhardtparameters/traditionalsteinhardtparameters.html>`_ based on which we could calculate the disorder variable.

.. code:: python

    fcc.find_neighbors(method='cutoff', cutoff='adaptive')

Once the neighbors are found, we can calculate the Steinhardt parameter
value. In this example :math:`q=6` will be used.

.. code:: python

    fcc.calculate_q(6)

Finally, disorder parameter can be calculated.

.. code:: python

    fcc.calculate_disorder()

The calculated disorder value can be accessed for each atom using the
:attr:`~pyscal.catom.disorder` variable.

.. code:: python

    atoms = fcc.atoms

.. code:: python

    disorder = [atom.disorder for atom in atoms]

.. code:: python

    np.mean(disorder)




.. parsed-literal::

    -1.041556887034408e-16



As expected, for a perfect fcc structure, we can see that the disorder
is zero. The variation of disorder variable on a distorted lattice can
be explored now. We will once again use the ``noise`` keyword along with
:func:`~pyscal.crystal_structures.make_crystal` to create a distorted
lattice.

.. code:: python

    fcc_atoms_d1, fcc_box_d1 = pcs.make_crystal('fcc', lattice_constant=4, repetitions=[4,4,4], noise=0.01)
    fcc_d1 = pc.System()
    fcc_d1.atoms = fcc_atoms_d1
    fcc_d1.box = fcc_box_d1

Once again, find neighbors and then calculate disorder

.. code:: python

    fcc_d1.find_neighbors(method='cutoff', cutoff='adaptive')
    fcc_d1.calculate_q(6)
    fcc_d1.calculate_disorder()

Check the value of disorder

.. code:: python

    atoms_d1 = fcc_d1.atoms

.. code:: python

    disorder = [atom.disorder for atom in atoms_d1]

.. code:: python

    np.mean(disorder)




.. parsed-literal::

    0.013889967380485688



The value of average disorder for the system has increased with noise.
Finally trying with a high amount of noise.

.. code:: python

    fcc_atoms_d2, fcc_box_d2 = pcs.make_crystal('fcc', lattice_constant=4, repetitions=[4,4,4], noise=0.1)
    fcc_d2 = pc.System()
    fcc_d2.atoms = fcc_atoms_d2
    fcc_d2.box = fcc_box_d2

.. code:: python

    fcc_d2.find_neighbors(method='cutoff', cutoff='adaptive')
    fcc_d2.calculate_q(6)
    fcc_d2.calculate_disorder()

.. code:: python

    atoms_d2 = fcc_d2.atoms

.. code:: python

    disorder = [atom.disorder for atom in atoms_d2]
    np.mean(disorder)




.. parsed-literal::

    1.8469165876016702



The value of disorder parameter shows an increase with the amount of
lattice distortion. An averaged version of disorder parameter, averaged
over the neighbors for each atom can also be calculated as shown below.

.. code:: python

    fcc_d2.calculate_disorder(averaged=True)

.. code:: python

    atoms_d2 = fcc_d2.atoms
    disorder = [atom.avg_disorder for atom in atoms_d2]
    np.mean(disorder)




.. parsed-literal::

    1.850630088115515



The disorder parameter can also be calculated for values of Steinhardt
parameter other than 6. For example,

.. code:: python

    fcc_d2.find_neighbors(method='cutoff', cutoff='adaptive')
    fcc_d2.calculate_q([4, 6])
    fcc_d2.calculate_disorder(q=4, averaged=True)

.. code:: python

    atoms_d2 = fcc_d2.atoms
    disorder = [atom.disorder for atom in atoms_d2]
    np.mean(disorder)




.. parsed-literal::

    1.880741277448693



:math:`q=4`, for example, can be useful when measuring disorder in bcc
crystals


Getting started with pybop
--------------------------

This example illustrates basic functionality of pyscal python library by
setting up a system and reading in atoms.

.. code:: ipython3

    import pyscal.core as pc
    import numpy as np

The ``System`` class
~~~~~~~~~~~~~~~~~~~~

``System`` is the basic class of pyscal and is required to be setup in
order to perform any calculations. It can be set up as easily as-

.. code:: ipython3

    sys = pc.System()

``sys`` is a ``System`` object. But at this point, it is completely
empty. We have to populate the system with two major information- \* the
simulation box dimensions \* the information regarding individual atoms.

| Let us try setting up a small system, which is the bcc unitcell of
  lattice constant 1. The simulation box dimensions of such a unit cell
  would be [[0.0, 1.0], [0.0, 1.0], [0.0, 1.0]] where the first set
  correspond to the x axis, second to y axis and so on.
| The unitcell has 2 atoms and their positions are [0,0,0] and [0.5,
  0.5, 0.5].

.. code:: ipython3

    sys.set_box([[0.0, 1.0], [0.0, 1.0], [0.0, 1.0]])

We can easily check if everything worked by getting the box dimensions

.. code:: ipython3

    box = sys.get_box()
    box




.. parsed-literal::

    [[0.0, 1.0], [0.0, 1.0], [0.0, 1.0]]



The ``Atom`` class
~~~~~~~~~~~~~~~~~~

| The next part is assigning the atoms. This can be done using the
  ``Atom`` class. Here, we will only look at the basic properties of
  ``Atom`` class. For a more detailed description, check the examples.
| Now lets create two atoms.

.. code:: ipython3

    atom1 = pc.Atom()
    atom2 = pc.Atom()

Now two empty atom objects are created. The major poperties of an atom
are its positions and id. Lets set this up.

.. code:: ipython3

    atom1.set_pos([0., 0., 0.])
    atom1.set_id(0)
    atom2.set_pos([0.5, 0.5, 0.5])
    atom2.set_id(1)

Alternatively, atom objects can also be set up as

.. code:: ipython3

    atom1 = pc.Atom(pos=[0., 0., 0.], id=0)
    atom2 = pc.Atom(pos=[0.5, 0.5, 0.5], id=1)

We can check the details of the atom by querying it

.. code:: ipython3

    x1 = atom1.get_x()
    x1




.. parsed-literal::

    [0.0, 0.0, 0.0]



Combining ``System`` and ``Atom``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Now that we have created the atoms, we can assign them to the system. We
can also assign the same box we created before.

.. code:: ipython3

    sys = pc.System()
    sys.assign_atoms([atom1, atom2], box)

That sets up the system completely. It has both of it's constituents -
atoms and the simulation box. We can check if everything works
correctly.

.. code:: ipython3

    atoms = sys.get_atoms()


.. parsed-literal::

    /home/users/menonsqr/anaconda3/envs/ml/lib/python3.7/site-packages/pyscal-1.0.1-py3.7-linux-x86_64.egg/pyscal/core.py:585: UserWarning: If the loc of atom is changed and set to system, it will overwrite the existing data, if any.
      warnings.warn("If the loc of atom is changed and set to system, it will overwrite the existing data, if any.")


This returns all the atoms of the system. Once you have all the atoms,
you can modify any one and set it back to the system. The following
statement will set the type of the first atom to 2.

.. code:: ipython3

    atom = atoms[0]
    atom.set_type(2)

Lets verify if it was done properly

.. code:: ipython3

    atom.get_type()




.. parsed-literal::

    2



Now we can push the atom back to the system

.. code:: ipython3

    sys.set_atom(atom)

We can also get individual atoms from the system instead of getting all
of them

.. code:: ipython3

    atom = sys.get_atom(0)

the above statement will return the atom at position 0

Reading in an input file
~~~~~~~~~~~~~~~~~~~~~~~~

| We are all set! The ``System`` is ready for calculations. However, in
  most realistic simulation situations, we have many atoms and it can be
  difficult to set each of them
| individually. In this situation we can read in input file directly. An
  example input file containing 500 atoms in a simulation box can be
  read in automatically. The file we use for this example is a file of
  the `lammps-dump <https://lammps.sandia.gov/doc/dump.html>`__ format.
  ``pyscal`` can also read in POSCAR files. In principle, ``pyscal``
  only needs the atom positions and simulation box size, so you can
  write a python function to process the input file, extract the details
  and pass to ``pyscal``.

.. code:: ipython3

    sys = pc.System()
    sys.read_inputfile('conf.dump')

Once again, lets check if the box dimensions are read in correctly

.. code:: ipython3

    box = sys.get_box()
    box




.. parsed-literal::

    [[-7.66608, 11.1901], [-7.66915, 11.1931], [-7.74357, 11.2676]]



Now we can get all atoms that belong to this system

.. code:: ipython3

    atoms = sys.get_atoms()
    len(atoms)




.. parsed-literal::

    500



| We can see that all the atoms are read in correctly and there are 500
  atoms in total. Once again, individual atom properties can be
| accessed as before.

.. code:: ipython3

    atoms[0].get_x()




.. parsed-literal::

    [-5.66782, -6.06781, -6.58151]



Thats it! Now we are ready for some calculations. You can find more in
the examples section of the documentation.

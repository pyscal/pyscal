
Getting started with pybop
--------------------------

| This example illustrates basic functionality of pybop python library
  by setting up a system and reading in atoms. An interactive version of
  this example would be made available soon.
| We start by importing the module

.. code:: ipython2

    import pybop.core as pc
    import numpy as np

The System class
~~~~~~~~~~~~~~~~

``System`` is the basic class of pydoc and is required to be setup in
order to perform any calculations. It can be set up as easily as-

.. code:: ipython2

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

.. code:: ipython2

    sys.set_box([[0.0, 1.0], [0.0, 1.0], [0.0, 1.0]])

We can easily check if everything worked by getting the box dimensions

.. code:: ipython2

    box = sys.get_box()
    box




.. parsed-literal::

    [0.0, 1.0, 0.0, 1.0, 0.0, 1.0]



The Atom class
~~~~~~~~~~~~~~

| The next part is assigning the atoms. This can be done using the
  ``Atom`` class. Here, we will only look at the basic properties of
  ``Atom`` class. For a more detailed description, check the examples.
| Now lets create two atoms.

.. code:: ipython2

    atom1 = pc.Atom()
    atom2 = pc.Atom()

Now two empty atom objects are created. The major poperties of an atom
are its positions and id. Lets set this up.

.. code:: ipython2

    atom1.set_x([0, 0, 0])
    atom1.set_id(0)
    atom2.set_x([0.5, 0.5, 0.5])
    atom2.set_id(1)

We can check the details of the atom by querying it

.. code:: ipython2

    x1 = atom1.get_x()
    x1




.. parsed-literal::

    [0.0, 0.0, 0.0]



Reading in an input file
~~~~~~~~~~~~~~~~~~~~~~~~

| We are all set! The ``System`` is ready for calculations. However, in
  most realistic simulation situations, we have many atoms and it can be
  difficult to set each of them
| individually. In this situation we can read in input file directly. An
  example input file containing 500 atoms in a simulation box can be
  read in automatically.

.. code:: ipython2

    sys = pc.System()
    sys.read_inputfile('conf.dump')

Once again, lets check if the box dimensions are read in correctly

.. code:: ipython2

    box = sys.get_box()
    box




.. parsed-literal::

    [-7.66608, 11.1901, -7.66915, 11.1931, -7.74357, 11.2676]



Now we can get all atoms that belong to this system

.. code:: ipython2

    atoms = sys.get_allatoms()
    len(atoms)




.. parsed-literal::

    500



| We can see that all the atoms are read in correctly and there are 500
  atoms in total. Once again, individual atom properties can be
| accessed as before.

.. code:: ipython2

    atoms[0].get_x()




.. parsed-literal::

    [-5.66782, -6.06781, -6.58151]



Thats it! Now we are ready for some calculations. You can find more in
the examples section of the documentation.

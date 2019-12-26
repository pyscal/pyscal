Getting started with pyscal
---------------------------

This example illustrates basic functionality of pyscal python library by
setting up a system and the atoms.

.. code:: python

    import pyscal.core as pc
    import numpy as np

The :class:`~pyscal.core.System` class
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

 :class:`~pyscal.core.System` is the basic class of pyscal and is required to be setup in
order to perform any calculations. It can be set up as-

.. code:: python

    sys = pc.System()

``sys`` is a :class:`~pyscal.core.System` object. But at this point, it is completely
empty. We have to provide the system with the following information-

* the simulation box dimensions
* the positions of individual atoms.

Let us try to set up a small system, which is the bcc unitcell of
lattice constant 1. The simulation box dimensions of such a unit cell
would be [[0.0, 1.0], [0.0, 1.0], [0.0, 1.0]] where the first set
correspond to the x axis, second to y axis and so on.
The unitcell has 2 atoms and their positions are [0,0,0] and [0.5, 0.5, 0.5].

.. code:: python

    sys.box = [[0.0, 1.0], [0.0, 1.0], [0.0, 1.0]]

We can easily check if everything worked by getting the box dimensions

.. code:: python

    sys.box




.. parsed-literal::

    [[0.0, 1.0], [0.0, 1.0], [0.0, 1.0]]



The :class:`~pyscal.catom.Atom` class
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The next part is assigning the atoms. This can be done using the
:class:`~pyscal.catom.Atom` class. Here, we will only look at the basic properties of
:class:`~pyscal.catom.Atom` class. For a more detailed description, check the
`documentation <https://pyscal.readthedocs.io/en/latest/pyscal.html#module-pyscal.catom>`__.
Now let us create two atoms.

.. code:: python

    atom1 = pc.Atom()
    atom2 = pc.Atom()

Now two empty atom objects are created. The basic poperties of an atom
are its positions and id. There are various other properties which can
be set here. A detailed description can be found
`here <https://pyscal.readthedocs.io/en/latest/pyscal.html#module-pyscal.catom>`__.

.. code:: python

    atom1.pos = [0., 0., 0.]
    atom1.id = 0
    atom2.pos = [0.5, 0.5, 0.5]
    atom2.id = 1

Alternatively, atom objects can also be set up as

.. code:: python

    atom1 = pc.Atom(pos=[0., 0., 0.], id=0)
    atom2 = pc.Atom(pos=[0.5, 0.5, 0.5], id=1)

We can check the details of the atom by querying it

.. code:: python

    atom1.pos




.. parsed-literal::

    [0.0, 0.0, 0.0]



Combining :class:`~pyscal.core.System` and :class:`~pyscal.catom.Atom`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Now that we have created the atoms, we can assign them to the system. We
can also assign the same box we created before.

.. code:: python

    sys = pc.System()
    sys.atoms = [atom1, atom2]
    sys.box = [[0.0, 1.0], [0.0, 1.0], [0.0, 1.0]]

That sets up the system completely. It has both of it's constituents -
atoms and the simulation box. We can check if everything works
correctly.

.. code:: python

    sys.atoms




.. parsed-literal::

    [<pyscal.catom.Atom at 0x7fb343025830>, <pyscal.catom.Atom at 0x7fb343025b30>]



This returns all the atoms of the system. Alternatively a single atom
can be accessed by,

.. code:: python

    atom = sys.get_atom(1)

The above call will fetch the atom at position 1 in the list of all
atoms in the system.

Once you have all the atoms, you can modify any one and add it back to
the list of all atoms in the system. The following statement will set
the type of the first atom to 2.

.. code:: python

    atom = sys.atoms[0]
    atom.type = 2

Lets verify if it was done properly

.. code:: python

    atom.type




.. parsed-literal::

    2



Now we can push the atom back to the system with the new type

.. code:: python

    sys.set_atom(atom)

.. warning::

      Due to Atom being a completely C++ class, it is necessary to use :func:`~pyscal.core.System.get_atom`
      and :func:`~pyscal.core.System.set_atom` to access individual atoms and set them back into the system
      object after modification. A list of all atoms however can be accessed directly by :attr:`~pyscal.core.System.atoms`

Reading in an input file
~~~~~~~~~~~~~~~~~~~~~~~~

We are all set! The :class:`~pyscal.core.System` is ready for calculations. However, in
most realistic simulation situations, we have many atoms and it can be
difficult to set each of them
individually. In this situation we can read in input file directly. An
example input file containing 500 atoms in a simulation box can be
read in automatically. The file we use for this example is a file of
the `lammps-dump <https://lammps.sandia.gov/doc/dump.html>`__ format.
pyscal can also read in POSCAR files. In principle, pyscal
only needs the atom positions and simulation box size, so you can
write a python function to process the input file, extract the details
and pass to pyscal.

.. code:: python

    sys = pc.System()
    sys.read_inputfile('conf.dump')

Once again, lets check if the box dimensions are read in correctly

.. code:: python

    sys.box




.. parsed-literal::

    [[-7.66608, 11.1901], [-7.66915, 11.1931], [-7.74357, 11.2676]]



Now we can get all atoms that belong to this system

.. code:: python

    len(sys.atoms)




.. parsed-literal::

    500



We can see that all the atoms are read in correctly and there are 500
atoms in total. Once again, individual atom properties can be
accessed as before.

.. code:: python

    sys.atoms[0].pos




.. parsed-literal::

    [-5.66782, -6.06781, -6.58151]



Thats it! Now we are ready for some calculations. You can find more in
the examples section of the documentation.

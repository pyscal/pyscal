
Calculating coordination numbers
--------------------------------

| In this example, we will read in a snapshot from an MD simulation and
  then calculate the coordination number distribution.
| This example assumes that you read the basic example.

.. code:: ipython3

    import pybop.core as pc
    import numpy as np
    import matplotlib.pyplot as plt

Read in a file
~~~~~~~~~~~~~~

The first step is setting up a system. We can create atoms and
simulation box using the ``pybop.crystal_structures`` module. Let’s
start by importing the module.

.. code:: ipython3

    import pybop.crystal_structures as pcs

.. code:: ipython3

    atoms, box = pcs.make_crystal('bcc', lattice_constant= 4.00, repetitions=[6,6,6])

The above function creates an bcc crystal of 6x6x6 unit cells with a
lattice constant of 4.00 along with a simulation box that encloses the
particles. We can then create a ``System`` and assign the atoms and box
to it.

.. code:: ipython3

    sys = pc.System()
    sys.assign_atoms(atoms, box)

Calculating neighbors
~~~~~~~~~~~~~~~~~~~~~

We start by calculating the neighbors of each atom in the system. There
are two ways to do this, using a ``cutoff`` method and using a
``voronoi`` polyhedra method. We will try with both of them. First we
try with cutoff system.

Cutoff method
^^^^^^^^^^^^^

Cutoff method takes cutoff distance value and finds all atoms within the
cutoff distance of the host atom.

.. code:: ipython3

    sys.get_neighbors(method='cutoff', cutoff=4.1)

Now lets get all the atoms.

.. code:: ipython3

    atoms = sys.get_atoms()

lets try accessing the coordination number of an atom

.. code:: ipython3

    atoms[0].get_coordination()




.. parsed-literal::

    14



As we would expect for a bcc type lattice, we see that the atom has 14
neighbors (8 in the first shell and 6 in the second). Lets try a more
interesting example by reading in a bcc system with thermal vibrations.
Thermal vibrations lead to distortion in atomic positions, and hence
there will be a distribution of coordination numbers.

.. code:: ipython3

    sys = pc.System()
    sys.read_inputfile('conf.dump')
    sys.get_neighbors(method='cutoff', cutoff=3.6)
    atoms = sys.get_atoms()

We can loop over all atoms and create a histogram of the results

.. code:: ipython3

    coord = [atom.get_coordination() for atom in atoms]

Now lets plot and see the results

.. code:: ipython3

    nos, counts = np.unique(coord, return_counts=True)
    plt.bar(nos, counts, color="#AD1457")
    plt.ylabel("density")
    plt.xlabel("coordination number")
    plt.title("Cutoff method")




.. parsed-literal::

    Text(0.5, 1.0, 'Cutoff method')




.. image:: output_23_1.png


Voronoi method
~~~~~~~~~~~~~~

Voronoi method calculates the voronoi polyhedra of all atoms. Any atom
that shares a voronoi face area with the host atom are considered
neighbors. Voronoi polyhedra is calculated using the Voro++ code.
However, you dont need to install this specifically as it is linked to
pybop.

.. code:: ipython3

    sys.get_neighbors(method='voronoi')

Once again, lets get all atoms and find their coordination

.. code:: ipython3

    atoms = sys.get_allatoms()
    coord = [atom.get_coordination() for atom in atoms]

And visualise the results

.. code:: ipython3

    nos, counts = np.unique(coord, return_counts=True)
    plt.bar(nos, counts, color="#AD1457")
    plt.ylabel("density")
    plt.xlabel("coordination number")
    plt.title("Voronoi method")




.. parsed-literal::

    Text(0.5, 1.0, 'Voronoi method')




.. image:: output_30_1.png


Finally..
~~~~~~~~~

Both methods find the coordination number, and the results are
comparable. Cutoff method is very sensitive to the choice of cutoff
radius, but voronoi method can slightly overestimate the neighbors due
to thermal vibrations.

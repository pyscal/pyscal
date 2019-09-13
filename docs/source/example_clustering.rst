Distinction of solid liquid atoms and clustering
------------------------------------------------

In this example, we will take one snapshot from a molecular dynamics
simulation which has a solid cluster in liquid. The task is to identify
solid atoms and cluster them. More details about the method can be found
`here <https://pyscal.readthedocs.io/en/latest/solidliquid.html>`_.

The first step is, of course, importing all the necessary module. For
visualisation, we will use `Ovito <https://www.ovito.org/>`__ [1]_.

.. figure:: system1.png
   :alt: original system

   liquid system with a solid cluster

Importing modules,

.. code:: python

    import pyscal.core as pc

Now we need to read this into a system, and calculate neighbors. Here we
will use a cutoff method to find neighbors. More details about finding
neighbors can be found
`here <https://pyscal.readthedocs.io/en/latest/nearestneighbormethods.html#>`__.

.. code:: python

    sys = pc.System()
    sys.read_inputfile('cluster.dump')
    sys.find_neighbors(method='cutoff', cutoff=3.63)

Now, the neighbors are found. Next step is finding solid atoms. This can
be done using ``System.find_solids`` method. There are few parameters
that can be set, which can be found in detail `here <https://pyscal.readthedocs.io/en/latest/solidliquid.html>`_.

.. code:: python

    sys.find_solids(bonds=6, threshold=0.5, avgthreshold=0.6)




.. parsed-literal::

    176



The above statement found all the solid atoms, also clustered them and
gave the size of the largest cluster of solid atoms as 176 atoms. If you
do not want to cluster atoms, clustering can be turned off by passing
``cluster=False`` to the above function. In order to visualise the solid
atoms, we need to first get the `issolid` value from each atom. `issolid`
gives a value of 1 if the atom is solid, 0 otherwise. `Atom.get_solid` method
can be used to access this value.

.. code:: python

    atoms = sys.get_atoms()
    solids = [atom.get_solid() for atom in atoms]

In order to visualise in Ovito, we need to first write it out
to a trajectory file. This can be done easily with the help of
``pyscal.traj_process`` module. We will save an extra column to the dump file 
with the keyword 'solid' and provide the calculated solidity array as the values of each
atom.

.. code:: python

    import pyscal.traj_process as ptp

.. code:: python

    ptp.write_structure(sys, 'sys.solid.dat', customkey='solid', customvals=solids)

We can now visualise this file in Ovito. After opening the file in
Ovito, the modifier `compute
property <https://ovito.org/manual/particles.modifiers.compute_property.html>`__
can be selected. The ``Output property`` should be ``selection`` and in
the expression field, ``solid==0`` can be selected to select all the non
solid atoms. Applying a modifier `delete selected
particles <https://ovito.org/manual/particles.modifiers.delete_selected_particles.html>`__
can be applied to delete all the non solid particles. The system after
removing all the liquid atoms is shown below.

.. figure:: system2.png
   :alt: system with only solid

   system with only solid atoms

You can see that there are some solid atoms that are not part of the
large cluster in the centre. Here is where the
clustering functions that pyscal offers comes into play. If you used
``find_clusters`` with ``cluster=True``, the clustering is already done.
We can just get the necessary information out. `Atom.get_cluster <https://pyscal.readthedocs.io/en/latest/pyscal.html#pyscal.core.Atom.get_cluster>`_ can
be used for this. The function gives an array of four integer values. The first value is
1 if the atom is solid, the second value if 1 if the atom is solid and has liquid neighbors,
the third value is 1 if the atom belongs to the largest cluster of solid atoms, and the fourth value gives the number of the cluster that the atom belongs to. We are
interested in the third value which is 1 if the atom belongs to the
largest cluster, 0 otherwise.

.. code:: python

    largest_cluster = [atom.get_cluster()[2] for atom in atoms]

Once again we will save this information to a file and visualise it in
Ovito.

.. code:: python

    ptp.write_structure(sys, 'sys.cluster.dat', customkey='cluster', customvals=largest_cluster)

The system visualised in Ovito following similar steps as above is shown
below.

.. figure:: system3.png
   :alt: system with only largest solid cluster

   system with only the largest cluster of solid atoms

Now we can see that all the solid atoms that do not belong to the largest cluster are removed and only the largest solid cluster is identified. Clustering can be done over any
property. Not just solid atoms. The following example with the same
system will illustrate this.

Clustering based on a custom property
-------------------------------------

The find the clusters based on a custom property, the
`System.clusters_atoms <https://pyscal.readthedocs.io/en/latest/pyscal.html#pyscal.core.System.cluster_atoms>`_ method has to be used. The simulation box
shown above has the centre roughly at (25, 25, 25). For the custom
clustering, we will cluster all atoms within a distance of 10 from the
the rough centre of the box at (25, 25, 25). This is not a really hard
task, but it would give a glimpse into the clustering method. Lets
define a function with checks the above condition.

.. code:: python

    def check_distance(atom):
        #get position of atom
        pos = atom.get_pos()
        #calculate distance from (25, 25, 25)
        dist = ((pos[0]-25)**2 + (pos[1]-25)**2 + (pos[2]-25)**2)**0.5
        #check if dist < 10
        return (dist <= 10)

The above function would return True or False depending on a condition
and takes the Atom as an argument. These are the two important
conditions to be satisfied. Now we can pass this function to cluster.
But first, set up system and find neighbors.

.. code:: python

    sys = pc.System()
    sys.read_inputfile('cluster.dump')
    sys.find_neighbors(method='cutoff', cutoff=3.63)

Now cluster

.. code:: python

    sys.cluster_atoms(check_distance)




.. parsed-literal::

    242



There are 242 atoms in the cluster! Once again we can check this, save
to a file and visualise in Ovito.

.. code:: python

    atoms = sys.get_atoms()
    largest_cluster = [atom.get_cluster()[2] for atom in atoms]

.. code:: python

    ptp.write_structure(sys, 'sys.dist.dat', customkey='cluster', customvals=largest_cluster)

.. figure:: system4.png
   :alt: custom clustering

   atoms at a distance of 10 from (25, 25, 25)

It looks like everything worked. Any atom property, or any property can
be used to cluster the atoms!


.. [1] A. Stukowski, Modelling Simul. Mater. Sci. Eng. 18, 2010
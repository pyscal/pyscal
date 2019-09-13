Classification of atoms as solid or liquid
------------------------------------------

pyscal can also be used to distinguish solid and liquid atoms. The
classification is based on Steinhardt's parameters, specifically the
:math:`q_6` value. The method defines two neighboring atoms :math:`i`
and :math:`j` as having solid bonds if a parameter [1]_,

.. math::  s_{ij} = \sum_{m=-6}^6 q_{6m}(i) q_{6m}^*(j) \geq \mathrm{threshold} 

Additionally, a second order parameter is used to improve the
distinction in solid-liquid boundaries [2]_. This is defined by the criteria,

.. math::  \langle s_{ij} \rangle > \mathrm{avgthreshold} 

If a particle has :math:`n` number of bonds with
:math:`s_{ij} \geq \mathrm{threshold}` and the above condition is also
satisfied, it is considered as a solid. The solid atoms can be clustered
to find the largest solid cluster of atoms. Clustering based on a
different criteria, any criteria is straightforward. Please check the
examples on how to do this.

Finding solid atoms in liquid start with reading in a file and
calculation of neighbors.

.. code:: python

    import pyscal.core as pc
    sys = pc.System()
    sys.read_inputfile('conf.dump')
    sys.get_neighbors(method='cutoff', cutoff=4)

Once again, there are various methods for finding neighbors. Please
check
`here <https://pyscal.readthedocs.io/en/latest/nearestneighbormethods.html#>`__
for details on neighbor calculation methods. Once the neighbors are
calculated, solid atoms can be found directly by,

.. code:: python

    sys.find_solids(bonds=6, threshold=0.5, avgthreshold=0.6, cluster=True)

``bonds`` set the number of minimum bonds a particle should have (as
defined above), ``threshold`` and ``avgthreshold`` are the same
quantities that appear in the equations above. Setting the keyword
``cluster`` to True returns the size of the largest solid cluster. It is
also possible to check if each atom is solid or not.

.. code:: python

    atoms = sys.get_atoms()
    solids = [atom.get_solid() for atom in atoms]

.. [1] Auer, S, Frenkel, D. Adv Polym Sci 173, 2005
.. [2] Bokeloh, J, Rozas, RE, Horbach, J, Wilde, G, Phys. Rev. Lett. 107, 2011
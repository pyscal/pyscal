
Voronoi tessellation to identify local structures
-------------------------------------------------

Voronoi tessellation can be used for identification of local structure by counting the number of faces
of the Voronoi polyhedra of an atom [1]_ [2]_. For each atom a vector :math:` \langle n_3, n_4, n_5, n_6 \rangle` can be calculated
where :math:`n_3` is the number of Voronoi faces of the associated Voronoi polyhedron with three vertices,
:math:`n_4` is with four vertices and so on. Each perfect crystal structure such as a signature vector, for example,
bcc can be identified by :math:`\langle 0 6 0 8 \rangle` and fcc can be identified using :math:`\langle 0 12 0 0 \rangle`.
It is also a useful tool for identifying icosahedral structure which has the fingerprint :math:`\langle 0 0 12 0 \rangle`.
In pyscal, the voronoi vector can be calculated using,

.. code:: python

    import pyscal.core as pc
    sys = pc.System()
    sys.read_inputfile('conf.dump')
    sys.find_neighbors(method='voronoi')
    sys.calculate_vorovector()

The vector for each atom can be accessed using :attr:`~pyscal.catom.Atom.vorovector`. Furthermore, the associated Voronoi
volume of the polyhedron, which may be indicative of the local structure, is also automatically calculated when finding
neighbors using :func:`~pyscal.core.System.find_neighbors`. This value for each atom can be accessed by
:attr:`~pyscal.catom.Atom.volume`. An averaged version of the volume, which is averaged over the neighbors of an atom
can be accessed using :attr:`~pyscal.catom.Atom.avg_volume`.

.. [1] Finney, J. L. PRS 319 ,1970
.. [2] Tanemura et al., PTP 58, 1977

.. note::

    Associated methods-
    :func:`~pyscal.core.System.find_neighbors`
    :func:`~pyscal.core.System.calculate_vorovector`
    :attr:`~pyscal.catom.Atom.vorovector`    
    :attr:`~pyscal.catom.Atom.volume`
    :attr:`~pyscal.catom.Atom.avg_volume`

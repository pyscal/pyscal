# Voronoi tessellation to identify local structures

Voronoi tessellation can be used for identification of local structure
by counting the number of faces of the Voronoi polyhedra of an
atom {cite}`Finney1970,Tanemura1977`. For each atom a vector $\langle n_3~n_4~n_5~n_6 \rangle$
can be calculated where $n_3$ is the number of Voronoi faces of the
associated Voronoi polyhedron with three vertices, $n_4$ is with four
vertices and so on. Each perfect crystal structure such as a signature
vector, for example, bcc can be identified by $\langle 0~6~0~8 \rangle$
and fcc can be identified using $\langle 0~12~0~0 \rangle$. It is also a
useful tool for identifying icosahedral structure which has the
fingerprint $\langle 0~0~12~0 \rangle$. In pyscal, the voronoi vector
can be calculated using,

``` python
import pyscal.core as pc
sys = pc.System()
sys.read_inputfile('conf.dump')
sys.find_neighbors(method='voronoi')
sys.calculate_vorovector()
```

The vector for each atom can be accessed using
[Atom.vorovector](https://docs.pyscal.org/en/latest/pyscal.html#pyscal.catom.Atom.vorovector).
Furthermore, the associated Voronoi volume of the polyhedron, which may
be indicative of the local structure, is also automatically calculated
when finding neighbors using
[System.find_neighbors](https://docs.pyscal.org/en/latest/pyscal.html#pyscal.core.System.find_neighbors).
This value for each atom can be accessed by
[Atom.volume](https://docs.pyscal.org/en/latest/pyscal.html#pyscal.catom.Atom.volume). An averaged
version of the volume, which is averaged over the neighbors of an atom
can be accessed using [Atom.avg_volume](https://docs.pyscal.org/en/latest/pyscal.html#pyscal.catom.Atom.avg_volume).


## References

```{bibliography} ../references.bib
:filter: docname in docnames
:style: unsrt
```
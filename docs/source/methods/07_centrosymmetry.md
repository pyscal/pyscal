# Centrosymmetry parameter

Centrosymmetry parameter (CSP) was introduced by Kelchner et
al. {cite}`Kelchner1998` to identify defects in crystals. The parameter measures the loss of local symmetry. For an atom with $N$ nearest neighbors, the parameter is given
by,

> $$\mathrm{CSP} = \sum_{i=1}^{N/2} \big | \textbf{r}_i + \textbf{r}_{i+N/2} \big |^2$$

$\textbf{r}_i$ and $\textbf{r}_{i+N/2}$ are vectors from the central
atom to two opposite pairs of neighbors. There are two main methods to
identify the opposite pairs of neighbors as described in [this
publication](https://arxiv.org/abs/2003.08879). The first of the
approaches is called Greedy Edge
Selection (GES) {cite}`Stukowski2012`
and is implemented in [LAMMPS](https://lammps.sandia.gov/) and
[Ovito](https://www.ovito.org/). GES algorithm calculates a weight
$w_{ij} = |\textbf{r}_i + \textbf{r}_j|$ for all combinations of
neighbors around an atom and calculates CSP over the smallest $N/2$
weights.

A centrosymmetry parameter calculation using GES algorithm can be
carried out as follows-

``` python
import pyscal.core as pc
sys = pc.System()
sys.read_inputfile('conf.dump')
sys.find_neighbors(method='voronoi')
sys.calculate_centrosymmetry(nmax = 12)
```

`nmax` parameter specifies the number of nearest neighbors to be
considered for the calculation of CSP. The second algorithm is called
the Greedy Vertex
Matching {cite}`Bulatov2006` and is
implemented in [AtomEye](http://li.mit.edu/Archive/Graphics/A/) and
[Atomsk](https://atomsk.univ-lille.fr/). This algorithm orders the
neighbors atoms in order of increasing distance from the central atom.
From this list, the closest neighbor is paired with its lowest weight
partner and both atoms removed from the list. This process is continued
until no more atoms are remaining in the list. CSP calculation using
this algorithm can be carried out by,

``` python
import pyscal.core as pc
sys = pc.System()
sys.read_inputfile('conf.dump')
sys.find_neighbors(method='voronoi')
sys.calculate_centrosymmetry(nmax = 12, algorithm = "gvm")
```

## References

```{bibliography} ../references.bib
:filter: docname in docnames
:style: unsrt
```
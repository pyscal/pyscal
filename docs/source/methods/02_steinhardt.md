# Steinhardt's parameters

Steinhardt\'s bond orientational order parameters {cite}`Steinhardt1983` are a set of
parameters based on [spherical
harmonics](https://en.wikipedia.org/wiki/Spherical_harmonics) to explore
the local atomic environment. These parameters have been used
extensively for various uses such as distinction of crystal structures,
identification of solid and liquid atoms and identification of
defects {cite}`Steinhardt1983`.

These parameters, which are rotationally and translationally invariant
are defined by,

> $$ q_l (i) =  \Big(  \frac{4\pi}{2l+1}  \sum_{m=-l}^l | q_{lm}(i) |^2 \Big )^{\frac{1}{2}} $$

where,

> $$ q_{lm} (i) =  \frac{1}{N(i)} \sum_{j=1}^{N(i)} Y_{lm}(\pmb{r}_{ij}) $$

in which $Y_{lm}$ are the spherical harmonics and $N(i)$ is the number
of neighbours of particle $i$, $\pmb{r}_{ij}$ is the vector connecting
particles $i$ and $j$, and $l$ and $m$ are both intergers with
$m \in [-l,+l]$. Various parameters have found specific uses, such as
$q_2$ and $q_6$ for identification of crystallinity, $q_6$ for
identification of solidity, and $q_4$ and $q_6$ for distinction of
crystal structures {cite}`Mickel2013`. Commonly this method uses a cutoff radius to
identify the neighbors of an atom. The cutoff can be chosen based on
[different methods
available](../part2/intro.md).
Once the cutoff is chosen and neighbors are calculated, the calculation
of Steinhardt\'s parameters is straightforward.

``` python
sys.calculate_q([4, 6])
q = sys.get_qvals([4, 6])
```

## Averaged Steinhardt's parameters

At high temperatures, thermal vibrations affect the atomic positions.
This in turn leads to overlapping distributions of $q_l$ parameters,
which makes the identification of crystal structures difficult. To
address this problem, the averaged version $\bar{q}_l$ of Steinhardt\'s
parameters was introduced by Lechner and Dellago {cite}`Lechner2008`. $\bar{q}_l$ is
given by,

> $$\bar{q}_l (i) =  \Big(  \frac{4\pi}{2l+1}  \sum_{m=-l}^l \Big| \frac{1}{\tilde{N}(i)} \sum_{k=0}^{\tilde{N}(i)} q_{lm}(k) \Big|^2 \Big )^{\frac{1}{2}}$$

where the sum from $k=0$ to $\tilde{N}(i)$ is over all the neighbors and
the particle itself. The averaged parameters takes into account the
first neighbor shell and also information from the neighboring atoms and
thus reduces the overlap between the distributions. Commonly $\bar{q}_4$
and $\bar{q}_6$ are used in identification of crystal structures.
Averaged versions can be calculated by setting the keyword
`averaged=True` as follows.

``` python
sys.calculate_q([4, 6], averaged=True)
q = sys.get_qvals([4, 6], averaged=True)
```

## Voronoi weighted Steinhardt's parameters

In order to improve the resolution of crystal structures Mickel et
al {cite}`Mickel2013` proposed weighting the contribution of each neighbor to the
Steinhardt parameters by the ratio of the area of the Voronoi facet
shared between the neighbor and host atom. The weighted parameters are
given by,

> $$q_{lm} (i) =  \frac{1}{N(i)} \sum_{j=1}^{N(i)} \frac{A_{ij}}{A} Y_{lm}(\pmb{r}_{ij})$$

where $A_{ij}$ is the area of the Voronoi facet between atoms $i$ and
$j$ and $A$ is the sum of the face areas of atom $i$. In pyscal, the
area weights are already assigned during the neighbor calculation phase
when the Voronoi method is used to calculate neighbors in the
[System.find_neighbors](https://docs.pyscal.org/en/latest/pyscal.html#pyscal.core.System.find_neighbors). The
Voronoi weighted Steinhardt\'s parameters can be calculated as follows,

``` python
sys.find_neighbors(method='voronoi')
sys.calculate_q([4, 6])
q = sys.get_qvals([4, 6])
```

The weighted Steinhardt\'s parameters can also be averaged as described
above. Once again, the keyword `averaged=True` can be used for this
purpose.

``` python
sys.find_neighbors(method='voronoi')
sys.calculate_q([4, 6], averaged=True)
q = sys.get_qvals([4, 6], averaged=True)
```

It was also proposed that higher powers of the weight {cite}`Haeberle2019`
$\frac{A_{ij}^{\alpha}}{A(\alpha)}$ where $\alpha = 2, 3$ can also be
used, where $A(\alpha) = \sum_{j=1}^{N(i)} A_{ij}^{\alpha}$ The value of
this can be set using the keyword `voroexp` during the neighbor
calculation phase.

``` python
sys.find_neighbors(method='voronoi', voroexp=2)
```

If the value of `voroexp` is set to 0, the neighbors would be found
using Voronoi method, but the calculated Steinhardt\'s parameters will
not be weighted.


## References

```{bibliography} ../references.bib
:filter: docname in docnames
:style: unsrt
```
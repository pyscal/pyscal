# Classification of atoms as solid or liquid

pyscal can also be used to distinguish solid and liquid atoms. The
classification is based on [Steinhardt\'s
parameters](https://pyscal.readthedocs.io/en/latest/steinhardtparameters.html),
specifically $q_6$. The method defines two neighboring atoms $i$ and $j$
as having solid bonds if a parameter $s_{ij}$ {cite}`Auer2005`,

> $$s_{ij} = \sum_{m=-6}^6 q_{6m}(i) q_{6m}^*(j) \geq \mathrm{threshold}$$

Additionally, a second order parameter is used to improve the
distinction in solid-liquid boundaries {cite}`Bokeloh2014`. This is defined by the
criteria,

> $$\langle s_{ij} \rangle > \mathrm{avgthreshold}$$

If a particle has $n$ number of bonds with
$s_{ij} \geq \mathrm{threshold}$ and the above condition is also
satisfied, it is considered as a solid. The solid atoms can be clustered
to find the largest solid cluster of atoms. 

Finding solid atoms in liquid start with reading in a file and
calculation of neighbors.

``` python
import pyscal.core as pc
sys = pc.System()
sys.read_inputfile('conf.dump')
sys.find_neighbors(method='cutoff', cutoff=4)
```

Once again, there are various methods for finding neighbors. Please
check
[here](../part2/intro.md)
for details on neighbor calculation methods. Once the neighbors are
calculated, solid atoms can be found directly by,

``` python
sys.find_solids(bonds=6, threshold=0.5, avgthreshold=0.6, cluster=True)
```

`bonds` set the number of minimum bonds a particle should have (as
defined above), `threshold` and `avgthreshold` are the same quantities
that appear in the equations above. Setting the keyword `cluster` to
True returns the size of the largest solid cluster. It is also possible
to check if each atom is solid or not.

``` python
atoms = sys.atom
solids = [atom.solid for atom in atoms]
```

## References

```{bibliography} ../references.bib
:filter: docname in docnames
:style: unsrt
```
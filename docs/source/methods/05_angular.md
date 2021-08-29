# Angular parameters

## Angular criteria for identification of diamond structure

Angular parameter introduced by Uttormark et al {cite}`Uttormark1993` is used to measure
the tetrahedrality of local atomic structure. An atom belonging to
diamond structure has four nearest neighbors which gives rise to six
three body angles around the atom. The angular parameter $A$ is then
defined as,

> $$ A = \sum_{i=1}^6 (\cos(\theta_i)+\frac{1}{3})^2$$

An atom belonging to diamond structure would show the value of angular
params close to 0. Angular parameter can be calculated in pyscal using
the following method -

``` python
import pyscal.core as pc
sys = pc.System()
sys.read_inputfile('conf.dump')
sys.find_neighbors(method='cutoff', cutoff='adaptive')
sys.calculate_angularcriteria()
```

The calculated angular criteria value can be accessed for each atom
using [Atom.angular](https://docs.pyscal.org/en/latest/pyscal.html#pyscal.catom.Atom.angular).

## $\chi$ parameters for structural identification

$\chi$ parameters introduced by Ackland and Jones {cite}`Ackland2006` measures all local
angles created by an atom with its neighbors and creates a histogram of
these angles to produce vector which can be used to identify structures.
After finding the neighbors of an atom, $\cos \theta_{ijk}$ for atoms j
and k which are neighbors of i is calculated for all combinations of i,
j and k. The set of all calculated cosine values are then added to a
histogram with the following bins - \[-1.0, -0.945, -0.915, -0.755,
-0.705, -0.195, 0.195, 0.245, 0.795, 1.0\]. Compared to $\chi$
parameters from $\chi_0$ to $\chi_7$ in the associated publication, the
vector calculated in pyscal contains values from $\chi_0$ to $\chi_8$
which is due to an additional $\chi$ parameter which measures the number
of neighbors between cosines -0.705 to -0.195. The $\chi$ vector is
characteristic of the local atomic environment and can be used to
identify crystal structures, details of which can be found in the
publication[^2].

$\chi$ parameters can be calculated in pyscal using,

``` python
import pyscal.core as pc
sys = pc.System()
sys.read_inputfile('conf.dump')
sys.find_neighbors(method='cutoff', cutoff='adaptive')
sys.calculate_chiparams()
```

The calculated values for each atom can be accessed using
[Atom.chiparams](https://docs.pyscal.org/en/latest/pyscal.html#pyscal.catom.Atom.chiparams).

## References

```{bibliography} ../references.bib
:filter: docname in docnames
:style: unsrt
```
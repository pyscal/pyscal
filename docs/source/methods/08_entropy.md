# Entropy - Enthalpy parameters

## Entropy fingerprint

The entropy parameter was introduced by Piaggi et
al {cite}`Piaggi2017` for identification of defects and
distinction between solid and liquid. The entropy paramater $s_s^i$ is
defined as,

> $$s_s^i = -2\pi\rho k_B \int_0^{r_m} [g_m^i(r)\ln g_m^i(r) - g_m^i(r) + 1] r^2 dr$$

where $r_m$ is the upper bound of integration and $g_m^i$ is radial
distribution function centered on atom $i$,

> $$g_m^i(r) = \frac{1}{4\pi\rho r^2} \sum_j \frac{1}{\sqrt{2\pi\sigma^2}} \exp{-(r-r_{ij})^2/(2\sigma^2)}$$

$r_{ij}$ is the interatomic distance between atom $i$ and its neighbors
$j$ and $\sigma$ is a broadening parameter.

The averaged version of entropy parameters $\bar{s}_s^i$ can be
calculated in two ways, either using a simple averaging over the
neighbors given by,

> $$\bar{s}_s^i = \frac{\sum_j s_s^j + s_s^i}{N + 1}$$

or using a switching function as described below,

> $$\bar{s}_s^i = \frac{\sum_j s_s^i f(r_{ij}) + s_s^i}{\sum_j f(r_{ij}) + 1}$$

$f(r_{ij})$ is a switching parameter which depends on $r_a$ which is the
cutoff distance. The switching function shows a value of 1 for
$r_{ij} << r_a$ and 0 for $r_{ij} >> r_a$. The switching function is
given by,

> $$f(r_{ij}) = \frac{1-(r_{ij}/r_a)^N}{1-(r_{ij}/r_a)^M}$$

Entropy parameters can be calculated in pyscal using the following code,

``` python
import pyscal.core as pc
sys = pc.System()
sys.read_inputfile('conf.dump')
sys.find_neighbors(method="cutoff", cutoff=0)
lattice_constant=4.00
sys.calculate_entropy(1.4*lattice_constant, averaged=True)
atoms = sys.atoms
entropy = [atom.entropy for atom in atoms]
average_entropy = [atom.avg_entropy for atom in atoms]
```

The value of $r_m$ is provided in units of lattice constant. Further
parameters shown above, such as $\sigma$ can be specified using the
various keyword arguments. The above code does a simple averaging over
neighbors. The switching function can be used by,

``` python
sys.calculate_entropy(1.4*lattice_constant, ra=0.9*lattice_constant, switching_function=True, averaged=True)
```

In pyscal, a slightly different version of $s_s^i$ is calculated. This
is given by,

> $$s_s^i = -\rho \int_0^{r_m} [g_m^i(r)\ln g_m^i(r) - g_m^i(r) + 1] r^2 dr$$

The prefactor $2\pi k_B$ is dropped in the entropy values calculated in
pyscal.

## References

```{bibliography} ../references.bib
:filter: docname in docnames
:style: unsrt
```
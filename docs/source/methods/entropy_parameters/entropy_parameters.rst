Entropy fingerprint
===================

The entropy parameter was introduced by `Piaggi et
al <https://doi.org/10.1063/1.4998408>`__ for identification of defects
and distinction between solid and liquid. The entropy paramater
:math:`s_s^i` is defined as,

.. math::


   s_s^i = -2\pi\rho k_B \int_0^{r_m} [g_m^i(r)\ln g_m^i(r) - g_m^i(r) + 1] r^2 dr

where :math:`r_m` is the upper bound of integration and :math:`g_m^i` is
radial distribution function centered on atom :math:`i`,

.. math::


   g_m^i(r) = \frac{1}{4\pi\rho r^2} \sum_j \frac{1}{\sqrt{2\pi\sigma^2}} \exp{-(r-r_{ij})^2/(2\sigma^2)}

:math:`r_{ij}` is the interatomic distance between atom :math:`i` and
its neighbors :math:`j` and :math:`\sigma` is a broadening parameter.

The averaged version of entropy parameters :math:`\bar{s}_s^i` can also
be calculated using the following expression,

.. math::


   \bar{s}_s^i = \frac{\sum_j s_s^i f(r_{ij}) + s_s^i}{\sum_j f(r_{ij}) + 1}

:math:`f(r_{ij})` is a switching parameter which depends on :math:`r_a`
which is the cutoff distance. The switching function shows a value of 1
for :math:`r_{ij} << r_a` and 0 for :math:`r_{ij} >> r_a`. The switching
function is given by,

.. math::


   f(r_{ij}) = \frac{1-(r_{ij}/r_a)^N}{1-(r_{ij}/r_a)^M}

Entropy parameters can be calculated in pyscal using the following code,

.. code:: python

    import pyscal.core as pc
    sys = pc.System()
    sys.read_inputfile('conf.dump')
    sys.find_neighbors(method="cutoff", cutoff=0)
    lattice_constant=4.00
    sys.calculate_entropy(1.4*lattice_constant, ra=0.9*lattice_constant, averaged=True)
    atoms = sys.atoms
    entropy = [atom.entropy for atom in atoms]
    average_entropy = [atom.avg_entropy for atom in atoms]

The values of :math:`r_m` and :math:`r_a` are provided in units of
lattice constant. Further parameters shown above, such as :math:`\sigma`
can be specified using the various keyword arguments.

In pyscal, a slightly different version of :math:`s_s^i` is calculated.
This is given by,

.. math::


   s_s^i = -\rho \int_0^{r_m} [g_m^i(r)\ln g_m^i(r) - g_m^i(r) + 1] r^2 dr

The prefactor :math:`2\pi k_B` is dropped in the entropy values
calculated in pyscal.
